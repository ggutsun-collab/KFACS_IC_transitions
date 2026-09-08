###############################################################################
## 15_ic13_headtohead.R   (Phase 3.1-3.2 — novelty 방어 분석)   v260903
##
## 리뷰어의 핵심 반론: "IC 와 Fried 가 5개 중 4개 지표를 공유하므로 IC 가 robust 층을
## 층화한다는 것은 grip·gait 를 연속형으로 쓴 것과 다르지 않다."
## 이 스크립트는 그 반론에 두 가지로 답합니다.
##   A. IC-13 (HGS, GS_ms, EXH, loss_of_Bwt 제외) 로 18개 전이를 재추정
##      -> IC-17 과 나란히. 공유 지표 없이도 연관이 유지되는가.
##   B. Robust 층(wave 1) 5년 사망·복합악화 head-to-head
##      M0 base / M1 +grip+gait(연속) / M2 +IC-17 / M3 +IC-13 / M4 +grip+gait+IC-13 / M5 +grip+gait+IC-17
##      -> Harrell C, ΔC(부트스트랩), 5년 AUC(timeROC), IDI, IPCW 순이익(DCA)
##      전체 코호트(+frailty phenotype 를 base 에 포함)도 참고로 산출.
## 모든 추정은 m=20 MI 세트에서 반복 후 Rubin 풀링.
##
## 입력 : KFACS_mi_pmm_m20.rds (13_), pre-imp 원본 (IC-13 anchor 적합용)
## 출력 : ST_IC13_Transitions, ST_HeadToHead_Robust, ST_HeadToHead_All, KFACS_IC13_anchor_fiml.rds
## 실행 : source("15_ic13_headtohead.R", encoding = "UTF-8")   (B=200 기준 20-40분)
###############################################################################

H2H_VERSION <- "v260903"
message("\n=== 15_ic13_headtohead ", H2H_VERSION, " ===")

## ── 헬퍼 적재 ──────────────────────────────────────────────────────────
## 위치 탐색과 헬퍼 로드는 R/_bootstrap.R 한 곳에서만 정의합니다.
## 모든 경로(PROJECT_ROOT / IMP_DIR / FIG_DIR)는 R/00_setup.R 에서만 정합니다.
local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of)) { d <- dirname(of); break }
  }
  cand <- c(if (!is.null(d)) file.path(d, c("R", "../R", "../../R", ".", "..")),
            "R", "../R", "../../R", ".", "..")
  for (p in cand)
    if (file.exists(file.path(p, "_bootstrap.R"))) {
      source(file.path(p, "_bootstrap.R"), encoding = "UTF-8"); break }
})
kf_init()
.need(c("survival", "lavaan", "readxl", "openxlsx", "timeROC"))
HAS_TROC <- requireNamespace("timeROC", quietly = TRUE)
suppressPackageStartupMessages({ library(survival); library(lavaan) })

###############################################################################
## 1. 설정
###############################################################################
MI_STEM  <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
PRE_FILE <- get0("PRE_FILE", ifnotfound = file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260618_state.xlsx"))  # 00_setup.R 에서 정함
TAU      <- 5
B_BOOT   <- 200          # ΔC 부트스트랩 (세트당). 빠른 점검은 50
NB_THR   <- c(0.05, 0.10, 0.15, 0.20)
MIN_M_OK <- 0.8
SEED     <- 20260903
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
OVERLAP <- c("HGS", "GS_ms", "EXH", "loss_of_Bwt")
IND13 <- setdiff(IND17, OVERLAP)
BIF13 <- '
  g    =~ Balance+rev_CST+Appetite+rev_logMAR+rev_PTA+Orientation+Memory+
          Attention_Calculation+Language+Visuospatial+Negative_affect+Positive_affect+Motivation
  loco =~ Balance+rev_CST
  cogn =~ Orientation+Memory+Attention_Calculation+Language+Visuospatial
  psyc =~ Negative_affect+Positive_affect+Motivation
  g ~~ 0*loco + 0*cogn + 0*psyc '
BIF13_ALT <- sub("\n  loco =~ Balance\\+rev_CST", "", BIF13)   # loco 특정요인 불안정 시
BIF13_ALT <- sub("0\\*loco \\+ ", "", BIF13_ALT)

MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))); SETS <- names(MI); M <- length(SETS)
set.seed(SEED)

###############################################################################
## 2. IC-13 anchor (pre-imp 관측 자료, FIML) — 13_ 의 IC-17 anchor 와 같은 방식
###############################################################################
tmp <- file.path(tempdir(), "pre.xlsx"); file.copy(PRE_FILE, tmp, overwrite = TRUE)
pre <- as.data.frame(readxl::read_excel(tmp, sheet = "전체_long"))
pre$id <- as.character(pre$id); pre$wave <- as.integer(as_num(pre$wave))
pre <- pre[!pre$id %in% c("kf161224", "kf170602", "kf171295"), ]
dw <- tapply(pre$death_wave, pre$id, function(x) suppressWarnings(as.numeric(x[1])))
pre <- pre[is.na(dw[pre$id]) | pre$wave < dw[pre$id], ]
for (v in IND13) pre[[v]] <- as_num(pre[[v]])
pre <- pre[rowSums(!is.na(pre[, IND13])) > 0, ]
MU13 <- vapply(pre[IND13], mean, numeric(1), na.rm = TRUE); SD13 <- vapply(pre[IND13], sd, numeric(1), na.rm = TRUE)
Z13 <- pre[, IND13]; for (v in IND13) Z13[[v]] <- (pre[[v]] - MU13[[v]]) / SD13[[v]]
fit13 <- tryCatch(cfa(BIF13, data = Z13, std.lv = TRUE, estimator = "MLR", missing = "fiml"), error = function(e) NULL)
if (is.null(fit13) || !lavInspect(fit13, "converged") || any(diag(lavInspect(fit13, "est")$theta) < 0)) {
  message("  IC-13: loco 2지표 특정요인 불안정 -> general 에만 적재")
  fit13 <- cfa(BIF13_ALT, data = Z13, std.lv = TRUE, estimator = "MLR", missing = "fiml")
}
cat(sprintf("  IC-13 anchor: CFI=%.3f RMSEA=%.3f (n=%d person-waves)\n",
            fitMeasures(fit13, "cfi.robust"), fitMeasures(fit13, "rmsea.robust"), nrow(Z13)))
st <- standardizedSolution(fit13); st <- st[st$op == "=~" & st$lhs == "g", ]
SIGN13 <- if (mean(st$est.std) < 0) -1 else 1
saveRDS(list(fit = fit13, mu = MU13, sd = SD13, sign = SIGN13, indicators = IND13),
        file.path(IMP_DIR, "KFACS_IC13_anchor_fiml.rds"))

## ★ 회귀법 요인점수를 모수로 직접 계산 (lavPredict 의 FIML 케이스별 경로 회피 — 수학적으로 동일)
##   eta_hat = Psi Lambda' (Lambda Psi Lambda' + Theta)^-1 (y - nu),  std.lv 이므로 Psi = I
.fs_weights <- function(fit) {
  E <- lavInspect(fit, "est"); L <- E$lambda; Th <- E$theta; Ps <- E$psi
  nu <- as.numeric(E$nu); names(nu) <- rownames(L)
  S  <- L %*% Ps %*% t(L) + Th
  W  <- Ps %*% t(L) %*% solve(S)                     # 요인 x 지표
  list(W = W[rownames(W) == "g", , drop = FALSE], nu = nu, items = rownames(L))
}
FS13 <- .fs_weights(fit13)
score13 <- function(d) {
  ok <- stats::complete.cases(d[, IND13]); Z <- as.matrix(d[ok, FS13$items])
  for (v in FS13$items) Z[, v] <- (as_num(Z[, v]) - MU13[[v]]) / SD13[[v]]
  Zc <- sweep(Z, 2, FS13$nu[FS13$items])
  s <- rep(NA_real_, nrow(d)); s[ok] <- SIGN13 * as.numeric(Zc %*% t(FS13$W)); s
}
## 검증: lavPredict 와 일치하는지 표본 200행으로 확인 (느린 경로는 200행이면 충분히 빠름)
{
  ok <- stats::complete.cases(Z13); Zs <- Z13[which(ok)[seq_len(min(200, sum(ok)))], ]
  a <- as.numeric(lavPredict(fit13, newdata = Zs)[, "g"])
  b <- as.numeric(sweep(as.matrix(Zs), 2, FS13$nu[FS13$items]) %*% t(FS13$W))
  cat(sprintf("  요인점수 직접계산 vs lavPredict: r = %.6f, 최대 절대차 = %.2e\n", cor(a, b), max(abs(a - b))))
}
## within-sex W1 z (14_ 의 apply_scale 과 동일 규칙), 임의의 점수 열에 적용
ws_z <- function(d, col) {
  w1 <- d[!duplicated(d$id), ]; z <- d[[col]]
  for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
    z[k] <- (d[[col]][k] - mean(w1[[col]][k1], na.rm = TRUE)) / sd(w1[[col]][k1], na.rm = TRUE) }
  z
}

###############################################################################
## 3. 보조 함수 — 판별력 지표
###############################################################################
ADJ_WS <- setdiff(CFG$ADJ, "sex_f")         # within-sex 점수 모형: 성별은 척도에 흡수
risk5 <- function(cx, nd) {                 # 5년 예측 위험
  bh <- basehaz(cx, centered = FALSE); H <- bh$hazard[max(which(bh$time <= TAU))]
  1 - exp(-H * exp(predict(cx, newdata = nd, type = "lp", reference = "zero")))
}
cidx <- function(cx) { cc <- concordance(cx); c(C = cc$concordance, se = sqrt(cc$var)) }
ipcw_nb <- function(ft, ev, p, thr) {       # IPCW 순이익 at TAU (Vickers; 중도절단 가중)
  ## 중도절단 분포: 5년 도달(행정적 절단)은 중도절단이 아니므로 제외
  cens <- as.integer(ev == 0 & ft < TAU)
  G  <- survfit(Surv(ft, cens) ~ 1)
  Gf <- stepfun(G$time, c(1, G$surv))
  tt <- pmin(ft, TAU); d5 <- as.integer(ev == 1 & ft <= TAU)
  known <- (ev == 1 & ft <= TAU) | ft >= TAU
  g <- pmax(Gf(tt - 1e-8), 0.05)             # 직전 시점의 G; 극단 가중 방지
  w <- ifelse(known, 1 / g, 0); n <- length(ft)
  vapply(thr, function(pt) { pos <- p >= pt
    (sum(w * pos * d5) - sum(w * pos * (1 - d5)) * pt / (1 - pt)) / n }, numeric(1))
}
idi <- function(ft, ev, p) {                # 5년 사건/비사건 평균 위험 차 (판별 기울기)
  d5 <- ev == 1 & ft <= TAU; a5 <- ft >= TAU
  mean(p[d5]) - mean(p[a5])
}
auc5 <- function(ft, ev, lp) {
  if (!HAS_TROC) return(c(AUC = NA, se = NA))
  r <- try(timeROC::timeROC(T = ft, delta = ev, marker = lp, cause = 1, times = c(TAU - 1e-4, TAU), iid = TRUE), silent = TRUE)
  if (inherits(r, "try-error")) return(c(AUC = NA, se = NA))
  k <- length(r$AUC)
  c(AUC = unname(r$AUC[k]), se = unname(r$inference$vect_sd_1[k]))
}

## Table 3 형 사건시각 (14_ 와 동일)
mk_events <- function(dat) {
  dat <- dat[order(dat$id, dat$time), ]; sp <- split(seq_len(nrow(dat)), dat$id)
  do.call(rbind, lapply(sp, function(ix) {
    s <- dat$state[ix]; tm <- dat$time[ix]; ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
    n <- length(tm); if (n < 1 || !is.finite(s[1])) return(NULL)
    s0 <- s[1]; tmax <- max(tm)
    fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] > s0) else NA_real_
    kd <- which(is.finite(s) & s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
    data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax, t_wor = tw, t_dth = td)
  }))
}
mk_surv <- function(E, out) {
  inf <- function(x) ifelse(is.finite(x), x, Inf)
  if (out == "wor") { ft <- pmin(inf(E$t_wor), E$tmax, TAU); ev <- as.integer(is.finite(E$t_wor) & E$t_wor <= ft) }
  else { ft <- pmin(inf(E$t_dth), E$tmax, TAU); ev <- as.integer(is.finite(E$t_dth) & E$t_dth <= ft) }
  data.frame(id = E$id, s0 = E$s0, ft = pmax(ft, 1e-6), ev = ev)
}

## ★ robust 층은 정의상 Fried 기준 0개(chs_total = 0 상수)이므로 CHS 연속점수는 비교 대상이 될 수 없습니다.
##   "cut-point 까지의 거리" 반론의 실체는 grip·gait 의 연속값이므로 그것을 비교 모형으로 씁니다.
MODELS <- list(M0 = character(0), M1 = c("hgs_z", "gs_z"), M2 = "IC17_z", M3 = "IC13_z",
               M4 = c("hgs_z", "gs_z", "IC13_z"), M5 = c("hgs_z", "gs_z", "IC17_z"))
MLAB <- c(M0 = "Base", M1 = "+ grip strength + gait speed (continuous)", M2 = "+ IC-17",
          M3 = "+ IC-13 (no Fried items)", M4 = "+ grip + gait + IC-13", M5 = "+ grip + gait + IC-17")

## 한 세트·한 모집단·한 결과에 대해 5개 모형의 지표
h2h_block <- function(S, base_adj, pop, out, set) {
  fml <- function(x) as.formula(paste("Surv(ft, ev) ~", paste(c(base_adj, x), collapse = " + ")))
  fits <- lapply(MODELS, function(x) coxph(fml(x), data = S))
  C  <- t(vapply(fits, cidx, numeric(2)))
  P  <- lapply(fits, function(f) risk5(f, S)); LP <- lapply(fits, function(f) predict(f, type = "lp"))
  ## ΔC 부트스트랩 (참가자 재표집, 5 모형 동시)
  dC <- matrix(NA_real_, B_BOOT, length(MODELS)); n <- nrow(S)
  for (b in seq_len(B_BOOT)) {
    ix <- sample.int(n, n, replace = TRUE); Sb <- S[ix, ]
    cb <- vapply(MODELS, function(x) { f <- try(coxph(fml(x), data = Sb), silent = TRUE)
      if (inherits(f, "try-error")) NA_real_ else concordance(f)$concordance }, numeric(1))
    dC[b, ] <- cb - cb[1]
  }
  do.call(rbind, lapply(seq_along(MODELS), function(j) {
    a <- auc5(S$ft, S$ev, LP[[j]]); nb <- ipcw_nb(S$ft, S$ev, P[[j]], NB_THR)
    data.frame(set = set, pop = pop, outcome = out, model = names(MODELS)[j],
               n = n, events = sum(S$ev), C = C[j, 1], C_se = C[j, 2],
               dC = mean(dC[, j], na.rm = TRUE), dC_se = sd(dC[, j], na.rm = TRUE),
               AUC5 = a[1], AUC5_se = a[2], IDI = idi(S$ft, S$ev, P[[j]]) - idi(S$ft, S$ev, P[[1]]),
               NB05 = nb[1], NB10 = nb[2], NB15 = nb[3], NB20 = nb[4], stringsAsFactors = FALSE)
  }))
}

###############################################################################
## 4. 세트별 실행
###############################################################################
TR <- list(); HH <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, M))
  load_imputed(set, stem = MI_STEM); raw <- load_kfacs(); d <- prep_long(raw)
  d$IC13 <- score13(d)
  d$IC17_z <- ws_z(d, "gLIC"); d$IC13_z <- ws_z(d, "IC13")
  d$HGS_num <- as_num(d$HGS); d$GS_num <- as_num(d$GS_ms)
  d$hgs_z <- ws_z(d, "HGS_num"); d$gs_z <- ws_z(d, "GS_num")
  d$chs_z  <- { x <- as_num(d$chs_total); w1 <- d[!duplicated(d$id), ]
                (x - mean(as_num(w1$chs_total), na.rm = TRUE)) / sd(as_num(w1$chs_total), na.rm = TRUE) }
  ## ── A. 18 전이, IC-17 vs IC-13 (within-sex, per s.d.) ──
  for (ex in c("IC17_z", "IC13_z")) {
    dd <- d; dd$gLIC_z <- dd[[ex]]                    # build_intervals 는 gLIC_z 만 실어 나릅니다
    for (sv in c("state", "frailty_3cat")) {
      iv <- build_intervals(dd, sv, "rolling")
      r <- run_transition_table(iv, exposure_terms = "gLIC_z", adj = ADJ_WS)
      r$exposure <- ex; r$set <- set; r$system <- sv
      TR[[length(TR) + 1]] <- r[, c("set", "exposure", "system", "transition", "IRR", "lcl", "ucl", "p", "n_events", "pyears")]
    }
  }
  ## ── B. head-to-head ──
  b <- get_baseline(d); E <- mk_events(d)
  b$frail_f <- factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail"))
  for (out in c("dth", "wor")) {
    S0 <- merge(b, mk_surv(E, out), by = "id"); S0 <- S0[S0$s0 < 4, ]
    ## robust 층: base = age, comorbid, edu, income, area (+sex: CHS 는 성별 기준 컷이 이미 반영됨)
    Sr <- S0[S0$frail_f == "Robust", ]
    HH[[length(HH) + 1]] <- h2h_block(Sr, CFG$ADJ, "Robust stratum", out, set)
    ## 전체 코호트: base 에 frailty phenotype 포함
    HH[[length(HH) + 1]] <- h2h_block(S0, c(CFG$ADJ, "frail_f"), "Whole cohort", out, set)
  }
}
TR <- do.call(rbind, TR); HH <- do.call(rbind, HH)
save_vals(TR, "ST_IC13_Transitions_perset.csv", FIG_DIR); save_vals(HH, "ST_HeadToHead_perset.csv", FIG_DIR)

###############################################################################
## 5. Rubin 풀링
###############################################################################
pool <- function(b, se) { ok <- is.finite(b) & is.finite(se)
  if (sum(ok) < MIN_M_OK * M) return(c(est = NA, se = NA, lcl = NA, ucl = NA, p = NA, m = sum(ok)))
  rubin_pool(b[ok], se[ok]) }
fmt_irr <- function(e, l, u) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(u)), "NE")
fmt_num <- function(e, l, u, d = 3) ifelse(is.finite(e), sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), e, l, u), "-")
fmt_p <- function(p) ifelse(!is.finite(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## A
TR$b <- log(TR$IRR); TR$se <- (log(TR$ucl) - log(TR$lcl)) / (2 * qnorm(0.975))
TR$key <- paste(TR$exposure, TR$system, TR$transition, sep = "|")
PA <- do.call(rbind, lapply(split(TR, TR$key), function(z) { r <- pool(z$b, z$se)
  data.frame(exposure = z$exposure[1], system = z$system[1], transition = z$transition[1],
             events = round(median(z$n_events)), est = r["est"], lcl = r["lcl"], ucl = r["ucl"], p = r["p"]) }))
keys <- unique(PA[, c("system", "transition")])
TA <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  a <- PA[PA$exposure == "IC17_z" & PA$system == keys$system[k] & PA$transition == keys$transition[k], ]
  b <- PA[PA$exposure == "IC13_z" & PA$system == keys$system[k] & PA$transition == keys$transition[k], ]
  ret <- if (is.finite(a$est) && is.finite(b$est) && a$p < 0.05 && abs(a$est) > log(1.05)) 100 * b$est / a$est else NA
  data.frame(System = ifelse(keys$system[k] == "state", "Disability, 3 states", "Frailty phenotype"),
             Transition = keys$transition[k], Events = a$events,
             `IC-17, IRR per +1 s.d. (95% CI)` = fmt_irr(a$est, a$lcl, a$ucl), `P` = fmt_p(a$p),
             `IC-13, IRR per +1 s.d. (95% CI)` = fmt_irr(b$est, b$lcl, b$ucl), `P ` = fmt_p(b$p),
             `Effect retained, %` = ifelse(is.finite(ret), sprintf("%.0f", ret), "-"),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
TA <- TA[order(match(TA$System, c("Disability, 3 states", "Frailty phenotype"))), ]
save_table(TA, "ST_IC13_Transitions",
  title = "Supplementary Table | State-transition rates per +1 s.d. of intrinsic capacity with (IC-17) and without (IC-13) the four indicators shared with the frailty phenotype",
  footnotes = c(
    "IC-13 removes grip strength, gait speed, exhaustion and weight loss, the four indicators that also define the Fried phenotype; its general factor was estimated once on the observed data by full-information maximum likelihood and applied unchanged to every completed dataset. Both scores are standardized within sex on the wave-1 distribution.",
    sprintf("Incidence rate ratios from transition-specific Poisson models with a log(person-time) offset and participant-clustered standard errors, adjusted for age, baseline comorbidity count, education, household income and area of residence, fitted in each of %d completed datasets and combined by Rubin's rules.", M),
    "Effect retained is the IC-13 log rate ratio as a percentage of the IC-17 log rate ratio, shown where the IC-17 estimate was significant and non-negligible."))

## B
HH$key <- paste(HH$pop, HH$outcome, HH$model, sep = "|")
PB <- do.call(rbind, lapply(split(HH, HH$key), function(z) {
  c1 <- pool(z$C, z$C_se); dc <- pool(z$dC, z$dC_se); au <- pool(z$AUC5, z$AUC5_se)
  data.frame(pop = z$pop[1], outcome = z$outcome[1], model = z$model[1],
             n = round(median(z$n)), events = round(median(z$events)),
             C = c1["est"], C_l = c1["lcl"], C_u = c1["ucl"],
             dC = dc["est"], dC_l = dc["lcl"], dC_u = dc["ucl"], dC_p = dc["p"],
             AUC = au["est"], AUC_l = au["lcl"], AUC_u = au["ucl"],
             IDI = mean(z$IDI), NB05 = mean(z$NB05), NB10 = mean(z$NB10), NB15 = mean(z$NB15), NB20 = mean(z$NB20)) }))
OUTL <- c(dth = "5-year mortality", wor = "5-year composite worsening")
TB <- data.frame(Population = PB$pop, Outcome = OUTL[PB$outcome], Model = MLAB[PB$model], n = PB$n, Events = PB$events,
                 `Harrell C (95% CI)` = fmt_num(PB$C, PB$C_l, PB$C_u),
                 `Delta C vs base (95% CI)` = ifelse(PB$model == "M0", "reference", fmt_num(PB$dC, PB$dC_l, PB$dC_u)),
                 `P` = ifelse(PB$model == "M0", "-", fmt_p(PB$dC_p)),
                 `AUC at 5 years (95% CI)` = fmt_num(PB$AUC, PB$AUC_l, PB$AUC_u),
                 `IDI vs base` = ifelse(PB$model == "M0", "reference", sprintf("%+.4f", PB$IDI)),
                 `Net benefit, threshold 5%` = sprintf("%.4f", PB$NB05), `10%` = sprintf("%.4f", PB$NB10),
                 `15%` = sprintf("%.4f", PB$NB15), `20%` = sprintf("%.4f", PB$NB20),
                 check.names = FALSE, stringsAsFactors = FALSE)
TB <- TB[order(TB$Population, TB$Outcome, match(TB$Model, MLAB)), ]
for (pp in unique(TB$Population)) {
  save_table(TB[TB$Population == pp, -1], if (pp == "Robust stratum") "ST_HeadToHead_Robust" else "ST_HeadToHead_All",
    title = sprintf("Supplementary Table | Incremental prognostic value of intrinsic capacity beyond the frailty phenotype: %s", tolower(pp)),
    footnotes = c(
      if (pp == "Robust stratum") "Cox models among participants classified as robust (0 Fried criteria) at wave 1. The base model contains age, sex, baseline comorbidity count, education, household income and area of residence."
      else "Cox models in the whole analytic cohort at wave 1. The base model contains age, sex, baseline comorbidity count, education, household income, area of residence and the three-level frailty phenotype.",
      "Grip strength and gait speed are entered as continuous wave-1 values standardized within sex, so that the comparison tests whether intrinsic capacity adds information beyond the continuous form of the phenotype's own performance items (the distance to the frailty cut-points). Within the robust stratum the CHS count is zero by definition and cannot serve as a comparator. IC-17 and IC-13 are the general-factor scores standardized within sex; IC-13 omits the four indicators shared with the phenotype.",
      sprintf("Delta C is the change in Harrell's C relative to the base model, with confidence intervals from %d participant-level bootstrap resamples in each completed dataset combined by Rubin's rules across %d datasets. Time-dependent AUC at 5 years uses inverse-probability-of-censoring weighting (timeROC). IDI is the difference in discrimination slope (mean predicted 5-year risk in those who had the event minus those who did not) relative to the base model. Net benefit is computed at 5 years with inverse-probability-of-censoring weights; a higher value at a given risk threshold indicates greater clinical utility of acting on the model's predictions.", B_BOOT, M),
      "Within the robust stratum the frailty phenotype offers no gradation by construction; the comparison therefore tests whether intrinsic capacity, and intrinsic capacity purged of the phenotype's own items, adds prognostic information where a categorical phenotype cannot."))
}

cat("\n", strrep("=", 70), "\nPhase 3.1-3.2 요약\n", strrep("=", 70), "\n", sep = "")
ok <- TA$`Effect retained, %` != "-"
cat(sprintf("  IC-13 효과 유지율 (IC-17 유의 전이 %d개): 중앙 %s%% · 범위 %s-%s%%\n", sum(ok),
            median(as.numeric(TA$`Effect retained, %`[ok])), min(as.numeric(TA$`Effect retained, %`[ok])), max(as.numeric(TA$`Effect retained, %`[ok]))))
for (oo in c("dth", "wor")) { rb <- PB[PB$pop == "Robust stratum" & PB$outcome == oo, ]
  for (mm in c("M1", "M2", "M3", "M4", "M5")) { z <- rb[rb$model == mm, ]
    cat(sprintf("  robust %-4s %-40s C %.3f  ΔC %+.3f (%+.3f to %+.3f) P=%s  IDI %+.4f  NB10 %.4f\n", oo, MLAB[mm],
                z$C, z$dC, z$dC_l, z$dC_u, fmt_p(z$dC_p), z$IDI, z$NB10)) } }
cat("  -> M4(grip+gait+IC-13) 의 ΔC 가 M1(grip+gait) 보다 크면 IC 는 연속 grip·gait 너머의 정보를 더합니다.\n")
nb_ok <- all(is.finite(PB$NB10)) && all(PB$NB10 > -0.12) && all(PB$NB10 < 0.6)
cat(sprintf("  순이익 범위 점검: %s (NB10 %.3f ~ %.3f)\n", if (nb_ok) "정상" else "★ 비정상 — 가중치 확인", min(PB$NB10), max(PB$NB10)))
cat(strrep("=", 70), "\n=== 완료 ===\n")
