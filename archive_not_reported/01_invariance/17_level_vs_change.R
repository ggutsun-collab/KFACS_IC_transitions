###############################################################################
## 17_level_vs_change.R   (Phase 3.4 — 수준 대 변화, 결과별 비대칭)   v260903
##
## 원고의 "level > change" 섹션은 사망만 다루고, 상호보정 change 항의 역방향 계수(1.30)를
## 길게 변명하며, 2년 간격의 측정오차로 스스로를 무력화합니다. 재설계:
##   A. 이산시간 Poisson (MI m=20, within-sex): 다음 구간의 결과 ~ 현재 수준 + 직전 2년 변화
##      결과 4종: 사망 / ADL 악화(Normal->Mild|Severe) / Normal->Severe / Robust->Frail
##      -> 수준과 변화의 상대적 정보량을 '결과별로' 비교 (비대칭이 핵심 메시지)
##      + 현재 수준·이전 수준 재모수화 (원고 Fig 2e 와 동일 논리)
##   B. JMbayes2: 사망과 복합악화(첫 ADL 악화 또는 사망)에 대해
##      functional_forms = value + slope  -> 현재 값과 순간 기울기의 HR 을 동시에
##      (MI 세트 1개, 기본 MI01; JM 은 세트당 10-20분)
##
## 출력: ST_LevelChange_Discrete (A), ST_JM_ValueSlope (B), P3_jm_fits.rds
## 실행: source("17_level_vs_change.R", encoding = "UTF-8")   (A 5-10분, B 20-40분)
###############################################################################

LC_VERSION <- "v260903"
message("\n=== 17_level_vs_change ", LC_VERSION, " ===")

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
.need(c("survival", "nlme", "JMbayes2", "openxlsx"))
suppressPackageStartupMessages(library(survival))

###############################################################################
## 1. 설정
###############################################################################
MI_STEM  <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
JM_SETS  <- c("MI01")                  # JM 을 적합할 세트 (추가: c("MI01","MI02"))
JM_ITER  <- 12000L; JM_BURN <- 4000L; JM_CHAINS <- 3L; JM_THIN <- 5L
JM_CACHE <- file.path(FIG_DIR, "P3_jm_fits.rds")
MIN_M_OK <- 0.8
MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))); SETS <- names(MI); M <- length(SETS)
ADJ_WS <- setdiff(CFG$ADJ, "sex_f")

ws_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]; z <- d$gLIC
  for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
    z[k] <- (d$gLIC[k] - mean(w1$gLIC[k1], na.rm = TRUE)) / sd(w1$gLIC[k1], na.rm = TRUE) }
  d$gLIC_z <- z; d
}
rvc <- function(m, cl) { v <- try(sandwich::vcovCL(m, cluster = cl, type = "HC0"), silent = TRUE)
  if (inherits(v, "try-error")) vcov(m) else v }
## 구간 자료에 직전 수준·변화를 붙임 (구간 k 는 wave k 에서 시작; k>=2 만)
add_prev <- function(iv, d) {
  key <- paste(d$id, d$wave); z <- d$gLIC_z
  iv$prev_z <- z[match(paste(iv$id, iv$k - 1), key)]        # 직전 wave 의 수준
  iv$chg_z  <- iv$gLIC_z - iv$prev_z                        # 직전 2년 변화 (같은 s.d. 단위)
  iv
}
fit_lc <- function(dd, rhs) {
  f <- as.formula(paste("y ~", rhs, "+", paste(ADJ_WS, collapse = "+"), "+ from_lab + offset(log(dur))"))
  dd$from_lab <- droplevels(factor(dd$from_lab))
  if (nlevels(dd$from_lab) < 2) f <- as.formula(paste("y ~", rhs, "+", paste(ADJ_WS, collapse = "+"), "+ offset(log(dur))"))
  m <- try(glm(f, family = poisson(), data = dd), silent = TRUE)
  if (inherits(m, "try-error") || !m$converged) return(NULL)
  V <- rvc(m, dd$id); b <- coef(m); se <- sqrt(diag(V))[names(b)]
  data.frame(term = names(b), b = unname(b), se = unname(se))
}

###############################################################################
## 2. A — 이산시간 수준·변화 (MI m=20)
###############################################################################
OUTS <- list(
  list(id = "death",   lab = "Death (from any living state)",            sv = "state",        from = c("Normal","Mild","Severe"), to = "Death"),
  list(id = "adl_wor", lab = "Incident ADL disability (Normal -> Mild or Severe)", sv = "state", from = "Normal", to = c("Mild","Severe")),
  list(id = "severe",  lab = "Normal -> Severe ADL disability",         sv = "state",        from = "Normal", to = "Severe"),
  list(id = "r2f",     lab = "Robust -> Frail",                          sv = "frailty_3cat", from = "Robust", to = "Frail"),
  list(id = "pf2f",    lab = "Pre-frail -> Frail",                       sv = "frailty_3cat", from = "Pre-frail", to = "Frail"))
SPECS <- list(level_only = "gLIC_z", change_only = "chg_z", joint = "gLIC_z + chg_z", reparam = "gLIC_z + prev_z")
A <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]; msg(sprintf("[A %s] %d/%d", set, i, M))
  load_imputed(set, stem = MI_STEM); d <- ws_scale(prep_long(load_kfacs()))
  if (!"wave" %in% names(d)) d$wave <- round(d$time / CFG$WAVE_GAP) + 1
  d$wave <- as.integer(round(as_num(d$wave)))
  for (O in OUTS) {
    iv <- add_prev(build_intervals(d, O$sv, "rolling"), d)
    dd <- iv[iv$from_lab %in% O$from & is.finite(iv$prev_z), ]
    dd$y <- as.integer(dd$to_lab %in% O$to)
    for (sp in names(SPECS)) {
      r <- fit_lc(dd, SPECS[[sp]]); if (is.null(r)) next
      r <- r[r$term %in% c("gLIC_z", "chg_z", "prev_z"), ]
      r$set <- set; r$outcome <- O$id; r$spec <- sp; r$events <- sum(dd$y); r$n_int <- nrow(dd)
      A[[length(A) + 1]] <- r
    }
  }
}
A <- do.call(rbind, A); save_vals(A, "ST_LevelChange_perset.csv", FIG_DIR)
A$key <- paste(A$outcome, A$spec, A$term, sep = "|")
PA <- do.call(rbind, lapply(split(A, A$key), function(z) {
  ok <- is.finite(z$b) & is.finite(z$se); if (sum(ok) < MIN_M_OK * M) return(NULL)
  r <- rubin_pool(z$b[ok], z$se[ok])
  data.frame(outcome = z$outcome[1], spec = z$spec[1], term = z$term[1], events = round(median(z$events)), n_int = round(median(z$n_int)),
             est = r["est"], lcl = r["lcl"], ucl = r["ucl"], p = r["p"]) }))
rownames(PA) <- NULL
fi <- function(e, l, u) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(u)), "-")
g <- function(o, sp, tm) { z <- PA[PA$outcome == o & PA$spec == sp & PA$term == tm, ]; if (!nrow(z)) c(NA, NA, NA) else c(z$est, z$lcl, z$ucl) }
TA <- do.call(rbind, lapply(OUTS, function(O) {
  ev <- PA$events[PA$outcome == O$id][1]; ni <- PA$n_int[PA$outcome == O$id][1]
  l <- g(O$id, "level_only", "gLIC_z"); c1 <- g(O$id, "change_only", "chg_z")
  jl <- g(O$id, "joint", "gLIC_z"); jc <- g(O$id, "joint", "chg_z"); rc <- g(O$id, "reparam", "gLIC_z"); rp <- g(O$id, "reparam", "prev_z")
  data.frame(Outcome = O$lab, Intervals = ni, Events = ev,
             `Level alone` = fi(l[1], l[2], l[3]), `Change alone` = fi(c1[1], c1[2], c1[3]),
             `Level, mutually adjusted` = fi(jl[1], jl[2], jl[3]), `Change, mutually adjusted` = fi(jc[1], jc[2], jc[3]),
             `Current level (reparameterized)` = fi(rc[1], rc[2], rc[3]), `Previous level (reparameterized)` = fi(rp[1], rp[2], rp[3]),
             `Share of information in change, %` = ifelse(is.finite(jc[1]) & is.finite(jl[1]), sprintf("%.0f", 100 * abs(jc[1]) / (abs(jc[1]) + abs(jl[1]))), "-"),
             check.names = FALSE, stringsAsFactors = FALSE) }))
save_table(TA, "ST_LevelChange_Discrete",
  title = "Supplementary Table | Prognostic information in the current level of intrinsic capacity and in its change over the preceding two years, by outcome",
  footnotes = c(
    sprintf("Discrete-time Poisson models of the transition over the next between-visit interval, restricted to intervals preceded by a capacity measurement, fitted in %d completed datasets and combined by Rubin's rules. Rate ratios per +1 s.d. (within-sex scale); level and change are in the same units. Models are adjusted for age, baseline comorbidity count, education, household income, area of residence and, where more than one origin state contributes, the origin state; standard errors are clustered by participant.", M),
    "The reparameterized model replaces the change term by the level at the preceding visit; because change = current level - previous level, the two parameterizations are algebraically equivalent, and a mutually adjusted change coefficient above 1 with a strongly protective level coefficient indicates that, at a given current level, a larger recent gain implies a lower previous level.",
    "Share of information in change is |log RR(change)| / (|log RR(change)| + |log RR(level)|) from the mutually adjusted model, a descriptive index of how much of the joint prognostic information the recent change carries for each outcome.",
    "Change over a two-year interval is measured with more error than level and the interval is coarse relative to the time scale of functional decline; the comparison concerns the prognostic information available in each summary at a given encounter, not the biological importance of change."))
cat("\n[A] 결과별 수준·변화\n"); print(TA[, c(1, 3, 4, 5, 6, 7, 10)], row.names = FALSE)

###############################################################################
## 3. B — JMbayes2 value + slope (사망, 복합악화)
###############################################################################
if (!requireNamespace("JMbayes2", quietly = TRUE)) stop("JMbayes2 가 필요합니다.")
mk_jm_data <- function(d) {
  long <- d[is.finite(d$time) & is.finite(d$gLIC_z), c("id", "time", "gLIC_z", "age0", "sex_f", "comorbid_bl")]
  long <- long[order(long$id, long$time), ]
  dd <- d[order(d$id, d$time), ]
  sv <- do.call(rbind, lapply(split(dd, dd$id), function(p) {
    n <- nrow(p); s <- p$state; tm <- p$time; s0 <- s[1]
    kw <- which(is.finite(s[-1]) & s[-1] > s0); tw <- if (length(kw)) tm[-1][kw[1]] else NA
    data.frame(id = p$id[1], t_death = p$followup_years[n], ev_death = as.integer(p$death_event[n] == 1),
               t_wor = if (is.finite(tw)) tw else p$followup_years[n], ev_wor = as.integer(is.finite(tw) | p$death_event[n] == 1),
               age0 = p$age0[1], sex_f = p$sex_f[1], comorbid_bl = p$comorbid_bl[1], s0 = s0) }))
  sv <- sv[is.finite(sv$t_death) & sv$t_death > 0 & sv$s0 < 4, ]
  sv$t_wor <- pmin(sv$t_wor, sv$t_death); sv$t_wor <- pmax(sv$t_wor, 1e-3)
  long <- long[long$id %in% sv$id, ]
  long <- merge(long, sv[, c("id", "t_death")], by = "id"); long <- long[long$time <= long$t_death, ]
  long <- long[order(long$id, long$time), ]; sv <- sv[sv$id %in% unique(long$id), ]; sv <- sv[order(sv$id), ]
  list(long = long, sv = sv)
}
fit_jm <- function(J, outcome) {
  sv <- J$sv; long <- J$long
  if (outcome == "wor") { sv$t <- sv$t_wor; sv$ev <- sv$ev_wor; long <- long[long$time <= long$t_death, ] } else { sv$t <- sv$t_death; sv$ev <- sv$ev_death }
  ## 복합악화: 사건 이후의 종단 측정은 제외 (내생성)
  if (outcome == "wor") { long <- merge(long, sv[, c("id", "t")], by = "id"); long <- long[long$time <= long$t, ]; long <- long[order(long$id, long$time), ] }
  sv <- sv[sv$id %in% unique(long$id), ]
  lmeFit <- try(nlme::lme(gLIC_z ~ time + age0 + sex_f, random = ~ time | id, data = long,
                          control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)), silent = TRUE)
  if (inherits(lmeFit, "try-error")) lmeFit <- nlme::lme(gLIC_z ~ time + age0 + sex_f, random = ~ 1 | id, data = long,
                                                        control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
  coxFit <- coxph(Surv(t, ev) ~ age0 + sex_f + comorbid_bl, data = sv, x = TRUE, model = TRUE)
  jm <- JMbayes2::jm(coxFit, list(lmeFit), time_var = "time",
                     functional_forms = list("gLIC_z" = ~ value(gLIC_z) + slope(gLIC_z)),
                     n_iter = JM_ITER, n_burnin = JM_BURN, n_chains = JM_CHAINS, n_thin = JM_THIN, seed = CFG$SEED)
  jm
}
JMF <- if (file.exists(JM_CACHE)) readRDS(JM_CACHE) else list()
for (set in JM_SETS) {
  load_imputed(set, stem = MI_STEM); d <- ws_scale(prep_long(load_kfacs())); J <- mk_jm_data(d)
  msg(sprintf("[B %s] N=%d, deaths=%d, worsening events=%d, long obs=%d", set, nrow(J$sv), sum(J$sv$ev_death), sum(J$sv$ev_wor), nrow(J$long)))
  for (o in c("dth", "wor")) {
    kk <- paste(set, o)
    if (!is.null(JMF[[kk]])) { msg("  캐시 재사용: ", kk); next }
    t0 <- Sys.time(); JMF[[kk]] <- fit_jm(J, o)
    msg(sprintf("  [%s] JM 완료 %.1f분", kk, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    saveRDS(JMF, JM_CACHE)
  }
}
jm_row <- function(fit, set, o) {
  s <- summary(fit); S <- s$Survival
  rn <- rownames(S); vi <- grep("value\\(gLIC_z\\)", rn); si <- grep("slope\\(gLIC_z\\)", rn)
  rh <- tryCatch(max(unlist(fit$statistics$Rhat), na.rm = TRUE), error = function(e) NA)
  data.frame(set = set, outcome = o,
             value_hr = exp(S[vi, "Mean"]), value_lo = exp(S[vi, "2.5%"]), value_hi = exp(S[vi, "97.5%"]),
             slope_hr = exp(S[si, "Mean"]), slope_lo = exp(S[si, "2.5%"]), slope_hi = exp(S[si, "97.5%"]),
             rhat = rh, n = fit$model_data$n, events = sum(fit$model_data$event))
}
JR <- do.call(rbind, lapply(names(JMF), function(kk) { p <- strsplit(kk, " ")[[1]]; jm_row(JMF[[kk]], p[1], p[2]) }))
OL <- c(dth = "All-cause mortality", wor = "Composite worsening (first ADL worsening or death)")
TB <- data.frame(Outcome = OL[JR$outcome], Dataset = JR$set, n = JR$n, Events = JR$events,
                 `Current value, HR per +1 s.d. (95% CrI)` = sprintf("%.2f (%.2f-%.2f)", JR$value_hr, JR$value_lo, JR$value_hi),
                 `Instantaneous slope, HR per +1 s.d./year (95% CrI)` = sprintf("%.2f (%.2f-%.2f)", JR$slope_hr, JR$slope_lo, JR$slope_hi),
                 `Max Rhat` = sprintf("%.3f", JR$rhat), check.names = FALSE, stringsAsFactors = FALSE)
save_table(TB, "ST_JM_ValueSlope",
  title = "Supplementary Table | Joint longitudinal-survival models with the current value and the instantaneous slope of intrinsic capacity",
  footnotes = c(
    "Bayesian joint models (JMbayes2) linking a linear mixed model for the within-sex standardized capacity score (fixed effects of time, baseline age and sex; random intercept and slope) to a Cox model adjusted for age, sex and baseline comorbidity count, through both the current value and the current rate of change of the trajectory. Three chains of 12,000 iterations (4,000 burn-in), package-default priors.",
    "For composite worsening, longitudinal measurements after the event are excluded and the event time is the first visit at which a worse ADL state was recorded, or death.",
    "Fitted in the first completed dataset of the multiple-imputation series; the discrete-time comparison in the accompanying table uses all completed datasets."))
cat("\n[B] JM value + slope\n"); print(TB, row.names = FALSE)

cat("\n", strrep("=", 70), "\n수준 대 변화 — 판정\n", strrep("=", 70), "\n", sep = "")
sh <- TA$`Share of information in change, %`
cat(sprintf("  변화의 정보 비중 (상호보정): %s\n", paste(sprintf("%s %s%%", TA$Outcome, sh), collapse = " | ")))
cat("  -> 사망에서 낮고 기능 악화(Normal->Severe, Robust->Frail)에서 높으면 '결과별 비대칭' 메시지가 성립합니다.\n")
cat("     그렇지 않으면 이 섹션은 Supplementary 로 내리고 본문에서는 한 문장으로 처리하십시오.\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
