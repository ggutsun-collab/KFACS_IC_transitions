###############################################################################
## 14_mi_pool_tables.R   (Phase 2 — m=20 완성 데이터에서 Table 2·3 재추정 + Rubin 풀링)
##                                                                        v260903
## 입력 : IMP_DIR/KFACS_mi_pmm_m20.rds (13_ 산출), 01_functions_common 의 rubin_pool()
## 노출 : "withinsex" (Phase 1 결정, primary) 와 "pooled" (기존 primary, ED 로 이동)
##
## 산출 (FIG_DIR):
##   Table2_MI            — 18 전이 IRR per +1 s.d. 및 T2/T3 vs T1, 두 척도, m 풀링
##   Table3_MI            — wave-1 frailty × within-sex IC 삼분위: 5년 악화·사망 (KM 위험, Cox HR)
##                          + ADL 회복(원인별 Cox), robust 층 per-s.d. HR
##   MI_FMI               — 전이별 대체분율(FMI), 단일대체 대비 CI 폭 변화
##
## 실행: source("14_mi_pool_tables.R", encoding = "UTF-8")   (10-25분)
###############################################################################

POOL_VERSION <- "v260903"
message("\n=== 14_mi_pool_tables ", POOL_VERSION, " ===")

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
.need(c("survival", "openxlsx"))
suppressPackageStartupMessages(library(survival))

###############################################################################
## 1. 설정
###############################################################################
MI_STEM   <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
SCALES    <- c("withinsex", "pooled")
TAU       <- 5
MIN_EV_HR <- 5; MIN_N_HR <- 30; MIN_N_RISK <- 10
MIN_M_OK  <- 0.8      # 전체 m 중 이 비율 이상 추정된 항만 풀링

MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
SETS <- names(MI); M <- length(SETS)
msg(sprintf("MI 세트 %d개", M))

###############################################################################
## 2. 척도 적용 — prep_long 결과의 gLIC_z / gLIC_tert_W1cut 를 교체
##    withinsex: wave-1 성별 평균·SD 로 z, 성별 삼분위 절단점
###############################################################################
apply_scale <- function(d, scale) {
  if (scale == "pooled") return(d)
  w1 <- d[!duplicated(d$id), ]
  z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE))
  }
  d$gLIC_z <- z
  d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS)
  d
}
## within-sex 척도에서는 sex 가 표준화에 흡수되므로 보정군에서 뺍니다 (ED Table 1 정의)
adj_for <- function(scale) if (scale == "withinsex") setdiff(CFG$ADJ, "sex_f") else CFG$ADJ

###############################################################################
## 3. Table 3 용 개인별 사건시각 (22_table3 의 mk_events / mk_surv 동형)
###############################################################################
mk_events <- function(dat) {
  dat <- dat[order(dat$id, dat$time), ]
  sp <- split(seq_len(nrow(dat)), dat$id)
  do.call(rbind, lapply(sp, function(ix) {
    s <- dat$state[ix]; tm <- dat$time[ix]; ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
    n <- length(tm); if (n < 1 || !is.finite(s[1])) return(NULL)
    s0 <- s[1]; tmax <- max(tm)
    fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tr <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] < s0 & s[-1] >= 1) else NA_real_
    tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] > s0) else NA_real_
    kd <- which(is.finite(s) & s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
    data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax, t_rec = tr, t_wor = tw, t_dth = td)
  }))
}
mk_surv <- function(E, out) {
  inf <- function(x) ifelse(is.finite(x), x, Inf)
  if (out == "rec") { ft <- pmin(inf(E$t_rec), inf(E$t_dth), E$tmax, TAU)
    ev <- rep(0L, nrow(E)); ev[is.finite(E$t_rec) & E$t_rec <= ft] <- 1L
    ev[ev == 0L & is.finite(E$t_dth) & E$t_dth <= ft] <- 2L
  } else if (out == "wor") { ft <- pmin(inf(E$t_wor), E$tmax, TAU); ev <- as.integer(is.finite(E$t_wor) & E$t_wor <= ft)
  } else { ft <- pmin(inf(E$t_dth), E$tmax, TAU); ev <- as.integer(is.finite(E$t_dth) & E$t_dth <= ft) }
  data.frame(id = E$id, s0 = E$s0, ft = pmax(ft, 1e-6), ev = ev)
}
## 5년 위험 + logit 척도 SE (Rubin 풀링용)
risk5 <- function(ft, ev, competing) {
  if (competing) {
    if (!any(ev == 1L)) return(c(p = 0, se = NA))
    ef <- factor(ev, levels = 0:2, labels = c("censor", "event", "compete"))
    f <- survfit(Surv(ft, ef) ~ 1); s <- summary(f, times = TAU, extend = TRUE)
    j <- match("event", f$states); p <- s$pstate[1, j]; sd <- s$std.err[1, j]
  } else {
    if (!any(ev > 0)) return(c(p = 0, se = NA))
    f <- survfit(Surv(ft, ev) ~ 1); s <- summary(f, times = TAU, extend = TRUE)
    p <- 1 - s$surv[1]; sd <- s$std.err[1] * s$surv[1]    # std.err(surv) 는 log 척도 -> 확률 척도
  }
  p <- min(max(p, 1e-6), 1 - 1e-6)
  c(p = p, se = sd / (p * (1 - p)))                        # delta -> logit 척도
}

###############################################################################
## 4. 세트별 적합
###############################################################################
T2 <- list(); T3 <- list(); RS <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, M))
  load_imputed(set, stem = MI_STEM)
  d0 <- prep_long(load_kfacs())
  for (sc in SCALES) {
    d <- apply_scale(d0, sc); adj <- adj_for(sc)
    ## ── Table 2: 전이 IRR (per-s.d. 와 삼분위) ──
    for (sv in c("state", "frailty_3cat")) {
      iv <- build_intervals(d, sv, "rolling")
      r1 <- run_transition_table(iv, exposure_terms = "gLIC_z", adj = adj)
      r2 <- run_transition_table(iv, exposure_terms = "gLIC_tert", adj = adj)
      r  <- rbind(r1, r2); r$scale <- sc; r$set <- set
      T2[[length(T2) + 1]] <- r[, c("set","scale","system","transition","term","IRR","lcl","ucl","p","n_events","pyears","note")]
    }
    ## ── Table 3: wave-1 셀 ──
    b <- get_baseline(d)
    b$tert <- factor(b$gLIC_tert_W1cut, levels = CFG$TERT_LABELS)
    b$fr   <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))
    b <- b[!is.na(b$tert) & !is.na(b$fr), ]
    E <- mk_events(d)
    ADJ <- intersect(adj, names(b))
    for (out in c("wor", "dth", "rec")) {
      S <- merge(b[, c("id", "fr", "tert", ADJ)], mk_surv(E, out), by = "id")
      S <- if (out == "rec") S[S$s0 %in% c(2, 3), ] else S[S$s0 < 4, ]
      S$cell <- factor(paste(S$fr, S$tert, sep = "|"),
                       levels = as.vector(t(outer(levels(S$fr), CFG$TERT_LABELS, paste, sep = "|"))))
      S$cell <- droplevels(S$cell); S$cell <- relevel(S$cell, ref = "Robust|T3")
      cx <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ cell +", paste(ADJ, collapse = "+"))), data = S), silent = TRUE)
      for (cl in levels(S$cell)) {
        z <- S[S$cell == cl, ]; ne <- sum(z$ev == 1L)
        rk <- risk5(z$ft, z$ev, competing = out == "rec")
        nm <- paste0("cell", cl); bhr <- NA; sehr <- NA
        if (cl != "Robust|T3" && !inherits(cx, "try-error") && nm %in% names(coef(cx)) &&
            ne >= MIN_EV_HR && nrow(z) >= MIN_N_HR) {
          bhr <- coef(cx)[nm]; sehr <- sqrt(vcov(cx)[nm, nm]) }
        T3[[length(T3) + 1]] <- data.frame(set = set, scale = sc, outcome = out, cell = cl,
                                           n = nrow(z), events = ne, p = rk["p"], se_logit = rk["se"],
                                           b_hr = bhr, se_hr = sehr, stringsAsFactors = FALSE)
      }
      ## robust 층 내 per-s.d. HR (사망·악화)
      if (out != "rec") {
        R <- merge(b[b$fr == "Robust", c("id", "gLIC_z", ADJ)], mk_surv(E, out), by = "id")
        cr <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ gLIC_z +", paste(ADJ, collapse = "+"))), data = R), silent = TRUE)
        if (!inherits(cr, "try-error"))
          RS[[length(RS) + 1]] <- data.frame(set = set, scale = sc, outcome = out, n = nrow(R), events = sum(R$ev == 1L),
                                             b = coef(cr)["gLIC_z"], se = sqrt(vcov(cr)["gLIC_z", "gLIC_z"]))
      }
    }
  }
}
T2 <- do.call(rbind, T2); T3 <- do.call(rbind, T3); RS <- do.call(rbind, RS)
save_vals(T2, "MI_Table2_perset.csv", FIG_DIR); save_vals(T3, "MI_Table3_perset.csv", FIG_DIR)

###############################################################################
## 5. Rubin 풀링
###############################################################################
pool_rows <- function(b, se) {
  ok <- is.finite(b) & is.finite(se)
  if (sum(ok) < MIN_M_OK * M) return(c(est = NA, se = NA, lcl = NA, ucl = NA, p = NA, m = sum(ok), fmi = NA))
  r <- rubin_pool(b[ok], se[ok])
  bvar <- stats::var(b[ok]); ubar <- mean(se[ok]^2)
  fmi <- (1 + 1/sum(ok)) * bvar / (ubar + (1 + 1/sum(ok)) * bvar)
  c(r, fmi = fmi)
}
fmt_ci <- function(e, l, u) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(u)), "NE")
fmt_p  <- function(p) ifelse(!is.finite(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## Table 2
T2$b <- log(T2$IRR); T2$se <- (log(T2$ucl) - log(T2$lcl)) / (2 * qnorm(0.975))
T2$key <- paste(T2$scale, T2$system, T2$transition, T2$term, sep = "|")
P2 <- do.call(rbind, lapply(split(T2, T2$key), function(z) {
  r <- pool_rows(z$b, z$se)
  data.frame(scale = z$scale[1], system = z$system[1], transition = z$transition[1], term = z$term[1],
             n_events = round(median(z$n_events)), pyears = round(median(z$pyears)),
             est = r["est"], se = r["se"], lcl = r["lcl"], ucl = r["ucl"], p = r["p"], m = r["m"], fmi = r["fmi"],
             stringsAsFactors = FALSE)
}))
rownames(P2) <- NULL
T2TAB <- do.call(rbind, lapply(SCALES, function(sc) {
  z <- P2[P2$scale == sc, ]
  keys <- unique(z[, c("system", "transition")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    w <- z[z$system == keys$system[k] & z$transition == keys$transition[k], ]
    g <- function(term) { x <- w[w$term == term, ]; if (!nrow(x)) return(c(NA, NA, NA, NA)); c(x$est, x$lcl, x$ucl, x$p) }
    a <- g("gLIC_z"); t2 <- g("gLIC_tertT2"); t3 <- g("gLIC_tertT3")
    data.frame(Scale = sc, System = ifelse(keys$system[k] == "state", "Disability, 3 states", "Frailty phenotype"),
               Transition = keys$transition[k], Events = w$n_events[1], `Person-years` = w$pyears[1],
               `IRR per +1 s.d. (95% CI)` = fmt_ci(a[1], a[2], a[3]), P = fmt_p(a[4]),
               `T2 vs T1 (95% CI)` = fmt_ci(t2[1], t2[2], t2[3]), `P ` = fmt_p(t2[4]),
               `T3 vs T1 (95% CI)` = fmt_ci(t3[1], t3[2], t3[3]), `P  ` = fmt_p(t3[4]),
               FMI = ifelse(is.finite(w$fmi[w$term == "gLIC_z"]), sprintf("%.2f", w$fmi[w$term == "gLIC_z"]), "-"),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}))
save_table(T2TAB, "Table2_MI",
  title = "Table 2 | Intrinsic-capacity level and state-transition rates (multiple imputation, m = 20, Rubin's rules)",
  footnotes = c(
    sprintf("Incidence rate ratios from transition-specific Poisson models with a log(person-time) offset and participant-clustered standard errors, fitted in each of %d completed datasets generated by multiple imputation with predictive mean matching and combined by Rubin's rules. Events and person-years are medians across completed datasets.", M),
    "Within-sex scale: the general-factor score standardized within sex on the wave-1 distribution, with sex-specific tertile cut-points; models are adjusted for age, baseline comorbidity count, education, household income and area of residence (sex is absorbed by the standardization). Pooled scale: standardized on the whole wave-1 distribution; models additionally adjusted for sex.",
    "FMI, fraction of missing information for the per-s.d. estimate. NE, not estimable: the tertile contrast could not be estimated in at least 80% of completed datasets (a tertile contributed fewer than 5 events).",
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

## Table 3
T3$key <- paste(T3$scale, T3$outcome, T3$cell, sep = "|")
P3 <- do.call(rbind, lapply(split(T3, T3$key), function(z) {
  rr <- pool_rows(stats::qlogis(z$p), z$se_logit)
  rh <- pool_rows(z$b_hr, z$se_hr)
  data.frame(scale = z$scale[1], outcome = z$outcome[1], cell = z$cell[1],
             n = round(median(z$n)), events = round(median(z$events)),
             risk = if (is.finite(rr["est"])) 100 * stats::plogis(rr["est"]) else 100 * mean(z$p),
             risk_lo = 100 * stats::plogis(rr["lcl"]), risk_hi = 100 * stats::plogis(rr["ucl"]),
             hr = rh["est"], hr_lo = rh["lcl"], hr_hi = rh["ucl"], hr_p = rh["p"], stringsAsFactors = FALSE)
}))
rownames(P3) <- NULL
P3$risk[P3$n < MIN_N_RISK] <- NA
fmt_r <- function(e, l, u) ifelse(!is.finite(e), "-", ifelse(is.finite(l), sprintf("%.1f (%.1f-%.1f)", e, l, u), sprintf("%.1f", e)))
fmt_h <- function(cell, e, l, u) ifelse(cell == "Robust|T3", "1.00 (reference)", fmt_ci(e, l, u))
OUTLAB <- c(wor = "Composite worsening", dth = "Death (all cause)", rec = "Recovery from ADL disability")
T3TAB <- do.call(rbind, lapply(SCALES, function(sc) do.call(rbind, lapply(c("wor", "dth", "rec"), function(o) {
  z <- P3[P3$scale == sc & P3$outcome == o, ]
  z <- z[order(match(sub("\\|.*", "", z$cell), c("Robust", "Pre-frail", "Frail")), sub(".*\\|", "", z$cell)), ]
  data.frame(Scale = sc, Outcome = OUTLAB[o], `Frailty phenotype` = sub("\\|.*", "", z$cell),
             `IC tertile` = sub(".*\\|", "", z$cell), n = z$n, Events = z$events,
             `5-yr risk, % (95% CI)` = fmt_r(z$risk, z$risk_lo, z$risk_hi),
             `HR (95% CI)` = fmt_h(z$cell, z$hr, z$hr_lo, z$hr_hi), P = fmt_p(z$hr_p),
             check.names = FALSE, stringsAsFactors = FALSE)
}))))
save_table(T3TAB, "Table3_MI",
  title = "Table 3 | Five-year worsening, death and recovery according to frailty phenotype and intrinsic-capacity tertile at wave 1 (multiple imputation, m = 20)",
  footnotes = c(
    "Strata are defined at wave 1 by frailty phenotype and by tertile of the intrinsic-capacity general factor; on the within-sex scale tertile cut-points are sex-specific so that each tertile contains the same proportion of women. Follow-up is censored at 5 years.",
    "Composite worsening is the first transition to a worse ADL disability state or death; absolute risks are one minus the Kaplan-Meier estimate. Recovery is defined among wave-1 ADL disability carriers with death as a competing event; absolute risks are Aalen-Johansen cumulative incidences and hazard ratios are cause-specific.",
    sprintf("Estimates were obtained in each of %d completed datasets and combined by Rubin's rules (hazard ratios on the log scale, absolute risks on the logit scale). Hazard ratios are omitted for cells with fewer than %d events or %d participants; absolute risks are omitted for cells with fewer than %d participants.", M, MIN_EV_HR, MIN_N_HR, MIN_N_RISK),
    "Reference cell is the robust stratum, highest capacity tertile. Models are adjusted for age, baseline comorbidity count, education, household income and area of residence, plus sex on the pooled scale."))

## robust 층 per-s.d.
RS$key <- paste(RS$scale, RS$outcome)
PRS <- do.call(rbind, lapply(split(RS, RS$key), function(z) { r <- pool_rows(z$b, z$se)
  data.frame(Scale = z$scale[1], Outcome = OUTLAB[z$outcome[1]], n = round(median(z$n)), Events = round(median(z$events)),
             `HR per +1 s.d. (95% CI)` = fmt_ci(r["est"], r["lcl"], r["ucl"]), P = fmt_p(r["p"]),
             FMI = sprintf("%.2f", r["fmi"]), check.names = FALSE) }))
save_table(PRS, "Table3_MI_robust_perSD",
  title = "Supplementary Table | Intrinsic capacity and five-year outcomes within the robust stratum, per +1 s.d. (multiple imputation)",
  footnotes = "Cox models restricted to participants classified as robust at wave 1, adjusted as in Table 3; estimates combined by Rubin's rules across 20 completed datasets.")

## FMI 요약
fm <- P2[P2$term == "gLIC_z" & is.finite(P2$fmi), ]
cat("\n", strrep("=", 70), "\nMI 풀링 요약\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  per-s.d. 추정치 %d개 | FMI 중앙 %.2f · 최대 %.2f (%s)\n", nrow(fm),
            median(fm$fmi), max(fm$fmi), fm$transition[which.max(fm$fmi)]))
for (sc in SCALES) { z <- fm[fm$scale == sc, ]
  cat(sprintf("  [%s] 유의(P<0.05) %d/%d\n", sc, sum(z$p < 0.05), nrow(z))) }
cat("  -> Table2_MI / Table3_MI 가 원고 Table 2·3 의 새 원천입니다 (기존 kNN 단일대체는 ED 민감도).\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
