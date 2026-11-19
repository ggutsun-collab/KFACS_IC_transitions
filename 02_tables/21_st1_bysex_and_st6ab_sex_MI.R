###############################################################################
## 21_st1_bysex_and_st6ab_sex_MI.R                                      v261117
## (from 19_coauthor_additions.R v260912c; the horizon analyses of the former
##  Supplementary Table 18 were removed with that table)
##
## Outputs (FIG_DIR):
##   ST1_Baseline_bySex_MI          Supplementary Table 1: baseline characteristics
##                                  by sex and within-sex tertile of IC at wave 1
##   ST6a_SexCovariate_Table2_MI    Supplementary Table 6a: transition rate ratios
##                                  with sex added as a covariate
##   ST6b_SexCovariate_Table3_MI    Supplementary Table 6b: within-state hazard
##                                  ratios with sex added as a covariate
## Same scale, adjustment and pooling rules as 14_. Covariates use the imputed
## education, income and area (R/04_cov_imputation.R).
###############################################################################
message("\n=== 21_st1_bysex_and_st6ab_sex_MI v261117 ===")

## -- helpers (R/_bootstrap.R -> kf_init(); all paths come from R/00_setup.R) --
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
## 1. 설정 (14_ 와 동일)
###############################################################################
MI_STEM   <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")
TAU_MAIN  <- 5                       # 본분석 지평
MIN_EV_HR <- 5; MIN_N_HR <- 30; MIN_N_RISK <- 10; MIN_M_OK <- 0.8
LV <- CFG$TERT_LABELS

MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
SETS <- names(MI); M <- length(SETS); rm(MI); invisible(gc())
msg(sprintf("MI 세트 %d개", M))

## within-sex 척도 (14_ apply_scale 동형)
apply_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]
  z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = LV, right = TRUE))
  }
  d$gLIC_z <- z; d$gLIC_tert_W1cut <- factor(tt, levels = LV); d
}
ADJ_WS  <- setdiff(CFG$ADJ, "sex_f")     # 본분석 (within-sex, 성별 미포함)
ADJ_SEX <- CFG$ADJ                        # ST19: 성별을 공변량으로 되돌림

## 개인별 사건시각 (14_ mk_events 동형; 상태변수 선택 가능)
mk_events <- function(dat, var = "state") {
  dat <- dat[order(dat$id, dat$time), ]
  sp <- split(seq_len(nrow(dat)), dat$id)
  do.call(rbind, lapply(sp, function(ix) {
    s <- as_num(dat[[var]][ix]); tm <- dat$time[ix]; ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
    n <- length(tm); if (n < 1 || !is.finite(s[1])) return(NULL)
    s0 <- s[1]; tmax <- max(tm)
    fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tr <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] < s0) else NA_real_               # 첫 개선 (한 단계 이상; 상태 코드는 1..3, 사망 4)
    tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] > s0) else NA_real_               # 첫 악화 (사망 포함)
    kd <- which(is.finite(s) & s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
    data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax, t_rec = tr, t_wor = tw, t_dth = td)
  }))
}
mk_surv <- function(E, out, tau) {
  inf <- function(x) ifelse(is.finite(x), x, Inf)
  if (out == "rec") { ft <- pmin(inf(E$t_rec), inf(E$t_dth), E$tmax, tau)
    ev <- rep(0L, nrow(E)); ev[is.finite(E$t_rec) & E$t_rec <= ft] <- 1L
    ev[ev == 0L & is.finite(E$t_dth) & E$t_dth <= ft] <- 2L
  } else if (out == "wor") { ft <- pmin(inf(E$t_wor), E$tmax, tau); ev <- as.integer(is.finite(E$t_wor) & E$t_wor <= ft)
  } else { ft <- pmin(inf(E$t_dth), E$tmax, tau); ev <- as.integer(is.finite(E$t_dth) & E$t_dth <= ft) }
  data.frame(id = E$id, s0 = E$s0, ft = pmax(ft, 1e-6), ev = ev)
}
risk_at <- function(ft, ev, tau, competing) {
  if (competing) {
    if (!any(ev == 1L)) return(c(p = 0, se = NA))
    ef <- factor(ev, levels = 0:2, labels = c("censor", "event", "compete"))
    f <- survfit(Surv(ft, ef) ~ 1); s <- summary(f, times = tau, extend = TRUE)
    j <- match("event", f$states); p <- s$pstate[1, j]; sd <- s$std.err[1, j]
  } else {
    if (!any(ev > 0)) return(c(p = 0, se = NA))
    f <- survfit(Surv(ft, ev) ~ 1); s <- summary(f, times = tau, extend = TRUE)
    p <- 1 - s$surv[1]; sd <- s$std.err[1] * s$surv[1]
  }
  p <- min(max(p, 1e-6), 1 - 1e-6); c(p = p, se = sd / (p * (1 - p)))
}
pool_rows <- function(b, se) {
  ok <- is.finite(b) & is.finite(se)
  if (sum(ok) < MIN_M_OK * M) return(c(est = NA, se = NA, lcl = NA, ucl = NA, p = NA, m = sum(ok), fmi = NA))
  r <- rubin_pool(b[ok], se[ok]); bvar <- stats::var(b[ok]); ubar <- mean(se[ok]^2)
  c(r, fmi = (1 + 1/sum(ok)) * bvar / (ubar + (1 + 1/sum(ok)) * bvar))
}
fmt_ci <- function(e, l, u) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(u)), "NE")
fmt_r  <- function(e, l, u) ifelse(!is.finite(e), "-", sprintf("%.1f (%.1f-%.1f)", e, l, u))
fmt_p  <- function(p) ifelse(!is.finite(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
n_switch <- function(v) { v <- v[!is.na(v)]; if (length(v) < 2) 0L else sum(v[-1] != v[-length(v)]) }

###############################################################################
## 2. 세트별 적합
###############################################################################
B <- vector("list", M); EXTRA <- NULL
T2S <- list(); T3S <- list(); RSS <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, M))
  load_imputed(set, stem = MI_STEM)
  d <- apply_scale(prep_long(load_kfacs()))
  b <- get_baseline(d)
  b$gLIC_z <- d$gLIC_z[match(b$id, d$id)]
  b$tert   <- factor(d$gLIC_tert_W1cut[match(b$id, d$id)], levels = LV)
  b$fr     <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))

  ## ── (A) ST2 기저표용 세트별 데이터 ──
  if (is.null(EXTRA)) {
    be <- try(attach_extras(b, verbose = FALSE), silent = TRUE)
    EXTRA <- if (!inherits(be, "try-error"))
      be[, intersect(c("id", "alcohol_heavy", "social_activity"), names(be)), drop = FALSE] else data.frame(id = b$id)
  }
  bb <- merge(b, EXTRA, by = "id", all.x = TRUE, suffixes = c("", ".x"))
  bb$educ_high   <- as.integer(as_num(bb$edu_high_bl) >= 5)      # 8-level code, imputed (13b_); >= middle school
  bb$income_high <- as.integer(as_num(bb$income_high) == 1)      # imputed (13b_)
  bb$married     <- as.integer(as_num(bb$marri_status) == 1)
  bb$female      <- as.integer(as_num(bb$sex_unified) == 2)
  bb$state_chr   <- as.character(bb$state_lab); bb$frailty_chr <- as.character(bb$frailty_lab)
  bb$area_chr    <- as.character(as_num(bb$area_3cat))
  tr <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p)
    data.frame(id = p$id[1], n_tr_adl = n_switch(as_num(p$state)), n_tr_frailty = n_switch(as_num(p$frailty_3cat)))))
  bb <- merge(bb, tr, by = "id", all.x = TRUE)
  bb$grp <- factor(paste(ifelse(bb$female == 1, "Women", "Men"), bb$tert, sep = " "),
                   levels = c(paste("Men", LV), paste("Women", LV)))
  B[[i]] <- bb

  E_adl <- mk_events(d, "state")
  bx <- b[!is.na(b$tert) & !is.na(b$fr), ]

  ## ── (D) ST19: 성별 공변량 추가 (Table 2 · Table 3 · 층 내 per-s.d.) ──
  for (sv in c("state", "frailty_3cat")) {
    iv <- build_intervals(d, sv, "rolling")
    for (mod in c("primary", "sex")) {                    # primary = 14_ 와 동일(성별 미포함), sex = 성별 추가
      adj <- if (mod == "sex") ADJ_SEX else ADJ_WS
      r1 <- run_transition_table(iv, exposure_terms = "gLIC_z",    adj = adj)
      r2 <- run_transition_table(iv, exposure_terms = "gLIC_tert", adj = adj)
      r <- rbind(r1, r2); r$set <- set; r$model <- mod
      T2S[[length(T2S) + 1]] <- r[, c("set", "model", "system", "transition", "term", "IRR", "lcl", "ucl", "p", "n_events", "pyears")]
    }
  }
  for (out in c("wor", "dth")) {
    S <- merge(bx[, c("id", "fr", "tert", "gLIC_z", ADJ_SEX)], mk_surv(E_adl, out, TAU_MAIN), by = "id"); S <- S[S$s0 < 4, ]
    S$cell <- factor(paste(S$fr, S$tert, sep = "|"), levels = as.vector(t(outer(levels(S$fr), LV, paste, sep = "|"))))
    S$cell <- droplevels(S$cell); S$cell <- relevel(S$cell, ref = "Robust|T3")
    cx <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ cell +", paste(ADJ_SEX, collapse = "+"))), data = S), silent = TRUE)
    for (cl in setdiff(levels(S$cell), "Robust|T3")) {
      z <- S[S$cell == cl, ]; ne <- sum(z$ev == 1L); nm <- paste0("cell", cl)
      if (!inherits(cx, "try-error") && nm %in% names(coef(cx)) && ne >= MIN_EV_HR && nrow(z) >= MIN_N_HR)
        T3S[[length(T3S) + 1]] <- data.frame(set = set, outcome = out, cell = cl, n = nrow(z), events = ne,
                                             b_hr = coef(cx)[nm], se_hr = sqrt(vcov(cx)[nm, nm]))
    }
    for (st in c("Robust", "Pre-frail")) {
      R <- S[S$fr == st, ]
      cr <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ gLIC_z +", paste(ADJ_SEX, collapse = "+"))), data = R), silent = TRUE)
      if (!inherits(cr, "try-error"))
        RSS[[length(RSS) + 1]] <- data.frame(set = set, stratum = st, outcome = out, n = nrow(R), events = sum(R$ev == 1L),
                                             b = coef(cr)["gLIC_z"], se = sqrt(vcov(cr)["gLIC_z", "gLIC_z"]))
    }
  }
}
T2S <- do.call(rbind, T2S); T3S <- do.call(rbind, T3S); RSS <- do.call(rbind, RSS)
save_vals(T2S, "ST6a_perset_table2.csv", FIG_DIR)

###############################################################################
## 3. ST2 — 성별 × 삼분위 기저 특성 (20_ 의 행 생성기를 6개 열로 확장)
###############################################################################
GL <- levels(B[[1]]$grp)
.g <- function(b) b$grp
.blank <- function(nm) { r <- data.frame(Characteristic = nm, stringsAsFactors = FALSE); for (l in GL) r[[l]] <- ""; r$P_men <- ""; r$P_women <- ""; r }
.medp <- function(ps) { ps <- ps[is.finite(ps)]; if (length(ps)) stats::median(ps) else NA }
p_within <- function(b, x, sexlab, type) {
  k <- grepl(paste0("^", sexlab), as.character(b$grp)) & !is.na(b$grp)
  g <- droplevels(factor(b$tert[k])); xx <- x[k]; ok <- !is.na(xx) & !is.na(g)
  if (sum(ok) < 10 || nlevels(g) < 2) return(NA)
  p <- if (type == "cont") try(stats::kruskal.test(xx[ok] ~ g[ok])$p.value, silent = TRUE) else {
    tb <- try(table(xx[ok], g[ok]), silent = TRUE)
    if (!inherits(tb, "try-error") && all(dim(tb) >= c(2, 2))) try(suppressWarnings(stats::chisq.test(tb)$p.value), silent = TRUE) else NA }
  if (inherits(p, "try-error")) NA else p
}
row_cont <- function(var, nm, dg = 1) {
  st <- lapply(B, function(b) { x <- as_num(b[[var]]); g <- .g(b)
    c(unlist(lapply(GL, function(l) c(mean(x[g == l], na.rm = TRUE), stats::sd(x[g == l], na.rm = TRUE)))),
      p_within(b, x, "Men", "cont"), p_within(b, x, "Women", "cont")) })
  A <- colMeans(do.call(rbind, st)); r <- data.frame(Characteristic = nm, stringsAsFactors = FALSE)
  for (k in seq_along(GL)) r[[GL[k]]] <- sprintf("%.*f (%.*f)", dg, A[2*k - 1], dg, A[2*k])
  r$P_men <- fmt_p(.medp(vapply(st, function(z) z[length(z) - 1], 0))); r$P_women <- fmt_p(.medp(vapply(st, function(z) z[length(z)], 0))); r
}
row_bin <- function(var_or_x, nm, show_p = TRUE) {
  st <- lapply(B, function(b) {
    x <- if (is.character(var_or_x)) as.integer(as_num(b[[var_or_x]])) else as.integer(var_or_x(b)); g <- .g(b)
    c(unlist(lapply(GL, function(l) { xi <- x[g == l]; c(sum(xi == 1, na.rm = TRUE), sum(!is.na(xi))) })),
      if (show_p) c(p_within(b, x, "Men", "bin"), p_within(b, x, "Women", "bin")) else c(NA, NA)) })
  A <- colMeans(do.call(rbind, st)); f <- function(n, N) if (N > 0) sprintf("%d (%.1f)", round(n), 100 * n / N) else "-"
  r <- data.frame(Characteristic = nm, stringsAsFactors = FALSE)
  for (k in seq_along(GL)) r[[GL[k]]] <- f(A[2*k - 1], A[2*k])
  r$P_men <- if (show_p) fmt_p(.medp(vapply(st, function(z) z[length(z) - 1], 0))) else ""
  r$P_women <- if (show_p) fmt_p(.medp(vapply(st, function(z) z[length(z)], 0))) else ""; r
}
row_multi <- function(var, nm, lev, lab = lev) {
  h <- .blank(nm)
  pm <- vapply(B, function(b) { x <- as.character(b[[var]]); x[!x %in% as.character(lev)] <- NA; p_within(b, x, "Men", "bin") }, 0)
  pw <- vapply(B, function(b) { x <- as.character(b[[var]]); x[!x %in% as.character(lev)] <- NA; p_within(b, x, "Women", "bin") }, 0)
  h$P_men <- fmt_p(.medp(pm)); h$P_women <- fmt_p(.medp(pw))
  rs <- lapply(seq_along(lev), function(k) row_bin(function(b) { x <- as.character(b[[var]]); as.integer(x == as.character(lev[k])) },
                                                   paste0("    ", lab[k]), show_p = FALSE))
  do.call(rbind, c(list(h), rs))
}
has <- function(v) v %in% names(B[[1]]) && any(is.finite(as_num(B[[1]][[v]])))
opt_bin <- function(v, nm) if (has(v)) row_bin(v, nm) else NULL
R2 <- Filter(Negate(is.null), list(
  row_cont("age_unified", "Age, years"),
  row_cont("phys_bmi", "Body mass index, kg/m2"),
  row_bin ("educ_high", "Education >= middle school"),
  row_bin ("married", "Married or with partner"),
  row_bin ("living_alone", "Living alone"),
  row_bin ("income_high", "Household income, highest tertile"),
  row_multi("area_chr", "Area of residence", as.character(1:3), c("Metropolitan", "Urban", "Rural")),
  row_cont("comorbid_count", "Comorbidity count"),
  row_bin ("polypharmacy", "Polypharmacy (>=5 medications)"),
  row_bin ("smoking_current", "Current smoker"),
  opt_bin ("alcohol_heavy", "Alcohol >=2-3 times per week, among current drinkers"),
  row_cont("pa_met_min_wk", "Physical activity, MET-min per week", 0),
  row_bin ("strength_2plus", "Resistance exercise >=2 times per week"),
  opt_bin ("social_activity", "Social participation"),
  row_cont("gLIC_z", "Intrinsic capacity, within-sex z-score", 2),
  row_cont("gLIC", "Intrinsic capacity, factor score", 2),
  row_multi("state_chr", "ADL disability", CFG$STATE_LAB[1:3]),
  row_multi("frailty_chr", "Frailty phenotype", CFG$FRAILTY_LAB[1:3]),
  row_cont("followup_years", "Follow-up, years"),
  row_cont("n_tr_adl", "ADL state transitions per participant"),
  row_cont("n_tr_frailty", "Frailty state transitions per participant"),
  row_bin ("death_event", "Deaths during follow-up")))
TAB2 <- do.call(rbind, R2)
TAB2 <- TAB2[, c("Characteristic", paste("Men", LV), "P_men", paste("Women", LV), "P_women")]
nn <- round(rowMeans(sapply(B, function(b) as.integer(table(.g(b))))))
names(TAB2) <- c("Characteristic", sprintf("%s (n = %d)", GL[1:3], nn[1:3]), "P (men)", sprintf("%s (n = %d)", GL[4:6], nn[4:6]), "P (women)")
save_table(TAB2, "ST1_Baseline_bySex_MI",
  title = "Supplementary Table 1 | Baseline characteristics and follow-up status by sex and sex-specific tertiles of intrinsic capacity at wave 1",
  footnotes = c(
    sprintf("Values are mean (s.d.) or n (%%), combined across %d multiply imputed datasets; counts are rounded means. Tertiles are defined within sex on the wave-1 distribution of the general-factor score (Methods), so that T1-T3 denote relative positions within each sex.", M),
    "P values compare tertiles within men and within women separately (Kruskal-Wallis or chi-square; median across imputed datasets).",
    "Education, household income and area of residence include values imputed within each completed dataset (Methods); alcohol intake and social participation are observed wave-1 values (percentages use the number of participants with a non-missing value as the denominator).",
    "The unstandardized factor score shows that the same within-sex tertile corresponds to different absolute capacity in women and men.",
    "ADL, activities of daily living; MET, metabolic equivalent of task; s.d., standard deviation."))

OUTLAB <- c(wor = "Composite worsening", dth = "Death")
###############################################################################
## 5. ST19 — 성별을 공변량으로 추가한 민감도
###############################################################################
T2S$b <- log(T2S$IRR); T2S$se <- (log(T2S$ucl) - log(T2S$lcl)) / (2 * qnorm(0.975))
T2S$key <- paste(T2S$model, T2S$system, T2S$transition, T2S$term, sep = "|")
P2S <- do.call(rbind, lapply(split(T2S, T2S$key), function(z) { r <- pool_rows(z$b, z$se)
  data.frame(model = z$model[1], system = z$system[1], transition = z$transition[1], term = z$term[1], n_events = round(median(z$n_events)),
             est = r["est"], lcl = r["lcl"], ucl = r["ucl"], p = r["p"], stringsAsFactors = FALSE) }))
keys <- unique(P2S[P2S$model == "sex", c("system", "transition")])
TAB19a <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  g <- function(mod, term) { x <- P2S[P2S$model == mod & P2S$system == keys$system[k] & P2S$transition == keys$transition[k] & P2S$term == term, ]
    if (!nrow(x)) c(NA, NA, NA, NA) else c(x$est, x$lcl, x$ucl, x$p) }
  a <- g("sex", "gLIC_z"); a0 <- g("primary", "gLIC_z"); t3 <- g("sex", "gLIC_tertT3"); t30 <- g("primary", "gLIC_tertT3")
  ev <- P2S$n_events[P2S$model == "sex" & P2S$system == keys$system[k] & P2S$transition == keys$transition[k] & P2S$term == "gLIC_z"]
  data.frame(System = ifelse(keys$system[k] == "state", "ADL disability", "Frailty phenotype"), Transition = keys$transition[k], Events = ev[1],
             `IRR per +1 s.d., primary (95% CI)` = fmt_ci(a0[1], a0[2], a0[3]), `IRR per +1 s.d., + sex (95% CI)` = fmt_ci(a[1], a[2], a[3]), P = fmt_p(a[4]),
             `T3 vs T1, primary (95% CI)` = fmt_ci(t30[1], t30[2], t30[3]), `T3 vs T1, + sex (95% CI)` = fmt_ci(t3[1], t3[2], t3[3]),
             `Change in log IRR, %` = ifelse(is.finite(a[1]) & is.finite(a0[1]) & abs(a0[1]) > 0.05, sprintf("%+.0f", 100 * (a[1] - a0[1]) / abs(a0[1])), "-"),
             check.names = FALSE, stringsAsFactors = FALSE) }))
TAB19a <- TAB19a[order(match(TAB19a$System, c("ADL disability", "Frailty phenotype"))), ]
T3S$key <- paste(T3S$outcome, T3S$cell); RSS$key <- paste(RSS$stratum, RSS$outcome)
TAB19b <- rbind(
  do.call(rbind, lapply(split(T3S, T3S$key), function(z) { r <- pool_rows(z$b_hr, z$se_hr)
    data.frame(Analysis = "Table 3 cell versus robust T3", Outcome = OUTLAB[z$outcome[1]], Contrast = z$cell[1], n = round(median(z$n)), Events = round(median(z$events)),
               `HR + sex (95% CI)` = fmt_ci(r["est"], r["lcl"], r["ucl"]), P = fmt_p(r["p"]), check.names = FALSE) })),
  do.call(rbind, lapply(split(RSS, RSS$key), function(z) { r <- pool_rows(z$b, z$se)
    data.frame(Analysis = "Within-category, per +1 s.d.", Outcome = OUTLAB[z$outcome[1]], Contrast = z$stratum[1], n = round(median(z$n)), Events = round(median(z$events)),
               `HR + sex (95% CI)` = fmt_ci(r["est"], r["lcl"], r["ucl"]), P = fmt_p(r["p"]), check.names = FALSE) })))
save_table(TAB19a, "ST6a_SexCovariate_Table2_MI",
  title = "Supplementary Table 6a | Transition rates per +1 s.d. of within-sex standardized intrinsic capacity with and without sex as an additional covariate",
  footnotes = c(
    "Primary: the models of Table 2 (within-sex standardization; adjusted for age, baseline comorbidity count, education, household income and area of residence). + sex: the same models with sex added as a covariate.",
    "Within-sex standardization equalizes the distribution of the exposure between men and women at wave 1 but does not remove any independent association of sex with the outcome; adding sex tests whether such an association changes the estimate for intrinsic capacity.",
    sprintf("Estimates combined across %d completed datasets by Rubin's rules.", M)))
save_table(TAB19b, "ST6b_SexCovariate_Table3_MI",
  title = "Supplementary Table 6b | Five-year hazard ratios within frailty states with sex added as a covariate",
  footnotes = "Cox models as in Table 3 with sex added to the covariates; compare with the primary estimates in Table 3.")

###############################################################################
## 6. 판정 출력
###############################################################################
cat("\nST19a — 성별 공변량 추가 시 log IRR 변화 (%):\n"); print(TAB19a[, c("Transition", "IRR per +1 s.d., primary (95% CI)", "IRR per +1 s.d., + sex (95% CI)", "Change in log IRR, %")], row.names = FALSE)
chg <- suppressWarnings(as.numeric(TAB19a[["Change in log IRR, %"]])); cat(sprintf("  |변화| 중앙값 %.1f%%, 최대 %.1f%%\n", median(abs(chg), na.rm = TRUE), max(abs(chg), na.rm = TRUE)))
cat(strrep("=", 70), "\n=== 완료 ===\n")
