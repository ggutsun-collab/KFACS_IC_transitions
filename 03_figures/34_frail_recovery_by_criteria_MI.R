###############################################################################
## 34_frail_recovery_by_criteria_MI.R                                   v261117
##
## Table 4 (frail row) and Supplementary Table 5. Among participants frail at
## wave 1, five-year improvement (first visit in the pre-frail or robust state)
## above versus below the within-frail median of IC, overall and by the number
## of frailty criteria met at wave 1 (3 versus 4-5).
##   - cumulative incidence: Aalen-Johansen, death as a competing event,
##     pooled over the 20 completed datasets on the complementary log-log scale
##   - hazard ratio: cause-specific Cox, adjusted for the primary covariates
## The within-frail median is taken over all frail participants in each
## completed dataset, so the subgroups share one split.
## Reproduction check (missing-category covariates, 2026-10-05 run):
##   all frail 73.8 vs 44.4%, HR 2.05 (1.38-3.05); 3 criteria 73.5 vs 47.0%,
##   HR 1.82 (1.09-3.02); 4-5 criteria 76.0 vs 40.8%, HR 2.50 (1.07-5.83).
###############################################################################

message("\n=== 34_frail_recovery_by_criteria_MI v261117b ===")
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
suppressPackageStartupMessages(library(survival))

TAU  <- 5
SETS <- names(readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))))
GRPS <- c("All frail", "3 criteria", "4-5 criteria")

## first improvement from frailty; death competing; censored at last observation.
## Same definition as 31_ (Table 4, frail row): follow-up is not truncated, the
## cumulative incidence is read at TAU, and the last observation includes death.
frail_events <- function(d, ids) {
  d$fst <- as.integer(as_num(d$frail_state))
  x <- d[d$id %in% ids & is.finite(d$time) & is.finite(d$fst), ]; x <- x[order(x$id, x$time), ]
  do.call(rbind, lapply(split(x, x$id), function(p) {
    tt <- p$time - p$time[1]; s <- p$fst
    k_imp <- which(seq_along(s) > 1 & s %in% 1:2)[1]
    k_dth <- which(s %in% 4)[1]
    t_imp <- if (!is.na(k_imp)) tt[k_imp] else Inf
    t_dth <- if (!is.na(k_dth)) tt[k_dth] else Inf
    ft <- min(t_imp, t_dth, max(tt))
    ev <- if (t_imp <= ft) 1L else if (t_dth <= ft) 2L else 0L
    data.frame(id = p$id[1], fti = max(ft, 1e-6), evi = ev)
  }))
}
cif_at <- function(ft, ev) {
  f <- survfit(Surv(ft, factor(ev, 0:2, c("censor", "improve", "death"))) ~ 1)
  s <- summary(f, times = TAU, extend = TRUE); j <- match("improve", f$states)
  p <- min(max(s$pstate[1, j], 1e-6), 1 - 1e-6); se <- s$std.err[1, j]
  u <- -log(1 - p)
  c(g = log(u), se = se / ((1 - p) * u))      # complementary log-log scale
}
icll <- function(g) 100 * (1 - exp(-exp(g)))

R <- list()
for (set in SETS) {
  load_imputed(set, stem = MI_STEM)
  d <- ws_scale(prep_long(load_kfacs()))
  K <- sm_person_frame(d); F <- K[K$fr == "Frail", ]
  w1 <- d[!duplicated(d$id), c("id", "chs_total")]
  F$ncrit <- as_num(w1$chs_total[match(F$id, w1$id)])
  F$above <- as.integer(F$z > stats::median(F$z))
  F <- merge(F, frail_events(d, F$id), by = "id")
  for (g in GRPS) {
    x <- switch(g, "All frail" = F, "3 criteria" = F[F$ncrit %in% 3, ], "4-5 criteria" = F[F$ncrit >= 4, ])
    for (a in 0:1) {
      z <- x[x$above == a, ]; c5 <- cif_at(z$fti, z$evi)
      R[[length(R) + 1]] <- data.frame(set = set, grp = g, above = a, n = nrow(z), events = sum(z$evi == 1),
                                       deaths = sum(z$evi == 2), g = c5["g"], se = c5["se"], b = NA, v = NA)
    }
    m <- try(coxph(as.formula(paste("Surv(fti, evi == 1) ~ above +", paste(SM_COV_PRIMARY, collapse = " + "))),
                   data = droplevels(x)), silent = TRUE)
    if (!inherits(m, "try-error"))
      R[[length(R) + 1]] <- data.frame(set = set, grp = g, above = 9, n = nrow(x), events = NA, deaths = NA,
                                       g = NA, se = NA, b = coef(m)["above"], v = vcov(m)["above", "above"])
  }
  msg("[", set, "] frail n = ", nrow(F))
}
R <- do.call(rbind, R)
OUT <- do.call(rbind, lapply(GRPS, function(g) {
  r <- function(a) { z <- R[R$grp == g & R$above == a, ]; p <- sm_rubin(z$g, z$se^2)
    c(n = median(z$n), ev = median(z$events), dt = median(z$deaths),
      cif = unname(icll(p["est"])), lo = unname(icll(p["lo"])), hi = unname(icll(p["hi"]))) }
  b0 <- r(0); b1 <- r(1); h <- R[R$grp == g & R$above == 9, ]; ph <- sm_rubin(h$b, h$v)
  data.frame(Group = g,
             `Below median: n` = b0["n"], `Below: improved, % (95% CI)` = sprintf("%.1f (%.1f-%.1f)", b0["cif"], b0["lo"], b0["hi"]),
             `Above median: n` = b1["n"], `Above: improved, % (95% CI)` = sprintf("%.1f (%.1f-%.1f)", b1["cif"], b1["lo"], b1["hi"]),
             `HR above vs below (95% CI)` = sm_fmt(ph), P = fmt_p(ph["p"]),
             `Events below/above` = paste(b0["ev"], b1["ev"], sep = "/"), `Deaths below/above` = paste(b0["dt"], b1["dt"], sep = "/"),
             check.names = FALSE, row.names = NULL)
}))
cat("\n=== Supplementary Table 5 (and Table 4, frail row = 'All frail') ===\n"); print(OUT, row.names = FALSE)
save_table(OUT, "ST5_FrailRecovery_byCriteria_MI",
  title = "Supplementary Table 5 | Five-year improvement from frailty above versus below the within-frail median of intrinsic capacity, by number of frailty criteria at wave 1",
  footnotes = c("Improvement is the first visit in the pre-frail or robust state. Cumulative incidence by the Aalen-Johansen estimator with death as a competing event, combined across 20 completed datasets on the complementary log-log scale; hazard ratios are cause-specific, adjusted for age, baseline comorbidity count, education, household income and area of residence.",
                "The median of intrinsic capacity is taken over all frail participants in each completed dataset; n, events and deaths are medians across datasets."))
