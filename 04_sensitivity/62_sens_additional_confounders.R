###############################################################################
## 62_sens_additional_confounders.R                                     v261117
##
## Supplementary Table 6c,d. Adds living alone, polypharmacy, smoking and
## physical activity (wave 1; missing as a separate category; physical
## activity in tertiles of MET-min/week) to the primary models:
##   c. 18 transition-specific Poisson models, rate ratio per +1 s.d. of IC
##   d. within-state 5-year Cox models (robust, pre-frail), HR per +1 s.d.
## Primary covariates: age, wave-1 comorbidity count, education, household
## income, area of residence (imputed; R/04_cov_imputation.R).
## Attenuation = 100 x (log primary - log extended) / log primary.
###############################################################################

message("\n=== 62_sens_additional_confounders v261117 ===")
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
for (p in c("nnet", "MASS")) if (!requireNamespace(p, quietly = TRUE)) stop("Package missing: ", p)

SETS <- names(readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))))
RA <- list(); RB <- list(); NIV <- integer(0)
for (set in SETS) {
  load_imputed(set, stem = MI_STEM)
  d  <- ws_scale(prep_long(load_kfacs()))
  K  <- sm_person_frame(d); iv <- sm_intervals(d, K); NIV <- c(NIV, nrow(iv))
  for (j in seq_len(nrow(SM_TR))) {
    a <- sm_fit_transition(iv, j, SM_COV_PRIMARY); b <- sm_fit_transition(iv, j, SM_COV_EXTRA)
    RA[[length(RA) + 1]] <- data.frame(set = set, j = j, ev = a["ev"], b0 = a["b"], v0 = a["v"], b1 = b["b"], v1 = b["v"])
  }
  for (cl in SM_WITHIN) {
    a <- sm_fit_within(K, cl[1], cl[2], SM_COV_PRIMARY); b <- sm_fit_within(K, cl[1], cl[2], SM_COV_EXTRA)
    RB[[length(RB) + 1]] <- data.frame(set = set, key = paste(cl, collapse = "|"), n = a["n"], ev = a["ev"],
                                       b0 = a["b"], v0 = a["v"], b1 = b["b"], v1 = b["v"])
  }
  msg(sprintf("[%s] %d intervals", set, nrow(iv)))
}
pool <- function(z, lab) {
  p0 <- sm_rubin(z$b0, z$v0); p1 <- sm_rubin(z$b1, z$v1)
  data.frame(Estimate = lab, Events = median(z$ev), Primary = sm_fmt(p0), P_primary = fmt_p(p0["p"]),
             Extended = sm_fmt(p1), P_extended = fmt_p(p1["p"]),
             Attenuation_pct = round(100 * (p0["est"] - p1["est"]) / p0["est"], 1),
             sig0 = p0["p"] < 0.05, sig1 = p1["p"] < 0.05, row.names = NULL, stringsAsFactors = FALSE)
}
A  <- do.call(rbind, RA); B <- do.call(rbind, RB)
TA <- do.call(rbind, lapply(seq_len(nrow(SM_TR)), function(j) pool(A[A$j == j, ], SM_TR$label[j])))
TB <- do.call(rbind, lapply(unique(B$key), function(k) pool(B[B$key == k, ], k)))
cat("\nIntervals per set:", paste(range(NIV), collapse = "-"), " (Table 2: 11,571)\n")
cat("\n=== Supplementary Table 6c: transitions ===\n"); print(TA[, 1:7], row.names = FALSE)
cat("\n=== Supplementary Table 6d: within-state ===\n"); print(TB[, 1:7], row.names = FALSE)
s <- TA$sig0
cat(sprintf("\nSignificant in primary: %d; still significant: %d; attenuation (significant only) max %.1f%%; within-state %.1f-%.1f%%\n",
            sum(s, na.rm = TRUE), sum(s & TA$sig1, na.rm = TRUE), max(TA$Attenuation_pct[s], na.rm = TRUE),
            min(TB$Attenuation_pct), max(TB$Attenuation_pct)))
utils::write.csv(TA[, 1:7], file.path(FIG_DIR, "ST6c_confounders_transitions.csv"), row.names = FALSE)
utils::write.csv(TB[, 1:7], file.path(FIG_DIR, "ST6d_confounders_withinstate.csv"), row.names = FALSE)
msg("Saved ST6c/ST6d to ", FIG_DIR)
