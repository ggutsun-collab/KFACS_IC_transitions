###############################################################################
## 64_sens_covariate_missing_category.R                                 v261117
##
## Supplementary Table 6e. Compares the primary analysis (education, household
## income and area of residence imputed; 13b_) with the earlier coding in which
## a missing value is an explicit "Unknown" category, for
##   - the 18 transition rate ratios per +1 s.d. of IC (Table 2)
##   - the within-state 5-year HRs per +1 s.d. (Table 3)
## Change on the log scale = 100 x (log missing-category - log imputed) / |log imputed|.
###############################################################################

message("\n=== 64_sens_covariate_missing_category v261117 ===")
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
RA <- list(); RB <- list()
fit_all <- function(imputed) {
  options(kf.cov_impute = imputed)
  load_imputed(set, stem = MI_STEM)
  d <- ws_scale(prep_long(load_kfacs())); K <- sm_person_frame(d); iv <- sm_intervals(d, K)
  list(tr = lapply(seq_len(nrow(SM_TR)), function(j) sm_fit_transition(iv, j, SM_COV_PRIMARY)),
       wi = lapply(SM_WITHIN, function(cl) sm_fit_within(K, cl[1], cl[2], SM_COV_PRIMARY)))
}
for (set in SETS) {
  I <- fit_all(TRUE); M <- fit_all(FALSE)
  for (j in seq_len(nrow(SM_TR)))
    RA[[length(RA) + 1]] <- data.frame(j = j, ev = I$tr[[j]]["ev"], b0 = I$tr[[j]]["b"], v0 = I$tr[[j]]["v"],
                                       b1 = M$tr[[j]]["b"], v1 = M$tr[[j]]["v"])
  for (k in seq_along(SM_WITHIN))
    RB[[length(RB) + 1]] <- data.frame(k = k, ev = I$wi[[k]]["ev"], b0 = I$wi[[k]]["b"], v0 = I$wi[[k]]["v"],
                                       b1 = M$wi[[k]]["b"], v1 = M$wi[[k]]["v"])
  msg("[", set, "] done")
}
options(kf.cov_impute = NULL)
pool <- function(z, lab) {
  p0 <- sm_rubin(z$b0, z$v0); p1 <- sm_rubin(z$b1, z$v1)
  data.frame(Estimate = lab, Events = median(z$ev), Imputed_primary = sm_fmt(p0), P_imputed = fmt_p(p0["p"]),
             Missing_category = sm_fmt(p1), P_missing_category = fmt_p(p1["p"]),
             Change_logscale_pct = round(100 * (p1["est"] - p0["est"]) / abs(p0["est"]), 1),
             row.names = NULL, stringsAsFactors = FALSE)
}
A <- do.call(rbind, RA); B <- do.call(rbind, RB)
T6e <- rbind(do.call(rbind, lapply(seq_len(nrow(SM_TR)), function(j) pool(A[A$j == j, ], SM_TR$label[j]))),
             do.call(rbind, lapply(seq_along(SM_WITHIN), function(k)
               pool(B[B$k == k, ], paste("Within-state:", paste(SM_WITHIN[[k]], collapse = " "))))))
cat("\n=== Supplementary Table 6e ===\n"); print(T6e, row.names = FALSE)
cat(sprintf("\nmax |change| on log scale: %.1f%%\n", max(abs(T6e$Change_logscale_pct), na.rm = TRUE)))
utils::write.csv(T6e, file.path(FIG_DIR, "ST6e_covariate_missing_category.csv"), row.names = FALSE)
msg("Saved ST6e to ", FIG_DIR)
options(kf.cov_impute = NULL)   # back to the default (imputed covariates)
