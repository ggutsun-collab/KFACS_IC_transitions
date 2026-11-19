###############################################################################
## 04_cov_imputation.R - imputed baseline covariates, applied at load time
##
## Education (edu_high_bl), household income (income_high) and area of
## residence (area_3cat) were missing at wave 1 for 2, 192 and 16 participants.
## 01_invariance/13b_impute_baseline_covariates.R imputes them once within each
## of the 20 completed datasets and writes IMP_DIR/KFACS_cov_imputed_m20.rds.
##
## load_imputed() (R/02_theme_and_output.R) passes every multiply imputed set
## through apply_cov_imputation(), so prep_long() derives edu_bl_f,
## income_bl_f, area_bl_f and ses_index from the imputed values. The kNN
## single-imputation sets (IMP_STEM) are not touched.
##
## options(kf.cov_impute = FALSE) restores the earlier coding, in which a
## missing value is an explicit "Unknown" level. 13b_ uses it to impute from
## the original data, and 04_sensitivity/64_ uses it for the comparison
## reported in Supplementary Table 6e.
###############################################################################

COVIMP_VERSION <- "2026-11-17"
COVIMP_FILE    <- file.path(IMP_DIR, "KFACS_cov_imputed_m20.rds")
COVIMP_VARS    <- c("edu_high_bl", "income_high", "area_3cat")
.covimp_cache  <- new.env()

apply_cov_imputation <- function(d, set, stem = MI_STEM) {
  if (isFALSE(getOption("kf.cov_impute", TRUE))) return(d)
  if (!identical(stem, MI_STEM)) return(d)
  if (!file.exists(COVIMP_FILE)) {
    message("[cov-imputation] ", basename(COVIMP_FILE),
            " not found: missing-category coding. Run 01_invariance/13b_ first.")
    return(d)
  }
  mt <- file.mtime(COVIMP_FILE)
  if (is.null(.covimp_cache$x) || !identical(.covimp_cache$mtime, mt)) {
    .covimp_cache$x <- readRDS(COVIMP_FILE); .covimp_cache$mtime <- mt
  }
  CI  <- .covimp_cache$x
  if (!all(COVIMP_VARS %in% names(CI$sets[[1]])))
    stop("[cov-imputation] ", COVIMP_FILE, " is not in the v1.1.0 format (expected columns ",
         paste(COVIMP_VARS, collapse = ", "), "). It was probably written by an earlier ",
         "script. Delete it and run 01_invariance/13b_impute_baseline_covariates.R.")
  key <- if (is.numeric(set)) names(CI$sets)[set] else as.character(set)
  imp <- CI$sets[[key]]
  if (is.null(imp)) {
    message("[cov-imputation] set ", key, " not in ", basename(COVIMP_FILE), ": unchanged")
    return(d)
  }
  idx <- match(as.character(d$id), imp$id)
  n <- setNames(integer(length(COVIMP_VARS)), COVIMP_VARS)
  for (v in COVIMP_VARS) {
    if (!v %in% names(d) || !v %in% names(imp)) next
    miss <- is.na(d[[v]]) & !is.na(idx)
    d[[v]][miss] <- imp[[v]][idx[miss]]
    n[v] <- length(unique(d$id[miss]))
  }
  message(sprintf("[cov-imputation] %s: imputed for %d / %d / %d participants (education / income / area)",
                  key, n[1], n[2], n[3]))
  d
}
