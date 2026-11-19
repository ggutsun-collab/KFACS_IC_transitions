###############################################################################
## 13b_impute_baseline_covariates.R                                     v261117
##
## Imputes the three baseline covariates that 13_ left missing (education
## edu_high_bl, household income income_high, area of residence area_3cat;
## 2, 192 and 16 participants) ONCE within each of the 20 completed datasets.
## IC, ADL and frailty states and deaths in the completed datasets are not
## changed.
##
## Method: chained logistic (binary) and multinomial (area) regression, one
## draw per completed dataset after 10 iterations, with the regression
## coefficients drawn from their sampling distribution (proper imputation).
## Predictors: age, sex, wave-1 comorbidity count, IC, ADL and frailty state,
## living alone, polypharmacy, smoking, physical activity, the death
## indicator and the Nelson-Aalen cumulative hazard.
##
## Output: IMP_DIR/KFACS_cov_imputed_m20.rds, read by load_imputed() through
## R/04_cov_imputation.R in every later script. Run once, after 13_.
###############################################################################

COVIMP_SCRIPT_VERSION <- "v261117"
message("\n=== 13b_impute_baseline_covariates ", COVIMP_SCRIPT_VERSION, " ===")

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
options(kf.cov_impute = FALSE)        # impute from the original (missing) values

SEED <- 20261117
MAXIT <- 10

## Nelson-Aalen cumulative hazard at each participant's own time
.na_cumhaz <- function(t, e) {
  f <- survfit(Surv(t, e) ~ 1); f$cumhaz[pmax(findInterval(t, f$time), 1)]
}
## proper single draw for one categorical variable
.draw_cat <- function(y, X) {
  obs <- !is.na(y); if (all(obs)) return(y)
  lev <- levels(droplevels(y[obs]))
  dfo <- data.frame(.y = droplevels(y[obs]), X[obs, , drop = FALSE])
  mm  <- model.matrix(~ ., X)
  if (length(lev) == 2) {
    m <- glm(.y ~ ., family = binomial, data = dfo)
    b <- MASS::mvrnorm(1, coef(m), vcov(m)); b[is.na(b)] <- 0
    pr <- 1 / (1 + exp(-(mm[!obs, names(b), drop = FALSE] %*% b)))
    y[!obs] <- lev[1 + (runif(sum(!obs)) < pr)]
  } else {
    m <- nnet::multinom(.y ~ ., data = dfo, trace = FALSE, maxit = 500)
    cf <- coef(m)
    B <- matrix(MASS::mvrnorm(1, as.vector(t(cf)), vcov(m)), nrow = nrow(cf), byrow = TRUE,
                dimnames = dimnames(cf))
    eta <- cbind(0, mm[!obs, colnames(B), drop = FALSE] %*% t(B))
    pr <- exp(eta - apply(eta, 1, max)); pr <- pr / rowSums(pr)
    y[!obs] <- lev[apply(pr, 1, function(p) sample.int(length(p), 1, prob = p))]
  }
  y
}
.impute_chain <- function(J, targets, maxit = MAXIT) {
  targets <- targets[vapply(targets, function(v) anyNA(J[[v]]), TRUE)]
  if (!length(targets)) return(J)
  miss <- lapply(setNames(targets, targets), function(v) is.na(J[[v]]))
  for (v in targets) { o <- J[[v]][!miss[[v]]]; J[[v]][miss[[v]]] <- sample(o, sum(miss[[v]]), TRUE) }
  for (it in seq_len(maxit)) for (v in targets) {
    y <- J[[v]]; y[miss[[v]]] <- NA
    J[[v]] <- .draw_cat(y, J[, setdiff(names(J), v), drop = FALSE])
  }
  J
}

SETS <- names(readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))))
STORE <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]
  raw <- load_imputed(set, stem = MI_STEM)
  d   <- ws_scale(prep_long(load_kfacs()))
  K   <- sm_person_frame(d)
  w1r <- raw[!duplicated(as.character(raw$id)), ]          # first row per participant
  rv  <- w1r[match(K$id, as.character(w1r$id)), COVIMP_VARS]
  J <- data.frame(edu_high_bl = factor(rv$edu_high_bl), income_high = factor(rv$income_high),
                  area_3cat = factor(rv$area_3cat),
                  age_c = K$age_c, sex_f = K$sex_f, comorbid_bl = K$comorbid_bl, z = K$z,
                  adl0 = factor(K$adl0), fr = K$fr, alone = K$alone, poly = K$poly,
                  smoke = K$smoke, pa = K$pa, ev = K$ev, H = .na_cumhaz(K$t_fu, K$ev))
  J <- J[, vapply(J, function(x) length(unique(na.omit(x))) > 1, TRUE)]
  set.seed(SEED + i)
  Jc <- .impute_chain(J, intersect(COVIMP_VARS, names(J)))
  out <- data.frame(id = K$id, stringsAsFactors = FALSE)
  for (v in COVIMP_VARS) {
    out[[v]] <- if (v %in% names(Jc)) as.numeric(as.character(Jc[[v]])) else as.numeric(rv[[v]])
    out[[paste0(v, "_was_missing")]] <- is.na(rv[[v]])
  }
  STORE[[set]] <- out
  msg(sprintf("[%s] imputed: education %d, income %d, area %d", set,
              sum(is.na(rv$edu_high_bl)), sum(is.na(rv$income_high)), sum(is.na(rv$area_3cat))))
}

saveRDS(list(sets = STORE, vars = COVIMP_VARS, seed = SEED, maxit = MAXIT,
             method = "chained logistic/multinomial regression, 1 draw per completed dataset",
             created = format(Sys.time()), stem = MI_STEM), COVIMP_FILE)
cat("\nObserved vs imputed distribution, set ", SETS[1], "\n", sep = "")
x <- STORE[[1]]
for (v in COVIMP_VARS) {
  m <- x[[paste0(v, "_was_missing")]]
  print(rbind(observed = table(factor(x[[v]][!m], levels = sort(unique(x[[v]])))),
              imputed  = table(factor(x[[v]][m],  levels = sort(unique(x[[v]]))))))
}
msg("Saved: ", COVIMP_FILE)
options(kf.cov_impute = NULL)   # back to the default (imputed covariates)
