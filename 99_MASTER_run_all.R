###############################################################################
## 99_MASTER_run_all.R — 전체 파이프라인 실행 (이 파일 하나만 Source)
##
## KFACS intrinsic-capacity longitudinal study — one-click reproduction.
## 실행 전 확인할 것은 R/00_setup.R 의 경로 세 줄뿐입니다.
##
##  · 실패해도 나머지는 계속 진행하고 마지막에 요약을 냅니다.
##  · 콘솔 전체가 RUN_log_<날짜시각>.txt 로 저장됩니다.
##  · 스크립트마다 "이번에 만들어진 파일 수" 를 기록합니다 — 오류 없이 끝났는데
##    산출물이 0개면 요약표에 '의심' 으로 표시됩니다 (조용한 실패 방지).
##  · Expected run time: about 1-2 hours (12_/13_/13b_ are skipped when their outputs exist)
###############################################################################
options(warn = 1)

## ── 스크립트 위치 결정 ───────────────────────────────────────────────────
## 위치 탐색과 헬퍼 로드는 R/_bootstrap.R 한 곳에서만 정의합니다.
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
## ── 실행 목록 (원고 순서) ───────────────────────────────────────────────
## 빼고 싶은 것이 있으면 SKIP 에 파일명을 적으십시오.
##   예: SKIP <- c("03_figures/32_figure4_JM_MI.R")   # JM 이 오래 걸릴 때
SKIP <- c()

## 한 번만 만들면 되는 산출물이 이미 있으면 그 단계는 자동으로 건너뜁니다.
## 다시 만들려면 해당 파일을 지우거나 FORCE_REBUILD <- TRUE 로 두십시오.
FORCE_REBUILD <- FALSE
ONCE <- list(
  "01_invariance/12_sex_invariance_and_scores.R"     = file.path(FIG_DIR, "P1_scores.rds"),
  "01_invariance/13_mi_pmm_m20.R"                    = file.path(IMP_DIR, paste0(MI_STEM, ".rds")),
  "01_invariance/13b_impute_baseline_covariates.R"   = file.path(IMP_DIR, "KFACS_cov_imputed_m20.rds"),
  "01_invariance/18b_longitudinal_cfa_pairs.R"       = file.path(FIG_DIR, "LongCFA_pairs_fits.rds"))

## v1.1.0 (2026-11-17): run list = display items of the revised manuscript
## ("Intrinsic capacity as a marker of resilience across the frailty spectrum").
## Analyses no longer reported are kept, unchanged, in archive_not_reported/.
SCRIPTS <- c(
  ## -- A. Measurement model and imputation (inputs to everything below) ------
  "01_invariance/12_sex_invariance_and_scores.R",    # Supplementary Tables 2-3 (bifactor model, sex invariance); within-sex scale
  "01_invariance/13_mi_pmm_m20.R",                   # mice PMM m = 20 -> IMP_DIR/<MI_STEM>.rds  [30-90 min]
  "01_invariance/13b_impute_baseline_covariates.R",  # education, income, area imputed within each completed dataset
  "01_invariance/18b_longitudinal_cfa_pairs.R",      # adjacent-wave longitudinal CFA
  "01_invariance/18c_longitudinal_cfa_tables.R",     #   -> Supplementary Table 4
  ## -- B. Main tables ----------------------------------------------------------
  "02_tables/20_table1_baseline_MI.R",               # Table 1
  "01_invariance/14_mi_pool_tables.R",               # Table 2, Table 3 (cells)
  "02_tables/21_st1_bysex_and_st6ab_sex_MI.R",       # Supplementary Tables 1, 6a, 6b
  ## -- C. Figures and recovery ------------------------------------------------
  "03_figures/30_figure1_transition_rates_MI.R",     # Figure 1
  "03_figures/31_figure2_recovery_MI.R",             # Table 4 (pre-frail and mild-disability rows); Supplementary Fig. 1
  "03_figures/34_frail_recovery_by_criteria_MI.R",   # Table 4 (frail row); Supplementary Table 5
  "03_figures/35_edfig1_hospitalization_fall_MI.R",  # Extended Data Fig. 1; Supplementary Table 7
  "03_figures/43_edfig1_subgroups_MI.R",             # Supplementary Fig. 2 (subgroups, set MI01)
  "04_sensitivity/57_edfig2_wave_stability_MI.R",    # Supplementary Fig. 3 (waves, set MI01)
  ## -- D. Sensitivity ----------------------------------------------------------
  "04_sensitivity/62_sens_additional_confounders.R", # Supplementary Tables 6c, 6d; Table 3 per-s.d. rows (primary column)
  "04_sensitivity/64_sens_covariate_missing_category.R", # Supplementary Table 6e
  ## -- E. Reporting ------------------------------------------------------------
  "05_reporting/90_collect_captions.R")              # collect captions -- always last

## ── 콘솔 로그 ───────────────────────────────────────────────────────────
STAMP    <- format(Sys.time(), "%y%m%d_%H%M")
LOG_PATH <- file.path(FIG_DIR, paste0("RUN_log_", STAMP, ".txt"))
.logcon  <- file(LOG_PATH, open = "wt", encoding = "UTF-8")
sink(.logcon, split = TRUE)
sink(.logcon, type = "message", append = TRUE)
on.exit({ sink(type = "message"); sink(); try(close(.logcon), silent = TRUE) }, add = TRUE)
cat("KFACS 파이프라인 시작 ", format(Sys.time()), "\n\n", sep = "")

## ── 사전 점검 ───────────────────────────────────────────────────────────
miss <- SCRIPTS[!file.exists(file.path(KF_REPO, SCRIPTS))]
if (length(miss))
  stop("다음 스크립트가 없습니다. 패키지가 완전한지 확인하십시오:\n",
       paste0("    - ", miss, collapse = "\n"))
## 13b_ output must be in the v1.1.0 format; an older file (written by 65_) is removed
.cf <- file.path(IMP_DIR, "KFACS_cov_imputed_m20.rds")
if (file.exists(.cf)) {
  .ok <- try(all(c("edu_high_bl", "income_high", "area_3cat") %in% names(readRDS(.cf)$sets[[1]])), silent = TRUE)
  if (!isTRUE(.ok)) {
    file.rename(.cf, sub("\\.rds$", "_old_format.rds", .cf))
    cat("Covariate imputation file was in an old format -> renamed; 13b_ will run.\n")
  }
}
for (s in names(ONCE))
  if (!FORCE_REBUILD && s %in% SCRIPTS && file.exists(ONCE[[s]]) && !(s %in% SKIP)) {
    cat("이미 있음 -> 건너뜀: ", s, "  (", basename(ONCE[[s]]), ")\n", sep = "")
    SKIP <- c(SKIP, s)
  }
need_pkg <- c("ggplot2", "patchwork", "scales", "sandwich", "lmtest", "survival",
              "lavaan", "mice", "nnet", "MASS", "openxlsx", "readxl", "officer", "flextable")
have <- vapply(need_pkg, requireNamespace, logical(1), quietly = TRUE)
if (any(!have)) {
  cat("*** 설치되지 않은 패키지: ", paste(need_pkg[!have], collapse = ", "), "\n")
  try(utils::install.packages(need_pkg[!have], repos = "https://cloud.r-project.org"),
      silent = TRUE)
}

## ── 실행 ────────────────────────────────────────────────────────────────
.outputs_since <- function(t0) {
  f <- list.files(FIG_DIR, full.names = TRUE, recursive = TRUE)
  f[file.info(f)$mtime >= t0 & !grepl("RUN_log_", f)]
}
T_ALL <- Sys.time()
res <- data.frame(script = SCRIPTS, status = NA_character_, minutes = NA_real_,
                  outputs = NA_integer_, detail = "", stringsAsFactors = FALSE)
for (i in seq_along(SCRIPTS)) {
  s <- SCRIPTS[i]
  cat(strrep("=", 70), "\n[", i, "/", length(SCRIPTS), "] ", s,
      "   (", format(Sys.time(), "%H:%M:%S"), ")\n", strrep("=", 70), "\n", sep = "")
  if (s %in% SKIP) { res$status[i] <- if (s %in% names(ONCE)) "SKIP(산출물 있음)" else "SKIP(사용자)"; next }
  t0 <- Sys.time()
  r  <- try(source(file.path(KF_REPO, s), encoding = "UTF-8", local = new.env()),
            silent = TRUE)
  res$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  nf <- .outputs_since(t0)
  res$outputs[i] <- length(nf)
  if (inherits(r, "try-error")) {
    res$status[i] <- "실패"
    res$detail[i] <- gsub("[\r\n]+", " ", conditionMessage(attr(r, "condition")))
    cat("\n*** 실패: ", res$detail[i], "\n\n", sep = "")
  } else if (!length(nf) && !(s %in% c("01_invariance/13_mi_pmm_m20.R", "01_invariance/13b_impute_baseline_covariates.R"))) {
    res$status[i] <- "성공(산출물 0 — 의심)"
  } else res$status[i] <- "성공"
}

## ── 요약 ────────────────────────────────────────────────────────────────
tot <- round(as.numeric(difftime(Sys.time(), T_ALL, units = "mins")), 1)
cat("\n", strrep("=", 70), "\n요약   (총 ", tot, "분, 종료 ",
    format(Sys.time(), "%m-%d %H:%M"), ")\n", strrep("=", 70), "\n", sep = "")
print(res[, c("script", "status", "minutes", "outputs")], row.names = FALSE)
if (any(res$status == "실패")) {
  cat("\n실패 상세\n")
  for (i in which(res$status == "실패"))
    cat(sprintf("  %s\n    %s\n", res$script[i], res$detail[i]))
}
try(utils::write.csv(res, file.path(FIG_DIR, "RUN_summary.csv"),
                     row.names = FALSE, fileEncoding = "UTF-8"), silent = TRUE)

out <- list.files(FIG_DIR, pattern = "\\.(pdf|tiff|png|xlsx|csv|md|docx)$",
                  full.names = TRUE)
out <- out[file.info(out)$mtime >= T_ALL]
cat("\n이번 실행에서 만들어진 파일 ", length(out), "개  (", FIG_DIR, ")\n", sep = "")
if (length(out)) {
  o <- data.frame(file = basename(out), KB = round(file.info(out)$size / 1024, 1))
  print(o[order(o$file), ], row.names = FALSE)
}
cat("\n  콘솔 로그 : ", LOG_PATH, "\n", sep = "")
cat("  요약표    : ", file.path(FIG_DIR, "RUN_summary.csv"), "\n", sep = "")

sink(type = "message"); sink(); try(close(.logcon), silent = TRUE)
message("완료. 로그: ", LOG_PATH)
