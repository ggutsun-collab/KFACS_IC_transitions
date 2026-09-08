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
##  · 예상 소요: 2-4시간 (13_ MI, 15_ 부트스트랩, 16_ msm, 32_ JM 포함 시; 12_/13_ 은 산출물이 있으면 건너뜀)
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
  "01_invariance/12_sex_invariance_and_scores.R" = file.path(FIG_DIR, "P1_scores.rds"),
  "01_invariance/13_mi_pmm_m20.R"                = file.path(IMP_DIR, paste0(MI_STEM, ".rds")))

SCRIPTS <- c(
  ## ── A. 측정모형과 대체 (주 분석의 입력을 만듭니다) ──────────────────────
  "01_invariance/12_sex_invariance_and_scores.R",  # 성별 다집단 bifactor, 잠재평균 차, within-sex 척도 근거 (kNN 세트 입력)
  "01_invariance/13_mi_pmm_m20.R",                 # mice PMM m=20 (wide) -> IMP_DIR/<MI_STEM>.rds  [30-90분]
  "01_invariance/10_measurement_invariance_v3.R",  # wave 다집단 불변성 (관측자료 + FIML, 4 특정요인)
  "01_invariance/11_invariance_sensitivity.R",     #   불변모형 점수로 연관 재추정
  "01_invariance/18b_longitudinal_cfa_pairs.R",    # 인접 wave 쌍 종단 CFA (within-person dependence)
  "01_invariance/18c_longitudinal_cfa_tables.R",   #   -> Supplementary Table (longitudinal invariance)
  ## ── B. 주 분석 (MI 20세트, within-sex 척도) ─────────────────────────────
  "01_invariance/14_mi_pool_tables.R",             # Table 2, Table 3, ED Table 1 (Rubin 풀링)
  "02_tables/20_table1_baseline_MI.R",             # Table 1
  "01_invariance/15_ic13_headtohead.R",            # IC-13 유지율, robust 층 ΔC/IDI 대 grip+gait  [20-40분]
  "01_invariance/16_msm_pmatrix.R",                # 연속시간 msm 2년 전이확률 (N_SETS 세트)  [30-120분]
  "01_invariance/17_level_vs_change.R",            # 수준 대 변화 (A 이산; B JM value+slope 는 비식별 보고)
  ## ── C. 그림 ─────────────────────────────────────────────────────────────
  "03_figures/30_figure1_transition_rates_MI.R",   # Figure 1
  "03_figures/31_figure2_recovery_MI.R",           # Figure 2 (+ Supplementary Fig 1 frail median split)
  "03_figures/33_figure3_strata_MI.R",             # Figure 3 (14_ 의 MI_Table3_perset.csv 사용)
  "03_figures/32_figure4_JM_MI.R",                 # Figure 4 JM 동적예측 (MI01)  [가장 오래 걸림]
  "03_figures/43_edfig1_subgroups_MI.R",           # ED Fig 1 하위군 (MI01)
  "04_sensitivity/57_edfig2_wave_stability_MI.R",  # ED Fig 2 파동 안정성 (MI01)
  ## ── D. 민감도 (원 kNN 단일대체 세트; 원고에 명시) ───────────────────────
  "04_sensitivity/50_sens_lead_time.R",            # 역인과 (lead time)
  "04_sensitivity/51_sens_evalue.R",               # E-value
  "04_sensitivity/52_sens_exposure_definition.R",  # 노출 정의
  "04_sensitivity/53_sens_discrimination_ci.R",    # ΔAUC/ΔC 부트스트랩
  "04_sensitivity/54_sens_msm_compare.R",          # 연속시간 msm 대 이산시간
  "04_sensitivity/55_sens_complete_case.R",        # 완전자료
  "04_sensitivity/56_sens_adjustment_reduced.R",   # 축소 보정군
  "04_sensitivity/58_sens_mi_pooling.R",           # 4-세트 풀링
  "04_sensitivity/60_sens_deathtime_ph.R",         # 사망시점 배정 + PH 진단
  "04_sensitivity/61_sens_observed_ic.R",          # 관측 지표 한정
  "04_sensitivity/59_sens_summary.R",              # 민감도 판정 요약 (50-61 뒤에!)
  ## ── E. 보고 ─────────────────────────────────────────────────────────────
  "05_reporting/90_collect_captions.R")            # 각주 모으기 — 반드시 마지막

## ★ 2026-09-08 변경 (투고본 rev9 파이프라인)
##   - 주 분석이 kNN 단일대체·pooled 척도에서 mice PMM m=20·within-sex 척도로 바뀌었습니다.
##     이전 표·그림 스크립트(20-23, 30-33, 43, 57 의 옛 버전)는 legacy_single_imputation/ 로
##     옮겼고 실행 목록에서 뺐습니다. 민감도 50-61 은 결정대로 kNN 세트를 그대로 씁니다.
##   - 12_/13_ 은 산출물이 있으면 자동으로 건너뜁니다 (ONCE).
##   - 모든 스크립트가 같은 헤더(R/_bootstrap.R -> kf_init())를 쓰고, 경로는 R/00_setup.R 에서만 정합니다.

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
for (s in names(ONCE))
  if (!FORCE_REBUILD && s %in% SCRIPTS && file.exists(ONCE[[s]]) && !(s %in% SKIP)) {
    cat("이미 있음 -> 건너뜀: ", s, "  (", basename(ONCE[[s]]), ")\n", sep = "")
    SKIP <- c(SKIP, s)
  }
need_pkg <- c("ggplot2", "patchwork", "scales", "sandwich", "lmtest", "survival",
              "lavaan", "mice", "msm", "JMbayes2", "nlme", "openxlsx", "readxl",
              "timeROC", "officer", "flextable")
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
  } else if (!length(nf) && !(s %in% c("01_invariance/11_invariance_sensitivity.R", "01_invariance/13_mi_pmm_m20.R"))) {
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
