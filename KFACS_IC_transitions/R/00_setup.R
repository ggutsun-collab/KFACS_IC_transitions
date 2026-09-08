###############################################################################
## 00_setup.R — 경로와 실행 옵션의 단일 원천 (이 파일만 편집하면 됩니다)
##
## KFACS intrinsic-capacity longitudinal study — analysis code
## All paths used by every script in this package are defined HERE and only
## here. Edit the three paths below to match your system, then run
## 99_MASTER_run_all.R (or any individual numbered script).
###############################################################################

## ── 1. 사용자 편집 구간 ──────────────────────────────────────────────────
## 경로는 이 세 줄에서만 정합니다. 저장소에는 개인 폴더 경로를 두지 않습니다.
##
## 우선순위:
##   1) 환경변수 KFACS_ROOT 가 있으면 그것을 씁니다
##        Sys.setenv(KFACS_ROOT = "D:/work/Nat_Aging")   # 또는 .Renviron 에 기록
##   2) 없으면 현재 작업디렉터리의 상위에서 data/ 와 result/ 를 찾습니다
##   3) 그래도 없으면 아래 PROJECT_ROOT_FALLBACK 을 직접 고치십시오
PROJECT_ROOT_FALLBACK <- "."          ## 예: "D:/work/Nat_Aging"

.kf_root <- function() {
  e <- Sys.getenv("KFACS_ROOT", unset = "")
  if (nzchar(e) && dir.exists(e)) return(normalizePath(e, winslash = "/"))
  for (cand in c(".", "..", "../..")) {
    if (dir.exists(file.path(cand, "data"))) return(normalizePath(cand, winslash = "/"))
  }
  normalizePath(PROJECT_ROOT_FALLBACK, winslash = "/", mustWork = FALSE)
}
PROJECT_ROOT <- .kf_root()

## 입력 1: 대체 데이터 폴더 — 아래 두 파일이 모두 이 폴더에 있어야 합니다
IMP_DIR  <- Sys.getenv("KFACS_IMP_DIR",
                       unset = file.path(PROJECT_ROOT, "data/260730_imputation_longitudinal"))
## 입력 1a: 단일대체(kNN) 4-세트 — 12_(성별 불변성·기준 점수)와 04_sensitivity/50-61 이 읽습니다
IMP_STEM <- Sys.getenv("KFACS_IMP_STEM",
                       unset = "KFACS_imputed_4sets_LONGITUDINAL_260730")
## 입력 1b: 다중대체(mice PMM, m = 20) — 주 분석(13_ 이후의 모든 표·그림)이 읽습니다
##          13_mi_pmm_m20.R 이 만들며, 파일명은 IMP_DIR/<MI_STEM>.rds 입니다
MI_STEM  <- Sys.getenv("KFACS_MI_STEM", unset = "KFACS_mi_pmm_m20")
## 입력 2: 대체 전 종단 master (00b 가 만듦) — 13_ / 15_ / 18_ / 18b_ / 10_ 이 읽습니다
PRE_FILE <- Sys.getenv("KFACS_PRE_FILE",
                       unset = file.path(PROJECT_ROOT, "data",
                                         "KFACS_master_FINAL_preimput_260618_state.xlsx"))

## 출력: 모든 표·그림·근거수치가 이 폴더 하나로만 저장됩니다
FIG_DIR  <- Sys.getenv("KFACS_OUT_DIR",
                       unset = file.path(PROJECT_ROOT, "result/260812_submission"))

## ── 2. 검증 (여기서부터는 편집하지 않습니다) ─────────────────────────────
## 입력 데이터가 실제로 있는지 - 없으면 시작 전에 즉시 멈춥니다.
## (kNN 세트 또는 MI 세트 중 하나라도 있으면 진행. 어느 쪽이 없는지는 메시지로 알립니다)
.kf_knn_ok <- any(file.exists(file.path(IMP_DIR, paste0(IMP_STEM, c(".rds", ".xlsx", ".csv")))))
.kf_mi_ok  <- file.exists(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
if (!.kf_knn_ok && !.kf_mi_ok)
  stop("대체 데이터를 찾을 수 없습니다.\n",
       "  IMP_DIR  = ", IMP_DIR, "\n",
       "  IMP_STEM = ", IMP_STEM, "   (kNN 4-세트)\n",
       "  MI_STEM  = ", MI_STEM,  "   (mice PMM m=20)\n",
       "  -> 00_setup.R 의 경로를 실제 위치로 수정하거나 00_measurement/ 와 13_ 을 먼저 실행하십시오.")
if (!.kf_knn_ok) message("주의: kNN 세트(", IMP_STEM, ")가 없습니다 — 12_ 과 04_sensitivity/50-61 은 실행되지 않습니다.")
if (!.kf_mi_ok)  message("주의: MI 세트(", MI_STEM, ".rds)가 없습니다 — 13_mi_pmm_m20.R 을 먼저 실행하십시오.")
if (!file.exists(PRE_FILE))
  message("주의: 대체 전 master(", basename(PRE_FILE), ")가 없습니다 — 10_/13_/15_/18_/18b_ 는 실행되지 않습니다.")

## 출력 폴더는 절대 원자료/과거 결과 폴더가 될 수 없습니다.
## (과거에 출력 경로가 옛 결과 폴더를 가리켜 결과가 섞이는 사고가 있었습니다)
if (normalizePath(FIG_DIR, mustWork = FALSE) ==
    normalizePath(IMP_DIR, mustWork = FALSE) ||
    grepl("result/260730", gsub("\\\\", "/", FIG_DIR), fixed = TRUE))
  stop("FIG_DIR 이 보호된 폴더를 가리키고 있습니다: ", FIG_DIR)

if (!dir.exists(FIG_DIR)) {
  ok <- dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!ok && !dir.exists(FIG_DIR)) stop("출력 폴더를 만들 수 없습니다: ", FIG_DIR)
}

KFACS_SETUP_VERSION <- "2026-09-08"
message("00_setup.R loaded  [", KFACS_SETUP_VERSION, "]")
message("  PROJECT_ROOT = ", PROJECT_ROOT)
message("  입력 IMP_DIR = ", IMP_DIR, "   [kNN: ", IMP_STEM, " | MI: ", MI_STEM, "]")
message("  입력 PRE_FILE = ", PRE_FILE)
message("  출력 FIG_DIR = ", FIG_DIR)
