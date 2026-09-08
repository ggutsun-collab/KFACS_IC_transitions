###############################################################################
## 91_sf_cleanup.R — Supplementary 산출물 정리 (1회 실행, 재실행 안전)
##
##  · 폐지된 보충그림/보충표 산출물을 _retired_260813 폴더로 이동 (삭제 아님)
##  · SupplFig2_WaveStability 로 번호 정리 (구명 파일이 있으면 개명)
##  · 남는 정식 보충그림: SupplFig1_Subgroups · SupplFig2_WaveStability
###############################################################################

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
RET <- file.path(FIG_DIR, "_retired_260813")
dir.create(RET, showWarnings = FALSE)

## ── 1. 폐지 산출물 이동 ─────────────────────────────────────────────────
retire_pat <- c(
  "^SupplFig1_TransitionIRR_level",   # 구 SF1 (Table 2 중복)
  "^SupplFig2_Frailty_x_ICtertile",   # 구 SF2 (Table 3 중복)
  "^SupplFig3_ClinicalUtility",       # 구 SF3 (삭제 결정)
  "^SupplFig4_Subgroups",             # 구 SF4 (신 SF1 로 대체)
  "^SupplFig_Subgroups",              # 과도기 이름
  "^SF1_transition_IRR_level",        # 구 SF1 근거수치
  "^SF2_risk5", "^SF2_HR",            # 구 SF2 근거수치
  "^SF3_",                            # 구 SF3 근거수치 (판별력 boot 는 아래서 개명)
  "^SF4_subgroup_values",             # 구 SF4 근거수치 (신판은 SF1_subgroup_values)
  "^SupplTable2_TransitionsAdj",      # 폐지 ST (23 과 중복)
  "^Table4_Attenuation", "^Table4_values")  # 구명 (신명 ST_FrailtyAdjustment)
moved <- 0L
for (pt in retire_pat) {
  f <- list.files(FIG_DIR, pattern = pt, full.names = TRUE)
  f <- f[!file.info(f)$isdir]
  for (x in f) {
    if (file.rename(x, file.path(RET, basename(x)))) {
      cat("  이동:", basename(x), "\n"); moved <- moved + 1L
    } else cat("  !! 이동 실패(파일 열림?):", basename(x), "\n")
  }
}
cat(sprintf("  폐지 산출물 %d개 -> %s\n", moved, RET))

## ── 2. 번호 정리 ────────────────────────────────────────────────────────
ren <- rbind(
  c("SupplFig_WaveStability.pdf",          "SupplFig2_WaveStability.pdf"),
  c("SupplFig_WaveStability.tiff",         "SupplFig2_WaveStability.tiff"),
  c("SupplFig_WaveStability.png",          "SupplFig2_WaveStability.png"),
  c("SupplFig_WaveStability_caption.txt",  "SupplFig2_WaveStability_caption.txt"),
  c("SF3_discrimination_boot.csv",         "ST_Discrimination_boot.csv"))
for (i in seq_len(nrow(ren))) {
  a <- file.path(FIG_DIR, ren[i, 1]); b <- file.path(FIG_DIR, ren[i, 2])
  if (file.exists(a) && !file.exists(b)) {
    file.rename(a, b); cat("  개명:", ren[i, 1], "->", ren[i, 2], "\n")
  }
}

## ── 3. 현황 ─────────────────────────────────────────────────────────────
cat("\n  정식 보충그림 현황:\n")
for (f in c("SupplFig1_Subgroups.pdf", "SupplFig2_WaveStability.pdf"))
  cat(sprintf("   %-32s %s\n", f,
              if (file.exists(file.path(FIG_DIR, f))) "[있음]"
              else "[없음 - 해당 스크립트를 실행하십시오 (SF1=43, SF2=57)]"))
cat("\n=== 정리 완료 ===\n")
