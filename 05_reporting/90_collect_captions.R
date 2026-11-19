###############################################################################
## KF_CollectCaptions_v260801.R
## 그림 각주 모으기 — 각 그림이 만든 *_caption.txt 를 Figure_captions.md 로 합칩니다.
##
##   그림 파일(pdf/tiff)에는 제목도 각주도 들어가지 않습니다(투고 규격).
##   각주는 그림을 만들 때마다 <그림파일명>_caption.txt 로 따로 저장되고,
##   이 스크립트가 그것들을 Figure 1 -> Figure 2 -> Supplementary 순으로
##   하나의 문서로 묶습니다. 본문 파일에 그대로 붙여넣으시면 됩니다.
##
##   그림을 전부 다시 돌린 뒤 이 파일 하나만 Source 하십시오. 1초.
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
if (!exists("collect_captions"))
  stop("02_theme_and_output.R 이 구버전입니다(THEME_VERSION 2026-08-01 필요).")

fp <- collect_captions(FIG_DIR)

cat("\n=== 그림 각주 모음 완료 ===\n")
if (!is.null(fp)) {
  cat("  파일: ", fp, "\n", sep = "")
  cat("\n---- 미리보기 ----\n")
  writeLines(utils::head(readLines(fp, encoding = "UTF-8", warn = FALSE), 30))
} else {
  cat("  *_caption.txt 가 하나도 없습니다.\n")
  cat("  -> 그림 스크립트를 최신본(THEME_VERSION 2026-08-01)으로 다시 돌리십시오.\n")
}
