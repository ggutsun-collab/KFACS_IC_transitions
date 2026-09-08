###############################################################################
## 02_theme_and_output.R  (figure standards, palettes, save helpers)
## Nature Aging 그림 규격 공통 테마 / 팔레트 / 저장 헬퍼
##
## 사용: 각 Figure 스크립트 맨 위에서
##   source("02_theme_and_output.R")
##
## ── Nature 계열 규격 요약 ────────────────────────────────────────────────
##  폭   : 단일단 88 mm · 1.5단 120 mm · 양단 180 mm  (초과 금지)
##  글꼴 : sans-serif (Arial/Helvetica). 최종 5–7 pt.
##  패널 : 소문자 볼드 a, b, c ... (좌상단, 축 바깥)
##  선   : 0.25–1 pt. 격자는 최소화, 3D·그림자·그라데이션 금지.
##  색   : 색맹 안전. 색만으로 정보를 전달하지 말 것(모양/라벨 병기).
##  파일 : 벡터 PDF 우선, 래스터는 300–600 dpi TIFF.
###############################################################################

###############################################################################
## 경로는 00_setup.R 한 곳에서만 정의합니다. 이 파일은 그것을 물려받습니다.
###############################################################################
if (!exists("IMP_DIR") || !exists("IMP_STEM") || !exists("FIG_DIR"))
  stop("IMP_DIR / IMP_STEM / FIG_DIR 이 정의되지 않았습니다.\n",
       "  -> 00_setup.R 을 먼저 source 하십시오 (99_MASTER_run_all.R 은 자동으로 합니다).")

## ── 00_common 버전 확인 ──────────────────────────────────────────────────
## 구버전을 섞어 쓰면 default_transitions() 인자 오류, N=3,014 등이 발생합니다.
if (!exists("COMMON_VERSION") || COMMON_VERSION < "2026-07-30") {
  stop("00_common 이 구버전입니다.\n",
       "  현재: ", if (exists("COMMON_VERSION")) COMMON_VERSION else "(버전 표시 없음 = 매우 오래된 버전)",
       "\n  필요: 2026-07-30 이상",
       "\n  -> 최신 00_common_KFACS_IC_260728.R 로 교체한 뒤 다시 source 하십시오.",
       "\n  (증상: default_transitions() 사용되지 않은 인자 (iv) / N 이 3,011 이 아니라 3,014)")
}

THEME_VERSION <- "2026-08-01"

###############################################################################
## ★ 보조 원본 파일 — 대체 데이터셋에 실리지 않은 Table 1 기술변수용
##   (음주 / 사회활동은 분석 공변량이 아니어서 대체 대상에서 빠졌습니다.
##    재대체는 필요 없고, 표를 만들 때 원본에서 붙이면 됩니다.)
##   파일이 아래 이름으로 data 폴더 어디에 있든 자동으로 찾습니다.
###############################################################################
DATA_DIR    <- dirname(IMP_DIR)
SOCIAL_FILE <- ""   # 비워두면 자동 탐색. 직접 지정하려면 전체 경로를 넣으십시오.
PREIMP_FILE <- ""   # 〃

.need <- function(p) {
  m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(m)) install.packages(m, repos = "https://cloud.r-project.org")
}
.need(c("ggplot2", "patchwork", "scales"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

## ── 규격 상수 ────────────────────────────────────────────────────────────
NA_WIDTH  <- list(single = 88, onehalf = 120, double = 180)   # mm
BASE_PT   <- 7      # 본문 글자
LAB_PT    <- 8      # 패널 라벨
LINE_PT   <- 0.35   # 선 굵기 (mm 환산은 ggplot linewidth)

###############################################################################
## ★★ 색 규칙 — 전 그림(Fig 1, Fig 2, Suppl Fig 1, Suppl Fig 2)의 단일 원천 ★★
##
##  원칙: "같은 색 = 같은 뜻". 그림마다 다른 팔레트를 쓰지 않습니다.
##        각 그림 스크립트 안에서 색코드를 직접 쓰는 것을 금지합니다.
##
##   파랑   #0072B2   최상 / 회복  : Normal, Robust, T3 (high IC), Recovery
##   진회색 #595959   중간         : Mild,   Pre-frail, T2 (mid)
##   주황   #E69F00   최악(비사망) : Severe, Frail,     T1 (low IC), Worsening
##   주홍   #D55E00   사망         : Death, Mortality  ← 모든 그림에서 동일
##   연회색 #A6A6A6   참조선/보조
##
##  세 등급이 파랑 - 중립 - 주황 의 대칭 발산형이라 "하나만 튀는" 느낌이 없고,
##  따뜻함=나쁨 / 차가움=좋음 이라는 방향은 그대로 읽힙니다.
##
##  즉 disability(Normal-Mild-Severe), frailty(Robust-Pre-frail-Frail),
##  IC tertile(T3-T2-T1) 세 가지가 모두 "좋음→나쁨" 동일 색 순서를 갖고,
##  사망은 어느 그림에서든 주홍 한 가지입니다.
##  (Okabe–Ito 기반 → 색맹 안전. 색 외에 모양/라벨도 항상 병기합니다.)
###############################################################################
OI <- c(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7",
        grey="#999999")

###############################################################################
## ★ 중간 등급의 색 (2026-07-31 수정)
##
##  이전 판: 중간 = 하늘색(#56B4E9). 그러면 T2·T3 가 둘 다 파랑 계열이라
##  T1(주황) 혼자 튀어 보입니다 -- "순서 있는 척도"가 아니라 "하나만 다른 것"
##  으로 읽힙니다(2026-07-31 지적).
##
##  수정: 중간을 **중립 회색**으로 두어 주황 - 회색 - 파랑 의 대칭적인
##  발산형(diverging) 척도로 만듭니다. 세 색의 색상(hue)이 모두 달라지므로
##  튀는 느낌이 사라지고, "따뜻함=나쁨 / 차가움=좋음" 규칙은 그대로 유지됩니다.
##
##  MID_SCHEME 로 바꿀 수 있습니다 (이 한 줄만 고치면 4개 그림에 모두 반영):
##    "neutral"  주황 - 진회색 - 파랑   (기본, 권장)
##    "cool"     주황 - 하늘색 - 파랑   (이전 판)
###############################################################################
MID_SCHEME <- "neutral"

## 의미축(semantic anchor) — 아래 팔레트는 전부 이 5색에서 파생됩니다
SEV_COL <- c(good  = OI[["blue"]],                                  # 최상 / 회복
             mid   = if (MID_SCHEME == "cool") OI[["skyblue"]] else "#595959",
             bad   = OI[["orange"]],                                # 최악(생존) / 악화
             death = OI[["vermillion"]],                            # 사망 — 항상 고정
             aux   = "#A6A6A6")                                     # 참조선 등(연회색)

## 1) 상태(state) --------------------------------------------------------
PAL_STATE   <- c("Normal"    = SEV_COL[["good"]],  "Mild"  = SEV_COL[["mid"]],
                 "Severe"    = SEV_COL[["bad"]],   "Death" = SEV_COL[["death"]])
PAL_FRAILTY <- c("Robust"    = SEV_COL[["good"]],  "Pre-frail" = SEV_COL[["mid"]],
                 "Frail"     = SEV_COL[["bad"]],   "Death" = SEV_COL[["death"]])

## 2) IC 삼분위 ---------------------------------------------------------
## T3(고IC)=좋음=파랑, T1(저IC)=나쁨=주황. frailty/state 와 방향이 같습니다.
LAB_TERTILE     <- c(T1 = "T1 (low IC)", T2 = "T2 (mid)", T3 = "T3 (high IC)")
PAL_TERTILE     <- c("T1" = SEV_COL[["bad"]], "T2" = SEV_COL[["mid"]],
                     "T3" = SEV_COL[["good"]])
PAL_TERTILE_LAB <- stats::setNames(unname(PAL_TERTILE[c("T1","T2","T3")]),
                                   unname(LAB_TERTILE[c("T1","T2","T3")]))

## 3) 전이 종류 / 결과(outcome) -----------------------------------------
PAL_KIND    <- c("Worsening" = SEV_COL[["bad"]],  "Recovery"  = SEV_COL[["good"]],
                 "Death"     = SEV_COL[["death"]])
PAL_OUTCOME <- c("Worsening" = SEV_COL[["bad"]],  "Incident frailty" = SEV_COL[["bad"]],
                 "Recovery"  = SEV_COL[["good"]], "Mortality" = SEV_COL[["death"]],
                 "Death"     = SEV_COL[["death"]])

## 4) 하위군(성별·연령 등) ----------------------------------------------
## ★ 하위군은 '좋음/나쁨' 축이 아니므로 중증도 팔레트를 쓰면 안 됩니다.
##   (남성=주황 이면 "남성이 나쁨"으로 읽힙니다)
##   중증도 4색과 겹치지 않는 검정/보라 2색을 별도로 둡니다.
SUB_COL <- c(OI[["black"]], OI[["purple"]], OI[["green"]], OI[["yellow"]])
pal_sub <- function(lv) stats::setNames(SUB_COL[seq_along(lv)], lv)

## 어떤 라벨이 와도 규칙에 맞는 색을 돌려주는 조회 함수 (지역 팔레트 금지용)
pal_of <- function(x) {
  tab <- c(PAL_STATE, PAL_FRAILTY, PAL_TERTILE, PAL_TERTILE_LAB, PAL_KIND, PAL_OUTCOME)
  tab <- tab[!duplicated(names(tab))]
  out <- unname(tab[as.character(x)])
  out[is.na(out)] <- SEV_COL[["aux"]]
  stats::setNames(out, as.character(x))
}

## 모든 그림 각주 끝에 붙일 색 설명 (심사자가 그림 간 색을 대조할 수 있게)
## ★ MID_SCHEME 를 바꾸면 이 문장도 자동으로 따라갑니다.
##   (하드코딩해 두면 색만 바뀌고 각주는 옛 색을 말하는 사고가 납니다 — 실제로 발생)
COLOUR_NOTE <- paste0(
  "Colours follow one scheme across all main and supplementary figures: ",
  "blue, best state (normal, robust, highest IC tertile) or recovery; ",
  if (MID_SCHEME == "cool") "light blue, intermediate state; "
  else "neutral grey, intermediate state; ",
  "orange, worst non-fatal state (severe, frail, lowest IC tertile) or worsening; ",
  "vermillion, death. The palette is colour-vision-deficiency safe and each ",
  "category is also identified by its axis label.")

## 색 규칙 확인용 (콘솔에서 pal_card() 실행)
pal_card <- function() {
  z <- c(PAL_STATE, PAL_FRAILTY, PAL_TERTILE_LAB, PAL_KIND)
  z <- z[!duplicated(paste(names(z), z))]
  print(data.frame(label = names(z), colour = unname(z), row.names = NULL))
  invisible(z)
}

## ── 테마 ─────────────────────────────────────────────────────────────────
theme_na <- function(base = BASE_PT) {
  theme_classic(base_size = base, base_family = "") +
    theme(
      text             = element_text(colour = "black"),
      axis.text        = element_text(size = base, colour = "black"),
      axis.title       = element_text(size = base),
      axis.line        = element_line(linewidth = LINE_PT, colour = "black"),
      axis.ticks       = element_line(linewidth = LINE_PT, colour = "black"),
      axis.ticks.length= unit(1.2, "mm"),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(size = base, face = "bold", hjust = 0),
      legend.key.size  = unit(3, "mm"),
      legend.text      = element_text(size = base),
      legend.title     = element_text(size = base),
      legend.position  = "top",
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      plot.title       = element_text(size = base, face = "plain", hjust = 0),
      plot.tag         = element_text(size = LAB_PT, face = "bold"),
      plot.tag.position= c(0, 1),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

## 패널 라벨 (a, b, c ...) — patchwork 조합용
tag_na <- function() plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = LAB_PT, face = "bold"))

###############################################################################
## 그림 제목/각주 서식 통일 — 네 그림이 같은 글자크기·정렬을 갖도록 이것만 씁니다.
##   fig <- fig + annot_na(title = "...", caption = paste0("...", COLOUR_NOTE))
##   색 설명(COLOUR_NOTE)은 자동으로 각주 끝에 붙습니다(중복되면 붙이지 않음).
###############################################################################
## ★ 각주 줄바꿈: ggplot 은 caption 을 자동으로 접지 않아 오른쪽이 잘립니다.
##   그림 폭(mm)에 맞춰 글자수를 계산해 직접 접습니다.
wrap_cap <- function(x, width_mm = NA_WIDTH$double, pt = BASE_PT - 1) {
  if (is.null(x) || !nzchar(x)) return(x)
  ## 대략 1 글자 폭 ~= 0.5 * pt (pt -> mm 환산 0.3528)
  n <- max(60, floor(width_mm / (0.5 * pt * 0.3528)))
  paste(strwrap(x, width = n), collapse = "\n")
}

###############################################################################
## ★ 제목·각주를 그림에 넣지 않습니다 (2026-08-01)
##  학술지 투고 규격에서 figure legend 는 본문 파일에 텍스트로 들어가고,
##  그림 파일에는 패널 문자와 축만 있어야 합니다. 그림 안에 각주를 박으면
##  (1) 편집부가 다시 떼어내야 하고 (2) 각주가 세로 공간을 잡아먹어 패널이
##  눌리며 (3) 축 제목과 겹치는 사고가 납니다. 실제로 Figure 1 에서 났습니다.
##
##  EMBED_CAPTION <- TRUE 로 두면 예전처럼 그림 안에 넣습니다(검토용).
##  기본값 FALSE 에서는 save_na() 가 각주를 <파일명>_caption.txt 로 따로 씁니다.
###############################################################################
EMBED_CAPTION <- FALSE
.KF_CAP <- new.env(parent = emptyenv())
.KF_CAP$pending <- NULL

annot_na <- function(title = NULL, caption = NULL, tag = FALSE, colour_note = TRUE,
                     width_mm = NA_WIDTH$double) {
  if (!is.null(caption) && colour_note &&
      !grepl("Colours follow one scheme", caption, fixed = TRUE))
    caption <- paste0(caption, " ", COLOUR_NOTE)
  ## 원문(줄바꿈 없는 상태)을 따로 보관 — 본문에 붙여넣을 때는 접혀 있으면 안 됩니다.
  .KF_CAP$pending <- list(title = title, caption = caption)
  if (!isTRUE(EMBED_CAPTION))
    return(plot_annotation(
      tag_levels = if (isTRUE(tag)) "a" else NULL,
      theme = theme(plot.tag = element_text(size = LAB_PT, face = "bold"))))
  plot_annotation(
    title = title, caption = wrap_cap(caption, width_mm),
    tag_levels = if (isTRUE(tag)) "a" else NULL,
    theme = theme(
      plot.title   = element_text(size = LAB_PT,   face = "bold",  hjust = 0),
      plot.caption = element_text(size = BASE_PT - 1, face = "plain", hjust = 0,
                                  colour = "grey15", lineheight = 1.15),
      plot.tag     = element_text(size = LAB_PT, face = "bold")))
}

## 모아둔 각주를 한 파일로 합칩니다. RUN_ALL 마지막에 부르면 됩니다.
collect_captions <- function(dir = FIG_DIR, out = "Figure_captions.md") {
  fs <- sort(list.files(dir, pattern = "_caption\\.txt$", full.names = TRUE))
  if (!length(fs)) { message("각주 파일이 없습니다: ", dir); return(invisible(NULL)) }
  ord <- c("Figure1", "Figure2", "Figure2ALT", "SupplFig1", "SupplFig2",
           "SupplFig3", "SupplFig4")
  key <- basename(fs)
  rk  <- vapply(key, function(k) {
    i <- which(vapply(ord, function(o) startsWith(k, o), logical(1)))
    if (length(i)) min(i) else length(ord) + 1L }, integer(1))
  fs <- fs[order(rk, key)]
  txt <- unlist(lapply(fs, function(f)
    c(readLines(f, encoding = "UTF-8", warn = FALSE), "")), use.names = FALSE)
  fp <- file.path(dir, out)
  writeLines(txt, fp, useBytes = TRUE)
  message("captions: ", fp, "  (", length(fs), "개)")
  invisible(fp)
}

## ── 저장 ─────────────────────────────────────────────────────────────────
## PDF(벡터) + TIFF(600 dpi, LZW) 동시 저장. 폭은 mm 로 지정.
save_na <- function(p, file, width_mm = NA_WIDTH$double, height_mm = 120, dir = ".") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  stem <- file.path(dir, sub("\\.(pdf|tif|tiff|png)$", "", file))
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::cairo_pdf)
  ok <- try(ggplot2::ggsave(paste0(stem, ".tiff"), p, width = width_mm, height = height_mm,
                            units = "mm", dpi = 600, compression = "lzw"), silent = TRUE)
  if (inherits(ok, "try-error"))
    ggplot2::ggsave(paste0(stem, ".png"), p, width = width_mm, height = height_mm,
                    units = "mm", dpi = 600)
  message("saved: ", stem, ".pdf / .tiff")
  ## ── 각주를 별도 텍스트로 (그림 파일에는 넣지 않습니다) ──────────────
  cp <- .KF_CAP$pending
  if (!is.null(cp)) {
    ttl <- basename(sub("\\.(pdf|tif|tiff|png)$", "", file))
    ## Figure1_TransitionRates -> "Figure 1" 처럼 사람이 읽을 라벨
    lab <- sub("^SupplFig(\\d+).*$", "Supplementary Figure \\1", ttl)
    lab <- sub("^Figure(\\d+)ALT.*$", "Figure \\1 (alternative)", lab)
    lab <- sub("^Figure(\\d+)[^0-9].*$", "Figure \\1", lab)
    body <- c(paste0("## ", lab,
                     if (!is.null(cp$title) && nzchar(cp$title)) paste0(" | ", cp$title) else ""),
              "", paste0("**File:** ", basename(stem), ".pdf / .tiff"), "",
              if (!is.null(cp$caption)) cp$caption else "")
    writeLines(body, paste0(stem, "_caption.txt"), useBytes = TRUE)
    message("caption: ", stem, "_caption.txt")
    .KF_CAP$pending <- NULL
  }
  invisible(stem)
}

## 값 표 저장 (그림의 근거 숫자를 항상 함께 남긴다)
save_vals <- function(x, file, dir = ".") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fp <- file.path(dir, file)
  utils::write.csv(x, fp, row.names = FALSE)
  message("values: ", fp)
  invisible(fp)
}

## log 축 눈금 (HR 그림용)
## 그림용 라벨: 파일에는 ASCII "->" 로 저장하고, 그림에서만 화살표로 표시
arrow_lab <- function(x) gsub("->", "\u2192", as.character(x), fixed = TRUE)

hr_axis <- function(breaks = c(0.25, 0.5, 1, 2, 4)) {
  list(scale_x_log10(breaks = breaks,
                     labels = format(breaks, drop0trailing = TRUE)),
       geom_vline(xintercept = 1, linetype = "22", linewidth = LINE_PT, colour = "grey40"))
}

###############################################################################
## 표 저장 — Nature 계열 서식
##   제목 1행 + 헤더 + 본문 + 빈 줄 + 각주.  csv 와 xlsx 를 함께 남깁니다.
##   (Nature 는 표를 이미지가 아니라 편집가능 파일로 요구합니다)
###############################################################################
save_table <- function(x, file, title = NULL, footnotes = NULL, dir = FIG_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  stem <- file.path(dir, sub("\\.(csv|xlsx)$", "", file))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  ## 표 파일은 ASCII 로 통일합니다. en-dash/dagger 등은 로케일에 따라 깨지고
  ## 필드가 잘리는 사고가 납니다(조판에서 어차피 다시 처리됨).
  .ascii <- function(z) {
    z <- as.character(z); z[is.na(z)] <- ""
    z <- gsub("\u2013|\u2014|\u2212", "-", z)
    z <- gsub("\u2020", "*", z); z <- gsub("\u2021", "**", z)
    z <- gsub("\u2265", ">=", z); z <- gsub("\u2264", "<=", z)
    z <- gsub("\u2003|\u00a0", " ", z)
    z <- gsub("\u2192", "->", z)
    out <- iconv(z, to = "ASCII//TRANSLIT", sub = "?")
    ifelse(is.na(out), gsub("[^\x20-\x7E]", "?", z), out)
  }
  x[]  <- lapply(x, .ascii)
  names(x) <- .ascii(names(x))
  if (!is.null(title))     title     <- .ascii(title)
  if (!is.null(footnotes)) footnotes <- .ascii(footnotes)

  pad <- function(v) c(v, rep("", ncol(x) - 1))
  blk <- rbind(if (!is.null(title)) pad(title),
               names(x), as.matrix(x), pad(""),
               if (!is.null(footnotes)) t(vapply(footnotes, pad, character(ncol(x)))))
  ## ★ fileEncoding 을 쓰면 native -> UTF-8 재인코딩이 일어나 en-dash(\u2013),
  ##   dagger(\u2020) 등이 깨지고 문자열이 잘립니다. 바이트로 직접 씁니다.
  esc <- function(v) paste0("\"", gsub("\"", "\"\"", v), "\"")
  lines <- apply(blk, 1, function(r) paste(esc(as.character(r)), collapse = ","))
  con <- file(paste0(stem, ".csv"), open = "wb")
  writeBin(charToRaw("\ufeff"), con)                     # Excel 용 UTF-8 BOM
  writeLines(enc2utf8(lines), con, useBytes = TRUE)
  close(con)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook(); sh <- substr(basename(stem), 1, 31)
    openxlsx::addWorksheet(wb, sh)
    r <- 1
    if (!is.null(title)) {
      openxlsx::writeData(wb, sh, title, startRow = r, startCol = 1)
      openxlsx::addStyle(wb, sh, openxlsx::createStyle(textDecoration = "bold"), rows = r, cols = 1)
      r <- r + 2
    }
    openxlsx::writeData(wb, sh, x, startRow = r, headerStyle =
      openxlsx::createStyle(textDecoration = "bold", border = "TopBottom"))
    openxlsx::addStyle(wb, sh, openxlsx::createStyle(border = "bottom"),
                       rows = r + nrow(x), cols = seq_len(ncol(x)), gridExpand = TRUE, stack = TRUE)
    if (!is.null(footnotes)) {
      fr <- r + nrow(x) + 2
      for (k in seq_along(footnotes))
        openxlsx::writeData(wb, sh, footnotes[k], startRow = fr + k - 1, startCol = 1)
      openxlsx::addStyle(wb, sh, openxlsx::createStyle(fontSize = 9),
                         rows = fr:(fr + length(footnotes) - 1), cols = 1, gridExpand = TRUE)
    }
    openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(x)), widths = "auto")
    openxlsx::saveWorkbook(wb, paste0(stem, ".xlsx"), overwrite = TRUE)
  }
  message("table: ", stem, ".csv / .xlsx")
  invisible(stem)
}

###############################################################################
## 대체 데이터 로더
##   rds 우선(빠르고 자료형 보존), 없으면 xlsx.
##   한글 경로에서 readxl 이 실패하는 경우가 있어 ASCII 임시경로로 복사해 읽습니다.
##   반환: data.frame. 부수효과로 CFG$DATA_FILE 을 임시 csv 로 설정합니다.
###############################################################################
load_imputed <- function(set = "MAIN", dir = IMP_DIR, stem = IMP_STEM) {
  rds <- file.path(dir, paste0(stem, ".rds"))
  xls <- file.path(dir, paste0(stem, ".xlsx"))
  csv <- file.path(dir, paste0(stem, ".csv"))
  d <- NULL

  if (file.exists(rds)) {
    x <- readRDS(rds)
    d <- if (is.data.frame(x)) x else x[[set]]
    message("[load_imputed] rds: ", basename(rds), " [", set, "]")
  } else if (file.exists(xls)) {
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    tmp <- file.path(tempdir(), "imp_in.xlsx")          # 한글경로 회피
    file.copy(xls, tmp, overwrite = TRUE)
    sh <- readxl::excel_sheets(tmp)
    if (!set %in% sh) { warning("시트 '", set, "' 없음 -> 첫 시트 사용: ", sh[1]); set <- sh[1] }
    d <- as.data.frame(readxl::read_excel(tmp, sheet = set))
    message("[load_imputed] xlsx: ", basename(xls), " [", set, "]")
  } else if (file.exists(csv)) {
    d <- utils::read.csv(csv, stringsAsFactors = FALSE)
    message("[load_imputed] csv: ", basename(csv))
  } else {
    stop("대체 데이터를 찾을 수 없습니다.\n  찾은 위치: ", dir,
         "\n  기대 파일: ", paste0(stem, c(".rds", ".xlsx", ".csv"), collapse = " / "),
         "\n  -> 00_setup.R 의 IMP_DIR / IMP_STEM 을 확인하십시오.",
         "\n  현재 폴더 내용: ",
         if (dir.exists(dir)) paste(list.files(dir), collapse = ", ") else "(폴더 자체가 없음)")
  }

  ## 00_common 의 load_kfacs() 가 읽을 수 있도록 ASCII 임시 csv 로 넘김
  tmpc <- file.path(tempdir(), "imp_main.csv")
  utils::write.csv(d, tmpc, row.names = FALSE)
  if (exists("CFG", envir = globalenv())) {
    CFG <- get("CFG", envir = globalenv()); CFG$DATA_FILE <- tmpc
    assign("CFG", CFG, envir = globalenv())
  }
  message(sprintf("[load_imputed] %d행 x %d열", nrow(d), ncol(d)))
  invisible(d)
}

###############################################################################
## 보조 원본에서 Table 1 기술변수 붙이기
##   attach_extras(b)  ->  b 에 social_activity / alcohol_heavy 열을 추가
##
## 왜 이렇게 하나:
##   음주(alcohol_3cat)와 사회활동(social_activity)은 종단 재대체의 대상이
##   아니었습니다(분석 공변량 = age, sex, comorbidity). 따라서 다시 대체를
##   돌릴 필요가 없고, 표를 만들 때 wave-1 관측값을 그대로 붙이는 것이 맞습니다.
##   음주는 wave 1 에서 결측이 큽니다(약 28%) -> 관측된 사람만 분모로 씁니다.
###############################################################################
.kf_find <- function(patterns, hints = c(DATA_DIR, IMP_DIR, dirname(DATA_DIR), FIG_DIR, getwd())) {
  hints <- unique(hints[nzchar(hints)])
  for (h in hints) {
    if (!dir.exists(h)) next
    ff <- try(list.files(h, recursive = TRUE, full.names = TRUE, include.dirs = FALSE),
              silent = TRUE)
    if (inherits(ff, "try-error") || !length(ff)) next
    bn <- basename(ff)
    for (pt in patterns) {
      k <- which(grepl(pt, bn, ignore.case = TRUE))
      if (length(k)) return(ff[k[1]])
    }
  }
  NA_character_
}

.kf_read_any <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
    tmp <- file.path(tempdir(), paste0("kfx_", ext, ".", ext))   # 한글경로 회피
    file.copy(path, tmp, overwrite = TRUE)
    x <- try(as.data.frame(readxl::read_excel(tmp)), silent = TRUE)
  } else if (ext == "rds") {
    x <- try(readRDS(path), silent = TRUE)
    if (!inherits(x, "try-error") && !is.data.frame(x)) x <- x[[1]]
  } else {
    x <- try(utils::read.csv(path, stringsAsFactors = FALSE), silent = TRUE)
  }
  if (inherits(x, "try-error")) NULL else x
}

attach_extras <- function(b, verbose = TRUE) {
  attr(b, "extras_found")   <- character(0)
  attr(b, "extras_missing") <- character(0)
  found <- missing_ <- character(0)

  ## ── 1) 사회활동 ────────────────────────────────────────────────────────
  if (!"social_activity" %in% names(b)) {
    fp <- if (nzchar(SOCIAL_FILE) && file.exists(SOCIAL_FILE)) SOCIAL_FILE else
          .kf_find(c("^KFACS_social_activity.*\\.(xlsx|xls|csv)$",
                     "social.*activ.*\\.(xlsx|xls|csv)$"))
    s <- .kf_read_any(fp)
    if (!is.null(s) && all(c("id", "social_activity") %in% names(s))) {
      s <- s[!duplicated(s$id), c("id", "social_activity")]
      b <- merge(b, s, by = "id", all.x = TRUE)
      found <- c(found, sprintf("social_activity (%s, 매칭 %d/%d)",
                                basename(fp), sum(b$id %in% s$id), nrow(b)))
    } else {
      missing_ <- c(missing_, "social_activity")
      if (verbose) message("[attach_extras] 사회활동 파일을 찾지 못했습니다.\n",
        "  기대 이름: KFACS_social_activity_*.xlsx (id, social_activity 두 열)\n",
        "  찾은 곳  : ", DATA_DIR, " 이하\n",
        "  -> KF_theme 의 SOCIAL_FILE 에 전체 경로를 직접 넣으면 됩니다.")
    }
  }

  ## ── 2) 음주 (pre-imputation 원본의 wave 1) ─────────────────────────────
  if (!any(c("alcohol_heavy", "alcohol_3cat") %in% names(b))) {
    fp <- if (nzchar(PREIMP_FILE) && file.exists(PREIMP_FILE)) PREIMP_FILE else
          .kf_find(c("^KFACS_master.*preimput.*\\.(xlsx|xls|csv|rds)$",
                     "preimput.*\\.(xlsx|xls|csv|rds)$"))
    p <- .kf_read_any(fp)
    if (!is.null(p) && all(c("id", "wave") %in% names(p)) &&
        any(c("alcohol_heavy", "alcohol_3cat") %in% names(p))) {
      add <- intersect(c("alcohol_heavy", "alcohol_3cat"), names(p))
      w1  <- p[as_num(p$wave) == 1, c("id", add), drop = FALSE]
      w1  <- w1[!duplicated(w1$id), , drop = FALSE]
      b   <- merge(b, w1, by = "id", all.x = TRUE)
      ## alcohol_heavy 가 없으면 3분류에서 만듭니다 (3 = 주 2-3회 이상)
      if (!"alcohol_heavy" %in% names(b) && "alcohol_3cat" %in% names(b))
        b$alcohol_heavy <- as.integer(as_num(b$alcohol_3cat) == 3)
      found <- c(found, sprintf("alcohol_heavy (%s, 관측 %d/%d)", basename(fp),
                                sum(is.finite(as_num(b$alcohol_heavy))), nrow(b)))
    } else {
      missing_ <- c(missing_, "alcohol_heavy")
      if (verbose) message("[attach_extras] pre-imputation 원본을 찾지 못했습니다.\n",
        "  기대 이름: KFACS_master_FINAL_preimput_*.xlsx (id, wave, alcohol_3cat 포함)\n",
        "  찾은 곳  : ", DATA_DIR, " 이하\n",
        "  -> KF_theme 의 PREIMP_FILE 에 전체 경로를 직접 넣으면 됩니다.")
    }
  }

  if (verbose && length(found)) message("[attach_extras] 붙임: ", paste(found, collapse = " | "))
  attr(b, "extras_found")   <- found
  attr(b, "extras_missing") <- missing_
  b
}

if (!dir.exists(FIG_DIR)) {
  ok <- dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!ok && !dir.exists(FIG_DIR))
    warning("출력 폴더를 만들 수 없습니다: ", FIG_DIR, "\n  경로/권한을 확인하십시오.")
}
message("02_theme_and_output loaded  [THEME_VERSION = ", THEME_VERSION, "]")
message("  색 규칙: 파랑=최상/회복 · 하늘=중간 · 주황=최악(생존)/악화 · 주홍=사망")
message("           (확인: pal_card() )")
message("  입력 IMP_DIR = ", IMP_DIR, if (dir.exists(IMP_DIR)) "  [OK]" else "  [폴더 없음]")
message("  출력 FIG_DIR = ", FIG_DIR, if (dir.exists(FIG_DIR)) "  [OK]" else "  [생성 실패]")
