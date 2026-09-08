###############################################################################
## _bootstrap.R — 스크립트 위치 탐색과 헬퍼 로드 (한 곳에서만 정의)
##
## 이전에는 아래 로직이 25개 스크립트에 그대로 복제되어 있었습니다
## (스크립트당 약 35줄, 합계 약 900줄). 한 곳으로 모았습니다.
##
## 각 분석 스크립트는 맨 위에서 이 세 줄만 실행합니다:
##
##   for (.p in c("_bootstrap.R","R/_bootstrap.R","../R/_bootstrap.R","../_bootstrap.R"))
##     if (file.exists(.p)) { source(.p, encoding = "UTF-8"); break }
##   kf_init()
##
## 정의되는 것:
##   KF_HOME  헬퍼가 있는 폴더 (= 저장소의 R/)
##   KF_REPO  저장소 최상위 (= KF_HOME 의 부모)
##   kf_init()  00_setup / 01_functions_common / 02_theme_and_output /
##              03_functions_incremental 을 순서대로 적재. 이미 있으면 건너뜀.
##
## 스크립트가 R/ 이 아닌 하위 폴더(02_tables/ 등)에 있어도 동작합니다.
## 후보 위치를 하나씩 시험하되, 01_functions_common.R 이 실제로 있는 폴더만
## 채택합니다 — 스크립트 자신의 폴더를 헬퍼 폴더로 오인하지 않습니다.
###############################################################################

.kf_ok <- function(d)
  !is.null(d) && nzchar(d) && file.exists(file.path(d, "01_functions_common.R"))

## 후보 폴더에서 시작해 위로 네 단계까지, 각 단계에서 그 폴더와 그 폴더의 R/ 을 확인
.kf_resolve <- function(d) {
  if (is.null(d) || !nzchar(d)) return(NULL)
  cur <- normalizePath(d, winslash = "/", mustWork = FALSE)
  for (i in 0:4) {
    if (.kf_ok(cur)) return(cur)
    r <- file.path(cur, "R")
    if (.kf_ok(r)) return(normalizePath(r, winslash = "/", mustWork = FALSE))
    up <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(up, cur)) break
    cur <- up
  }
  NULL
}

.kf_home <- function() {
  cand <- character(0)
  ## 1) source() 로 실행된 경우 — 호출 스택에 있는 모든 파일 위치
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of))
      cand <- c(cand, dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
  }
  ## 2) Rscript --file=
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a))
    cand <- c(cand, dirname(normalizePath(sub("^--file=", "", a[1]),
                                          winslash = "/", mustWork = FALSE)))
  ## 3) RStudio 편집기의 Source 버튼
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(try(rstudioapi::isAvailable(), silent = TRUE))) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && !is.null(p) && nzchar(p)) cand <- c(cand, dirname(p))
  }
  ## 4) 작업디렉터리
  cand <- c(cand, normalizePath(".", winslash = "/"))

  for (d in unique(cand)) { h <- .kf_resolve(d); if (!is.null(h)) return(h) }
  stop("헬퍼 폴더(R/)를 찾지 못했습니다. 저장소 최상위를 작업디렉터리로 지정하십시오.")
}

KF_HOME <- .kf_home()
KF_REPO <- if (identical(basename(KF_HOME), "R"))
  normalizePath(file.path(KF_HOME, ".."), winslash = "/", mustWork = FALSE) else KF_HOME

.kf_src <- function(f) {
  p <- file.path(KF_HOME, f)
  if (!file.exists(p)) p <- file.path(normalizePath(".", winslash = "/"), f)
  if (!file.exists(p))
    stop("헬퍼 파일을 찾을 수 없습니다: ", f,
         "\n  KF_HOME = ", KF_HOME,
         "\n  -> 저장소 최상위를 작업디렉터리로 지정하십시오.")
  source(p, encoding = "UTF-8")   # local = FALSE -> 전역환경에 적재
  invisible(TRUE)
}

.kf_has <- function(...) all(vapply(c(...), exists, logical(1),
                                    envir = globalenv(), inherits = FALSE))

## 헬퍼 파일마다 '이미 적재되었는가' 를 판정할 표지를 둡니다.
## 표지에 반드시 점(.)으로 시작하는 이름을 하나 이상 넣습니다.
##   이유: save.image() / RStudio 의 .RData 복원은 ls() 기준이라 점으로 시작하는
##   객체를 저장하지 않습니다. 그래서 예전 세션을 복원하면 COMMON_VERSION 같은
##   보통 객체는 살아 있는데 .need() 같은 내부 함수는 사라진 상태가 됩니다.
##   버전 문자열만 보고 판정하면 이 상태를 '적재됨' 으로 오인해
##   ".need 를 찾을 수 없습니다" 로 터집니다.
.KF_HELPERS <- list(
  list(file = "00_setup.R",               marks = c("IMP_DIR", "FIG_DIR")),
  list(file = "01_functions_common.R",    marks = c("COMMON_VERSION", ".need", ".robust_ct")),
  list(file = "02_theme_and_output.R",    marks = c("THEME_VERSION", ".kf_find", ".KF_CAP")),
  list(file = "03_functions_incremental.R", marks = c("ADJ_VERSION")))
  ## 23 / 43 / 50 / 51 / 52 는 ADJ_VERSION 을 요구하면서도 03 을 스스로 불러오지
  ## 않았습니다. 여기서 함께 적재합니다.

kf_init <- function(quiet = FALSE, force = FALSE) {
  if (!quiet) message("KF_HOME = ", KF_HOME)
  for (h in .KF_HELPERS) {
    if (force || !.kf_has(h$marks)) {
      if (!force && .kf_has(h$marks[1]) && !quiet)
        message("  ", h$file, " 재적재 (내부 함수가 없는 상태였습니다 — ",
                ".RData 복원 세션에서 흔합니다)")
      .kf_src(h$file)
    }
  }
  bad <- unlist(lapply(.KF_HELPERS, function(h) h$marks[!vapply(h$marks, exists,
                logical(1), envir = globalenv(), inherits = FALSE)]))
  if (length(bad))
    stop("헬퍼 적재 후에도 없는 객체가 있습니다: ", paste(bad, collapse = ", "),
         "\n  -> kf_init(force = TRUE) 로 다시 시도하십시오.")
  invisible(KF_HOME)
}
