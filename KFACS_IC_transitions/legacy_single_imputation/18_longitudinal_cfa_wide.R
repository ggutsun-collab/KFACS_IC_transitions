###############################################################################
## 18_longitudinal_cfa_wide.R   (v260907)
## 종단 측정불변성 — wide-format longitudinal CFA (within-person dependence 반영)
##
## 왜 필요한가
##   Supplementary Table 12 는 wave 를 group 으로 두는 multi-group CFA 였습니다.
##   같은 사람이 다섯 group 에 모두 들어가므로 group 간 독립 가정이 깨집니다.
##   표준 대안은 다섯 wave 의 지표를 한 사람의 85개 변수로 펼쳐(wide) 한 모형에
##   넣고, (i) 같은 지표의 잔차를 wave 간에 상관시키고, (ii) 요인을 wave 간에
##   상관시키는 것입니다. 적재값·절편 동일성 제약은 그 안에서 겁니다.
##
## 하는 일
##   pre-imp 관측 자료(사망 후 행 제외) -> 17개 지표를 person-wave 풀 분포로 표준화
##   -> wide -> configural / metric / scalar 세 모형 (MLR, FIML)
##   -> Supplementary Table 12 형식의 적합도 표 + 스칼라 모형의 g 잠재평균(wave 별)
##
## 출력 (FIG_DIR): ST_LongInvariance_wide.csv, ST_LongInvariance_wide_latentmeans.csv,
##                 LongCFA_wide_fits.rds
## 실행 : source("18_longitudinal_cfa_wide.R", encoding = "UTF-8")
##        85개 관측변수 × n=3,011, FIML. 모형당 수 분 ~ 수십 분. 순서대로 3개 적합.
###############################################################################

LCFA_VERSION <- "v260907"
message("\n=== 18_longitudinal_cfa_wide ", LCFA_VERSION, " ===")

KF_HOME_FALLBACK <- ""
.kf_home <- function() {
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of)) {
      dd <- dirname(normalizePath(of, winslash = "/", mustWork = FALSE))
      if (file.exists(file.path(dd, "00_setup.R"))) return(dd)
    }
  }
  if (nzchar(KF_HOME_FALLBACK) && file.exists(file.path(KF_HOME_FALLBACK, "00_setup.R"))) return(KF_HOME_FALLBACK)
  if (file.exists("00_setup.R")) return(normalizePath(".", winslash = "/"))
  roots <- unique(c("C:/Users/user/OneDrive/바탕 화면/Nat_Aging", getwd(), dirname(getwd())))
  for (r in roots) {
    if (!dir.exists(r)) next
    hit <- list.files(r, pattern = "^00_setup\\.R$", recursive = TRUE, full.names = TRUE)
    hit <- hit[!grepl("backup|_retired|old", hit, ignore.case = TRUE)]
    if (length(hit)) return(dirname(hit[order(file.info(hit)$mtime, decreasing = TRUE)][1]))
  }
  stop("00_setup.R 이 있는 폴더를 찾지 못했습니다.")
}
KF_HOME <- .kf_home()
.kf_src <- function(f) { p <- file.path(KF_HOME, f)
  if (!file.exists(p)) stop("헬퍼 파일을 찾을 수 없습니다: ", p); source(p, encoding = "UTF-8") }
.kf_src("00_setup.R"); .kf_src("01_functions_common.R"); .kf_src("02_theme_and_output.R")
.need(c("lavaan", "readxl"))
suppressPackageStartupMessages(library(lavaan))

###############################################################################
## 1. 설정
###############################################################################
PRE_FILE <- file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260618_state.xlsx")
WAVES    <- 1:5
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
SPEC <- list(loco = c("Balance","GS_ms","rev_CST"),
             vita = c("loss_of_Bwt","Appetite","EXH","HGS"),
             cogn = c("Orientation","Memory","Attention_Calculation","Language","Visuospatial"),
             psyc = c("Negative_affect","Positive_affect","Motivation"))
## 잔차 상관: 같은 지표를 모든 wave 쌍에서 상관 (lag 무관). 수렴이 어려우면 ADJACENT_ONLY <- TRUE
ADJACENT_ONLY <- FALSE
EXCL_IDS <- c("kf161224", "kf170602", "kf171295")

###############################################################################
## 2. 자료: 관측 person-wave, 사망 후 행 제외, 풀 표준화, wide
###############################################################################
tmp <- file.path(tempdir(), "pre.xlsx"); file.copy(PRE_FILE, tmp, overwrite = TRUE)
pre <- as.data.frame(readxl::read_excel(tmp, sheet = "전체_long"))
pre$id <- as.character(pre$id); pre$wave <- as.integer(as_num(pre$wave))
pre <- pre[!pre$id %in% EXCL_IDS & pre$wave %in% WAVES, ]
dw <- tapply(pre$death_wave, pre$id, function(x) suppressWarnings(as.numeric(x[1])))
pre <- pre[is.na(dw[pre$id]) | pre$wave < dw[pre$id], ]
for (v in IND17) pre[[v]] <- as_num(pre[[v]])
pre <- pre[rowSums(!is.na(pre[, IND17])) > 0, ]
MU <- vapply(pre[IND17], mean, numeric(1), na.rm = TRUE); SD <- vapply(pre[IND17], sd, numeric(1), na.rm = TRUE)
for (v in IND17) pre[[v]] <- (pre[[v]] - MU[[v]]) / SD[[v]]       # 12_ / 13_ 와 같은 풀 표준화
msg(sprintf("관측 person-wave %d행, 참가자 %d명", nrow(pre), length(unique(pre$id))))

ids <- sort(unique(pre$id))
Wd <- data.frame(id = ids, stringsAsFactors = FALSE)
for (w in WAVES) {
  z <- pre[pre$wave == w, c("id", IND17)]; z <- z[!duplicated(z$id), ]
  names(z)[-1] <- paste0(IND17, "_w", w)
  Wd <- merge(Wd, z, by = "id", all.x = TRUE)
}
cat(sprintf("  wide: %d명 × %d변수 | wave별 관측: %s\n", nrow(Wd), ncol(Wd) - 1,
            paste(vapply(WAVES, function(w) sum(!is.na(Wd[[paste0("Balance_w", w)]])), integer(1)), collapse = "/")))

###############################################################################
## 3. lavaan 모형 생성기
##    level = "configural" | "metric" | "scalar"
##    - g_w, loco_w, vita_w, cogn_w, psyc_w ; 감각 2지표는 g 에만
##    - wave 내: bifactor 직교 (g ~~ 0*spec, spec ~~ 0*spec)
##    - wave 간: g_w ~~ g_v 자유, 같은 specific 끼리 자유, 그 외 0
##    - 같은 지표 잔차 wave 간 상관 자유
##    - 식별: wave 1 요인분산 1, 잠재평균 0 ; metric 이상에서 wave 2-5 분산 자유,
##            scalar 에서 wave 2-5 g 잠재평균 자유 (specific 평균은 0 고정)
###############################################################################
build_model <- function(level = c("configural", "metric", "scalar")) {
  level <- match.arg(level); L <- character(0)
  lab <- function(prefix, item, w) if (level == "configural") "" else paste0(prefix, "_", item, "*")
  for (w in WAVES) {
    iw <- function(x) paste0(x, "_w", w)
    L <- c(L, sprintf("g_w%d =~ %s", w, paste0(lab("lg", IND17, w), iw(IND17), collapse = " + ")))
    for (s in names(SPEC))
      L <- c(L, sprintf("%s_w%d =~ %s", s, w, paste0(lab(paste0("l", s), SPEC[[s]], w), iw(SPEC[[s]]), collapse = " + ")))
    ## 절편
    ilab <- if (level == "scalar") paste0("nu_", IND17, "*") else ""
    L <- c(L, sprintf("%s ~ %s1", iw(IND17), ilab))
    ## 요인분산 / 평균
    facs <- c(paste0("g_w", w), paste0(names(SPEC), "_w", w))
    if (w == 1 || level == "configural") {
      L <- c(L, sprintf("%s ~~ 1*%s", facs, facs))
    } else {
      L <- c(L, sprintf("%s ~~ NA*%s", facs, facs))
    }
    L <- c(L, sprintf("%s ~ %s1", facs, ifelse(facs == paste0("g_w", w) & w > 1 & level == "scalar", "NA*", "0*")))
    ## wave 내 직교
    for (a in seq_along(facs)) for (b in seq_along(facs)) if (a < b)
      L <- c(L, sprintf("%s ~~ 0*%s", facs[a], facs[b]))
  }
  ## wave 간 요인 공분산
  for (a in WAVES) for (b in WAVES) if (a < b) {
    L <- c(L, sprintf("g_w%d ~~ g_w%d", a, b))
    for (s in names(SPEC)) L <- c(L, sprintf("%s_w%d ~~ %s_w%d", s, a, s, b))
    for (s in names(SPEC)) { L <- c(L, sprintf("g_w%d ~~ 0*%s_w%d", a, s, b), sprintf("%s_w%d ~~ 0*g_w%d", s, a, b)) }
    for (s in names(SPEC)) for (t in names(SPEC)) if (s != t) L <- c(L, sprintf("%s_w%d ~~ 0*%s_w%d", s, a, t, b))
  }
  ## 같은 지표 잔차 상관
  for (v in IND17) for (a in WAVES) for (b in WAVES) if (a < b && (!ADJACENT_ONLY || b == a + 1))
    L <- c(L, sprintf("%s_w%d ~~ %s_w%d", v, a, v, b))
  paste(L, collapse = "\n")
}

###############################################################################
## 4. 적합
###############################################################################
fit_one <- function(level) {
  msg(sprintf("[%s] 적합 시작 ...", level)); t0 <- Sys.time()
  f <- try(cfa(build_model(level), data = Wd, estimator = "MLR", missing = "fiml",
               meanstructure = TRUE, std.lv = FALSE, auto.fix.first = FALSE,
               control = list(iter.max = 2000)), silent = TRUE)
  if (inherits(f, "try-error")) { message("  !! ", level, " 실패: ", conditionMessage(attr(f, "condition"))); return(NULL) }
  cat(sprintf("  [%s] converged=%s, %.1f분\n", level, lavInspect(f, "converged"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  f
}
FITS <- list(configural = fit_one("configural"), metric = fit_one("metric"), scalar = fit_one("scalar"))
saveRDS(FITS, file.path(FIG_DIR, "LongCFA_wide_fits.rds"))

###############################################################################
## 5. 적합도 표 (Supplementary Table 12 와 같은 열)
###############################################################################
fm <- function(f) if (is.null(f)) rep(NA_real_, 6) else {
  m <- fitMeasures(f, c("chisq.scaled", "df.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
  if (any(!is.finite(m[3:5]))) m[3:5] <- fitMeasures(f, c("cfi.scaled", "tli.scaled", "rmsea.scaled"))
  as.numeric(m) }
TAB <- data.frame(Model = c("Configural", "Metric (loadings)", "Scalar (+ intercepts)"),
                  do.call(rbind, lapply(FITS, fm)), stringsAsFactors = FALSE)
names(TAB)[-1] <- c("Chi-square", "df", "CFI", "TLI", "RMSEA", "SRMR")
TAB$`Delta CFI`   <- c(NA, diff(TAB$CFI)); TAB$`Delta RMSEA` <- c(NA, diff(TAB$RMSEA)); TAB$`Delta SRMR` <- c(NA, diff(TAB$SRMR))
crit <- function(dc, dr, ds, step) if (!is.finite(dc)) "-" else
  if (dc >= -0.010 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "supported" else
  if (dc >= -0.015 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "borderline (CFI criterion only)" else "not supported"
TAB$Verdict <- c("-", crit(TAB$`Delta CFI`[2], TAB$`Delta RMSEA`[2], TAB$`Delta SRMR`[2], "metric"),
                 crit(TAB$`Delta CFI`[3], TAB$`Delta RMSEA`[3], TAB$`Delta SRMR`[3], "scalar"))
## 척도화 카이제곱 차이검정 (참고용; n 이 크면 거의 항상 유의)
lrt <- try(lavTestLRT(FITS$configural, FITS$metric, FITS$scalar, method = "satorra.bentler.2001"), silent = TRUE)
if (!inherits(lrt, "try-error")) { TAB$`Scaled chi-square difference P` <- c(NA, lrt[["Pr(>Chisq)"]][2:3]) }
save_vals(TAB, "ST_LongInvariance_wide.csv", FIG_DIR)
print(TAB, row.names = FALSE)

## 스칼라 모형의 wave 별 g 잠재평균 (wave 1 = 0 기준) — 코호트 평균 IC 의 변화
if (!is.null(FITS$scalar)) {
  pe <- parameterEstimates(FITS$scalar)
  gm <- pe[pe$op == "~1" & grepl("^g_w", pe$lhs), c("lhs", "est", "se", "ci.lower", "ci.upper")]
  gv <- pe[pe$op == "~~" & grepl("^g_w", pe$lhs) & pe$lhs == pe$rhs, c("lhs", "est")]
  names(gv)[2] <- "variance"
  LM <- merge(gm, gv, by = "lhs"); LM <- LM[order(LM$lhs), ]
  save_vals(LM, "ST_LongInvariance_wide_latentmeans.csv", FIG_DIR)
  cat("\n  스칼라 모형 g 잠재평균 (wave 1 = 0):\n"); print(LM, row.names = FALSE)
}

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat("  이 표가 Supplementary Table 12 를 대체합니다 (multi-group -> wide longitudinal CFA).\n")
cat("  metric 지지 + scalar 판정을 본문 Methods '종단 측정불변성' 문단에 그대로 옮기십시오.\n")
cat("  scalar 가 'not supported' 면: 17_ 의 level/change 모형에 wave 고정효과가 들어 있는지 확인하고,\n")
cat("  partial-scalar 점수로 17_ 을 재실행한 결과를 Supplementary Table 10 각주에 추가하십시오.\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
