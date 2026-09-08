###############################################################################
## 18b_longitudinal_cfa_pairs.R   (v260908)
## 종단 측정불변성 — 인접 wave 쌍별 longitudinal CFA (within-person dependence 반영)
##
##  왜 쌍(pair)인가
##   85변수 5-wave 단일 모형(18_)은 모형당 90분 이상 걸리고 Heywood case 가 났습니다.
##   인접 wave 쌍(W1-W2, W2-W3, W3-W4, W4-W5)마다 34변수 모형을 적합하면
##   (i) 같은 사람의 두 측정 사이 의존성은 요인 상관 + 같은 지표 잔차 상관으로 처리되고,
##   (ii) 2년 변화 분석(17_)이 쓰는 것이 정확히 인접 두 wave 의 수준이므로,
##       인접 wave 간 scalar invariance 가 그 분석이 요구하는 바로 그 검정입니다.
##
##  출력 (FIG_DIR): ST_LongInvariance_pairs.csv       (쌍 × configural/metric/scalar 적합도)
##                  ST_LongInvariance_pairs_means.csv (scalar 모형의 g 잠재평균 차, 뒤 wave - 앞 wave)
##                  LongCFA_pairs_fits.rds
##  실행 : source("18b_longitudinal_cfa_pairs.R", encoding = "UTF-8")   (쌍당 1-3분 예상)
##  v260908c: 두 wave 모두 참석자로 한정, bounds 제거, 앞 모형 추정치를 시작값으로 사용, 쌍마다 중간 저장,
##            표 생성부 수정(nobs), Heywood 표 저장
##  v260908e: bounds 완전 제거(제약 솔버 회피). Heywood 항목의 잔차분산을 FIX_THETA 값으로 고정하고
##            새로 음수가 나오는 항목은 자동으로 고정 목록에 추가해 재적합
##  v260908f: 잔차분산을 고정한 항목은 wave 간 잔차 공분산도 0 으로 고정 (theta 양정치 보장)
###############################################################################

LCFA_VERSION <- "v260908f"
message("\n=== 18b_longitudinal_cfa_pairs ", LCFA_VERSION, " ===")

## ── 헬퍼 적재 ──────────────────────────────────────────────────────────
## 위치 탐색과 헬퍼 로드는 R/_bootstrap.R 한 곳에서만 정의합니다.
## 모든 경로(PROJECT_ROOT / IMP_DIR / FIG_DIR)는 R/00_setup.R 에서만 정합니다.
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
.need(c("lavaan", "readxl"))
suppressPackageStartupMessages(library(lavaan))

###############################################################################
## 1. 설정
###############################################################################
PRE_FILE <- get0("PRE_FILE", ifnotfound = file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260618_state.xlsx"))  # 00_setup.R 에서 정함
WAVES    <- 1:5
PAIRS    <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 5))
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
SPEC <- list(loco = c("Balance","GS_ms","rev_CST"),
             vita = c("loss_of_Bwt","Appetite","EXH","HGS"),
             cogn = c("Orientation","Memory","Attention_Calculation","Language","Visuospatial"),
             psyc = c("Negative_affect","Positive_affect","Motivation"))
EXCL_IDS <- c("kf161224", "kf170602", "kf171295")
VERBOSE  <- TRUE           # 반복 로그 출력 (진행 확인용). 끄려면 FALSE
REFIT    <- TRUE           # FALSE 면 저장된 LongCFA_pairs_fits.rds 를 읽어 표만 다시 만듭니다
FIX_THETA <- c(Appetite = 0.05)   # Heywood 항목의 잔차분산 고정값 (표준화 지표 기준; 신뢰도 0.95 에 해당)
MAX_REFIT <- 3                    # 새 Heywood 항목이 나오면 고정 목록에 추가해 최대 이 횟수까지 재적합

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
for (v in IND17) pre[[v]] <- (pre[[v]] - MU[[v]]) / SD[[v]]       # 12_/13_ 와 같은 풀 표준화
msg(sprintf("관측 person-wave %d행, 참가자 %d명", nrow(pre), length(unique(pre$id))))

ids <- sort(unique(pre$id))
Wd <- data.frame(id = ids, stringsAsFactors = FALSE)
for (w in WAVES) {
  z <- pre[pre$wave == w, c("id", IND17)]; z <- z[!duplicated(z$id), ]
  names(z)[-1] <- paste0(IND17, "_w", w)
  Wd <- merge(Wd, z, by = "id", all.x = TRUE)
}

###############################################################################
## 3. 쌍별 lavaan 모형 생성기 (wave a, b)
##    - g_w, loco_w, vita_w, cogn_w, psyc_w ; 감각 2지표는 g 에만
##    - wave 내 bifactor 직교 ; wave 간 g~~g, 같은 specific~~specific 자유, 그 외 0
##    - 같은 지표 잔차 wave 간 상관 자유
##    - 식별: wave a 요인분산 1·잠재평균 0 ; metric 이상에서 wave b 분산 자유,
##            scalar 에서 wave b 의 g 잠재평균 자유 (specific 평균 0 고정)
###############################################################################
build_pair <- function(a, b, level = c("configural", "metric", "scalar"), fix = FIX_THETA) {
  level <- match.arg(level); L <- character(0)
  lab <- function(prefix, items) if (level == "configural") "" else paste0(prefix, "_", items, "*")
  for (w in c(a, b)) {
    iw <- function(x) paste0(x, "_w", w)
    L <- c(L, sprintf("g_w%d =~ %s", w, paste0(lab("lg", IND17), iw(IND17), collapse = " + ")))
    for (s in names(SPEC))
      L <- c(L, sprintf("%s_w%d =~ %s", s, w, paste0(lab(paste0("l", s), SPEC[[s]]), iw(SPEC[[s]]), collapse = " + ")))
    ilab <- if (level == "scalar") paste0("nu_", IND17, "*") else ""
    L <- c(L, sprintf("%s ~ %s1", iw(IND17), ilab))
    facs <- c(paste0("g_w", w), paste0(names(SPEC), "_w", w))
    L <- c(L, if (w == a || level == "configural") sprintf("%s ~~ 1*%s", facs, facs) else sprintf("%s ~~ NA*%s", facs, facs))
    L <- c(L, sprintf("%s ~ %s1", facs, ifelse(facs == paste0("g_w", w) & w == b & level == "scalar", "NA*", "0*")))
    for (i in seq_along(facs)) for (j in seq_along(facs)) if (i < j) L <- c(L, sprintf("%s ~~ 0*%s", facs[i], facs[j]))
  }
  L <- c(L, sprintf("g_w%d ~~ g_w%d", a, b))
  for (s in names(SPEC)) L <- c(L, sprintf("%s_w%d ~~ %s_w%d", s, a, s, b),
                                sprintf("g_w%d ~~ 0*%s_w%d", a, s, b), sprintf("%s_w%d ~~ 0*g_w%d", s, a, b))
  for (s in names(SPEC)) for (t in names(SPEC)) if (s != t) L <- c(L, sprintf("%s_w%d ~~ 0*%s_w%d", s, a, t, b))
  for (v in IND17) L <- c(L, if (v %in% names(fix)) sprintf("%s_w%d ~~ 0*%s_w%d", v, a, v, b) else sprintf("%s_w%d ~~ %s_w%d", v, a, v, b))
  ## Heywood 항목: 잔차분산을 작은 양수로 고정 (bounds 대신 — 제약 솔버를 피함)
  for (v in names(fix)) for (w in c(a, b)) L <- c(L, sprintf("%s_w%d ~~ %s*%s_w%d", v, w, format(fix[[v]]), v, w))
  paste(L, collapse = "\n")
}

fit_pair <- function(a, b, level, start_fit = NULL, fix = FIX_THETA) {
  ## 두 wave 모두 참석한 사람으로 한정: 결측 패턴이 항목 수준 소수로 줄어 FIML 이 빠릅니다.
  va <- paste0(IND17, "_w", a); vb <- paste0(IND17, "_w", b)
  both <- rowSums(!is.na(Wd[, va])) > 0 & rowSums(!is.na(Wd[, vb])) > 0
  dat  <- Wd[both, c(va, vb)]
  for (round in seq_len(MAX_REFIT)) {
    t0 <- Sys.time()
    f <- try(cfa(build_pair(a, b, level, fix), data = dat, estimator = "MLR", missing = "fiml",
                 meanstructure = TRUE, std.lv = FALSE, auto.fix.first = FALSE,
                 bounds = "none",
                 start = if (is.null(start_fit)) "default" else start_fit,   # 앞 모형 추정치에서 출발
                 verbose = VERBOSE,
                 control = list(iter.max = 3000)), silent = TRUE)
    if (inherits(f, "try-error")) { message(sprintf("  !! W%d-W%d %s 실패: %s", a, b, level, conditionMessage(attr(f, "condition")))); return(NULL) }
    cat(sprintf("  W%d-W%d %-10s n=%d converged=%s %.1f분 | 고정: %s\n", a, b, level, nrow(dat),
                lavInspect(f, "converged"), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                paste(names(fix), collapse = ",")))
    th <- diag(lavInspect(f, "est")$theta)
    neg <- unique(sub("_w\\d+$", "", names(th)[th < 0]))
    neg <- setdiff(neg, names(fix))
    if (!length(neg)) break
    cat(sprintf("     * 음수 잔차분산 발생: %s -> 고정 목록에 추가하고 재적합\n", paste(neg, collapse = ", ")))
    fix <- c(fix, setNames(rep(FIX_THETA[[1]], length(neg)), neg))
    start_fit <- NULL                                        # 구조가 바뀌었으므로 기본 시작값
  }
  attr(f, "fixed_theta") <- fix
  f
}

###############################################################################
## 4. 적합
###############################################################################
FITS <- list()
if (!REFIT && file.exists(file.path(FIG_DIR, "LongCFA_pairs_fits.rds"))) {
  FITS <- readRDS(file.path(FIG_DIR, "LongCFA_pairs_fits.rds")); msg("저장된 적합 결과를 읽었습니다 (REFIT = FALSE)")
} else for (p in PAIRS) {
  key <- sprintf("W%d-W%d", p[1], p[2]); msg(sprintf("[%s]", key))
  f0 <- fit_pair(p[1], p[2], "configural")
  fx <- if (is.null(f0)) FIX_THETA else attr(f0, "fixed_theta")   # configural 에서 확정된 고정 목록을 이어서 사용
  f1 <- fit_pair(p[1], p[2], "metric", start_fit = f0, fix = fx)
  fx <- if (is.null(f1)) fx else attr(f1, "fixed_theta")
  f2 <- fit_pair(p[1], p[2], "scalar", start_fit = f1, fix = fx)
  FITS[[key]] <- list(configural = f0, metric = f1, scalar = f2)
  saveRDS(FITS, file.path(FIG_DIR, "LongCFA_pairs_fits.rds"))     # 쌍마다 중간 저장
}
saveRDS(FITS, file.path(FIG_DIR, "LongCFA_pairs_fits.rds"))

###############################################################################
## 5. 표  (적합이 끝난 뒤에는 이 절만 따로 실행해도 됩니다: rds 에서 읽습니다)
###############################################################################
if (!exists("FITS")) FITS <- readRDS(file.path(FIG_DIR, "LongCFA_pairs_fits.rds"))
fm <- function(f) {
  if (is.null(f)) return(rep(NA_real_, 7))
  g <- function(x) { v <- try(unname(fitMeasures(f, x)), silent = TRUE)
                     if (inherits(v, "try-error") || !length(v)) NA_real_ else as.numeric(v) }
  cfi <- g("cfi.robust"); tli <- g("tli.robust"); rm <- g("rmsea.robust")
  if (!is.finite(cfi)) cfi <- g("cfi.scaled")
  if (!is.finite(tli)) tli <- g("tli.scaled")
  if (!is.finite(rm))  rm  <- g("rmsea.scaled")
  c(g("chisq.scaled"), g("df.scaled"), cfi, tli, rm, g("srmr"), as.numeric(lavInspect(f, "nobs")))
}
crit <- function(dc, dr, ds, step) if (!is.finite(dc)) "-" else
  if (dc >= -0.010 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "supported" else
  if (dc >= -0.015 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "borderline (CFI criterion only)" else "not supported"
TAB <- do.call(rbind, lapply(names(FITS), function(key) {
  M <- do.call(rbind, lapply(FITS[[key]], fm))
  d <- data.frame(Pair = key, Model = c("Configural", "Metric (loadings)", "Scalar (+ intercepts)"),
                  n = M[, 7], `Chi-square` = M[, 1], df = M[, 2], CFI = M[, 3], TLI = M[, 4], RMSEA = M[, 5], SRMR = M[, 6],
                  check.names = FALSE, stringsAsFactors = FALSE)
  d$`Delta CFI` <- c(NA, diff(d$CFI)); d$`Delta RMSEA` <- c(NA, diff(d$RMSEA)); d$`Delta SRMR` <- c(NA, diff(d$SRMR))
  d$Verdict <- c("-", crit(d$`Delta CFI`[2], d$`Delta RMSEA`[2], d$`Delta SRMR`[2], "metric"),
                 crit(d$`Delta CFI`[3], d$`Delta RMSEA`[3], d$`Delta SRMR`[3], "scalar"))
  d
}))
save_vals(TAB, "ST_LongInvariance_pairs.csv", FIG_DIR)
print(TAB, row.names = FALSE, digits = 4)

## scalar 모형: 뒤 wave g 잠재평균 (앞 wave = 0, SD = 1 기준) = 코호트 평균 IC 의 2년 변화
MEANS <- do.call(rbind, lapply(names(FITS), function(key) {
  f <- FITS[[key]]$scalar; if (is.null(f)) return(NULL)
  pe <- parameterEstimates(f); b <- as.integer(sub("W\\d+-W", "", key))
  r <- pe[pe$op == "~1" & pe$lhs == paste0("g_w", b), ]
  v <- pe[pe$op == "~~" & pe$lhs == paste0("g_w", b) & pe$rhs == paste0("g_w", b), "est"]
  data.frame(Pair = key, `Latent mean difference (later - earlier), s.d. units` = r$est,
             `95% CI lower` = r$ci.lower, `95% CI upper` = r$ci.upper, `Later-wave variance` = v,
             check.names = FALSE) }))
save_vals(MEANS, "ST_LongInvariance_pairs_means.csv", FIG_DIR)
cat("\n  scalar 모형 g 잠재평균 차 (뒤 wave - 앞 wave):\n"); print(MEANS, row.names = FALSE, digits = 3)

## 고정된 잔차분산 목록 (각주용)
HEY <- do.call(rbind, lapply(names(FITS), function(key) do.call(rbind, lapply(names(FITS[[key]]), function(lv) {
  f <- FITS[[key]][[lv]]; if (is.null(f)) return(NULL)
  fx <- attr(f, "fixed_theta"); if (is.null(fx) || !length(fx)) return(NULL)
  data.frame(Pair = key, Model = lv, Indicator = names(fx), `Fixed residual variance` = as.numeric(fx), check.names = FALSE) }))))
if (!is.null(HEY)) { save_vals(HEY, "ST_LongInvariance_pairs_fixedtheta.csv", FIG_DIR)
  cat("\n  고정된 잔차분산 (Heywood 대응):\n"); print(HEY, row.names = FALSE) }
## 남은 음수 잔차분산이 있으면 경고
for (key in names(FITS)) for (lv in names(FITS[[key]])) { f <- FITS[[key]][[lv]]; if (is.null(f)) next
  th <- diag(lavInspect(f, "est")$theta); if (any(th < 0)) cat(sprintf("  ※ %s %s 에 여전히 음수 잔차분산: %s\n", key, lv, paste(names(th)[th < 0], collapse = ", "))) }

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat("  metric  판정:", paste(TAB$Verdict[TAB$Model == "Metric (loadings)"], collapse = " | "), "\n")
cat("  scalar  판정:", paste(TAB$Verdict[TAB$Model == "Scalar (+ intercepts)"], collapse = " | "), "\n")
cat("  -> 이 표가 Supplementary Table 12 를 대체합니다. 네 쌍 모두 scalar 지지면 2년 변화 분석의\n")
cat("     수준 비교 전제가 직접 확인된 것입니다. 일부 쌍이 borderline 이면 17_ 에 wave 고정효과를\n")
cat("     넣고(Methods 문단 참조) partial-scalar 점수 민감도를 ST 10 각주에 추가하십시오.\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
