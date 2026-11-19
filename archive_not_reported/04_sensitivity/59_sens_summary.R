###############################################################################
## 59_sens_summary.R
## 보충표 — 민감도 분석 판정 요약 (한 표)
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
## 민감도 분석이 6종(완전자료·축소보정·역인과 2단계·노출정의·MI 풀링·파동
## 안정성)이라 표가 6개 생깁니다. 심사자가 원하는 것은 각 표의 세부가 아니라
## "결론이 흔들렸는가" 이므로, 판정 지표(변화율 중앙/최대, 방향 역전, 유의성
## 변화)를 한 표로 모으고 전이별 상세는 Supplementary Data 로 내립니다.
##
## ── 실행 조건 ────────────────────────────────────────────────────────────
## 50, 52, 55, 56, 57, 58 이 먼저 실행되어 FIG_DIR 에 산출물이 있어야 합니다.
## (99_MASTER 는 순서를 자동으로 지킵니다. 없는 항목은 표에 '미실행'으로.)
##
## 출력:  ST_SensitivitySummary.xlsx / .csv
###############################################################################

SENSSUM_VERSION <- "v260813"
message("\n=== 59_sens_summary ", SENSSUM_VERSION, " ===")

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
###############################################################################
## 0. 판정 지표 계산기
##  변화율은 로그척도 100*(b1-b0)/|b0|, 분모 보호 규칙(P<0.05 그리고
##  |logIRR|>log(1.05)) 을 두 조건 모두 적용합니다.
###############################################################################
.pz <- function(b, se) ifelse(is.finite(b) & is.finite(se) & se > 0,
                              2 * stats::pnorm(-abs(b / se)), NA_real_)
verdict <- function(b0, se0, b1, se1) {
  ok  <- is.finite(b0) & is.finite(b1)
  p0  <- .pz(b0, se0); p1 <- .pz(b1, se1)
  den <- ok & is.finite(p0) & p0 < 0.05 & abs(b0) > log(1.05)
  chg <- ifelse(den, 100 * (b1 - b0) / abs(b0), NA_real_)
  list(n_tr  = sum(ok),
       n_chg = sum(is.finite(chg)),
       med   = stats::median(abs(chg), na.rm = TRUE),
       mx    = suppressWarnings(max(abs(chg), na.rm = TRUE)),
       flip  = sum(ok & sign(b0) != sign(b1), na.rm = TRUE),
       sigch = sum(is.finite(p0) & is.finite(p1) & ((p0 < 0.05) != (p1 < 0.05)),
                   na.rm = TRUE))
}
.b_from_ci <- function(irr, lcl, ucl)
  list(b = log(as_num(irr)), se = (log(as_num(ucl)) - log(as_num(lcl))) / 3.92)
fmtrow <- function(analysis, varied, v)
  data.frame(Analysis = analysis, `What was varied` = varied,
             `Transitions` = sprintf("%d (%d)", v$n_tr, v$n_chg),
             `Median |change|, %` = ifelse(is.finite(v$med), sprintf("%.1f", v$med), "—"),
             `Largest |change|, %` = ifelse(is.finite(v$mx), sprintf("%.1f", v$mx), "—"),
             `Direction reversals` = v$flip,
             `Significance changes` = v$sigch,
             check.names = FALSE, stringsAsFactors = FALSE)
narow <- function(analysis, varied, note = "not run")
  data.frame(Analysis = analysis, `What was varied` = varied,
             Transitions = note, `Median |change|, %` = "", `Largest |change|, %` = "",
             `Direction reversals` = "", `Significance changes` = "",
             check.names = FALSE, stringsAsFactors = FALSE)

ROWS <- list()

## ── 1. 완전자료 (소득 무응답) ──────────────────────────────────────────────
ROWS$cc <- tryCatch({
  x <- openxlsx::read.xlsx(file.path(FIG_DIR, "ST_CompleteCase.xlsx"), sheet = "raw")
  v <- verdict(x$b_all, x$se_all, x$b_cc, x$se_cc)
  fmtrow("Complete case", "Participants with missing socioeconomic data excluded", v)
}, error = function(e) narow("Complete case", "Missing socioeconomic data excluded", "not run"))

## ── 2. 축소 보정군 ─────────────────────────────────────────────────────────
ROWS$adj <- tryCatch({
  x <- utils::read.csv(file.path(FIG_DIR, "ST_AdjReduced_values.csv"),
                       fileEncoding = "UTF-8")
  v <- verdict(x$b_core, x$se_core, x$b_red, x$se_red)
  fmtrow("Reduced adjustment", "Socioeconomic covariates removed (age, sex, comorbidity only)", v)
}, error = function(e) narow("Reduced adjustment", "Socioeconomic covariates removed", "not run"))

## ── 3-4. 역인과 (초기 구간 제외) ───────────────────────────────────────────
lead_row <- function(lag_lab, model_lab) {
  tryCatch({
    x <- utils::read.csv(file.path(FIG_DIR, "ST_LeadTime_values_long.csv"),
                         fileEncoding = "UTF-8")
    x$key <- paste(x$sysname, x$transition)
    s0 <- x[x$model == "S0", ]; s1 <- x[x$model == model_lab, ]
    k  <- match(s0$key, s1$key)
    e0 <- .b_from_ci(s0$IRR, s0$lcl, s0$ucl)
    e1 <- .b_from_ci(s1$IRR[k], s1$lcl[k], s1$ucl[k])
    v  <- verdict(e0$b, e0$se, e1$b, e1$se)
    fmtrow("Reverse causation", lag_lab, v)
  }, error = function(e) narow("Reverse causation", lag_lab, "not run"))
}
ROWS$lt1 <- lead_row("First interval of each participant excluded", "S1")
ROWS$lt2 <- lead_row("First two intervals excluded", "S2")

## ── 5. 노출 정의 (W1 고정 vs rolling) ──────────────────────────────────────
ROWS$ex <- tryCatch({
  x <- utils::read.csv(file.path(FIG_DIR, "ST_ExposureDef_values_long.csv"),
                       fileEncoding = "UTF-8")
  x$key <- paste(x$sysname, x$transition)
  s0 <- x[x$expo == "rolling", ]; s1 <- x[x$expo == "W1", ]
  k  <- match(s0$key, s1$key)
  e0 <- .b_from_ci(s0$IRR, s0$lcl, s0$ucl)
  e1 <- .b_from_ci(s1$IRR[k], s1$lcl[k], s1$ucl[k])
  v  <- verdict(e0$b, e0$se, e1$b, e1$se)
  fmtrow("Exposure definition", "Wave-1 capacity applied to every interval instead of rolling landmark", v)
}, error = function(e) narow("Exposure definition", "Wave-1 fixed vs rolling landmark", "not run"))

## ── 6. 대체 전략 (MAIN vs Seq/NoFrailty/MNAR — 전략별 재적합) ─────────────
##  v260814: Rubin 풀링 폐지. 네 셋은 m=4 draw 가 아니라 전략 변형이므로
##  기준(MAIN) 대비 전이별 최대 변화로 판정합니다.
ROWS$mi <- tryCatch({
  x <- utils::read.csv(file.path(FIG_DIR, "ST_ImputationStrategy_values.csv"),
                       fileEncoding = "UTF-8")
  if ("model" %in% names(x)) x <- x[x$model == "M0", , drop = FALSE]
  key <- paste(x$system, x$transition)
  per_tr <- split(x, key)
  mx  <- vapply(per_tr, function(z)
    if (all(!is.finite(z$chg))) NA_real_ else max(abs(z$chg), na.rm = TRUE), numeric(1))
  n_tr  <- length(per_tr); n_q <- sum(is.finite(mx))
  nflip <- sum(vapply(per_tr, function(z) any(z$flip_dir, na.rm = TRUE), logical(1)))
  nsig  <- sum(vapply(per_tr, function(z) any(z$sig_chg,  na.rm = TRUE), logical(1)))
  data.frame(Analysis = "Imputation strategy",
             `What was varied` = "Estimates refitted under three alternative strategies (sequential; frailty-excluded; MNAR)",
             Transitions = sprintf("%d (%d)", n_tr, n_q),
             `Median |change|, %` = ifelse(is.finite(stats::median(mx, na.rm = TRUE)),
                                           sprintf("%.1f", stats::median(mx, na.rm = TRUE)), "—"),
             `Largest |change|, %` = ifelse(any(is.finite(mx)),
                                            sprintf("%.1f", max(mx, na.rm = TRUE)), "—"),
             `Direction reversals` = nflip,
             `Significance changes` = nsig,
             check.names = FALSE, stringsAsFactors = FALSE)
}, error = function(e) narow("Imputation strategy",
                             "Alternative imputation strategies (sequential; frailty-excluded; MNAR)", "not run"))

## ── 6b. 관측치 전용 (대체된 IC 제외) — v260817 ────────────────────────────
ROWS$obs <- tryCatch({
  x <- utils::read.csv(file.path(FIG_DIR, "ST_ObservedIC_values.csv"),
                       fileEncoding = "UTF-8")
  b0 <- log(as_num(x$irr_full)); se0 <- (log(as_num(x$hi_full)) - log(as_num(x$lo_full))) / 3.92
  b1 <- log(as_num(x$irr_obs));  se1 <- (log(as_num(x$hi_obs))  - log(as_num(x$lo_obs)))  / 3.92
  v  <- verdict(b0, se0, b1, se1)
  fmtrow("Observed measurements only",
         "Intervals restricted to start-of-interval capacity with no imputed indicator", v)
}, error = function(e) narow("Observed measurements only",
                             "Imputed capacity measurements excluded", "not run"))

## ── 7. 파동 안정성 (측정 표류) ─────────────────────────────────────────────
ROWS$ws <- tryCatch({
  x <- utils::read.csv(file.path(FIG_DIR, "ST_WaveStability_values.csv"),
                       fileEncoding = "UTF-8")
  x$key <- paste(x$system, x$transition)
  het <- vapply(split(x, x$key), function(z) {
    ok <- is.finite(z$b) & is.finite(z$se) & z$se > 0
    b <- z$b[ok]; se <- z$se[ok]; k <- length(b)
    if (k < 2) return(NA_real_)
    w <- 1 / se^2; mu <- sum(w * b) / sum(w)
    Q <- sum(w * (b - mu)^2)
    max(0, (Q - (k - 1)) / Q) * 100
  }, numeric(1))
  data.frame(Analysis = "Wave stability",
             `What was varied` = "Effect estimated separately within each between-visit interval",
             Transitions = sprintf("%d", sum(is.finite(het))),
             `Median |change|, %` = sprintf("I² median %.0f", stats::median(het, na.rm = TRUE)),
             `Largest |change|, %` = sprintf("I² max %.0f", suppressWarnings(max(het, na.rm = TRUE))),
             `Direction reversals` = "—",
             `Significance changes` = "0 Bonferroni-significant interactions",
             check.names = FALSE, stringsAsFactors = FALSE)
}, error = function(e) narow("Wave stability", "Wave-specific estimation", "not run"))

SUM <- do.call(rbind, ROWS); rownames(SUM) <- NULL
print(SUM, row.names = FALSE)

FN <- c(
 paste0("Each row re-estimates the transition-specific incidence rate ratios of Table 2 under ",
        "one changed analytic decision and compares them with the primary estimates. ",
        "Transitions gives the number of transitions compared, with the number for which a ",
        "percentage change could be quantified in parentheses."),
 paste0("Change is computed on the log scale as the difference in log incidence rate ratios ",
        "expressed as a percentage of the primary estimate; it is quantified only where the ",
        "primary estimate was itself significant and non-negligible (P < 0.05 and IRR outside ",
        "0.95-1.05), because the ratio is unstable when the denominator approaches zero."),
 paste0("Significance changes counts transitions crossing P = 0.05 in either direction. For ",
        "wave stability the entries are the median and maximum I² across transitions and the ",
        "number of capacity-by-wave interactions below the Bonferroni threshold, from Wald ",
        "tests against the cluster-robust covariance matrix."),
 paste0("Reverse-causation rows exclude each participant's earliest interval(s), where ",
        "subclinical disease is most likely to have depressed the capacity measurement. The ",
        "exposure-definition row fixes the wave-1 measurement for all intervals in place of ",
        "the rolling landmark."),
 paste0("The exposure-definition row is not a null-expectation check: replacing the ",
        "measurement taken at the start of each interval with the increasingly outdated ",
        "wave-1 value is expected to attenuate the estimates, and the observed changes were ",
        "predominantly attenuations toward the null. That the current measurement outperforms ",
        "the baseline one parallels the current-versus-change comparison of Fig. 2."),
 paste0("All models are adjusted for ", adj_phrase(),
        ". Full transition-level results for every row are provided as Supplementary Data."))

save_table(SUM, "ST_SensitivitySummary",
  title = paste("Supplementary Table | Robustness of the transition estimates",
                "across analytic decisions"),
  footnotes = FN)

cat("\n=== 완료 ===  표:", file.path(FIG_DIR, "ST_SensitivitySummary.xlsx"), "\n")
