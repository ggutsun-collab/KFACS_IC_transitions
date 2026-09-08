###############################################################################
## 18c_longitudinal_cfa_tables.R   (v260908)
## 18b 가 저장한 LongCFA_pairs_fits.rds 에서 표만 만듭니다 (재적합 없음).
##  - fitMeasures 는 모형당 한 번만 호출 (여러 번 호출하면 FIML baseline 을 매번 다시 계산해 매우 느림)
##  - 모형마다 진행 상황을 출력
##  출력: ST_LongInvariance_pairs.csv, ST_LongInvariance_pairs_means.csv,
##        ST_LongInvariance_pairs_fixedtheta.csv
##  실행: source("18c_longitudinal_cfa_tables.R", encoding = "UTF-8")   (2-5분)
###############################################################################
message("\n=== 18c_longitudinal_cfa_tables ===")

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
suppressPackageStartupMessages(library(lavaan))

RDS <- file.path(FIG_DIR, "LongCFA_pairs_fits.rds")
if (!file.exists(RDS)) stop("적합 결과가 없습니다: ", RDS, "  -> 18b 를 먼저 실행하십시오.")
FITS <- readRDS(RDS)
cat("  읽은 쌍:", paste(names(FITS), collapse = ", "), "\n")

## ── 1. 적합도: 모형당 fitMeasures 한 번 ────────────────────────────────────
WANT <- c("chisq.scaled", "df.scaled", "cfi.robust", "tli.robust", "rmsea.robust",
          "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
one_fm <- function(f, tag) {
  if (is.null(f)) { cat(sprintf("  %-22s (없음)\n", tag)); return(rep(NA_real_, 7)) }
  t0 <- Sys.time()
  m <- try(fitMeasures(f, WANT), silent = TRUE)
  if (inherits(m, "try-error")) m <- try(fitMeasures(f), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-22s fitMeasures 실패\n", tag)); return(rep(NA_real_, 7)) }
  g <- function(x) { v <- suppressWarnings(as.numeric(m[x])); if (length(v) && is.finite(v)) v else NA_real_ }
  cfi <- g("cfi.robust"); if (!is.finite(cfi)) cfi <- g("cfi.scaled")
  tli <- g("tli.robust"); if (!is.finite(tli)) tli <- g("tli.scaled")
  rms <- g("rmsea.robust"); if (!is.finite(rms)) rms <- g("rmsea.scaled")
  out <- c(g("chisq.scaled"), g("df.scaled"), cfi, tli, rms, g("srmr"), as.numeric(lavInspect(f, "nobs")))
  cat(sprintf("  %-22s CFI %.4f RMSEA %.4f SRMR %.4f  (%.1f초)\n", tag, cfi, rms, g("srmr"),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out
}
crit <- function(dc, dr, ds, step) if (!is.finite(dc)) "-" else
  if (dc >= -0.010 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "supported" else
  if (dc >= -0.015 && (dr <= 0.015 || ds <= ifelse(step == "metric", 0.030, 0.010))) "borderline (CFI criterion only)" else "not supported"

cat("\n[적합도 계산]\n")
TAB <- do.call(rbind, lapply(names(FITS), function(key) {
  M <- rbind(one_fm(FITS[[key]]$configural, paste(key, "configural")),
             one_fm(FITS[[key]]$metric,     paste(key, "metric")),
             one_fm(FITS[[key]]$scalar,     paste(key, "scalar")))
  d <- data.frame(Pair = key, Model = c("Configural", "Metric (loadings)", "Scalar (+ intercepts)"),
                  n = M[, 7], `Chi-square` = M[, 1], df = M[, 2], CFI = M[, 3], TLI = M[, 4], RMSEA = M[, 5], SRMR = M[, 6],
                  check.names = FALSE, stringsAsFactors = FALSE)
  d$`Delta CFI` <- c(NA, diff(d$CFI)); d$`Delta RMSEA` <- c(NA, diff(d$RMSEA)); d$`Delta SRMR` <- c(NA, diff(d$SRMR))
  d$Verdict <- c("-", crit(d$`Delta CFI`[2], d$`Delta RMSEA`[2], d$`Delta SRMR`[2], "metric"),
                 crit(d$`Delta CFI`[3], d$`Delta RMSEA`[3], d$`Delta SRMR`[3], "scalar"))
  d
}))
save_vals(TAB, "ST_LongInvariance_pairs.csv", FIG_DIR)
cat("\n"); print(TAB, row.names = FALSE, digits = 4)

## ── 2. scalar 모형의 g 잠재평균 차 (뒤 wave - 앞 wave; 앞 wave 평균 0, SD 1) ──
cat("\n[잠재평균]\n")
MEANS <- do.call(rbind, lapply(names(FITS), function(key) {
  f <- FITS[[key]]$scalar; if (is.null(f)) return(NULL)
  pe <- parameterEstimates(f); b <- as.integer(sub("W\\d+-W", "", key))
  r <- pe[pe$op == "~1" & pe$lhs == paste0("g_w", b), ]
  v <- pe[pe$op == "~~" & pe$lhs == paste0("g_w", b) & pe$rhs == paste0("g_w", b), "est"]
  data.frame(Pair = key, `Latent mean difference (later - earlier), s.d. units` = r$est,
             `95% CI lower` = r$ci.lower, `95% CI upper` = r$ci.upper, `Later-wave variance` = v,
             check.names = FALSE) }))
save_vals(MEANS, "ST_LongInvariance_pairs_means.csv", FIG_DIR)
print(MEANS, row.names = FALSE, digits = 3)

## ── 3. 고정된 잔차분산 / 남은 음수 잔차분산 ──────────────────────────────────
cat("\n[잔차분산 점검]\n")
HEY <- do.call(rbind, lapply(names(FITS), function(key) do.call(rbind, lapply(names(FITS[[key]]), function(lv) {
  f <- FITS[[key]][[lv]]; if (is.null(f)) return(NULL)
  fx <- attr(f, "fixed_theta"); if (is.null(fx) || !length(fx)) return(NULL)
  data.frame(Pair = key, Model = lv, Indicator = names(fx), `Fixed residual variance` = as.numeric(fx), check.names = FALSE) }))))
if (!is.null(HEY)) { save_vals(HEY, "ST_LongInvariance_pairs_fixedtheta.csv", FIG_DIR); print(HEY, row.names = FALSE) }
NEG <- NULL
for (key in names(FITS)) for (lv in names(FITS[[key]])) { f <- FITS[[key]][[lv]]; if (is.null(f)) next
  th <- diag(lavInspect(f, "est")$theta)
  if (any(th < 0)) NEG <- rbind(NEG, data.frame(Pair = key, Model = lv, Indicator = names(th)[th < 0], Value = th[th < 0])) }
if (is.null(NEG)) cat("  남은 음수 잔차분산 없음\n") else { cat("  ※ 남은 음수 잔차분산:\n"); print(NEG, row.names = FALSE, digits = 3) }

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat("  metric  판정:", paste(TAB$Verdict[TAB$Model == "Metric (loadings)"], collapse = " | "), "\n")
cat("  scalar  판정:", paste(TAB$Verdict[TAB$Model == "Scalar (+ intercepts)"], collapse = " | "), "\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
