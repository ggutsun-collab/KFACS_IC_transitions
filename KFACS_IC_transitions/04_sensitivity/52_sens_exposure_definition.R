###############################################################################
## KF_ST_ExposureDef_v260801.R
## Supplementary Table | 노출 정의 민감도 — rolling landmark vs wave-1 고정
##
## ── 무엇에 답하는가 ──────────────────────────────────────────────────────
##  주분석은 각 구간의 **시작 시점** IC 로 그 구간의 전이를 예측합니다
##  (rolling landmark). 대안은 모든 구간에 **wave-1** IC 를 쓰는 것입니다.
##
##    rolling : 최신 정보를 쓰므로 예측력이 높지만, 노출이 이미 진행 중인
##              변화를 반영할 수 있습니다(역인과 우려 -> ST_LeadTime 에서 검증).
##    W1 고정 : 노출이 확실히 선행하지만, 8년 전 값이라 회귀희석이 큽니다.
##
##  두 정의에서 방향과 크기가 유지되면, 결과가 노출 정의의 산물이 아닙니다.
##
## ── 해석 지침 ────────────────────────────────────────────────────────────
##  W1 고정에서 효과가 **약해지는 것은 정상**입니다 (측정시점이 멀수록
##  회귀희석). 중요한 것은 **방향이 유지되는가** 이지 크기가 같은가가 아닙니다.
###############################################################################

###############################################################################
## 이 스크립트는 작업디렉터리와 무관하게 동작합니다.
## 자기 위치를 찾아 R/_bootstrap.R 을 불러오고, kf_init() 이 헬퍼를 적재합니다.
## (source() / RStudio Source 버튼 / Rscript 모두 지원)
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
if (!exists("ADJ_VERSION")) stop("03_functions_incremental.R 을 R/ 폴더에 두십시오 (kf_init 이 적재합니다).")
load_imputed("MAIN")
d <- prep_long(load_kfacs())

blocks <- list(c("state",        "ADL disability, 3 states"),
               c("frailty_3cat", "Frailty phenotype"),
               c("adl_state",    "Any ADL disability (binary)"))

EXPO <- c(rolling = "Rolling landmark (primary)", W1 = "Wave-1 exposure")

LONG <- do.call(rbind, lapply(blocks, function(bk) {
  do.call(rbind, lapply(names(EXPO), function(e) {
    iv <- prep_iv_adj(d, bk[1], e)
    if (is.null(iv)) return(NULL)
    tl <- default_transitions(bk[1], iv)
    r  <- run_adj_set(iv, e, c(CFG$ADJ, "frail_f"), "gLIC_z", tl)
    if (is.null(r)) return(NULL)
    r$expo <- e; r$sysname <- bk[2]
    r
  }))
}))
if (is.null(LONG) || !nrow(LONG)) stop("추정 결과가 비어 있습니다.")
save_vals(LONG, "ST_ExposureDef_values_long.csv", FIG_DIR)

key <- function(x) paste(x$sysname, x$transition)
M <- split(LONG, LONG$expo)
W <- M[["rolling"]][, c("sysname", "ord", "transition", "n_events", "pyears")]
W$kind <- kind_of(sub(" ->.*$", "", W$transition), sub("^.*-> ", "", W$transition))
for (e in names(EXPO)) {
  src <- M[[e]]; i <- match(key(W), key(src))
  W[[paste0("irr_", e)]] <- src$IRR[i]; W[[paste0("lcl_", e)]] <- src$lcl[i]
  W[[paste0("ucl_", e)]] <- src$ucl[i]; W[[paste0("p_",   e)]] <- src$p[i]
}
sg <- function(l, u) is.finite(l) & is.finite(u) & (l > 1 | u < 1)
W$same_dir <- is.finite(W$irr_rolling) & is.finite(W$irr_W1) &
              (sign(log(W$irr_rolling)) == sign(log(W$irr_W1)))
W$ratio <- ifelse(is.finite(W$irr_rolling) & is.finite(W$irr_W1) &
                  abs(log(W$irr_rolling)) > 0.05,
                  log(W$irr_W1) / log(W$irr_rolling), NA_real_)

W$kind <- factor(W$kind, levels = c("Worsening", "Recovery", "Death"))
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$kind, W$ord), ]

TAB <- data.frame(
  System                    = as.character(W$sysname),
  Transition                = W$transition,
  Type                      = as.character(W$kind),
  Events                    = W$n_events,
  `Rolling landmark (95% CI)` = fmt_irr(W$irr_rolling, W$lcl_rolling, W$ucl_rolling),
  `P`                       = fmt_p(W$p_rolling),
  `Wave-1 exposure (95% CI)` = fmt_irr(W$irr_W1, W$lcl_W1, W$ucl_W1),
  `P `                      = fmt_p(W$p_W1),
  `Same direction`          = ifelse(is.na(W$same_dir), "-", ifelse(W$same_dir, "yes", "no")),
  `Effect retained`         = ifelse(is.finite(W$ratio), sprintf("%.0f%%", 100 * W$ratio), "-"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

n_est  <- sum(is.finite(W$irr_rolling) & is.finite(W$irr_W1))
n_same <- sum(W$same_dir, na.rm = TRUE)
s_r <- sum(sg(W$lcl_rolling, W$ucl_rolling)); s_w <- sum(sg(W$lcl_W1, W$ucl_W1))
md_ret <- stats::median(100 * W$ratio, na.rm = TRUE)

save_table(TAB, "ST_ExposureDef",
  title = "Supplementary Data | Sensitivity of the transition analyses to the definition of exposure",
  footnotes = c(
    "Models are as in the concurrent-frailty supplementary table (M1): adjusted for age, sex, baseline comorbidity count, education, household income, area of residence and the frailty phenotype at the start of the same interval, with a log(person-time) offset and standard errors clustered by participant.",
    "In the primary (rolling landmark) analysis each interval is predicted by the intrinsic capacity measured at its start; in the alternative every interval is predicted by the wave-1 measurement, which is up to eight years earlier.",
    "Attenuation under the wave-1 definition is expected, because a measurement taken further in advance is subject to greater regression dilution. The informative comparison is therefore whether the direction is preserved, not whether the magnitude is identical.",
    "Effect retained is the ratio of the two log incidence rate ratios, expressed as a percentage; it is left blank where the primary estimate is negligible (IRR between 0.95 and 1.05).",
    sprintf("Of %d transitions estimable under both definitions, %d had the same direction; %d reached P<0.05 under the rolling-landmark definition and %d under the wave-1 definition, and the median proportion of the effect retained was %.0f%%.",
            n_est, n_same, s_r, s_w, md_ret),
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

cat("\n=== 노출 정의 민감도 완료 ===\n")
cat(sprintf("  양쪽 추정 가능 %d개 | 방향 동일 %d개\n", n_est, n_same))
cat(sprintf("  유의: rolling %d개 -> W1 %d개 | 효과 유지율 중앙값 %.0f%%\n", s_r, s_w, md_ret))
if (n_est > 0 && n_same == n_est)
  cat("  ★ 노출 정의를 바꿔도 방향이 모두 유지됩니다.\n")
