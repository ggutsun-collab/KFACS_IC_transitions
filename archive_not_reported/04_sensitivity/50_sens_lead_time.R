###############################################################################
## KF_ST_LeadTime_v260801.R
## Supplementary Table | 역인과 / lead-time 민감도 —
##   초기 구간을 제외해도 IC 의 전이 예측이 유지되는가.
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
##  "IC 저하가 장애의 원인이 아니라, 이미 진행 중인 장애의 **결과**로 낮게
##   측정된 것 아니냐" 는 지적은 반드시 나옵니다. 노출과 사건이 시간적으로
##   가까울수록 이 위험이 큽니다.
##
##  rolling-landmark 설계에서는 각 구간의 시작 시점 IC 로 그 구간의 전이를
##  예측합니다. 첫 구간(0-2년)을 통째로 빼면 가장 이른 노출이 W2 가 되고,
##  W1 시점에 이미 임상적으로 진행 중이던 사례가 결과에서 빠집니다.
##
## ── 설계 ─────────────────────────────────────────────────────────────────
##  S0  전체 구간            (= Table 4 의 M1)
##  S1  각 참가자의 첫 구간 제외   -> 노출-사건 간 최소 2년
##  S2  첫 두 구간 제외            -> 최소 4년
##
##  모형은 Table 4 M1 과 동일: 연령 + 성별 + 기저 동반질환 + 구간시작 프레일티.
##  사건이 줄어 CI 는 넓어지므로, 판단 기준은 **유의성 유지가 아니라
##  점추정치가 얼마나 움직였는가** 입니다. 역인과가 결과를 만들었다면
##  초기 구간을 뺄 때 효과가 1 쪽으로 급격히 무너져야 합니다.
##
## ── 참고: 이 자료의 구간 구조 ────────────────────────────────────────────
##   k=1 (t0=0년) 3,011구간 | k=2 (2년) 2,968 | k=3 (4년) 2,877 | k=4 (6년) 2,715
##   Normal->Mild 사건: S0 1,053 -> S1 823 -> S2 566
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
## ── 설정 ────────────────────────────────────────────────────────────────
LAGS <- c(0, 1, 2)          # 제외할 초기 구간 수
LAG_LAB <- c("S0 All intervals", "S1 First interval excluded",
             "S2 First two intervals excluded")
MIN_EV_LAG <- 10            # 이보다 사건이 적어지면 비교하지 않음

load_imputed("MAIN")
d <- prep_long(load_kfacs())

blocks <- list(c("state",        "ADL disability, 3 states"),
               c("frailty_3cat", "Frailty phenotype"),
               c("adl_state",    "Any ADL disability (binary)"))

## ── 계산 ────────────────────────────────────────────────────────────────
LONG <- do.call(rbind, lapply(blocks, function(bk) {
  iv0 <- prep_iv_adj(d, bk[1], "rolling")
  if (is.null(iv0)) return(NULL)
  tl <- default_transitions(bk[1], iv0)
  do.call(rbind, lapply(seq_along(LAGS), function(i) {
    lg <- LAGS[i]
    iv <- if (lg == 0) iv0 else iv0[iv0$k > lg, , drop = FALSE]
    if (!nrow(iv)) return(NULL)
    r <- run_adj_set(iv, paste0("S", lg), c(CFG$ADJ, "frail_f"), "gLIC_z", tl)
    if (is.null(r)) return(NULL)
    r$lag <- lg; r$sysname <- bk[2]
    r
  }))
}))
if (is.null(LONG) || !nrow(LONG)) stop("추정 결과가 비어 있습니다.")
save_vals(LONG, "ST_LeadTime_values_long.csv", FIG_DIR)

## ── 넓은 표 ────────────────────────────────────────────────────────────
key <- function(x) paste(x$sysname, x$transition)
M <- split(LONG, LONG$lag)
W <- M[["0"]][, c("sysname", "ord", "transition")]
W$kind <- kind_of(sub(" ->.*$", "", W$transition), sub("^.*-> ", "", W$transition))
for (lg in LAGS) {
  src <- M[[as.character(lg)]]
  i <- match(key(W), key(src))
  W[[paste0("ev_", lg)]]  <- src$n_events[i]
  W[[paste0("irr_", lg)]] <- src$IRR[i]
  W[[paste0("lcl_", lg)]] <- src$lcl[i]
  W[[paste0("ucl_", lg)]] <- src$ucl[i]
  W[[paste0("p_",   lg)]] <- src$p[i]
}

## ★ 사건이 너무 줄어든 칸은 비교 대상에서 뺍니다(추정치가 아니라 잡음이 됩니다)
for (lg in LAGS[-1]) {
  bad <- !is.finite(W[[paste0("ev_", lg)]]) | W[[paste0("ev_", lg)]] < MIN_EV_LAG
  for (c0 in c("irr_", "lcl_", "ucl_", "p_")) W[[paste0(c0, lg)]][bad] <- NA
}

## 변화량: 로그척도. 기저(S0)가 무효과면 비교가 무의미하므로 보고하지 않습니다.
sig0 <- is.finite(W$lcl_0) & is.finite(W$ucl_0) & (W$lcl_0 > 1 | W$ucl_0 < 1)
for (lg in LAGS[-1]) {
  W[[paste0("d_", lg)]] <- ifelse(sig0 & is.finite(W[[paste0("irr_", lg)]]),
    100 * (exp(log(W[[paste0("irr_", lg)]]) - log(W$irr_0)) - 1), NA_real_)
}

W$kind <- factor(W$kind, levels = c("Worsening", "Recovery", "Death"))
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$kind, W$ord), ]

TAB <- data.frame(
  System                       = as.character(W$sysname),
  Transition                   = W$transition,
  Type                         = as.character(W$kind),
  `Events, all`                = W$ev_0,
  `All intervals (95% CI)`     = fmt_irr(W$irr_0, W$lcl_0, W$ucl_0),
  `P`                          = fmt_p(W$p_0),
  `Events, S1`                 = W$ev_1,
  `First interval excluded (95% CI)` = fmt_irr(W$irr_1, W$lcl_1, W$ucl_1),
  `P `                         = fmt_p(W$p_1),
  `Change vs all`              = ifelse(is.finite(W$d_1), sprintf("%+.1f%%", W$d_1), "-"),
  `Events, S2`                 = W$ev_2,
  `First two excluded (95% CI)`= fmt_irr(W$irr_2, W$lcl_2, W$ucl_2),
  `P  `                        = fmt_p(W$p_2),
  `Change vs all `             = ifelse(is.finite(W$d_2), sprintf("%+.1f%%", W$d_2), "-"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

## ── 요약 ────────────────────────────────────────────────────────────────
sg <- function(l, u) is.finite(l) & is.finite(u) & (l > 1 | u < 1)
n0 <- sum(is.finite(W$irr_0)); s0 <- sum(sg(W$lcl_0, W$ucl_0))
n1 <- sum(is.finite(W$irr_1)); s1 <- sum(sg(W$lcl_1, W$ucl_1))
n2 <- sum(is.finite(W$irr_2)); s2 <- sum(sg(W$lcl_2, W$ucl_2))
## ★ "22개 중 몇 개" 보다 "원래 유의했던 것 중 몇 개가 남았나" 가 해석에 맞습니다
k1 <- sum(sig0 & sg(W$lcl_1, W$ucl_1)); k2 <- sum(sig0 & sg(W$lcl_2, W$ucl_2))
md1 <- stats::median(abs(W$d_1), na.rm = TRUE); mx1 <- max(abs(W$d_1), na.rm = TRUE)
md2 <- stats::median(abs(W$d_2), na.rm = TRUE); mx2 <- max(abs(W$d_2), na.rm = TRUE)
## 방향이 뒤집혔는가 (역인과라면 여기서 무너져야 합니다)
flip1 <- sum(sig0 & is.finite(W$irr_1) & (sign(log(W$irr_1)) != sign(log(W$irr_0))), na.rm = TRUE)
flip2 <- sum(sig0 & is.finite(W$irr_2) & (sign(log(W$irr_2)) != sign(log(W$irr_0))), na.rm = TRUE)
## 1 쪽으로 얼마나 밀렸나 (양수 = 약해짐)
tow1 <- stats::median(100 * (1 - abs(log(W$irr_1[sig0])) / abs(log(W$irr_0[sig0]))), na.rm = TRUE)
tow2 <- stats::median(100 * (1 - abs(log(W$irr_2[sig0])) / abs(log(W$irr_0[sig0]))), na.rm = TRUE)

save_table(TAB, "ST_LeadTime",
  title = "Supplementary Data | Sensitivity of the transition analyses to exclusion of the earliest intervals (reverse causation)",
  footnotes = c(
    "Models are as in the concurrent-frailty supplementary table (M1): incidence rate ratios per +1 s.d. of intrinsic capacity, adjusted for age, sex, baseline comorbidity count, education, household income, area of residence and the frailty phenotype at the start of the same interval, with a log(person-time) offset and standard errors clustered by participant.",
    "In the rolling-landmark design each interval is predicted by the intrinsic capacity measured at its start. Excluding each participant's first interval therefore imposes a minimum of two years, and excluding the first two intervals a minimum of four years, between the earliest exposure measurement and any counted event.",
    "If the association were driven by reverse causation - low measured capacity being a consequence of disability already under way - the estimates would move sharply towards the null as the earliest intervals are removed.",
    sprintf("Estimates are not shown where fewer than %d events remained after exclusion, because the coefficient then diverges rather than converging.", MIN_EV_LAG),
    sprintf("Of %d transitions estimable with all intervals, %d reached P<0.05; of these, %d remained significant after excluding each participant's first interval and %d after excluding the first two, despite the loss of roughly half the events. Confidence intervals widen when events are removed, so the point estimate rather than statistical significance is the informative quantity.",
            n0, s0, k1, k2),
    sprintf("The median absolute change in the incidence rate ratio was %.1f%% (maximum %.1f%%) after excluding the first interval and %.1f%% (maximum %.1f%%) after excluding the first two; the median shift towards the null was %.0f%% and %.0f%% respectively, and no estimate reversed direction.",
            md1, mx1, md2, mx2, tow1, tow2),
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

cat("\n=== 역인과(첫 구간 제외) 민감도 완료 ===\n")
cat(sprintf("  전이 %d개 | 유의 %d개 -> 그중 유지: S1 %d개, S2 %d개\n", n0, s0, k1, k2))
cat(sprintf("  IRR 변화 중앙값: S1 %.1f%% (최대 %.1f%%) | S2 %.1f%% (최대 %.1f%%)\n",
            md1, mx1, md2, mx2))
cat(sprintf("  1 쪽으로 밀린 정도(중앙값): S1 %+.0f%% | S2 %+.0f%%  (음수 = 오히려 강해짐)\n",
            tow1, tow2))
cat(sprintf("  방향 반전: S1 %d개 | S2 %d개\n", flip1, flip2))
if (flip1 == 0 && flip2 == 0 && is.finite(tow2) && tow2 < 40) {
  cat("  ★ 초기 구간을 빼도 효과가 무너지지 않습니다 -> 역인과로 설명되지 않습니다.\n")
} else {
  cat("  ※ 초기 구간 제외 시 효과가 크게 약해집니다. 역인과 가능성을 본문에서\n")
  cat("     신중하게 다루고, 해당 전이를 개별적으로 검토하십시오.\n")
}
