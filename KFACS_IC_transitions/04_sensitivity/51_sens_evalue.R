###############################################################################
## KF_ST_Evalue_v260801.R
## Supplementary Table | 미측정 교란 민감도 (E-value)
##
## ── 무엇에 답하는가 ──────────────────────────────────────────────────────
##  "연령·성별·동반질환·프레일티를 보정했다지만, 측정 안 한 교란이 결과를
##   만든 것 아니냐" 는 지적에 숫자로 답합니다.
##
##  E-value = 관측된 연관을 완전히 설명해 없애려면, 미측정 교란요인이
##            노출(IC)과도 결과(전이)와도 **각각 최소 몇 배**의 연관을
##            가져야 하는가 (이미 보정한 공변량을 넘어서).
##            (VanderWeele & Ding, Ann Intern Med 2017)
##
##  두 개를 함께 봅니다.
##    · 점추정치의 E-value  : 연관을 0 으로 만들려면 필요한 강도
##    · 신뢰구간의 E-value  : 연관을 '비유의' 로 만들려면 필요한 강도
##      (귀무가설에 가까운 쪽 한계를 씁니다. CI 가 1 을 포함하면 E=1)
##
## ── 계산 ─────────────────────────────────────────────────────────────────
##  RR >= 1 : E = RR + sqrt(RR * (RR - 1))
##  RR <  1 : 먼저 1/RR 로 뒤집은 뒤 같은 식을 적용 (보호효과도 대칭)
##
##  ★ 이 식은 위험비/율비(rate ratio)에 직접 적용됩니다. 본 분석의 IRR 은
##    Poisson 율비이므로 변환 없이 그대로 씁니다. (위험비가 아니라 위험도
##    'hazard' 비였다면 흔한 결과에서 별도 근사가 필요합니다.)
##
## ── 왜 벤치마크가 필요한가 ───────────────────────────────────────────────
##  E-value 만 던지면 "3.2 가 큰가 작은가" 를 독자가 판단할 수 없습니다.
##  같은 모형에서 **이미 측정한 가장 강한 교란요인**(프레일티 표현형)의 연관
##  강도를 함께 보고합니다. 미측정 교란이 그보다 훨씬 강해야 한다면,
##  그런 요인이 남아 있을 가능성은 낮다고 논증할 수 있습니다.
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

## ── E-value ─────────────────────────────────────────────────────────────
evalue <- function(rr) {
  r <- ifelse(is.finite(rr) & rr > 0, ifelse(rr < 1, 1 / rr, rr), NA_real_)
  ifelse(is.na(r) | r < 1, NA_real_, r + sqrt(r * (r - 1)))
}
## 신뢰구간의 E-value: 귀무가설(1)에 가까운 쪽 한계. CI 가 1 을 포함하면 1.
evalue_ci <- function(rr, lcl, ucl) {
  lim <- ifelse(is.finite(rr) & rr > 1, lcl, ucl)
  cross <- (is.finite(rr) & rr > 1 & is.finite(lcl) & lcl <= 1) |
           (is.finite(rr) & rr < 1 & is.finite(ucl) & ucl >= 1)
  out <- evalue(lim)
  ifelse(!is.finite(lim), NA_real_, ifelse(cross, 1, out))
}

## ── 추정 ────────────────────────────────────────────────────────────────
## 주모형은 Table 4 의 M1 과 동일합니다.
IC <- do.call(rbind, lapply(blocks, function(bk) {
  iv <- prep_iv_adj(d, bk[1], "rolling")
  if (is.null(iv)) return(NULL)
  tl <- default_transitions(bk[1], iv)
  r <- run_adj_set(iv, "M1", c(CFG$ADJ, "frail_f"), "gLIC_z", tl)
  if (!is.null(r)) { r$sysname <- bk[2] }
  r
}))
if (is.null(IC) || !nrow(IC)) stop("추정 결과가 비어 있습니다.")

###############################################################################
## 벤치마크: 프레일티(Frail vs Robust)가 이 자료에서 실제로 얼마나 강한
## 교란요인인가.
##
## ★ 수정 (2026-07-31): 벤치마크 모형에 gLIC_z 를 넣으면 안 됩니다.
##   E-value 벤치마크가 묻는 것은 "이 상황에서 알려진 교란요인이 결과와
##   얼마나 강하게 연관되는가" 입니다. 노출(IC)까지 보정하면 IC 가 프레일티의
##   연관을 흡수해 벤치마크가 과소평가됩니다 -- 실행값에서 1.16 까지 내려갔고,
##   그러면 "미측정 교란은 프레일티보다 2배 강해야 한다" 는 논증이 실제보다
##   유리하게 부풀려집니다. 노출을 뺀 모형에서 뽑습니다.
##   (VanderWeele 의 benchmarking 은 다른 측정 공변량에는 조건화하되
##    노출에는 조건화하지 않습니다.)
###############################################################################
BM <- do.call(rbind, lapply(blocks, function(bk) {
  if (bk[1] == "frailty_3cat") return(NULL)   # 출발상태가 곧 프레일티 -> 불가
  iv <- prep_iv_adj(d, bk[1], "rolling")
  if (is.null(iv)) return(NULL)
  tl <- default_transitions(bk[1], iv)
  r <- run_adj_set(iv, "BM", CFG$ADJ, "frail_f", tl)   # ★ gLIC_z 제외
  if (!is.null(r)) { r$sysname <- bk[2] }
  r
}))

## ── 표 ──────────────────────────────────────────────────────────────────
W <- IC[, c("sysname", "ord", "transition", "n_events", "pyears",
            "IRR", "lcl", "ucl", "p")]
W$kind <- kind_of(sub(" ->.*$", "", W$transition), sub("^.*-> ", "", W$transition))
W$E_est <- evalue(W$IRR)
W$E_ci  <- evalue_ci(W$IRR, W$lcl, W$ucl)
W$sig   <- is.finite(W$lcl) & is.finite(W$ucl) & (W$lcl > 1 | W$ucl < 1)

## 벤치마크 붙이기 (Frail vs Robust 항만)
if (!is.null(BM) && nrow(BM)) {
  bmf <- BM[grepl("Frail$", BM$term) & !grepl("Pre", BM$term), ]
  if (!nrow(bmf)) bmf <- BM[grepl("frail_f", BM$term), ]
  i <- match(paste(W$sysname, W$transition), paste(bmf$sysname, bmf$transition))
  W$bm_irr <- bmf$IRR[i]; W$bm_lcl <- bmf$lcl[i]; W$bm_ucl <- bmf$ucl[i]
} else { W$bm_irr <- NA; W$bm_lcl <- NA; W$bm_ucl <- NA }

W$kind <- factor(W$kind, levels = c("Worsening", "Recovery", "Death"))
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$kind, W$ord), ]
save_vals(W, "ST_Evalue_values.csv", FIG_DIR)

TAB <- data.frame(
  System                  = as.character(W$sysname),
  Transition              = W$transition,
  Type                    = as.character(W$kind),
  Events                  = W$n_events,
  `IRR per +1 s.d. (95% CI)` = fmt_irr(W$IRR, W$lcl, W$ucl),
  `P`                     = fmt_p(W$p),
  `E-value, estimate`     = ifelse(is.finite(W$E_est), sprintf("%.2f", W$E_est), "-"),
  `E-value, CI limit`     = ifelse(is.finite(W$E_ci),  sprintf("%.2f", W$E_ci),  "-"),
  `Frailty (frail vs robust), same model` =
      ifelse(is.finite(W$bm_irr), fmt_irr(W$bm_irr, W$bm_lcl, W$bm_ucl), "-"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

## ── 요약 ────────────────────────────────────────────────────────────────
sw <- W[W$sig & is.finite(W$E_est), ]
md <- stats::median(sw$E_est, na.rm = TRUE)
mn <- min(sw$E_est, na.rm = TRUE); mx <- max(sw$E_est, na.rm = TRUE)
mdc <- stats::median(sw$E_ci, na.rm = TRUE)
bm_med <- stats::median(W$bm_irr[is.finite(W$bm_irr)], na.rm = TRUE)
top <- sw[order(-sw$E_est), ][seq_len(min(5, nrow(sw))), ]

save_table(TAB, "ST_Evalue",
  title = "Supplementary Table | Sensitivity of the transition analyses to unmeasured confounding (E-values)",
  footnotes = c(
    "Incidence rate ratios are from the primary models (the concurrent-frailty supplementary table, M1): adjusted for age, sex, baseline comorbidity count, education, household income, area of residence and the frailty phenotype at the start of the same interval.",
    "The E-value is the minimum strength of association, on the risk-ratio scale, that an unmeasured confounder would need to have with both intrinsic capacity and the transition, above and beyond the measured covariates, to explain away the observed association (VanderWeele and Ding, Ann Intern Med 2017). The E-value for the confidence limit is the strength needed to move the interval to include the null; it is 1 where the interval already includes the null.",
    "Protective associations are inverted before the E-value is computed, so the scale is symmetric around the null.",
    "The last column gives the association of the frailty phenotype (frail versus robust) with the same transition, adjusted for age, sex, comorbidity count, education, household income and area of residence but not for intrinsic capacity, as a benchmark for how strong a known confounder is in this setting; an unmeasured confounder would have to exceed this to explain the findings. Intrinsic capacity is deliberately excluded from the benchmark model, because conditioning on the exposure would absorb part of the confounder-outcome association and understate the benchmark. The benchmark is not available for the frailty system, where the origin state is the phenotype itself.",
    sprintf("Across the %d transitions that reached P<0.05, E-values for the point estimate ranged from %.2f to %.2f (median %.2f) and the median E-value for the confidence limit was %.2f.",
            nrow(sw), mn, mx, md, mdc),
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

cat("\n=== E-value (미측정 교란 민감도) 완료 ===\n")
cat(sprintf("  유의한 전이 %d개 | E-value 중앙값 %.2f (범위 %.2f-%.2f) | CI 기준 중앙값 %.2f\n",
            nrow(sw), md, mn, mx, mdc))
if (is.finite(bm_med))
  cat(sprintf("  벤치마크: 프레일티(frail vs robust) 연관 중앙값 IRR = %.2f  (IC 미보정)\n", bm_med))
cat("  E-value 상위 전이:\n")
for (i in seq_len(nrow(top)))
  cat(sprintf("    %-24s IRR %.2f -> E %.2f (CI 기준 %.2f)\n",
              top$transition[i], top$IRR[i], top$E_est[i], top$E_ci[i]))
cat("\n  해석: 위 값보다 강한 미측정 교란이 (이미 보정한 변수들을 넘어서)\n")
cat("        남아 있어야만 결과가 설명됩니다.\n")
