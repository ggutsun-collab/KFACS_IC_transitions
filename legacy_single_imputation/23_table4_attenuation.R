###############################################################################
## KF_Table4_Attenuation_v260801.R
## Table 4 | 전 전이에서 IC 효과가 프레일티 보정 후에도 남는가 (감쇠표)
##
##   M0  기본        : 연령 + 성별 + 기저 동반질환
##   M1  + 표현형    : M0 + interval 시작 시점 프레일티 표현형 (3범주)
##   M2  + CHS 연속  : M0 + interval 시작 시점 CHS 총점(0-5) 표준화값
##
## 왜 M2 가 필요한가
##   "3범주 표현형이 거칠어서 IC 효과가 남는 것 아니냐" 는 반론을 선제 차단합니다.
##   연속점수로 보정해도 남으면, 해상도 문제가 아니라는 근거가 됩니다.
##
## ── 주의: 순환성 ─────────────────────────────────────────────────────────
##   CHS 표현형 5항목(체중감소·소진·악력·보행속도·활동량) 중 4개가 gLIC 구성에도
##   쓰입니다. 따라서 M2 는 **보수적으로 과보정된** 모형입니다.
##   그래도 IC 효과가 남는다면 주장은 더 강해집니다(각주에 명시).
##   비순환 비교(질환 기반 결손지수)는 개별 질환 항목이 필요하므로 원자료가
##   확보되면 별도로 붙이십시오.
##
## ── 감쇠 계산 ────────────────────────────────────────────────────────────
##   (IRR_adj - 1)/(IRR_base - 1) 은 쓰지 않습니다. 기저 IRR 이 1 근처인 전이에서
##   분모가 0 에 가까워져 -1041% / +1948% 같은 값이 나옵니다(실측 확인).
##   로그척도로 계산하고, 기저효과가 무시할 만하면 공란으로 둡니다.
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
if (!exists("ADJ_VERSION") || ADJ_VERSION < "2026-07-31")
  stop("03_functions_incremental.R 이 없습니다/구버전입니다. R/ 폴더를 확인하십시오.")
load_imputed("MAIN")
d <- prep_long(load_kfacs())

blocks <- list(c("state",        "ADL disability, 3 states"),
               c("frailty_3cat", "Frailty phenotype"),
               c("adl_state",    "Any ADL disability (binary)"))

MODELS <- list(
  list(id = "M0", lab = "Base",             adj = CFG$ADJ),
  list(id = "M1", lab = "+ frailty phenotype", adj = c(CFG$ADJ, "frail_f")),
  list(id = "M2", lab = "+ CHS score (continuous)", adj = c(CFG$ADJ, "chs_z")))

LONG <- do.call(rbind, lapply(blocks, function(bk) {
  iv <- prep_iv_adj(d, bk[1], "rolling")
  if (is.null(iv)) return(NULL)
  tl <- default_transitions(bk[1], iv)
  r <- do.call(rbind, lapply(MODELS, function(m)
    run_adj_set(iv, m$id, m$adj, "gLIC_z", tl)))
  r$sysname <- bk[2]
  r
}))
if (is.null(LONG) || !nrow(LONG)) stop("추정 결과가 비어 있습니다.")
save_vals(LONG, "ST_FrailtyAdjustment_values_long.csv", FIG_DIR)

## ── 넓은 표로 ───────────────────────────────────────────────────────────
key <- function(x) paste(x$sysname, x$transition)
M   <- split(LONG, LONG$model)
W   <- M[["M0"]][, c("sysname", "ord", "transition", "n_events", "pyears")]
W$kind <- kind_of(sub(" ->.*$", "", W$transition), sub("^.*-> ", "", W$transition))

for (m in MODELS) {
  src <- M[[m$id]]
  i <- match(key(W), key(src))
  W[[paste0("irr_", m$id)]] <- src$IRR[i]
  W[[paste0("lcl_", m$id)]] <- src$lcl[i]
  W[[paste0("ucl_", m$id)]] <- src$ucl[i]
  W[[paste0("p_",   m$id)]] <- src$p[i]
}
W$att_M1 <- pct_attn(W$irr_M0, W$irr_M1, p_base = W$p_M0)
W$att_M2 <- pct_attn(W$irr_M0, W$irr_M2, p_base = W$p_M0)
## 프레일티 체계처럼 보정항이 자동으로 빠진 행은 중앙값 계산에서 제외
W$app_M1 <- attn_applied(W$irr_M0, W$irr_M1)
W$app_M2 <- attn_applied(W$irr_M0, W$irr_M2)

## 표시 순서: 체계 -> 악화 / 회복 / 사망 -> 중증도
W$kind <- factor(W$kind, levels = c("Worsening", "Recovery", "Death"))
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$kind, W$ord), ]

TAB <- data.frame(
  System                          = as.character(W$sysname),
  Transition                      = W$transition,
  Type                            = as.character(W$kind),
  Events                          = W$n_events,
  `Person-years`                  = round(W$pyears),
  `M0 Base (95% CI)`              = fmt_irr(W$irr_M0, W$lcl_M0, W$ucl_M0),
  `P`                             = fmt_p(W$p_M0),
  `M1 + frailty phenotype (95% CI)` = fmt_irr(W$irr_M1, W$lcl_M1, W$ucl_M1),
  `P `                            = fmt_p(W$p_M1),
  `Attenuation, M0 to M1`         = fmt_attn(W$att_M1),
  `M2 + CHS score (95% CI)`       = fmt_irr(W$irr_M2, W$lcl_M2, W$ucl_M2),
  `P  `                           = fmt_p(W$p_M2),
  `Attenuation, M0 to M2`         = fmt_attn(W$att_M2),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

## ── 본문에 쓸 요약 수치 ────────────────────────────────────────────────
sig  <- function(l, u) is.finite(l) & is.finite(u) & (l > 1 | u < 1)
n_t  <- sum(is.finite(W$irr_M0))
s0   <- sum(sig(W$lcl_M0, W$ucl_M0), na.rm = TRUE)
s1   <- sum(sig(W$lcl_M1, W$ucl_M1), na.rm = TRUE)
s2   <- sum(sig(W$lcl_M2, W$ucl_M2), na.rm = TRUE)
a1   <- stats::median(W$att_M1[W$app_M1], na.rm = TRUE)
a2   <- stats::median(W$att_M2[W$app_M2], na.rm = TRUE)
n1   <- sum(W$app_M1 & is.finite(W$att_M1)); n2 <- sum(W$app_M2 & is.finite(W$att_M2))
rec  <- W[W$kind == "Recovery" & sig(W$lcl_M1, W$ucl_M1), ]

save_table(TAB, "ST_FrailtyAdjustment",
  title = "Supplementary Table | Effect of intrinsic capacity on each state transition, before and after adjustment for concurrent frailty",
  footnotes = c(
    "Incidence rate ratios per +1 s.d. of wave-1 intrinsic capacity, from transition-specific Poisson models with a log(person-time) offset; standard errors are clustered by participant. IRR below 1 indicates a slower transition.",
    "M0 is adjusted for age, sex, baseline comorbidity count, education, household income and area of residence. M1 adds the frailty phenotype and M2 the continuous CHS frailty score, both recorded at the start of the same interval, so exposure and confounder are contemporaneous.",
    "M2 addresses the concern that a three-category phenotype is too coarse: if the effect of intrinsic capacity persists against a continuous frailty score, the residual association is not an artefact of categorisation.",
    "Four of the five CHS phenotype items (weight loss, exhaustion, grip strength, gait speed) also contribute to the intrinsic-capacity factor, so M2 is deliberately over-adjusted and its estimates are conservative.",
    "In the frailty block the origin state is the frailty phenotype itself, so the M1 term is constant within stratum and is dropped; M1 therefore equals M0 there by construction.",
    "Attenuation is computed on the log scale as 100 x (1 - log(IRR_adjusted) / log(IRR_base)) and is left blank when the base effect is negligible (IRR between 0.95 and 1.05), because the ratio is unstable there.",
    sprintf("Of %d estimable transitions, %d reached P<0.05 in M0, %d in M1 and %d in M2.",
            n_t, s0, s1, s2),
    sprintf("Median attenuation was %s in M1 (over the %d transitions in which the frailty term was actually retained) and %s in M2 (%d transitions). Attenuation is not reported where the unadjusted estimate was itself non-significant, because a shift between two null estimates is not interpretable.",
            if (is.finite(a1)) sprintf("%.0f%%", a1) else "not estimable", n1,
            if (is.finite(a2)) sprintf("%.0f%%", a2) else "not estimable", n2),
    sprintf("NE, not estimable: a stratum contributed fewer than %d events.", CFG$MIN_CELL_EVENTS),
    "CHS, Cardiovascular Health Study; CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

cat("\n=== Table 4 완료 ===\n")
cat(sprintf("  추정 가능한 전이 %d개\n", n_t))
cat(sprintf("  유의(P<0.05): M0 %d -> M1 %d -> M2 %d\n", s0, s1, s2))
cat(sprintf("  감쇠 중앙값 : M1 %s (보정 적용 %d전이) | M2 %s (%d전이)\n",
            if (is.finite(a1)) sprintf("%.0f%%", a1) else "-", n1,
            if (is.finite(a2)) sprintf("%.0f%%", a2) else "-", n2))
if (nrow(rec)) {
  cat("  ★ 프레일티 보정 후에도 유의한 '회복' 전이 (이 논문의 차별점):\n")
  for (i in seq_len(nrow(rec)))
    cat(sprintf("      %-22s %s  IRR %.2f (%.2f-%.2f)\n", rec$transition[i],
                rec$sysname[i], rec$irr_M1[i], rec$lcl_M1[i], rec$ucl_M1[i]))
} else cat("  ※ 프레일티 보정 후 유의한 회복 전이가 없습니다 -- 프레이밍을 재검토하십시오.\n")
