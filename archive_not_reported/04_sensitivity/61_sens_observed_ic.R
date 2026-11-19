###############################################################################
## 61_sens_observed_ic.R  (v260817)
## 관측치 전용 민감도 — 대체된 IC 가 결과를 만들었는가에 대한 직접 검증
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
##  추적 파동에서 지표 대체율이 상당합니다(칸 기준 W2 13% -> W5 35%).
##  MNAR 민감도는 "대체값을 악화시켜도 결론 불변"을 보이지만, 심사자는
##  "대체된 관측을 아예 빼도 같은가"를 물을 수 있습니다. 이 스크립트는
##  구간 시작시점의 IC 가 **17개 지표 전부 관측값**인 구간만으로 Table 2
##  전이모형(M0)을 재적합해 주분석과 비교합니다.
##
## ── 대체 마스크 복원 ─────────────────────────────────────────────────────
##  관측 칸은 네 전략(MAIN/Seq/NoFrailty/MNAR)에서 값이 동일하고, 대체 칸만
##  전략에 따라 달라집니다. 따라서 세 대안 전략 중 하나라도 MAIN 과 다른
##  칸 = 대체 칸입니다. (이산형 지표에서 우연히 같은 값이 뽑히면 소폭
##  과소집계될 수 있으나, 세 전략의 합집합이라 그 확률은 작습니다.)
##
##  출력: ST_ObservedIC.xlsx / _values.csv
###############################################################################

OBSIC_VERSION <- "v260817"
message("\n=== 61_sens_observed_ic ", OBSIC_VERSION, " ===")

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
## ── 1. 대체 마스크 (id x wave) ──────────────────────────────────────────
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
rds <- file.path(IMP_DIR, paste0(IMP_STEM, ".rds"))
xx  <- readRDS(rds)
if (is.data.frame(xx)) stop("RDS 에 전략 세트가 없습니다 -- 4-set 파일이 맞는지 확인하십시오.")
Am <- xx$MAIN
ALT3 <- xx[setdiff(names(xx), "MAIN")]
ind <- intersect(IND17, names(Am))
message(sprintf("  지표 %d/17개 확인 | 세트 %d개", length(ind), length(ALT3) + 1))
if (length(ind) < 17)
  warning("지표 ", 17 - length(ind), "개를 찾지 못했습니다: ",
          paste(setdiff(IND17, ind), collapse = ", "))
if (length(ind)) {
  imp_mat <- vapply(ind, function(v) {
    d <- rep(FALSE, nrow(Am))
    for (B in ALT3) if (v %in% names(B)) d <- d | (abs(B[[v]] - Am[[v]]) > 1e-9)
    d
  }, logical(nrow(Am)))
  imp_mat <- matrix(imp_mat, nrow = nrow(Am))
  n_imp <- rowSums(imp_mat, na.rm = TRUE)
} else {
  warning("대체 마스크를 만들 지표가 없습니다 -- 제한 없이 진행합니다 (결과는 주분석과 동일).")
  n_imp <- rep(0L, nrow(Am))
}
MASK <- data.frame(id = as.character(Am$id), wave = Am$wave, n_imp = n_imp)
MASK$obs_ok <- MASK$n_imp == 0
MASK$time <- (MASK$wave - 1) * 2
rm(xx, ALT3); invisible(gc())
message(sprintf("  전지표 관측 person-wave: %d / %d (%.1f%%)",
                sum(MASK$obs_ok), nrow(MASK), 100 * mean(MASK$obs_ok)))

## ── 2. 주분석 vs 관측치 전용 적합 ──────────────────────────────────────
load_imputed("MAIN")
d <- prep_long(load_kfacs())
OKKEY <- paste(MASK$id[MASK$obs_ok], MASK$time[MASK$obs_ok])

to_log <- function(r) { r$b <- log(r$IRR)
  r$se <- (log(r$ucl) - log(r$lcl)) / (2 * stats::qnorm(0.975)); r }

LONG <- NULL
for (sv in c("state", "frailty_3cat", "adl_state")) {
  iv <- try(build_intervals(d, sv, "rolling"), silent = TRUE)
  if (inherits(iv, "try-error") || is.null(iv)) next
  r0 <- to_log(run_transition_table(iv, exposure_terms = "gLIC_z"))
  ivr <- iv[paste(as.character(iv$id), iv$t0) %in% OKKEY, , drop = FALSE]
  message(sprintf("  [%s] 구간 %d -> 관측전용 %d (%.1f%%)",
                  sv, nrow(iv), nrow(ivr), 100 * nrow(ivr) / nrow(iv)))
  r1 <- try(run_transition_table(ivr, exposure_terms = "gLIC_z"), silent = TRUE)
  if (inherits(r1, "try-error") || is.null(r1)) { message("  !! ", sv, " 제한 적합 실패"); next }
  r1 <- to_log(r1)
  k  <- match(r0$transition, r1$transition)
  den <- is.finite(r0$p) & r0$p < 0.05 & abs(r0$b) > log(1.05)
  z <- data.frame(system = sv, transition = r0$transition,
                  ev_full = r0$n_events, irr_full = r0$IRR, lo_full = r0$lcl,
                  hi_full = r0$ucl, p_full = r0$p,
                  ev_obs = r1$n_events[k], irr_obs = r1$IRR[k], lo_obs = r1$lcl[k],
                  hi_obs = r1$ucl[k], p_obs = r1$p[k],
                  chg = ifelse(den & is.finite(r1$b[k]),
                               100 * (r1$b[k] - r0$b) / abs(r0$b), NA_real_),
                  flip_dir = sign(r0$b) != sign(r1$b[k]),
                  sig_chg  = (r0$p < 0.05) != (r1$p[k] < 0.05),
                  stringsAsFactors = FALSE)
  LONG <- rbind(LONG, z)
}
save_vals(LONG, "ST_ObservedIC_values.csv", FIG_DIR)

## ── 3. 표 ───────────────────────────────────────────────────────────────
fmt <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, h), "-")
syslab <- c(state = "Disability, 3 states", frailty_3cat = "Frailty phenotype",
            adl_state = "Any disability (binary)")
TAB <- data.frame(
  System = unname(syslab[LONG$system]),
  Transition = gsub("->", "→", LONG$transition),
  `Events, primary` = LONG$ev_full,
  `IRR, primary (95% CI)` = fmt(LONG$irr_full, LONG$lo_full, LONG$hi_full),
  `Events, observed only` = LONG$ev_obs,
  `IRR, observed only (95% CI)` = fmt(LONG$irr_obs, LONG$lo_obs, LONG$hi_obs),
  `Change in log IRR, %` = ifelse(is.finite(LONG$chg), sprintf("%.1f", LONG$chg), "-"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB <- TAB[order(match(TAB$System, syslab)), ]
TAB$System[duplicated(TAB$System)] <- ""

n_tr <- nrow(LONG); n_q <- sum(is.finite(LONG$chg))
med  <- stats::median(abs(LONG$chg), na.rm = TRUE)
mx   <- suppressWarnings(max(abs(LONG$chg), na.rm = TRUE))
nfl  <- sum(LONG$flip_dir & is.finite(LONG$chg), na.rm = TRUE)
nsig <- sum(LONG$sig_chg, na.rm = TRUE)

save_table(TAB, "ST_ObservedIC",
  title = "Supplementary Table | Transition estimates restricted to intervals with fully observed intrinsic-capacity indicators",
  footnotes = c(
    "Attendance at in-person assessments declined over follow-up, so an increasing share of indicator values was imputed. This analysis restricts each transition model to between-visit intervals whose start-of-interval capacity measurement was based entirely on observed indicators (no imputed value among the 17), and compares the estimates with the primary analysis.",
    "Imputed cells were identified as those at which at least one alternative imputation strategy produced a value different from the primary strategy; observed cells are identical across strategies by construction.",
    "Incidence rate ratios per +1 s.d. of intrinsic capacity from transition-specific Poisson models with a log(person-time) offset and participant-clustered robust standard errors, adjusted for age, sex, baseline comorbidity count, education, household income and area of residence.",
    "Change is computed on the log scale as a percentage of the primary estimate and is quantified only where the primary estimate was significant and non-negligible (P < 0.05 and IRR outside 0.95-1.05).",
    sprintf("Across %d transitions (%d with a quantifiable change) the median absolute change was %.1f%% and the maximum %.1f%%; there were %d direction reversals and %d significance changes.",
            n_tr, n_q, med, mx, nfl, nsig),
    "Outcome states at the destination visit may still incorporate imputed components for participants not assessed in person; the missing-not-at-random analysis bounds that channel separately.",
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  전이 %d개 (%d 산정) | 변화 중앙 %.1f%% · 최대 %.1f%% | 역전 %d | 유의성 변화 %d\n",
            n_tr, n_q, med, mx, nfl, nsig))
cat("=== 완료 ===  표:", file.path(FIG_DIR, "ST_ObservedIC.xlsx"), "\n")
