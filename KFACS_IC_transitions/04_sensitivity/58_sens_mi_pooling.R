###############################################################################
## 58_sens_mi_pooling.R  (v260814 — 전략별 비교판)
## 대체 "전략" 민감도 — MAIN / Seq / NoFrailty / MNAR 을 나란히 제시
##
## ── 왜 Rubin 통합이 아닌가 (v260814 변경 사유) ───────────────────────────
##  RDS 의 네 셋은 같은 대체모형에서 뽑은 m=4 draw 가 아니라 서로 다른
##  대체 전략입니다: MAIN(주분석) · Seq(순차 대체) · NoFrailty(쇠약 항목
##  제외 모형) · MNAR(비무작위결측 시나리오).
##  Rubin's rules 는 m 개 자료가 같은 사후예측분포에서 독립 추출되었다는
##  가정 위에 성립합니다. 전략이 다른 셋을 통합하면 대체 불확실성이 아니라
##  모형 간 차이를 분산으로 계산하게 되고 FMI 는 해석 불가능한 수가 됩니다
##  (구판에서 관찰된 max FMI 0.84 가 그 인공물입니다).
##  올바른 제시는 네 전략의 추정치를 나란히 놓고 결론이 바뀌는지 보이는
##  것입니다. MNAR 까지 견디면 그 자체가 강한 근거입니다.
##
## ── 출력 ─────────────────────────────────────────────────────────────────
##  ST_ImputationStrategy.xlsx / _values.csv
##  (구 ST_MIpooling* 은 실행 시 _retired_260813 으로 이동)
###############################################################################

MISTRAT_VERSION <- "v260814"
message("\n=== 58_sens_mi_pooling (전략별 비교) ", MISTRAT_VERSION, " ===")

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
## 1. 대체셋 목록
###############################################################################
imp_sets <- function(dir = IMP_DIR, stem = IMP_STEM) {
  rds <- file.path(dir, paste0(stem, ".rds"))
  if (file.exists(rds)) {
    x <- readRDS(rds)
    if (is.data.frame(x)) return(NULL)
    nm <- names(x); return(nm[vapply(x, is.data.frame, logical(1))])
  }
  NULL
}
SETS <- imp_sets()
if (is.null(SETS) || length(SETS) < 2)
  stop("대체셋을 찾지 못했습니다. IMP_DIR / IMP_STEM 을 확인하십시오.")
REF <- if ("MAIN" %in% SETS) "MAIN" else SETS[1]
ALT <- setdiff(SETS, REF)
msg(sprintf("대체셋 %d개: %s (기준 = %s)", length(SETS), paste(SETS, collapse = ", "), REF))

## 표시용 전략 이름 (영문)
STRAT_LAB <- c(MAIN = "Primary", Seq = "Sequential",
               NoFrailty = "Frailty-excluded", MNAR = "MNAR")
slab <- function(s) ifelse(s %in% names(STRAT_LAB), STRAT_LAB[s], s)

###############################################################################
## 2. 셋마다 전이모형 적합 (M0 = CORE 6변수)
###############################################################################
to_logscale <- function(r) {
  r$b  <- log(r$IRR)
  r$se <- (log(r$ucl) - log(r$lcl)) / (2 * stats::qnorm(0.975))
  r
}
one_set <- function(set) {
  msg(sprintf("[%s] 적합 중 ...", set))
  load_imputed(set)
  d <- prep_long(load_kfacs())
  out <- NULL
  for (sv in c("state", "frailty_3cat", "adl_state")) {
    iv <- try(build_intervals(d, sv, "rolling"), silent = TRUE)
    if (inherits(iv, "try-error") || is.null(iv)) next
    r0 <- to_logscale(run_transition_table(iv, exposure_terms = "gLIC_z"))
    r0$model <- "M0"; r0$system <- sv
    out <- rbind(out, r0[, c("system","model","transition","from","to","term",
                             "IRR","lcl","ucl","p","n_events","pyears","b","se")])
  }
  if (is.null(out)) stop('[', set, '] 전이모형을 하나도 적합하지 못했습니다.')
  w1 <- d[d$time == min(d$time, na.rm = TRUE), ]
  attr(out, "glic") <- stats::setNames(w1$gLIC_z, as.character(w1$id))
  out$set <- set
  out
}
RES <- lapply(SETS, one_set); names(RES) <- SETS

## gLIC 지문 — 전략마다 IC 점수가 실제로 다른지
G <- lapply(RES, attr, "glic")
ids <- Reduce(intersect, lapply(G, names))
GM <- do.call(cbind, lapply(G, function(v) v[ids]))
cr <- stats::cor(GM, use = "pairwise.complete.obs")
same <- all(abs(cr[lower.tri(cr)] - 1) < 1e-10)
cat("\ngLIC 의 셋 간 상관:\n"); print(round(cr, 6))

###############################################################################
## 3. 전략별 비교 (기준 = MAIN, 분모 보호 규칙 동일)
###############################################################################
ALL <- do.call(rbind, RES)
ALL$key <- paste(ALL$system, ALL$model, ALL$transition, ALL$term, sep = "|")
M0REF <- ALL[ALL$set == REF, ]

LONG <- do.call(rbind, lapply(split(ALL, ALL$key), function(z) {
  zr <- z[z$set == REF, ][1, ]
  do.call(rbind, lapply(ALT, function(s) {
    za <- z[z$set == s, ]
    if (!nrow(za)) return(NULL)
    za <- za[1, ]
    den <- is.finite(zr$p) && zr$p < 0.05 && abs(zr$b) > log(1.05)
    data.frame(system = zr$system, model = zr$model, transition = zr$transition,
               term = zr$term, n_events = zr$n_events, strategy = s,
               irr_main = zr$IRR, lo_main = zr$lcl, hi_main = zr$ucl, p_main = zr$p,
               irr_alt = za$IRR, lo_alt = za$lcl, hi_alt = za$ucl, p_alt = za$p,
               chg = if (den) 100 * (za$b - zr$b) / abs(zr$b) else NA_real_,
               flip_dir = sign(zr$b) != sign(za$b),
               sig_chg  = (zr$p < 0.05) != (za$p < 0.05),
               stringsAsFactors = FALSE)
  }))
}))
rownames(LONG) <- NULL
save_vals(LONG, "ST_ImputationStrategy_values.csv", FIG_DIR)

###############################################################################
## 4. 표 — 전이별로 네 전략의 IRR 을 나란히
###############################################################################
fmt <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, h), "-")
W <- M0REF
W$key2 <- paste(W$system, W$transition)
TAB <- data.frame(System = W$system,
                  Transition = gsub("->", "→", W$transition),
                  Events = W$n_events,
                  check.names = FALSE, stringsAsFactors = FALSE)
TAB[[paste0(slab(REF), " (95% CI)")]] <- fmt(W$IRR, W$lcl, W$ucl)
for (s in ALT) {
  za <- LONG[LONG$strategy == s, ]
  k  <- match(paste(W$system, W$transition), paste(za$system, za$transition))
  TAB[[paste0(slab(s), " (95% CI)")]] <- fmt(za$irr_alt[k], za$lo_alt[k], za$hi_alt[k])
}
mx_by <- vapply(seq_len(nrow(W)), function(i) {
  z <- LONG[LONG$system == W$system[i] & LONG$transition == W$transition[i], ]
  if (all(!is.finite(z$chg))) NA_real_ else max(abs(z$chg), na.rm = TRUE)
}, numeric(1))
sig_all <- vapply(seq_len(nrow(W)), function(i) {
  z <- LONG[LONG$system == W$system[i] & LONG$transition == W$transition[i], ]
  all(c(W$p[i], z$p_alt) < 0.05)
}, logical(1))
TAB[["Largest |change| vs primary, %"]] <-
  ifelse(is.finite(mx_by), sprintf("%.1f", mx_by), "-")
TAB[["Significant in all four"]] <- ifelse(sig_all, "yes",
  ifelse(W$p < 0.05, "no", "n.s. in primary"))
syslab <- c(state = "Disability, 3 states", frailty_3cat = "Frailty phenotype",
            adl_state = "Any disability (binary)")
TAB$System <- unname(syslab[TAB$System])
TAB <- TAB[order(match(TAB$System, syslab)), ]
TAB$System[duplicated(TAB$System)] <- ""

n_tr  <- nrow(W); n_q <- sum(is.finite(mx_by))
med_x <- stats::median(mx_by, na.rm = TRUE); max_x <- max(mx_by, na.rm = TRUE)
nflip <- sum(LONG$flip_dir, na.rm = TRUE); nsig <- sum(LONG$sig_chg, na.rm = TRUE)

save_table(TAB, "ST_ImputationStrategy",
  title = "Supplementary Data | Transition estimates under the primary and three alternative imputation strategies",
  footnotes = c(
    "The four completed longitudinal datasets are not multiple-imputation draws from a single model but deliberate strategy variants: the primary strategy used in all main analyses; a sequential-imputation variant; a variant excluding the frailty components from the imputation model; and a missing-not-at-random (MNAR) scenario. Rubin's rules do not apply across strategy variants, so estimates are presented side by side rather than pooled.",
    "Incidence rate ratios per +1 s.d. of intrinsic capacity from the primary transition models (adjusted for age, sex, baseline comorbidity count, education, household income and area of residence), refitted in full within each completed dataset.",
    "Change is computed on the log scale as the difference from the primary-strategy estimate expressed as a percentage, and is quantified only where the primary estimate was significant and non-negligible (P < 0.05 and IRR outside 0.95-1.05).",
    sprintf("Across %d transitions (%d with a quantifiable change) the median of the largest per-transition change was %.1f%% and the maximum %.1f%%; there were %d direction reversals and %d significance changes across all strategy comparisons.",
            n_tr, n_q, med_x, max_x, nflip, nsig),
    if (same) "The intrinsic-capacity score itself was identical across the four datasets; the strategies differ in the imputation of covariate and state information." else "The intrinsic-capacity score differed across datasets, so the comparison also reflects sensitivity of the score to the imputation strategy.",
    "CI, confidence interval; IRR, incidence rate ratio; MNAR, missing not at random; s.d., standard deviation."))

## 구판 산출물 정리
RET <- file.path(FIG_DIR, "_retired_260813"); dir.create(RET, showWarnings = FALSE)
for (f in list.files(FIG_DIR, pattern = "^ST_MIpooling", full.names = TRUE))
  if (file.rename(f, file.path(RET, basename(f)))) cat("  구판 이동:", basename(f), "\n")

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  전이 %d개 | 최대 변화 중앙 %.1f%% · 최대 %.1f%% | 방향 역전 %d | 유의성 변화 %d\n",
            n_tr, med_x, max_x, nflip, nsig))
cat("=== 완료 ===  표:", file.path(FIG_DIR, "ST_ImputationStrategy.xlsx"), "\n")
