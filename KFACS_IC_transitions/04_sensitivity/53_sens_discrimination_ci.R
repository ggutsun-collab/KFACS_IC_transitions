###############################################################################
## KF_ST_DiscriminationCI_v260801.R
## ST_Discrimination 을 '차이의 불확실성까지' 담은 판(ST_Discrimination) 으로 교체
##
## ── 왜 다시 만드는가 ─────────────────────────────────────────────────────
##  기존 표는 AUC 두 개와 그 차이만 있고 차이의 95% CI 도 P 도 없습니다.
##  Nature Aging 심사에서 ΔAUC 를 CI 없이 제시하면 거의 확실히 지적받습니다
##  (+0.039 와 +0.004 를 같은 근거로 읽을 수 없기 때문입니다).
##  또한 기존 표에는 '5-year incident ADL disability, 1년' 칸이 NA 로 남아
##  있는데, 이유를 각주로 밝히지 않으면 계산 실패로 읽힙니다.
##
## ── 방법 ────────────────────────────────────────────────────────────────
##  참가자 단위 부트스트랩(기본 500회). 각 반복에서 두 Cox 모형을 다시 적합해
##  ΔAUC(t) 와 ΔC-index 를 계산하고, 백분위수 95% CI 와 양측 부트스트랩 P 를
##  냅니다. 모형 적합 자체를 반복 안에 넣어야 '변수 추가로 인한 낙관 편향'이
##  CI 에 반영됩니다. C-index 차이의 CI 도 이 표에서 처음 제공됩니다.
##
## ── 실행 시간 ────────────────────────────────────────────────────────────
##  약 4-10분 (BOOT_B = 500 기준). 줄이려면 BOOT_B 를 200 으로 낮추십시오.
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
BOOT_B    <- 500                      # 부트스트랩 반복수
BOOT_SEED <- 20260801
AUC_TIMES <- c(1, 3, 5)

HAS_TROC <- requireNamespace("timeROC", quietly = TRUE)
if (!HAS_TROC) {
  try(install.packages("timeROC", repos = "https://cloud.r-project.org"), silent = TRUE)
  HAS_TROC <- requireNamespace("timeROC", quietly = TRUE)
}

load_imputed("MAIN")
d <- prep_long(load_kfacs()); b <- get_baseline(d)

## ── 결과변수 (SF3 과 동일 구성) ─────────────────────────────────────────
ev <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p) {
  n <- nrow(p); s0 <- p$state[1]; tw <- NA_real_
  if (n >= 2) { k <- which(p$state[-1] > s0 & p$state[-1] < 4); if (length(k)) tw <- p$time[k[1] + 1] }
  data.frame(id = p$id[1], dtime = p$followup_years[n], dev = as.integer(p$death_event[n] == 1),
             wtime = tw, cens = max(p$time[n], p$followup_years[n], na.rm = TRUE),
             s0 = s0, stringsAsFactors = FALSE)
}))
B <- merge(b[, c("id", "gLIC_z", "frailty_lab", "age_c", "sex_f", "comorbid_bl")], ev, by = "id")
B <- B[stats::complete.cases(B[, c("gLIC_z", "frailty_lab", "age_c", "sex_f", "comorbid_bl")]), ]
B$frailty_lab <- droplevels(factor(B$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))
B$wt <- ifelse(is.na(B$wtime), pmin(B$cens, B$dtime, na.rm = TRUE), B$wtime)
B$we <- as.integer(!is.na(B$wtime))

OUTCOMES <- list(
  list(key = "death", lab = "5-year mortality",
       dat = function(z = B) z, tv = "dtime", ev = "dev"),
  list(key = "disab", lab = "5-year incident ADL disability",
       dat = function(z = B) z[is.finite(z$s0) & z$s0 == 1, ], tv = "wt", ev = "we"))
RHS <- c(base = "age_c + sex_f + comorbid_bl + frailty_lab",
         ic   = "age_c + sex_f + comorbid_bl + frailty_lab + gLIC_z")

auc_at <- function(risk, tm, st, t) {
  if (HAS_TROC) {
    r <- try(timeROC::timeROC(T = tm, delta = st, marker = risk, cause = 1,
                              times = t, iid = FALSE), silent = TRUE)
    if (!inherits(r, "try-error")) {
      v <- as.numeric(r$AUC[length(r$AUC)])
      if (is.finite(v)) return(v)
    }
  }
  case <- st == 1 & tm <= t; ctrl <- tm >= t
  if (sum(case) < 5 || sum(ctrl) < 5) return(NA_real_)
  mean(outer(risk[case], risk[ctrl], function(a, b) (a > b) + 0.5 * (a == b)))
}

## 한 자료집합에서 (AUC1,AUC3,AUC5,C) x (base, ic) 를 계산
one_run <- function(z, tv, evv) {
  S <- survival::Surv(z[[tv]], z[[evv]])
  out <- list()
  for (k in names(RHS)) {
    f <- try(survival::coxph(stats::as.formula(paste("S ~", RHS[[k]])), data = z), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    lp <- as.numeric(stats::predict(f, newdata = z, type = "lp"))
    cc <- try(survival::concordance(f)$concordance, silent = TRUE)
    out[[k]] <- c(vapply(AUC_TIMES, function(t) auc_at(lp, z[[tv]], z[[evv]], t), numeric(1)),
                  if (inherits(cc, "try-error")) NA_real_ else cc)
  }
  names(out$base) <- names(out$ic) <- c(paste0("auc", AUC_TIMES), "cindex")
  out
}

set.seed(BOOT_SEED)
RES <- NULL
for (oc in OUTCOMES) {
  z <- oc$dat()
  if (nrow(z) < 100) next
  pt <- one_run(z, oc$tv, oc$ev)
  if (is.null(pt)) next
  ## 시점별 사건수 (NA 칸의 이유를 각주에 쓰기 위해)
  nev <- vapply(AUC_TIMES, function(t) sum(z[[oc$ev]] == 1 & z[[oc$tv]] <= t), numeric(1))

  msg(sprintf("[boot] %s: n=%d, 부트스트랩 %d회 ...", oc$lab, nrow(z), BOOT_B))
  ids <- unique(z$id)
  D <- matrix(NA_real_, BOOT_B, length(pt$base),
              dimnames = list(NULL, names(pt$base)))
  for (bb in seq_len(BOOT_B)) {
    sel <- sample(ids, length(ids), replace = TRUE)
    zb  <- z[match(sel, z$id), , drop = FALSE]
    zb$id <- paste0(zb$id, "_", seq_len(nrow(zb)))
    r <- try(one_run(zb, oc$tv, oc$ev), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) D[bb, ] <- r$ic - r$base
    if (bb %% 100 == 0) msg(sprintf("   %d / %d", bb, BOOT_B))
  }

  for (j in seq_along(pt$base)) {
    nmj <- names(pt$base)[j]
    v   <- D[, j]; v <- v[is.finite(v)]
    dhat <- pt$ic[j] - pt$base[j]
    if (length(v) < 0.5 * BOOT_B || !is.finite(dhat)) {
      lo <- hi <- pv <- NA_real_
    } else {
      qq <- stats::quantile(v, c(0.025, 0.975), names = FALSE)
      lo <- qq[1]; hi <- qq[2]
      pv <- 2 * min(mean(v <= 0), mean(v >= 0)); pv <- min(1, max(pv, 1 / (length(v) + 1)))
    }
    RES <- rbind(RES, data.frame(
      outcome = oc$lab, metric = nmj,
      time = if (nmj == "cindex") NA_real_ else as.numeric(sub("auc", "", nmj)),
      n = nrow(z),
      n_events = if (nmj == "cindex") sum(z[[oc$ev]] == 1) else nev[j],
      base = pt$base[j], ic = pt$ic[j], diff = dhat, lo = lo, hi = hi, p = pv,
      nboot_ok = length(v), stringsAsFactors = FALSE))
  }
}
if (is.null(RES)) stop("모형 적합이 모두 실패했습니다.")
save_vals(RES, "ST_Discrimination_boot.csv", FIG_DIR)

## ── 표 ──────────────────────────────────────────────────────────────────
f3 <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "\u2013")
fd <- function(e, l, h) ifelse(is.finite(e) & is.finite(l),
                               sprintf("%+.3f (%+.3f to %+.3f)", e, l, h),
                               ifelse(is.finite(e), sprintf("%+.3f", e), "\u2013"))
fp <- function(p) ifelse(!is.finite(p), "\u2013",
                  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

A <- RES[RES$metric != "cindex", ]
A <- A[order(match(A$outcome, vapply(OUTCOMES, `[[`, character(1), "lab")), A$time), ]
CX <- RES[RES$metric == "cindex", ]

TAUC <- data.frame(
  Outcome = A$outcome,
  `Time (years)` = sprintf("%g", A$time),
  `Events by t` = A$n_events,
  `AUC, base` = f3(A$base),
  `AUC, base + IC` = f3(A$ic),
  `Difference (95% CI)` = fd(A$diff, A$lo, A$hi),
  P = fp(A$p),
  check.names = FALSE, stringsAsFactors = FALSE)

## C-index 는 시점이 없으므로 결과별 1행으로 아래에 붙입니다.
TC <- data.frame(
  Outcome = CX$outcome, `Time (years)` = "Overall (C-index)",
  `Events by t` = CX$n_events,
  `AUC, base` = f3(CX$base), `AUC, base + IC` = f3(CX$ic),
  `Difference (95% CI)` = fd(CX$diff, CX$lo, CX$hi), P = fp(CX$p),
  check.names = FALSE, stringsAsFactors = FALSE)

TAB <- do.call(rbind, lapply(vapply(OUTCOMES, `[[`, character(1), "lab"), function(L)
  rbind(TAUC[TAUC$Outcome == L, ], TC[TC$Outcome == L, ])))
TAB <- TAB[!is.na(TAB$Outcome), ]
TAB$Outcome[duplicated(TAB$Outcome)] <- ""
rownames(TAB) <- NULL

na_rows <- A[!is.finite(A$base) | !is.finite(A$ic), ]
## 사건이 0건인 칸은 '검정력 부족' 이 아니라 '구조적으로 관측 불가' 입니다.
## 장애 발생은 2년 간격 방문에서만 확인되므로 1년 시점에는 사건이 있을 수 없습니다.
na_note <- NULL
if (nrow(na_rows)) {
  z0 <- na_rows[na_rows$n_events == 0, ]
  zf <- na_rows[na_rows$n_events > 0, ]
  bits <- character(0)
  if (nrow(z0))
    bits <- c(bits, paste0("At ", paste(unique(sprintf("%g year for %s", z0$time, z0$outcome)),
      collapse = " and "), " no events can be observed, because that outcome is ascertained only ",
      "at the biennial study visits; the cell is therefore structurally empty rather than imprecise."))
  if (nrow(zf))
    bits <- c(bits, paste0("The area under the curve is not estimable at ",
      paste(unique(sprintf("%g year for %s (%d events by that time)",
                           zf$time, zf$outcome, zf$n_events)), collapse = "; "),
      ", because too few events had accrued."))
  na_note <- paste(c(bits, "Empty cells are shown as a dash rather than as zero."), collapse = " ")
}

save_table(TAB, "ST_Discrimination",
  title = "Supplementary Table | Discrimination of models with and without intrinsic capacity, at the wave-1 landmark",
  footnotes = c(
    "Cox models fitted at the wave-1 landmark. The base model contains age, sex, baseline comorbidity count, education, household income, area of residence and the frailty phenotype; the second adds the intrinsic-capacity score.",
    if (HAS_TROC) "Time-dependent areas under the curve are cumulative/dynamic and account for censoring by inverse probability of censoring weighting (timeROC)."
    else "The timeROC package was unavailable, so time-dependent areas under the curve were computed among participants whose status was determined by the time point; this is a complete-case approximation.",
    sprintf(paste0("Differences and their confidence intervals come from %d participant-level bootstrap ",
                   "resamples in which both models were refitted within each resample, so the interval ",
                   "reflects the optimism of adding a predictor as well as sampling variability. ",
                   "Intervals are percentile intervals and P values are two-sided bootstrap P values."), BOOT_B),
    na_note,
    "The number of events accrued by each time point is given so that estimates based on few events can be recognised as such; the one-year estimates for mortality rest on 14 events and are correspondingly imprecise.",
    "Discrimination indices change little when a marker is added to a model that already contains strong predictors, and a small difference should not be read as evidence that the marker is uninformative; the calibration and decision-curve analyses and the transition-specific estimates in Table 4 are the more informative comparisons. Conversely, an association can be strong in the transition models while adding little to discrimination for a single time-to-event outcome, because discrimination summarises the ranking of all participants rather than the rate of a particular transition.",
    "AUC, area under the receiver operating characteristic curve; CI, confidence interval; IC, intrinsic capacity."))

cat("\n=== ST_Discrimination (차이의 CI 포함) 완료 ===\n")
print(TAB, row.names = FALSE)
cat("\n  원자료: ST_Discrimination_boot.csv\n")
