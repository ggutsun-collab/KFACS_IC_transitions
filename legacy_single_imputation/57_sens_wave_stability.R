###############################################################################
## KF_ST_WaveStability_v260811.R
## 보충표/보충그림 — 파동별 내재역량 효과의 안정성
##   (측정 불변성: 완전 스칼라 미지지에 대한 실증적 보강)
##
## ── 왜 이 분석이 필요한가 ────────────────────────────────────────────────
## ST1 에서 완전 스칼라 불변성이 지지되지 않았습니다. 지표 절편이 파동마다
## 표류하면, 참 능력이 같은 사람이 파동에 따라 다른 점수를 받습니다.
##
## ST3(부분 스칼라 점수로 재적합)이 있지만 그 논거에는 약점이 있습니다.
## 두 점수의 상관이 0.9967 이라 거의 같은 점수를 두 번 비교한 셈이고,
## "결과가 안 바뀐다"가 어느 정도 당연합니다. 리뷰어가 이 점을 짚습니다.
##
## 이 스크립트는 다른 각도에서 직접 확인합니다.
##   측정 표류가 연관을 만들어 내고 있다면, 표준편차당 IRR 이 파동에 따라
##   체계적으로 달라져야 합니다. 파동별로 따로 추정해서 흩어지는지 봅니다.
##   흩어지지 않으면, 표류가 있더라도 추정량을 움직이지 못한다는 뜻입니다.
##
## ── 무엇을 계산하는가 ────────────────────────────────────────────────────
##   1) 전이별 · 파동별 IRR (구간 시작 파동 k = 1..K 로 층화)
##   2) Cochran's Q 와 I^2  — 파동 간 이질성의 크기
##   3) IC x 파동 상호작용의 Wald 검정 (참가자 단위 강건 공분산 사용)
##      ★ LR 검정을 쓰지 않는 이유: 클러스터 강건 SE 를 쓰면 우도비 검정이
##        타당하지 않습니다. 강건 공분산에 대한 Wald 통계량이 맞습니다.
##
## 출력:  ST_WaveStability.xlsx / .csv  ·  SupplFig2_WaveStability.pdf/png
## 실행:  source("KF_ST_WaveStability_v260811.R", encoding = "UTF-8")
###############################################################################

WAVESTAB_VERSION <- "v260811"
message("\n=== KF_ST_WaveStability ", WAVESTAB_VERSION, " ===")

###############################################################################
## 0. 헬퍼 로드 (파일명이 정확히 일치하지 않아도 같은 계열의 최신본을 찾습니다)
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
###############################################################################
## 1. 설정
###############################################################################
MIN_EV_WAVE  <- 10   # 한 파동에서 이만큼 사건이 있어야 그 파동 추정치를 보고
MIN_EV_TOTAL <- 40   # 전체 사건이 이만큼은 되어야 전이를 표에 올림
MIN_WAVES    <- 3    # 추정 가능한 파동이 이만큼은 되어야 이질성을 계산
FIG_H        <- 170  # 그림 높이(mm)

###############################################################################
## 2. 보조 함수
###############################################################################
## 참가자 단위 강건 공분산 (Wald 검정에 행렬 전체가 필요합니다)
.rvcov <- function(m, cluster) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    v <- try(sandwich::vcovCL(m, cluster = cluster, type = "HC0"), silent = TRUE)
    if (!inherits(v, "try-error") && all(is.finite(v))) return(v)
  }
  stats::vcov(m)
}
## 다변량 Wald 검정
.wald <- function(b, V, sel) {
  ## ★ 이름으로 선택합니다. 위치 인덱스를 쓰면 별칭(aliased, NA) 계수가 있을 때
  ##   축소된 강건 공분산 행렬의 범위를 벗어나 subscript out of bounds 가 납니다.
  sel <- intersect(intersect(sel, names(b)), rownames(V))
  if (!length(sel)) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  bb <- b[sel]; VV <- V[sel, sel, drop = FALSE]
  ok <- is.finite(bb) & is.finite(diag(VV))
  if (!any(ok)) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  bb <- bb[ok]; VV <- VV[ok, ok, drop = FALSE]
  Vi <- try(solve(VV), silent = TRUE)
  if (inherits(Vi, "try-error")) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  st <- as.numeric(t(bb) %*% Vi %*% bb)
  c(stat = st, df = length(bb), p = stats::pchisq(st, length(bb), lower.tail = FALSE))
}
## Cochran's Q / I^2 (파동별 추정치의 이질성)
.het <- function(b, se) {
  ok <- is.finite(b) & is.finite(se) & se > 0
  b <- b[ok]; se <- se[ok]
  k <- length(b)
  if (k < 2) return(c(k = k, Q = NA_real_, df = NA_real_, p = NA_real_, I2 = NA_real_))
  w  <- 1 / se^2
  mu <- sum(w * b) / sum(w)
  Q  <- sum(w * (b - mu)^2)
  df <- k - 1
  c(k = k, Q = Q, df = df,
    p  = stats::pchisq(Q, df, lower.tail = FALSE),
    I2 = max(0, (Q - df) / Q) * 100)
}
fmt_ci <- function(e, l, u)
  ifelse(is.finite(e), sprintf("%.2f (%.2f\u2013%.2f)", e, l, u), "\u2013")
fmt_p <- function(p) ifelse(!is.finite(p), "\u2013",
                     ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

###############################################################################
## 3. 자료와 구간
###############################################################################
load_imputed("MAIN")
d <- prep_long(load_kfacs())
msg(sprintf("WaveStability: N=%d, person-waves=%d", length(unique(d$id)), nrow(d)))

SYS <- list(list(var = "state",        lab = "ADL disability, 3 states"),
            list(var = "frailty_3cat", lab = "Frailty phenotype"))

###############################################################################
## 4. 전이별 · 파동별 추정
###############################################################################
one_transition <- function(iv, from_lab, to_lab, sys_lab) {
  dsub <- iv[iv$from_lab == from_lab, , drop = FALSE]
  dsub$y <- as.integer(dsub$to_lab == to_lab)
  n_ev <- sum(dsub$y, na.rm = TRUE)
  if (n_ev < MIN_EV_TOTAL) return(NULL)

  ## 보정 공변량 중 실제로 변이가 있는 것만 (파동 내에서 상수가 되는 경우 대비)
  pick_rhs <- function(dd) {
    a <- intersect(CFG$ADJ, names(dd))
    a[vapply(a, function(x) length(unique(stats::na.omit(dd[[x]]))) > 1, logical(1))]
  }

  waves <- sort(unique(dsub$k))
  per <- do.call(rbind, lapply(waves, function(w) {
    dw <- dsub[dsub$k == w, , drop = FALSE]
    ev <- sum(dw$y, na.rm = TRUE)
    out <- data.frame(wave = w, events = ev, pyears = sum(dw$dur, na.rm = TRUE),
                      b = NA_real_, se = NA_real_, irr = NA_real_,
                      lo = NA_real_, hi = NA_real_, p = NA_real_)
    if (ev < MIN_EV_WAVE || length(unique(stats::na.omit(dw$gLIC_z))) < 5) return(out)
    rhs <- pick_rhs(dw)
    f <- stats::as.formula(paste0("y ~ gLIC_z",
           if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")) else "",
           " + offset(log(dur))"))
    m <- try(stats::glm(f, family = stats::poisson(), data = dw), silent = TRUE)
    if (inherits(m, "try-error") || !isTRUE(m$converged)) return(out)
    V <- .rvcov(m, dw$id)
    b <- stats::coef(m)["gLIC_z"]; s <- sqrt(V["gLIC_z", "gLIC_z"])
    if (!is.finite(b) || !is.finite(s) || s > 3 || abs(b) > 10) return(out)
    out$b <- b; out$se <- s; out$irr <- exp(b)
    out$lo <- exp(b - 1.96 * s); out$hi <- exp(b + 1.96 * s)
    out$p  <- 2 * stats::pnorm(-abs(b / s))
    out
  }))

  ## 전체(파동 통합) 추정 + 상호작용 검정
  rhs <- pick_rhs(dsub)
  f0 <- stats::as.formula(paste0("y ~ gLIC_z + factor(k)",
          if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")) else "",
          " + offset(log(dur))"))
  f1 <- stats::as.formula(paste0("y ~ gLIC_z * factor(k)",
          if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")) else "",
          " + offset(log(dur))"))
  m0 <- try(stats::glm(f0, family = stats::poisson(), data = dsub), silent = TRUE)
  m1 <- try(stats::glm(f1, family = stats::poisson(), data = dsub), silent = TRUE)
  pooled_irr <- NA_real_; pooled_lo <- NA_real_; pooled_hi <- NA_real_
  if (!inherits(m0, "try-error") && isTRUE(m0$converged)) {
    V0 <- .rvcov(m0, dsub$id)
    b0 <- stats::coef(m0)["gLIC_z"]; s0 <- sqrt(V0["gLIC_z", "gLIC_z"])
    if (is.finite(b0) && is.finite(s0) && s0 <= 3) {
      pooled_irr <- exp(b0); pooled_lo <- exp(b0 - 1.96 * s0); pooled_hi <- exp(b0 + 1.96 * s0)
    }
  }
  wald <- c(stat = NA_real_, df = NA_real_, p = NA_real_)
  if (!inherits(m1, "try-error") && isTRUE(m1$converged)) {
    b1 <- stats::coef(m1); V1 <- .rvcov(m1, dsub$id)
    sel  <- grep("^gLIC_z:factor\\(k\\)", names(b1), value = TRUE)
    wald <- .wald(b1, V1, sel)
  }
  h <- .het(per$b, per$se)

  list(per = cbind(system = sys_lab,
                   transition = paste(from_lab, "\u2192", to_lab), per),
       sum = data.frame(system = sys_lab,
                        transition = paste(from_lab, "\u2192", to_lab),
                        events = n_ev,
                        pooled = fmt_ci(pooled_irr, pooled_lo, pooled_hi),
                        n_waves = h["k"], Q = h["Q"], I2 = h["I2"],
                        p_het = h["p"], p_int = wald["p"],
                        stringsAsFactors = FALSE))
}

PER <- list(); SUM <- list()
for (S in SYS) {
  iv <- build_intervals(d, S$var, "rolling")
  if (is.null(iv)) next
  msg(sprintf("[WaveStab] %s: 구간 %d개, 파동 %s",
              S$lab, nrow(iv), paste(sort(unique(iv$k)), collapse = ",")))
  tl <- default_transitions(S$var, iv)
  for (tr in tl) {
    r <- one_transition(iv, tr[1], tr[2], S$lab)
    if (is.null(r)) next
    PER[[length(PER) + 1]] <- r$per
    SUM[[length(SUM) + 1]] <- r$sum
  }
}
PER <- do.call(rbind, PER)
SUM <- do.call(rbind, SUM)
rownames(SUM) <- NULL
if (is.null(SUM) || !nrow(SUM))
  stop("추정된 전이가 없습니다. MIN_EV_TOTAL 을 낮춰 보십시오.")

###############################################################################
## 5. 표
###############################################################################
WAVE_LAB <- function(w) sprintf("W%d\u2192W%d", w, w + 1)
wide <- reshape(PER[, c("system", "transition", "wave", "irr", "lo", "hi")],
                idvar = c("system", "transition"), timevar = "wave",
                direction = "wide")
wv <- sort(unique(PER$wave))
TAB <- SUM[, c("system", "transition", "events", "pooled")]
names(TAB) <- c("System", "Transition", "No. of events", "All waves pooled (95% CI)")
for (w in wv) {
  cn <- paste0(c("irr.", "lo.", "hi."), w)
  if (!all(cn %in% names(wide))) next
  k <- match(paste(TAB$System, TAB$Transition),
             paste(wide$system, wide$transition))
  TAB[[WAVE_LAB(w)]] <- fmt_ci(wide[[cn[1]]][k], wide[[cn[2]]][k], wide[[cn[3]]][k])
}
TAB[["I\u00b2, %"]]            <- ifelse(is.finite(SUM$I2), sprintf("%.0f", SUM$I2), "\u2013")
TAB[["P, heterogeneity"]]      <- fmt_p(SUM$p_het)
TAB[["P, IC \u00d7 wave"]]     <- fmt_p(SUM$p_int)
## System 은 첫 행에만
TAB$System[duplicated(TAB$System)] <- ""

n_tr    <- nrow(SUM)
n_int   <- sum(SUM$p_int < 0.05, na.rm = TRUE)
n_het   <- sum(SUM$p_het < 0.05, na.rm = TRUE)
i2_med  <- stats::median(SUM$I2, na.rm = TRUE)
i2_max  <- max(SUM$I2, na.rm = TRUE)
bonf    <- 0.05 / n_tr
n_bonf  <- sum(SUM$p_int < bonf, na.rm = TRUE)

save_table(TAB, "ST_WaveStability",
  title = paste("Supplementary Data | Stability of the effect of intrinsic capacity",
                "across study waves"),
  footnotes = c(
    paste("Incidence rate ratios per +1 s.d. of intrinsic capacity from transition-specific",
          "Poisson models with a log(person-time) offset and standard errors clustered by",
          "participant, adjusted for", adj_phrase(),
          "and fitted separately within each between-visit interval."),
    paste("Wave W\u2192W+1 denotes the interval that begins at that visit; the exposure is",
          "the intrinsic capacity measured at the start of the interval, so each column uses",
          "a different measurement occasion of the same instrument."),
    paste("This analysis addresses the departure from full scalar invariance reported in",
          "Supplementary Table 1. If drift in indicator intercepts across waves were",
          "generating the association, the effect estimated at one measurement occasion",
          "would differ systematically from the effect estimated at another. Stability across",
          "waves therefore bounds the practical consequence of the departure, independently",
          "of the invariance-based score used in Supplementary Table 3."),
    paste("I\u00b2 is the proportion of variation across waves not attributable to sampling",
          "error, computed from the wave-specific log rate ratios and their robust standard",
          "errors. P for heterogeneity is Cochran's Q. P for IC \u00d7 wave is a Wald test of the",
          "interaction terms in a single model fitted to all waves, evaluated against the",
          "cluster-robust covariance matrix; a likelihood-ratio test is not valid under",
          "cluster-robust estimation."),
    sprintf(paste("Estimates are shown for waves contributing at least %d events and for",
                  "transitions with at least %d events in total. Across %d transitions the",
                  "median I\u00b2 was %.0f%% (maximum %.0f%%); the interaction reached P<0.05 for",
                  "%d transition%s and P<%.4f (Bonferroni) for %d."),
            MIN_EV_WAVE, MIN_EV_TOTAL, n_tr, i2_med, i2_max,
            n_int, if (n_int == 1) "" else "s", bonf, n_bonf),
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))

save_vals(PER, "ST_WaveStability_values.csv", FIG_DIR)

###############################################################################
## 6. 그림 — 전이별 파동 IRR
###############################################################################
G <- PER[is.finite(PER$irr), ]
if (nrow(G)) {
  G$wave_lab <- factor(WAVE_LAB(G$wave), levels = WAVE_LAB(sort(unique(G$wave))))
  ## 전이 이름이 길어 패널 제목이 겹치므로 줄바꿈
  G$panel <- factor(G$transition, levels = unique(PER$transition))
  P_ALL <- ggplot(G, aes(x = irr, y = wave_lab)) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22, linewidth = 0.42,
                   colour = "#1C7293", na.rm = TRUE) +
    geom_point(size = 1.3, colour = "#1C7293") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c("0.25", "0.50", "1.00", "2.00", "4.00")) +
    facet_wrap(~ panel, ncol = 3, scales = "free_y") +
    labs(x = "Incidence rate ratio per +1 s.d. of intrinsic capacity (log scale)",
         y = "Interval beginning at wave") +
    theme_na() +
    theme(panel.spacing.x = unit(3.5, "mm"),
          panel.spacing.y = unit(2.6, "mm"),
          strip.text = element_text(face = "plain", size = 6.6))

  CAP <- paste0(
    "Incidence rate ratios per +1 s.d. of intrinsic capacity, estimated separately within ",
    "each between-visit interval. The exposure is the capacity measured at the start of the ",
    "interval, so each row of a panel uses a different measurement occasion of the same ",
    "instrument. Horizontal bars are 95% confidence intervals from standard errors clustered ",
    "by participant. Waves contributing fewer than ", MIN_EV_WAVE, " events are not shown. ",
    "If drift in indicator intercepts across waves were generating the association, estimates ",
    "would differ systematically between rows; across ", n_tr, " transitions the median I\u00b2 was ",
    sprintf("%.0f%%", i2_med), " and the interaction between capacity and wave reached P<0.05 for ",
    n_int, " transition", if (n_int == 1) "" else "s", ".")

  .af <- names(formals(annot_na))
  if (all(c("title", "caption") %in% .af)) {
    P_ALL <- P_ALL + annot_na(width_mm = NA_WIDTH$double,
      title = "Stability of the effect of intrinsic capacity across study waves",
      caption = CAP)
  } else {
    P_ALL <- P_ALL + annot_na(width_mm = NA_WIDTH$double)
    writeLines(c("Supplementary Figure | Stability of the effect of intrinsic capacity across study waves",
                 "", CAP),
               file.path(FIG_DIR, "SupplFig2_WaveStability_caption.txt"), useBytes = TRUE)
  }
  save_na(P_ALL, "SupplFig2_WaveStability", width_mm = NA_WIDTH$double,
          height_mm = FIG_H, dir = FIG_DIR)
}

###############################################################################
## 7. 판정
###############################################################################
cat("\n=== 파동 안정성 판정 =========================================\n")
cat(sprintf("  전이 %d개 · 파동별 추정치 %d개\n", n_tr, sum(is.finite(PER$irr))))
cat(sprintf("  I\u00b2 중앙값 %.0f%% (최대 %.0f%%)\n", i2_med, i2_max))
cat(sprintf("  이질성 P<0.05        : %d개\n", n_het))
cat(sprintf("  IC \u00d7 파동 상호작용 P<0.05 : %d개\n", n_int))
cat(sprintf("  Bonferroni (P<%.4f) 통과 : %d개\n", bonf, n_bonf))
cat("--------------------------------------------------------------\n")
if (is.finite(i2_med) && i2_med < 50 && n_bonf == 0) {
  cat("  판정: 파동 간 안정적입니다.\n")
  cat("  -> 완전 스칼라 불변성 미지지의 실질적 영향이 크지 않다는 실증 근거입니다.\n")
  cat("     Discussion 에서 ST3(부분 스칼라 점수) 보다 이 결과를 앞세우십시오.\n")
} else {
  cat("  판정: 파동 간 이질성이 있습니다. 아래를 확인하십시오.\n")
  bad <- SUM[order(SUM$p_int), c("transition", "events", "I2", "p_int")]
  print(utils::head(bad[is.finite(bad$p_int), ], 5), row.names = FALSE)
  cat("  -> 사건이 적은 전이에서 나온 이질성이면 우연으로 해석해도 됩니다.\n")
  cat("     사건이 많은 전이에서 나왔다면 측정 표류를 진지하게 다루어야 합니다.\n")
}
cat("==============================================================\n")
cat("\n=== 완료 ===  표:", file.path(FIG_DIR, "ST_WaveStability.xlsx"),
    "\n              그림:", file.path(FIG_DIR, "SupplFig2_WaveStability.pdf"), "\n")
