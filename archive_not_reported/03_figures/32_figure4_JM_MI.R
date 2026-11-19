###############################################################################
## 32_figure2_alternative.R  (v260813 — 미팅 반영판)
## Figure 2 본문 그림 — 4패널 동적예측 + 3가지 추가 (2026-08-12 미팅)
##
##  ① 각 랜드마크 패널에 '현재 위험집합 n' 표기 (current sample size)
##  ② 각 패널의 랜드마크 시점에 '현재 IC 수준(저/중/고)과 직전 2년 변화(Δ)'를
##     주석 한 줄로 명시 (별도 궤적 패널은 가상 프로파일의 구성 산물이라 제외)
##  ③ 새 패널 e: '현재 수준 vs 직전 변화' 직접 비교 —
##     다음 구간 사망에 대한 이산시간 Poisson 모형에서 두 항을 따로/동시에.
##     메시지: 변한 것보다 현재의 수준(point)이 결과를 가른다.
##
##  ── 실행 ────────────────────────────────────────────────────────────────
##  ★ 조인트모형을 다시 적합하지 않습니다. 31 이 저장한
##    F2_dynamic_predictions.csv / F2_observed_IC.csv 를 읽습니다.
##  ★ 단, ①과 ③에 원자료가 필요하므로 데이터를 로드합니다(약 1분).
###############################################################################

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

F_PRED <- file.path(FIG_DIR, "F2_dynamic_predictions_MI.csv")
F_OBS  <- file.path(FIG_DIR, "F2_observed_IC_MI.csv")
if (!file.exists(F_PRED))
  stop("F2_dynamic_predictions.csv 가 없습니다.\n  -> 31_figure2_dynamic_prediction_MI.R 을 먼저 한 번 실행하십시오.\n  찾은 곳: ", F_PRED)
PRED <- utils::read.csv(F_PRED, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
OBS  <- if (file.exists(F_OBS))
  utils::read.csv(F_OBS, stringsAsFactors = FALSE, fileEncoding = "UTF-8") else NULL

PRED <- PRED[is.finite(PRED$time) & is.finite(PRED$surv), ]
if (!nrow(PRED)) stop("예측자료가 비어 있습니다.")

## ── 라벨 정리 (기존 로직 유지) ──────────────────────────────────────────
strip_tag <- function(x) trimws(sub("^[a-z]\\s+", "", as.character(x)))
ord <- unique(as.character(PRED$profile))
ord <- ord[order(sub("^([a-z])\\s+.*$", "\\1", ord))]   # a -> b -> c
have <- strip_tag(ord)
PRED$prof <- factor(strip_tag(PRED$profile), levels = have)
SYNTH <- any(grepl("IC \\(|Median IC", have))
if (!SYNTH) message("[Fig2] 실제 참가자 프로파일 감지 -> 캡션을 그에 맞게 씁니다.")

LM     <- sort(unique(PRED$landmark))
GAP    <- if (exists("CFG") && is.finite(CFG$WAVE_GAP)) CFG$WAVE_GAP else 2
lab_lm <- function(l) sprintf("%s  W%d (%g y)", letters[match(l, LM)],
                              round(l / GAP) + 1, round(l, 1))
PRED$lmf <- factor(lab_lm(PRED$landmark), levels = lab_lm(LM))

COL_PROF <- stats::setNames(
  c(SEV_COL[["bad"]], SEV_COL[["mid"]], SEV_COL[["good"]])[seq_along(have)], have)

XMAX  <- max(PRED$time, na.rm = TRUE)
XBRK  <- seq(0, floor(XMAX), by = if (XMAX > 12) 4 else 2)
S_LO  <- max(0, floor(min(PRED$surv_lo, na.rm = TRUE) * 50) / 50)
S_BRK <- seq(ceiling(S_LO * 20) / 20, 1, by = 0.05)
msg(sprintf("Fig2: 랜드마크 %s | x 0-%.1f년 | 생존축 %.2f-1.00",
            paste(LM, collapse = ", "), XMAX, S_LO))

###############################################################################
## [1] 원자료 로드 — 위험집합 n(①)과 현재 vs 변화 모형(③)에 필요
###############################################################################
## ── v260904 MI 판: MI01 완성 데이터 + within-sex 척도 (14_ 와 동일 정의) ──
load_imputed("MI01", stem = get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20"))
d <- prep_long(load_kfacs())
local({
  w1 <- d[!duplicated(d$id), ]; z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE); z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE)) }
  d$gLIC_z <<- z; d$gLIC_tert_W1cut <<- factor(tt, levels = CFG$TERT_LABELS)
})
msg("[MI 판] MI01 + within-sex 척도 적용")

## 개인별 추적 종료(t)와 사망 여부
sv <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p) {
  n <- nrow(p)
  data.frame(id = p$id[1], t = p$followup_years[n],
             ev = as.integer(p$death_event[n] == 1), stringsAsFactors = FALSE)
}))
sv <- sv[is.finite(sv$t), ]

## ── ① 랜드마크별 위험집합 (그 시점에 생존해 추적 중인 인원) ─────────────
NR <- data.frame(landmark = LM,
                 n_risk  = vapply(LM, function(l) sum(sv$t > l), integer(1)),
                 n_death = vapply(LM, function(l) sum(sv$ev == 1 & sv$t <= l), integer(1)))
NR$lmf <- factor(lab_lm(NR$landmark), levels = lab_lm(LM))
msg(paste0("[Fig2] 위험집합: ",
           paste(sprintf("%gy n=%s", NR$landmark,
                         format(NR$n_risk, big.mark = ",")), collapse = " | ")))

###############################################################################
## [2] ③ 현재 수준 vs 직전 2년 변화 — 다음 구간 사망의 이산시간 모형
##   · 구간: rolling landmark (주 분석과 동일 구조)
##   · Δ = 구간 시작 시점 IC - 직전 wave IC  -> 첫 구간(직전 측정 없음)은 제외
##   · 출발상태(from)와 CORE 공변량 보정, 참가자 군집 강건 SE
##   · '따로' 모형과 '동시' 모형을 같은 부분표본에서 적합 (직접 비교 가능)
###############################################################################
iv <- build_intervals(d, "state", "rolling")
d2 <- d[order(d$id, d$time), ]
d2$k <- stats::ave(seq_len(nrow(d2)), d2$id, FUN = seq_along)
iv$z_prev <- d2$gLIC_z[match(paste(iv$id, iv$k - 1), paste(d2$id, d2$k))]
iv$dz     <- iv$gLIC_z - iv$z_prev

CC <- iv[is.finite(iv$dz) & is.finite(iv$gLIC_z) & iv$from_lab != "Death", ]
CC$y      <- as.integer(CC$to_lab == "Death")
CC$from_f <- droplevels(factor(CC$from_lab))
N_INT <- nrow(CC); N_DTH <- sum(CC$y)
msg(sprintf("[Fig2] 현재vs변화 모형: 구간 %s개 (직전 측정 보유), 사망 %d건",
            format(N_INT, big.mark = ","), N_DTH))

.rvcov <- function(m, cl) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    v <- try(sandwich::vcovCL(m, cluster = cl, type = "HC0"), silent = TRUE)
    if (!inherits(v, "try-error") && all(is.finite(v))) return(v)
  }
  stats::vcov(m)
}
fit_terms <- function(terms) {
  adj <- setdiff(CFG$ADJ, "sex_f")   # within-sex 척도: sex 는 표준화에 흡수
  adj <- adj[vapply(adj, function(v)
    v %in% names(CC) && length(unique(stats::na.omit(CC[[v]]))) > 1, logical(1))]
  f <- stats::as.formula(paste0("y ~ ", paste(c(terms, "from_f", adj), collapse = " + "),
                                " + offset(log(dur))"))
  m <- stats::glm(f, family = stats::poisson(), data = CC)
  V <- .rvcov(m, CC$id)
  out <- lapply(terms, function(tm) {
    b <- stats::coef(m)[tm]; s <- sqrt(V[tm, tm])
    data.frame(term = tm, irr = exp(b), lo = exp(b - 1.96 * s), hi = exp(b + 1.96 * s),
               p = 2 * stats::pnorm(-abs(b / s)), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
r_cur  <- fit_terms("gLIC_z")                 # 현재 수준만
r_dz   <- fit_terms("dz")                     # 변화만
r_both <- fit_terms(c("gLIC_z", "dz"))        # 동시
r_eq   <- fit_terms(c("gLIC_z", "z_prev"))    # 등가 재모수화: 현재 + 과거 수준

FOR <- rbind(
  cbind(r_cur,             model = "separate", what = "current"),
  cbind(r_both[r_both$term == "gLIC_z", ], model = "mutual", what = "current"),
  cbind(r_dz,              model = "separate", what = "change"),
  cbind(r_both[r_both$term == "dz", ],     model = "mutual", what = "change"))
FOR$label <- factor(c("Current level, separate model",
                      "Current level, adjusted for change",
                      "2-y change, separate model",
                      "2-y change, adjusted for current level"),
                    levels = rev(c("Current level, separate model",
                                   "Current level, adjusted for change",
                                   "2-y change, separate model",
                                   "2-y change, adjusted for current level")))
EQ <- cbind(rbind(r_eq), model = "levels", what = c("current", "previous"),
            label = c("Current level (level model)", "Previous level (level model)"))
save_vals(cbind(rbind(FOR, EQ), n_intervals = N_INT, n_deaths = N_DTH),
          "F2_current_vs_change_MI.csv", FIG_DIR)
print(FOR[, c("label", "irr", "lo", "hi", "p")], row.names = FALSE, digits = 3)

B_CUR <- FOR[FOR$what == "current" & FOR$model == "mutual", ]
B_DZ  <- FOR[FOR$what == "change"  & FOR$model == "mutual", ]
A_CUR <- FOR[FOR$what == "current" & FOR$model == "separate", ]
A_DZ  <- FOR[FOR$what == "change"  & FOR$model == "separate", ]
E_CUR <- r_eq[r_eq$term == "gLIC_z", ]; E_PRV <- r_eq[r_eq$term == "z_prev", ]

###############################################################################
## [3] 그림 — 상단: 생존 4패널(+n) / 중단: IC 궤적 띠(+Δ) / 하단: 비교 패널 f
###############################################################################
## ── 랜드마크 시점의 현재 IC(저/중/고)와 직전 2년 변화 — a–d 주석용 ──────
ICANN <- NULL
if (!is.null(OBS) && all(c("time", "gLIC_z", "landmark", "profile") %in% names(OBS))) {
  IC <- OBS[OBS$time <= OBS$landmark + 1e-9, ]
  IC$prof <- factor(strip_tag(IC$profile), levels = have)
  CUR <- do.call(rbind, lapply(split(IC, list(IC$landmark, IC$prof), drop = TRUE),
    function(z) { z <- z[order(z$time), ]
      data.frame(landmark = z$landmark[1], prof = z$prof[1],
                 z = z$gLIC_z[nrow(z)],
                 dz = if (nrow(z) >= 2) z$gLIC_z[nrow(z)] - z$gLIC_z[nrow(z) - 1] else NA_real_) }))
  ICANN <- do.call(rbind, lapply(split(CUR, CUR$landmark), function(z) {
    z <- z[order(match(z$prof, have)), ]
    data.frame(landmark = z$landmark[1],
               lab = sprintf("IC %s", paste(sprintf("%.1f", z$z), collapse = " / ")),
               stringsAsFactors = FALSE)
  }))
  ICANN$lmf <- factor(lab_lm(ICANN$landmark), levels = lab_lm(LM))
  DZ_PROF <- stats::median(CUR$dz, na.rm = TRUE)   # 프로파일 공통 기울기 (캡션용)
}
if (!exists("DZ_PROF")) DZ_PROF <- NA_real_

## ── 상단: 동적예측 (기존 구조 + 위험집합 n + IC 주석) ───────────────────
p_top <- ggplot(PRED, aes(time, surv, colour = prof, fill = prof)) +
  geom_ribbon(aes(ymin = surv_lo, ymax = surv_hi), colour = NA, alpha = 0.16) +
  geom_line(linewidth = 0.5) +
  geom_vline(aes(xintercept = landmark), data = unique(PRED[, c("landmark", "lmf")]),
             linetype = "22", linewidth = LINE_PT, colour = SEV_COL[["aux"]],
             inherit.aes = FALSE) +
  geom_text(data = NR, aes(x = 0.15, y = S_LO + 0.055 * (1 - S_LO),
                           label = sprintf("At risk n = %s", format(n_risk, big.mark = ","))),
            inherit.aes = FALSE, hjust = 0, vjust = 0, size = 2.25, colour = "black") +
  { if (!is.null(ICANN))
      geom_text(data = ICANN, aes(x = 0.15, y = S_LO + 0.005 * (1 - S_LO), label = lab),
                inherit.aes = FALSE, hjust = 0, vjust = 0, size = 2.0,
                lineheight = 0.95, colour = "grey25") } +
  facet_wrap(~ lmf, nrow = 1) +
  scale_colour_manual(values = COL_PROF, name = NULL, drop = FALSE) +
  scale_fill_manual(values = COL_PROF, name = NULL, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0, XMAX), breaks = XBRK,
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = S_BRK, labels = sprintf("%.2f", S_BRK),
                     expand = expansion(mult = 0.03)) +
  coord_cartesian(ylim = c(S_LO, 1)) +
  labs(x = "Years from baseline", y = "Predicted survival probability") +
  theme_na() +
  theme(strip.text      = element_text(face = "bold", hjust = 0),
        panel.spacing.x = unit(1.8, "mm"),
        legend.position = "top")

## ── 하단: 패널 e — 현재 수준 vs 직전 변화 ───────────────────────────────
FBRK <- c(0.25, 0.5, 1, 2)
p_bot <- ggplot(FOR, aes(irr, label)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = LINE_PT, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.45,
                 colour = "black") +
  geom_point(size = 1.6, colour = "black") +
  geom_text(aes(x = hi * 1.06, label = sprintf("%.2f (%.2f–%.2f)", irr, lo, hi)),
            hjust = 0, size = 2.1, colour = "black") +
  scale_x_continuous(trans = "log",
                     breaks = FBRK, labels = format(FBRK, drop0trailing = TRUE),
                     expand = expansion(mult = c(0.03, 0.34))) +
  labs(x = "IRR for death over the next interval, per +1 s.d. (log scale)",
       y = NULL, tag = "e") +
  theme_na() +
  theme(plot.tag = element_text(size = LAB_PT, face = "bold"),
        plot.tag.position = c(0, 1))

fig <- p_top / p_bot + patchwork::plot_layout(heights = c(2.35, 0.85))

## ── 각주 ────────────────────────────────────────────────────────────────
CAP <- paste0(
  "Panels a–d: each panel is a landmark. The joint model is updated with every ",
  "intrinsic-capacity measurement available up to the dashed line and then predicts ",
  "survival to the end of follow-up (", sprintf("%.1f", XMAX), " years); the prediction ",
  "window shortens and the credible intervals narrow as information accrues. ",
  "The number at risk is the number of participants alive and under follow-up at ",
  "the landmark (wave 1, n = ", format(NR$n_risk[1], big.mark = ","), "; wave ",
  length(LM), ", n = ", format(NR$n_risk[nrow(NR)], big.mark = ","), "). ",
  if (SYNTH) paste0(
    "The three curves are hypothetical participants who differ only in intrinsic capacity ",
    "(10th, 50th and 90th centile of the wave-1 distribution, thereafter following the ",
    "cohort-average trajectory), with age, sex and baseline comorbidity held at cohort ",
    "reference values, so the vertical distance between curves is attributable to intrinsic ",
    "capacity alone. ")
  else paste0(
    "The three curves are three real participants of the same sex, selected to lie near the ",
    "10th, 50th and 90th centile of wave-1 intrinsic capacity; because they are real people ",
    "they also differ in age and comorbidity, so the vertical distance between curves is not ",
    "attributable to intrinsic capacity alone. "),
  "Shaded bands are 95% credible intervals. The survival axis is truncated at ",
  sprintf("%.2f", S_LO), ". ",
  "The annotation in each panel gives the capacity supplied to the model at the landmark ",
  "(lowest / median / highest profile, wave-1 s.d. units); the three profiles decline in ",
  "parallel at the cohort mean rate (", sprintf("%+.2f", DZ_PROF), " s.d. per ", GAP,
  " years), so the separation of the survival curves is produced by the current level alone. ",
  "Panel e tests the comparison empirically, in the cohort rather than in profiles: ",
  "discrete-time Poisson models of death over the next between-visit interval (",
  format(N_INT, big.mark = ","), " intervals with a preceding measurement, ",
  N_DTH, " deaths), adjusted for origin state and for ", adj_phrase(setdiff(CFG$ADJ, "sex_f")),
  ", with standard errors clustered by participant; the change is expressed in the same ",
  "wave-1 s.d. units as the level. On its own, the change over the preceding ", GAP,
  " years did not predict death (IRR ", sprintf("%.2f, %.2f–%.2f", A_DZ$irr, A_DZ$lo, A_DZ$hi),
  "), whereas the current level did (IRR ",
  sprintf("%.2f, %.2f–%.2f", A_CUR$irr, A_CUR$lo, A_CUR$hi), "; ",
  sprintf("%.2f, %.2f–%.2f", B_CUR$irr, B_CUR$lo, B_CUR$hi),
  " with the change added). In the joint model the change term takes an IRR of ",
  sprintf("%.2f (%.2f–%.2f)", B_DZ$irr, B_DZ$lo, B_DZ$hi),
  "; because the current level is held fixed, a larger recent gain implies a lower level ",
  GAP, " years earlier, so a value above 1 reflects residual prognostic information in the ",
  "previous level, not a hazard of improvement — refitting with the current and previous ",
  "levels in place of the change gives IRRs of ",
  sprintf("%.2f (%.2f–%.2f)", E_CUR$irr, E_CUR$lo, E_CUR$hi), " and ",
  sprintf("%.2f (%.2f–%.2f)", E_PRV$irr, E_PRV$lo, E_PRV$hi),
  ". Prognosis therefore follows the level of intrinsic capacity, current and recent, ",
  "not the direction in which it has been moving. ",
  "Colours are those used for intrinsic-capacity tertiles elsewhere: orange, lowest ",
  "capacity; grey, intermediate; blue, highest. ",
  "The joint model uses the current value of intrinsic capacity only; it cannot ",
  "extrapolate beyond the observed follow-up.")

fig <- fig + annot_na(
  width_mm = NA_WIDTH$double, colour_note = FALSE,
  title = "Current intrinsic capacity, not its recent change, drives predicted survival",
  caption = CAP)

save_na(fig, "Figure4_JM_dynamic_predictions_MI", width_mm = NA_WIDTH$double,
        height_mm = 118, dir = FIG_DIR)

cat("\n=== Figure 2 (미팅 반영판) 완료 ===\n")
cat("  본문 그림 : Figure4_JM_dynamic_predictions_MI.pdf / .tiff\n")
cat("  각주      : Figure4_JM_dynamic_predictions_MI_caption.txt\n")
cat("  근거수치  : F2_current_vs_change.csv (패널 e), F2_dynamic_predictions.csv\n")
cat(sprintf("  핵심수치  : 변화 단독 %.2f (%.2f-%.2f) | 현재 단독 %.2f (%.2f-%.2f)\n",
            A_DZ$irr, A_DZ$lo, A_DZ$hi, A_CUR$irr, A_CUR$lo, A_CUR$hi))
cat(sprintf("              동시: 현재 %.2f, 변화 %.2f | 수준 재모수화: 현재 %.2f, 과거 %.2f\n",
            B_CUR$irr, B_DZ$irr, E_CUR$irr, E_PRV$irr))
