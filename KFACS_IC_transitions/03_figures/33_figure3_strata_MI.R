###############################################################################
## 33_figure3_strata_MI.R   (Phase 4 — Figure 3: robust·pre-frail 층 내 within-sex IC 삼분위별 5년 위험)
##                                                                        v260904
## 입력 : FIG_DIR/MI_Table3_perset.csv (14_ 산출). 없으면 프로젝트 내 탐색 -> 없으면 직접 계산(10분).
## 패널 : a 5년 사망, b 5년 복합 악화(첫 ADL 악화 또는 사망). x = within-sex 삼분위,
##        점 = 층(Robust, Pre-frail), 수직선 = 95% CI (logit 척도 Rubin 풀링).
##        frail 층은 최저 삼분위 한 칸(n≈240)만 성립하므로 점선 참조선으로만 표시.
## 제목 : "within the robust and pre-frail strata" — "within every frailty stratum" 문구 삭제.
## 산출 : Figure3_Strata_MI.pdf/.tiff, F3_strata_risk.csv
## 실행 : setwd("C:/Users/user/OneDrive/바탕 화면/Nat_Aging"); source("R/33_figure3_strata_MI.R", encoding="UTF-8")
###############################################################################
message("\n=== 33_figure3_strata_MI v260904 ===")

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
.need(c("ggplot2", "patchwork"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

###############################################################################
## 1. 세트별 결과 읽기 · 풀링
###############################################################################
F_IN <- file.path(FIG_DIR, "MI_Table3_perset.csv")
if (!file.exists(F_IN)) {                       # ① 프로젝트 어디든 있으면 찾아 쓴다
  hit <- list.files(PROJECT_ROOT, pattern = "^MI_Table3_perset\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(hit)) { F_IN <- hit[order(file.info(hit)$mtime, decreasing = TRUE)][1]; msg("세트별 결과 발견: ", F_IN) }
}
if (file.exists(F_IN)) {
  T3 <- utils::read.csv(F_IN, stringsAsFactors = FALSE)
} else {                                         # ② 없으면 14_ 의 Table 3 블록으로 직접 계산 (세트당 ~30초)
  msg("세트별 결과 파일이 없어 직접 계산합니다 (14_ Table 3 블록과 동일).")
  .need("survival"); suppressPackageStartupMessages(library(survival))
  MI_STEM <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함; TAU <- 5; MIN_EV_HR <- 5; MIN_N_HR <- 30
  MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))); SETS <- names(MI); rm(MI); invisible(gc())
  apply_scale <- function(d) {
    w1 <- d[!duplicated(d$id), ]; z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
    for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
      mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE); z[k] <- (d$gLIC[k] - mu) / sdv
      ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
      tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE)) }
    d$gLIC_z <- z; d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS); d }
  ADJ_WS <- setdiff(CFG$ADJ, "sex_f")
  mk_events <- function(dat) {
    dat <- dat[order(dat$id, dat$time), ]; sp <- split(seq_len(nrow(dat)), dat$id)
    do.call(rbind, lapply(sp, function(ix) {
      s <- dat$state[ix]; tm <- dat$time[ix]; ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
      n <- length(tm); if (n < 1 || !is.finite(s[1])) return(NULL)
      s0 <- s[1]; tmax <- max(tm)
      fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
      tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] > s0) else NA_real_
      kd <- which(is.finite(s) & s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
      data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax, t_wor = tw, t_dth = td) })) }
  mk_surv <- function(E, out) { inf <- function(x) ifelse(is.finite(x), x, Inf)
    if (out == "wor") { ft <- pmin(inf(E$t_wor), E$tmax, TAU); ev <- as.integer(is.finite(E$t_wor) & E$t_wor <= ft)
    } else { ft <- pmin(inf(E$t_dth), E$tmax, TAU); ev <- as.integer(is.finite(E$t_dth) & E$t_dth <= ft) }
    data.frame(id = E$id, s0 = E$s0, ft = pmax(ft, 1e-6), ev = ev) }
  risk5 <- function(ft, ev) { if (!any(ev > 0)) return(c(p = 0, se = NA))
    f <- survfit(Surv(ft, ev) ~ 1); s <- summary(f, times = TAU, extend = TRUE)
    p <- 1 - s$surv[1]; sd <- s$std.err[1] * s$surv[1]; p <- min(max(p, 1e-6), 1 - 1e-6); c(p = p, se = sd / (p * (1 - p))) }
  T3 <- list()
  for (i in seq_along(SETS)) { set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, length(SETS)))
    load_imputed(set, stem = MI_STEM); d <- apply_scale(prep_long(load_kfacs())); b <- get_baseline(d)
    b$tert <- d$gLIC_tert_W1cut[match(b$id, d$id)]
    b$fr <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail"))); b <- b[!is.na(b$tert) & !is.na(b$fr), ]
    E <- mk_events(d); ADJ <- intersect(ADJ_WS, names(b))
    for (out in c("wor", "dth")) {
      S <- merge(b[, c("id", "fr", "tert", ADJ)], mk_surv(E, out), by = "id"); S <- S[S$s0 < 4, ]
      S$cell <- factor(paste(S$fr, S$tert, sep = "|"), levels = as.vector(t(outer(levels(S$fr), CFG$TERT_LABELS, paste, sep = "|"))))
      S$cell <- droplevels(S$cell); S$cell <- relevel(S$cell, ref = "Robust|T3")
      cx <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ cell +", paste(ADJ, collapse = "+"))), data = S), silent = TRUE)
      for (cl in levels(S$cell)) { z <- S[S$cell == cl, ]; ne <- sum(z$ev == 1L); rk <- risk5(z$ft, z$ev)
        nm <- paste0("cell", cl); bhr <- NA; sehr <- NA
        if (cl != "Robust|T3" && !inherits(cx, "try-error") && nm %in% names(coef(cx)) && ne >= MIN_EV_HR && nrow(z) >= MIN_N_HR) {
          bhr <- coef(cx)[nm]; sehr <- sqrt(vcov(cx)[nm, nm]) }
        T3[[length(T3) + 1]] <- data.frame(set = set, scale = "withinsex", outcome = out, cell = cl, n = nrow(z), events = ne,
                                           p = rk[["p"]], se_logit = rk[["se"]], b_hr = bhr, se_hr = sehr, stringsAsFactors = FALSE) } } }
  T3 <- do.call(rbind, T3); save_vals(T3, "MI_Table3_perset.csv", FIG_DIR)
}
T3 <- T3[T3$scale == "withinsex" & T3$outcome %in% c("dth", "wor"), ]
M  <- length(unique(T3$set)); MIN_M_OK <- 0.8; MIN_N_RISK <- 10
msg(sprintf("세트 %d개, 셀 %d개", M, nrow(unique(T3[, c("outcome", "cell")]))))

pool <- function(x, se) {
  ok <- is.finite(x) & is.finite(se); if (sum(ok) < max(2, ceiling(MIN_M_OK * M))) return(c(NA, NA, NA, NA))
  r <- rubin_pool(x[ok], se[ok]); c(r["est"], r["lcl"], r["ucl"], r["p"])
}
T3$key <- paste(T3$outcome, T3$cell)
P <- do.call(rbind, lapply(split(T3, T3$key), function(z) {
  rr <- pool(stats::qlogis(z$p), z$se_logit); rh <- pool(z$b_hr, z$se_hr)
  data.frame(outcome = z$outcome[1], cell = z$cell[1],
             stratum = sub("\\|.*", "", z$cell[1]), tert = sub(".*\\|", "", z$cell[1]),
             n = round(median(z$n)), events = round(median(z$events)),
             risk = 100 * stats::plogis(rr[1]), lo = 100 * stats::plogis(rr[2]), hi = 100 * stats::plogis(rr[3]),
             hr = exp(rh[1]), hr_lo = exp(rh[2]), hr_hi = exp(rh[3]), hr_p = rh[4], stringsAsFactors = FALSE)
}))
rownames(P) <- NULL
P$risk[P$n < MIN_N_RISK] <- NA
P$outcome_lab <- factor(c(dth = "Death", wor = "Composite worsening")[P$outcome],
                        levels = c("Death", "Composite worsening"))
P$stratum <- factor(P$stratum, levels = c("Robust", "Pre-frail", "Frail"))
P$tert_lab <- factor(LAB_TERTILE[P$tert], levels = unname(LAB_TERTILE))
save_vals(P, "F3_strata_risk.csv", FIG_DIR)

MAIN <- P[P$stratum %in% c("Robust", "Pre-frail") & is.finite(P$risk), ]
REF  <- P[P$stratum == "Frail" & P$tert == "T1" & is.finite(P$risk), ]   # frail 층 참조선

###############################################################################
## 2. 그림 — 기존 Fig 3 디자인 (a,b 히트맵 9칸 + c 사망 forest) 을 MI 값으로 재현
###############################################################################
STR_LEV <- c("Robust", "Pre-frail", "Frail")
P$stratum <- factor(P$stratum, levels = STR_LEV)
P$tert    <- factor(P$tert, levels = CFG$TERT_LABELS)
P$lab_risk <- ifelse(is.finite(P$risk), sprintf("%.1f", P$risk), "NE")
P$lab_n    <- ifelse(is.finite(P$risk), sprintf("%d/%d", P$events, P$n), sprintf("n = %d", P$n))
P$small    <- P$n < 20
HEAT <- c(dth = SEV_COL[["death"]], wor = SEV_COL[["bad"]])

heat <- function(o, ttl) {
  z <- P[P$outcome == o, ]
  mx <- max(z$risk, na.rm = TRUE)
  ggplot(z, aes(tert, stratum)) +
    geom_tile(aes(fill = risk), colour = "white", linewidth = 1.2) +
    geom_text(aes(label = lab_risk, colour = small), size = (BASE_PT) / ggplot2::.pt, fontface = "bold", vjust = -0.15) +
    geom_text(aes(label = lab_n, colour = small), size = (BASE_PT - 3) / ggplot2::.pt, vjust = 1.6) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "grey45"), guide = "none") +
    scale_fill_gradient(low = "#FBF3EE", high = HEAT[[o]], limits = c(0, mx), na.value = "grey92",
                        name = "%", breaks = scales::pretty_breaks(4)) +
    scale_y_discrete(limits = rev(STR_LEV)) +
    scale_x_discrete(labels = c(T1 = "T1\nlowest", T2 = "T2", T3 = "T3\nhighest")) +
    coord_equal() +
    labs(title = ttl, x = "Intrinsic-capacity tertile (within sex)", y = NULL) +
    theme_na() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          legend.position = "right", legend.key.height = unit(5, "mm"), legend.key.width = unit(3, "mm"),
          plot.title = element_text(size = BASE_PT, face = "plain"))
}
pa <- heat("dth", "Five-year mortality"); pb <- heat("wor", "Five-year composite worsening")

## c: 사망 9칸 forest, 층별 그룹, robust T1 점선
zc <- P[P$outcome == "dth", ]
zc <- zc[order(zc$stratum, zc$tert), ]
zc$y <- (nrow(zc) + 3) - (seq_len(nrow(zc)) + as.integer(zc$stratum) - 1)   # 위에서 아래로 Robust T1..T3, 빈 칸, Pre-frail ..., Frail ...
zc$tert_lab <- factor(LAB_TERTILE[as.character(zc$tert)], levels = unname(LAB_TERTILE))
zc$lab <- ifelse(is.finite(zc$risk), sprintf("%.1f%%", zc$risk), sprintf("not estimable (n = %d)", zc$n))
ref <- zc$risk[zc$stratum == "Robust" & zc$tert == "T1"]
grp <- aggregate(y ~ stratum, zc, mean)
xmax <- max(zc$hi, na.rm = TRUE) * 1.15
pc <- ggplot(zc, aes(x = risk, y = y, colour = tert_lab)) +
  geom_vline(xintercept = ref, linetype = "22", colour = "grey35", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(shape = small), size = 2.3, fill = "white", na.rm = TRUE) +
  geom_text(aes(x = ifelse(is.finite(hi), hi, 0.5), label = lab), hjust = -0.15, size = (BASE_PT - 1) / ggplot2::.pt, colour = "grey20") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  scale_colour_manual(values = PAL_TERTILE_LAB, name = "Intrinsic-capacity tertile (within sex)") +
  scale_y_continuous(breaks = zc$y, labels = as.character(zc$tert), expand = expansion(add = 0.6),
                     sec.axis = dup_axis(breaks = grp$y, labels = as.character(grp$stratum), name = NULL)) +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Five-year mortality in all nine cells, grouped by frailty stratum", x = "Five-year mortality (%)", y = NULL) +
  theme_na() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.text.y.right = element_text(face = "bold", hjust = 0), axis.ticks.y = element_blank(),
        plot.title = element_text(size = BASE_PT, face = "plain"))

fig <- (pa | pb) / pc + plot_layout(heights = c(1.1, 0.75))

cell_txt <- function(o) { z <- P[P$outcome == o & is.finite(P$risk), ]
  paste(sprintf("%s %s %.1f%% (%.1f-%.1f)", z$stratum, z$tert, z$risk, z$lo, z$hi), collapse = "; ") }
CAP <- paste0(
  "a,b, Five-year risk of death from any cause (a) and of composite worsening, the first transition to a worse ADL ",
  "disability state or death (b), by frailty phenotype at wave 1 and within-sex tertile of intrinsic capacity at wave 1 ",
  "(n = 3,011; 500 deaths). Risks are one minus the Kaplan-Meier estimate at five years, estimated in each of ", M,
  " multiply imputed datasets and combined by Rubin's rules on the logit scale; each cell gives the risk (%) with the ",
  "number of events and participants beneath (medians across datasets); a and b use separate colour scales. Cells with ",
  "fewer than 20 participants (frail, middle tertile; n = 14) are shown for completeness only; NE, not estimable (frail, ",
  "highest tertile, n = 1). c, The same nine mortality cells as point estimates with 95% confidence intervals, grouped by ",
  "frailty stratum and ordered from the lowest to the highest tertile within each stratum; colour follows Figs. 1 and 2 ",
  "(orange, lowest; grey, middle; blue, highest) and the open circle marks the cell with 14 participants. The dashed line ",
  sprintf("marks the five-year mortality of robust participants in the lowest tertile (%.1f%%), which exceeds the point ", ref),
  "estimate for pre-frail participants in the middle tertile, although the intervals overlap. Numbers of participants are ",
  "100, 439 and 802 (robust), 664, 550 and 201 (pre-frail) and 240, 14 and 1 (frail) across the lowest, middle and highest ",
  "tertiles. Values: ", cell_txt("dth"), " (mortality); ", cell_txt("wor"), " (worsening). Adjusted hazard ratios are in Table 3.")
fig <- fig + annot_na(width_mm = NA_WIDTH$double, tag = TRUE,
  title = "Intrinsic capacity stratifies five-year risk within the robust and pre-frail strata", caption = CAP)
save_na(fig, "Figure3_Strata_MI", width_mm = NA_WIDTH$double, height_mm = 150, dir = FIG_DIR)

###############################################################################
## 3. 판정 출력
###############################################################################
cat("\n", strrep("=", 70), "\nFigure 3 — 층 × 삼분위 5년 위험 (%)\n", strrep("=", 70), "\n", sep = "")
show <- P[order(P$outcome, P$stratum, P$tert), c("outcome", "stratum", "tert", "n", "events", "risk", "lo", "hi", "hr", "hr_lo", "hr_hi")]
show[, 6:11] <- round(show[, 6:11], 2); print(show, row.names = FALSE)
r1 <- P$risk[P$outcome == "dth" & P$cell == "Robust|T1"]; p2 <- P$risk[P$outcome == "dth" & P$cell == "Pre-frail|T2"]
cat(sprintf("\n  핵심 문장 점검: robust 최저 삼분위 사망 %.1f%% vs pre-frail 중간 삼분위 %.1f%% -> %s\n",
            r1, p2, if (isTRUE(r1 > p2)) "성립" else "불성립 (본문 문장 수정 필요)"))
cat(strrep("=", 70), "\n=== 완료 ===\n")
