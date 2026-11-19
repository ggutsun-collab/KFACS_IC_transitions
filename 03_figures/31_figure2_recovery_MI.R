###############################################################################
## 31_figure2_recovery_MI.R   (Phase 4 — 새 본문 Figure 2: 회복 누적발생, within-sex 삼분위)
##                                                                        v260904
## 새 원고의 첫 결과가 recovery 이므로 본문 그림이 필요합니다. 기존 Fig 2(JM)는 ED 로.
##
## 패널 (wave-1 상태별, 사망 = 경쟁사건, Aalen-Johansen):
##   a  Pre-frail -> Robust           (wave-1 pre-frail)
##   b  Mild ADL disability -> Normal (wave-1 mild)
##   (frail 개선은 표에만: T3 n=0 이라 삼분위 곡선 불성립)
## 각 패널: 삼분위별 누적발생 곡선 + 95% CI 리본, m=20 세트에서 cloglog 척도 Rubin 풀링.
## 부속 표: 5년 누적발생(%), 원인특이 Cox HR (T1 기준, per +1 s.d.), 보정 = age, comorbid, SES.
##
## 산출 (FIG_DIR): Figure2_Recovery_MI.pdf/.tiff, F2_recovery_curves.csv, ST_Recovery_Fig2.xlsx
## 실행: setwd("C:/Users/user/OneDrive/바탕 화면/Nat_Aging"); source("R/31_figure2_recovery_MI.R", encoding="UTF-8")
##       (세트당 ~30초, 총 10분)
###############################################################################
message("\n=== 31_figure2_recovery_MI v260904 ===")

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
.need(c("survival", "ggplot2", "patchwork", "openxlsx"))
suppressPackageStartupMessages({ library(survival); library(ggplot2); library(patchwork) })

###############################################################################
## 1. 설정
###############################################################################
MI_STEM  <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
N_SETS   <- NA                 # NA = 전체. 빠른 점검 3.
TMAX     <- 8; GRID <- seq(0, TMAX, by = 0.1); TAU <- 5
MIN_M_OK <- 0.8; MIN_EV_HR <- 5
FIG_H    <- 85                 # mm, 1행 2패널

MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
SETS <- names(MI); if (is.finite(N_SETS)) SETS <- SETS[seq_len(min(N_SETS, length(SETS)))]
M <- length(SETS); rm(MI); invisible(gc())
msg(sprintf("MI 세트 %d개", M))

apply_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]
  z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE))
  }
  d$gLIC_z <- z; d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS); d
}
ADJ_WS <- setdiff(CFG$ADJ, "sex_f")

## 패널 정의: 상태변수, 출발 상태, 회복 판정
PANELS <- list(
  a = list(var = "frail_state", s0 = 2, hit = function(s) s == 1,
           lab = "Pre-frail to robust",              from = "Pre-frail at wave 1"),
  b = list(var = "state",       s0 = 2, hit = function(s) s == 1,
           lab = "Mild disability to normal",        from = "Mild ADL disability at wave 1"),
  ## frail 층: 성별 삼분위로는 나뉘지 않으므로(T3 n=0, T2 n=14) 층 내 중앙값으로 두 집단
  c = list(var = "frail_state", s0 = 3, hit = function(s) s %in% c(1, 2), grp = "half",
           lab = "Frail to pre-frail or robust",     from = "Frail at wave 1; split at the median capacity within the stratum"))
FIG_PANELS <- c("a", "b")          # 본문 Fig 2; c(frail 중앙값 분할, 사후)는 Supplementary Fig 1
GRP_LEV <- list(tert = CFG$TERT_LABELS, half = c("H1", "H2"))
GRP_LAB <- c(LAB_TERTILE, H1 = "Below median (within frail)", H2 = "Above median (within frail)")
GRP_PAL <- c(PAL_TERTILE_LAB, stats::setNames(c(SEV_COL[["bad"]], SEV_COL[["good"]]), GRP_LAB[c("H1", "H2")]))
grp_of <- function(pn) if (is.null(PANELS[[pn]]$grp)) "tert" else PANELS[[pn]]$grp

## 개인별 첫 회복 시각 / 사망 시각 (state==4 = 사망)
first_rec <- function(dat, var, s0, hit) {
  dat <- dat[order(dat$id, dat$time), ]
  sp <- split(seq_len(nrow(dat)), dat$id)
  do.call(rbind, lapply(sp, function(ix) {
    s <- dat[[var]][ix]; tm <- dat$time[ix]; ok <- is.finite(tm) & is.finite(s); s <- s[ok]; tm <- tm[ok]
    if (length(s) < 1 || s[1] != s0) return(NULL)
    tmax <- max(tm)
    kr <- which(hit(s[-1]) & s[-1] != 4); tr <- if (length(kr)) tm[-1][kr[1]] else NA_real_
    kd <- which(s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
    ft <- min(c(tr, td, tmax), na.rm = TRUE)
    ev <- if (is.finite(tr) && tr <= ft) 1L else if (is.finite(td) && td <= ft) 2L else 0L
    data.frame(id = dat$id[ix[1]], ft = max(ft, 1e-6), ev = ev)
  }))
}
## 누적발생 곡선 (Aalen-Johansen) — cloglog 척도 값과 SE
cif_curve <- function(ft, ev, grid) {
  ef <- factor(ev, levels = 0:2, labels = c("censor", "event", "compete"))
  f <- survfit(Surv(ft, ef) ~ 1); s <- summary(f, times = grid, extend = TRUE)
  j <- match("event", f$states); p <- s$pstate[, j]; se <- s$std.err[, j]
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  q <- log(-log(1 - p)); se_q <- se / ((1 - p) * abs(log(1 - p)))
  data.frame(t = grid, p = p, q = q, se_q = se_q)
}

###############################################################################
## 2. 세트별 계산
###############################################################################
CUR <- list(); HR <- list(); NS <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, M))
  load_imputed(set, stem = MI_STEM)
  d <- apply_scale(prep_long(load_kfacs()))
  b <- get_baseline(d)
  b$tert <- d$gLIC_tert_W1cut[match(b$id, d$id)]
  b$gLIC_z <- d$gLIC_z[match(b$id, d$id)]
  ADJ <- intersect(ADJ_WS, names(b))
  for (pn in names(PANELS)) {
    P <- PANELS[[pn]]
    E <- first_rec(d, P$var, P$s0, P$hit)
    if (is.null(E) || !nrow(E)) next
    S <- merge(b[, c("id", "tert", "gLIC_z", ADJ)], E, by = "id"); S <- S[!is.na(S$tert), ]
    if (grp_of(pn) == "half") {                       # 층 내 중앙값 분할 (세트별)
      md <- stats::median(S$gLIC_z, na.rm = TRUE)
      S$tert <- ifelse(S$gLIC_z > md, "H2", "H1")
    }
    LEV <- GRP_LEV[[grp_of(pn)]]
    for (tt in LEV) {
      z <- S[S$tert == tt, ]
      NS[[length(NS) + 1]] <- data.frame(set = set, panel = pn, tert = tt, n = nrow(z),
                                         events = sum(z$ev == 1L), deaths = sum(z$ev == 2L))
      if (nrow(z) < 2 || !any(z$ev == 1L)) next
      cc <- cif_curve(z$ft, z$ev, GRID); cc$set <- set; cc$panel <- pn; cc$tert <- tt
      CUR[[length(CUR) + 1]] <- cc
    }
    ## 원인특이 Cox (T1 기준 삼분위, per +1 s.d.)
    S$tert <- relevel(factor(S$tert, levels = LEV), ref = LEV[1])
    for (term in c("tert", "gLIC_z")) {
      cx <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~", term, "+", paste(ADJ, collapse = "+"))), data = S), silent = TRUE)
      if (inherits(cx, "try-error")) next
      cf <- coef(cx); V <- vcov(cx)
      for (nm in grep(paste0("^", term), names(cf), value = TRUE))
        HR[[length(HR) + 1]] <- data.frame(set = set, panel = pn, term = nm, b = cf[nm], se = sqrt(V[nm, nm]),
                                           n = nrow(S), events = sum(S$ev == 1L))
    }
  }
}
CUR <- do.call(rbind, CUR); HR <- do.call(rbind, HR); NS <- do.call(rbind, NS)

###############################################################################
## 3. 풀링
###############################################################################
pool <- function(x, se) {
  ok <- is.finite(x) & is.finite(se); if (sum(ok) < max(2, ceiling(MIN_M_OK * M))) return(c(NA, NA, NA, NA))
  r <- rubin_pool(x[ok], se[ok]); c(r["est"], r["lcl"], r["ucl"], r["p"])
}
## 곡선
CUR$key <- paste(CUR$panel, CUR$tert, CUR$t)
CV <- do.call(rbind, lapply(split(CUR, CUR$key), function(z) {
  r <- pool(z$q, z$se_q)
  data.frame(panel = z$panel[1], tert = z$tert[1], t = z$t[1],
             p = 100 * (1 - exp(-exp(r[1]))), lo = 100 * (1 - exp(-exp(r[2]))), hi = 100 * (1 - exp(-exp(r[3]))))
}))
CV <- CV[order(CV$panel, CV$tert, CV$t), ]; rownames(CV) <- NULL
## n / events (세트 중앙값)
NSm <- aggregate(cbind(n, events, deaths) ~ panel + tert, NS, function(x) round(median(x)))
## 5년 값
R5 <- CV[abs(CV$t - TAU) < 1e-8, ]
## HR
HR$key <- paste(HR$panel, HR$term)
PH <- do.call(rbind, lapply(split(HR, HR$key), function(z) { r <- pool(z$b, z$se)
  data.frame(panel = z$panel[1], term = z$term[1], n = round(median(z$n)), events = round(median(z$events)),
             hr = exp(r[1]), lo = exp(r[2]), hi = exp(r[3]), p = r[4]) }))
save_vals(CV, "F2_recovery_curves.csv", FIG_DIR)

## 부속 표
fmt_ci <- function(e, l, u) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, u), "NE")
fmt_p  <- function(p) ifelse(!is.finite(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
TAB <- do.call(rbind, lapply(names(PANELS), function(pn) do.call(rbind, lapply(GRP_LEV[[grp_of(pn)]], function(tt) {
  ns <- NSm[NSm$panel == pn & NSm$tert == tt, ]; r5 <- R5[R5$panel == pn & R5$tert == tt, ]
  h <- PH[PH$panel == pn & PH$term == paste0("tert", tt), ]
  data.frame(Transition = PANELS[[pn]]$lab, `Capacity group` = unname(GRP_LAB[tt]),
             n = if (nrow(ns)) ns$n else NA, `Recovery events` = if (nrow(ns)) ns$events else NA,
             `Deaths before recovery` = if (nrow(ns)) ns$deaths else NA,
             `5-yr cumulative incidence, % (95% CI)` = if (nrow(r5)) sprintf("%.1f (%.1f-%.1f)", r5$p, r5$lo, r5$hi) else "-",
             `Cause-specific HR vs lowest group (95% CI)` = if (tt %in% c("T1", "H1")) "1.00 (reference)" else if (nrow(h)) fmt_ci(h$hr, h$lo, h$hi) else "NE",
             P = if (tt %in% c("T1", "H1")) "" else if (nrow(h)) fmt_p(h$p) else "-",
             check.names = FALSE, stringsAsFactors = FALSE) }))))
SD <- PH[PH$term == "gLIC_z", ]
TAB2 <- data.frame(Transition = vapply(SD$panel, function(pn) PANELS[[pn]]$lab, ""), n = SD$n, `Recovery events` = SD$events,
                   `HR per +1 s.d. (95% CI)` = fmt_ci(SD$hr, SD$lo, SD$hi), P = fmt_p(SD$p), check.names = FALSE)
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "By_tertile"); openxlsx::writeData(wb, "By_tertile", TAB)
openxlsx::addWorksheet(wb, "Per_SD");     openxlsx::writeData(wb, "Per_SD", TAB2)
openxlsx::addWorksheet(wb, "Notes"); openxlsx::writeData(wb, "Notes", data.frame(Note = c(
  sprintf("Cumulative incidence of first recovery from the wave-1 state, with death as a competing event (Aalen-Johansen), estimated in each of %d completed datasets and combined by Rubin's rules on the complementary log-log scale.", M),
  "Cause-specific Cox models adjusted for age, baseline comorbidity count, education, household income and area of residence; capacity standardized within sex (sex absorbed).",
  "For the frail stratum the within-sex tertiles do not partition the stratum (240, 14 and 0 participants), so participants were divided at the median capacity within the stratum in each completed dataset. The frail-to-robust per-s.d. estimate is in Table 2.",
  "n, events and deaths are medians across completed datasets.")))
openxlsx::saveWorkbook(wb, file.path(FIG_DIR, "ST_Recovery_Fig2.xlsx"), overwrite = TRUE)

###############################################################################
## 4. 그림
###############################################################################
CV$tert_lab <- factor(GRP_LAB[CV$tert], levels = unname(GRP_LAB))
mk_panel <- function(pn, leg = FALSE) {
  z <- CV[CV$panel == pn, ]
  ns <- NSm[NSm$panel == pn, ]
  sub <- paste(sprintf("%s: n=%d, %d recovered", ns$tert, ns$n, ns$events), collapse = "; ")
  ggplot(z, aes(t, p, colour = tert_lab, fill = tert_lab)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_step(linewidth = 0.6) +
    scale_colour_manual(values = GRP_PAL, name = NULL, breaks = unname(GRP_LAB[GRP_LEV[[grp_of(pn)]]])) +
    scale_fill_manual(values = GRP_PAL, name = NULL, breaks = unname(GRP_LAB[GRP_LEV[[grp_of(pn)]]])) +
    scale_x_continuous(breaks = seq(0, TMAX, 2), limits = c(0, TMAX)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(title = PANELS[[pn]]$lab, subtitle = PANELS[[pn]]$from,
         x = "Years from wave 1", y = "Cumulative incidence of recovery (%)") +
    theme_na() +
    theme(plot.subtitle = element_text(size = BASE_PT - 2, colour = "grey30"),
          legend.position = "bottom", legend.direction = "horizontal")
}
pa <- mk_panel("a"); pb <- mk_panel("b"); pc <- mk_panel("c")
fig <- (pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ns_txt <- paste(vapply(FIG_PANELS, function(pn) { z <- NSm[NSm$panel == pn, ]
  sprintf("%s (%s)", PANELS[[pn]]$lab, paste(sprintf("%s n=%d", z$tert, z$n), collapse = ", ")) }, ""), collapse = "; ")
CAP <- paste0(
  "Cumulative incidence of the first recovery from the wave-1 state, by within-sex tertile of intrinsic capacity ",
  "at wave 1, with death treated as a competing event (Aalen-Johansen estimator). Shaded bands are 95% confidence ",
  sprintf("intervals. Curves were estimated in each of %d multiply imputed datasets and combined by Rubin's rules on ", M),
  "the complementary log-log scale. Tertiles are defined within sex from the wave-1 distribution of the general ",
  "capacity factor. The frail stratum is not shown because the within-sex tertiles do not partition it (240, 14 ",
  "and 0 participants); a post hoc division of frail participants at the median capacity within the stratum is ",
  "shown in Supplementary Fig. 1. Participants ",
  "at wave 1 and recovery events (medians across imputed datasets): ", ns_txt, ". Adjusted cause-specific hazard ",
  "ratios are given in Supplementary Table (ST_Recovery_Fig2).")
fig <- fig + annot_na(width_mm = NA_WIDTH$onehalf, tag = TRUE,
  title = "Recovery from pre-frailty and from mild disability by intrinsic-capacity tertile", caption = CAP)
save_na(fig, "Figure2_Recovery_MI", width_mm = NA_WIDTH$onehalf, height_mm = FIG_H, dir = FIG_DIR)

## ── Supplementary Fig 1: frail 층 사후 중앙값 분할 ────────────────────────
nsc <- NSm[NSm$panel == "c", ]
CAPC <- paste0(
  "Cumulative incidence of the first improvement to pre-frail or robust among participants who were frail at wave 1, ",
  "with death as a competing event (Aalen-Johansen estimator). Because the within-sex tertiles used in Fig. 2 do not ",
  "partition the frail stratum (240, 14 and 0 participants in the lowest, middle and highest tertiles), participants ",
  "were divided post hoc at the median capacity within the stratum in each imputed dataset (",
  paste(sprintf("%s: n = %d, %d improved", GRP_LAB[nsc$tert], nsc$n, nsc$events), collapse = "; "),
  sprintf("). Shaded bands are 95%% confidence intervals; curves were estimated in each of %d multiply imputed datasets ", M),
  "and combined by Rubin's rules on the complementary log-log scale. The corresponding per-s.d. estimate ",
  "(cause-specific hazard ratio for improvement per +1 s.d. of capacity) is given in the text and in Supplementary Table 11.")
figc <- pc + annot_na(width_mm = NA_WIDTH$single,
  title = "Improvement from frailty by capacity within the frail stratum (post hoc median split)", caption = CAPC)
save_na(figc, "SupplFig1_Recovery_Frail_MI", width_mm = NA_WIDTH$single, height_mm = 90, dir = FIG_DIR)

###############################################################################
## 5. 판정 출력
###############################################################################
cat("\n", strrep("=", 70), "\nFigure 2 (회복) — 5년 누적발생 % (T1 / T2 / T3) 및 HR\n", strrep("=", 70), "\n", sep = "")
print(TAB, row.names = FALSE); cat("\n"); print(TAB2, row.names = FALSE)
cat("\n  -> 원고 Results 2절의 [pre-frail T1/T2/T3 %] 자리는 By_tertile 시트 첫 세 행에서 채웁니다.\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
