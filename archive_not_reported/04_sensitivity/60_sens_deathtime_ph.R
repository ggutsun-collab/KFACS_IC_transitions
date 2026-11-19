###############################################################################
## 60_sens_deathtime_ph.R  (v260814)
## Table 3 보완 민감도 2종 — 심사자 선제 방어용
##
##  (A) 사망 시점 배정 민감도
##      본분석은 모든 사건 시각을 격자(2년 방문) 위에 둡니다. 사망은 "사망이
##      보고된 방문 시각"에 배정되므로 늦은쪽 극단입니다. 대안 배정으로
##      "마지막 생존 확인 방문과 보고 방문의 중간점"을 사용해 Table 3 의
##      복합 악화·사망·회복 Cox HR 이 달라지는지 봅니다.
##  (B) 비례위험(Schoenfeld) 진단
##      Table 3 의 세 Cox 모형(복합 악화 · 사망 · 회복 원인별)에 대해
##      cox.zph 전역 검정과 셀(노출) 항 검정을 보고합니다.
##
##  출력: ST_DeathTiming.xlsx / _values.csv  (PH 결과 포함)
###############################################################################

DTPH_VERSION <- "v260814"
message("\n=== 60_sens_deathtime_ph ", DTPH_VERSION, " ===")

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
if (!requireNamespace("survival", quietly = TRUE)) stop("survival 패키지가 필요합니다.")

TAU <- 5

## ── 1. 자료 (22와 동일한 기저 구성) ─────────────────────────────────────────
load_imputed("MAIN")
d <- prep_long(load_kfacs())
b <- get_baseline(d)
b$tert <- factor(b$gLIC_tert_W1cut, levels = c("T1", "T2", "T3"))
b$fr   <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))
b      <- b[!is.na(b$tert) & !is.na(b$fr), , drop = FALSE]
ADJ <- intersect(CFG$ADJ, names(b))

## ── 2. 사건 시각 (22의 mk_events + 마지막 생존 확인 시각 t_la) ─────────────
mk_events <- function(dat, sysvar) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  sp  <- split(seq_len(nrow(dat)), dat$id)
  res <- lapply(sp, function(ix) {
    s  <- dat[[sysvar]][ix]; tm <- dat$time[ix]
    ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
    n <- length(tm); if (n < 1) return(NULL)
    s0 <- s[1]; if (!is.finite(s0)) return(NULL)
    tmax <- max(tm)
    fst  <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tr <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] <  s0 & s[-1] >= 1) else NA_real_
    tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] >  s0)              else NA_real_
    kd <- which(is.finite(s) & s == 4)
    td <- if (length(kd)) tm[kd[1]] else NA_real_
    ka <- which(is.finite(s) & s < 4)
    ta <- if (length(ka)) max(tm[ka]) else NA_real_    # 마지막 생존 확인
    data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax,
               t_rec = tr, t_wor = tw, t_dth = td, t_la = ta)
  })
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

mk_surv <- function(E, which_out, death_time = c("grid", "midpoint")) {
  death_time <- match.arg(death_time)
  inf <- function(x) ifelse(is.finite(x), x, Inf)
  td  <- E$t_dth
  if (death_time == "midpoint") {
    mid <- (pmin(E$t_la, E$t_dth, na.rm = TRUE) + E$t_dth) / 2
    td  <- ifelse(is.finite(E$t_dth) & is.finite(E$t_la) & E$t_la < E$t_dth,
                  mid, E$t_dth)
  }
  tmax_eff <- ifelse(is.finite(td), pmax(td, E$t_la, na.rm = TRUE), E$tmax)
  if (which_out == "rec") {
    ft <- pmin(inf(E$t_rec), inf(td), tmax_eff, TAU)
    ev <- rep(0L, nrow(E))
    ev[is.finite(E$t_rec) & E$t_rec <= ft] <- 1L
    ev[ev == 0L & is.finite(td) & td <= ft] <- 2L
  } else if (which_out == "wor") {
    twd <- pmin(inf(E$t_wor), inf(td))                 # 복합: 악화 또는 사망 중 이른 것
    ft  <- pmin(twd, tmax_eff, TAU)
    ev  <- as.integer(is.finite(twd) & twd <= ft)
  } else {
    ft <- pmin(inf(td), tmax_eff, TAU)
    ev <- as.integer(is.finite(td) & td <= ft)
  }
  ft <- pmax(ft, 1e-6)
  data.frame(id = E$id, s0 = E$s0, ft = ft, ev = ev)
}

## ── 3. Cox 적합기 (전역 기준칸 Robust|T3) ──────────────────────────────────
fit_cells <- function(SV, base, restrict_ids = NULL) {
  M <- merge(SV, base[, c("id", "fr", "tert", ADJ)], by = "id")
  if (!is.null(restrict_ids)) M <- M[M$id %in% restrict_ids, , drop = FALSE]
  M$cell <- interaction(M$fr, M$tert, sep = "|")
  M$cell <- stats::relevel(factor(as.character(M$cell)), ref = "Robust|T3")
  fml <- stats::as.formula(paste("survival::Surv(ft, ev == 1L) ~ cell",
                                 if (length(ADJ)) paste("+", paste(ADJ, collapse = " + ")) else ""))
  cx <- try(survival::coxph(fml, data = M), silent = TRUE)
  list(model = cx, data = M)
}
cells_of <- function(fit) {
  if (inherits(fit$model, "try-error")) return(NULL)
  ci <- summary(fit$model)$conf.int; cf <- summary(fit$model)$coefficients
  rn <- rownames(ci)[startsWith(rownames(ci), "cell")]
  data.frame(cell = sub("^cell", "", rn),
             HR = ci[rn, "exp(coef)"], lo = ci[rn, "lower .95"], hi = ci[rn, "upper .95"],
             p = cf[rn, ncol(cf)], stringsAsFactors = FALSE)
}

## ── 4. 두 배정에서 세 결과 적합 ─────────────────────────────────────────────
E_adl <- mk_events(d, "state")                        # ADL 3-state 체계
b_use <- b[b$id %in% E_adl$id, , drop = FALSE]
rec_ids <- E_adl$id[E_adl$s0 %in% c(2, 3)]            # wave-1 장애 보유자 (회복 분석)

RES <- list(); PH <- NULL
for (out in c("wor", "dth", "rec")) {
  rid <- if (out == "rec") rec_ids else NULL
  f0 <- fit_cells(mk_surv(E_adl, out, "grid"),     b_use, rid)
  f1 <- fit_cells(mk_surv(E_adl, out, "midpoint"), b_use, rid)
  c0 <- cells_of(f0); c1 <- cells_of(f1)
  if (is.null(c0) || is.null(c1)) { message("  !! ", out, " 적합 실패"); next }
  m  <- merge(c0, c1, by = "cell", suffixes = c("_grid", "_mid"))
  m$outcome <- c(wor = "Composite worsening", dth = "Death", rec = "Recovery")[out]
  m$ratio_pct <- 100 * (log(m$HR_mid) - log(m$HR_grid)) / abs(log(pmax(m$HR_grid, 1e-12)))
  m$ratio_pct[!is.finite(m$ratio_pct) | abs(log(m$HR_grid)) < log(1.05)] <- NA
  ev0 <- tapply(f0$data$ev == 1L, f0$data$cell, sum)[m$cell]
  ev1 <- tapply(f1$data$ev == 1L, f1$data$cell, sum)[m$cell]
  m$ev_grid <- as.integer(ev0); m$ev_mid <- as.integer(ev1)
  RES[[out]] <- m
  ## PH 진단 (본배정 모형) — 특이행렬이면 단계적으로 단순화해 재검정
  ##  1단계: 전체 모형 그대로
  ##  2단계: 10명 미만 셀 제외
  ##  3단계: + 보정변수를 age/sex/comorbid 로 축소 (희소 Unknown 수준 제거)
  ##  4단계: + 보정변수를 age/sex 로 축소
  .zph_try <- function(Mz, adj) {
    Mz <- droplevels(Mz)
    if (nlevels(Mz$cell) < 2) return(NULL)
    refz <- if ("Robust|T3" %in% levels(Mz$cell)) "Robust|T3" else
            names(sort(table(Mz$cell), decreasing = TRUE))[1]
    Mz$cell <- stats::relevel(Mz$cell, ref = refz)
    adj <- intersect(adj, names(Mz))
    fmlz <- stats::as.formula(paste("survival::Surv(ft, ev == 1L) ~ cell",
                                    if (length(adj)) paste("+", paste(adj, collapse = " + ")) else ""))
    cz <- try(survival::coxph(fmlz, data = Mz), silent = TRUE)
    if (inherits(cz, "try-error")) return(NULL)
    z <- try(survival::cox.zph(cz), silent = TRUE)
    if (inherits(z, "try-error")) NULL else z
  }
  Mfull <- f0$data
  Msub  <- Mfull[Mfull$cell %in% names(which(table(Mfull$cell) >= 10)), , drop = FALSE]
  steps <- list(
    list(M = Mfull, adj = ADJ,                                A = ""),
    list(M = Msub,  adj = ADJ,                                A = "cells with fewer than 10 participants excluded"),
    list(M = Msub,  adj = c("age_c", "sex_f", "comorbid_bl"), A = "sparse cells excluded; adjustment reduced to age, sex and comorbidity"),
    list(M = Msub,  adj = c("age_c", "sex_f"),                A = "sparse cells excluded; adjustment reduced to age and sex"))
  z <- NULL; ph_note <- ""
  for (st in steps) {
    z <- .zph_try(st$M, st$adj)
    if (!is.null(z)) { ph_note <- st$A; break }
  }
  if (!is.null(z)) {
    if (nzchar(ph_note)) message("  * cox.zph(", out, "): ", ph_note)
    tb <- z$table
    PH <- rbind(PH, data.frame(
      outcome = m$outcome[1],
      p_cell   = if ("cell" %in% rownames(tb)) tb["cell", "p"] else NA_real_,
      p_global = tb["GLOBAL", "p"], note = ph_note, stringsAsFactors = FALSE))
  } else message("  * cox.zph 실패(", out, "): 모든 단계에서 특이행렬")
}
LONG <- do.call(rbind, RES); rownames(LONG) <- NULL
save_vals(LONG, "ST_DeathTiming_values.csv", FIG_DIR)

## ── 5. 표 ───────────────────────────────────────────────────────────────────
fmt <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, h), "-")
TAB <- data.frame(
  Outcome = LONG$outcome,
  Cell    = gsub("\\|", " / ", LONG$cell),
  `Events, visit-time` = LONG$ev_grid,
  `HR, visit-time (95% CI)` = fmt(LONG$HR_grid, LONG$lo_grid, LONG$hi_grid),
  `Events, midpoint` = LONG$ev_mid,
  `HR, midpoint (95% CI)` = fmt(LONG$HR_mid, LONG$lo_mid, LONG$hi_mid),
  `Change in log HR, %` = ifelse(is.finite(LONG$ratio_pct), sprintf("%.1f", LONG$ratio_pct), "-"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$Outcome[duplicated(TAB$Outcome)] <- ""

PHTAB <- data.frame(
  Outcome = PH$outcome,
  `P, exposure term (Schoenfeld)` = sprintf("%.3f", PH$p_cell),
  `P, global (Schoenfeld)`        = sprintf("%.3f", PH$p_global),
  Note = PH$note,
  check.names = FALSE, stringsAsFactors = FALSE)

mx  <- suppressWarnings(max(abs(LONG$ratio_pct), na.rm = TRUE))
med <- stats::median(abs(LONG$ratio_pct), na.rm = TRUE)
nfl <- sum(sign(log(LONG$HR_grid)) != sign(log(LONG$HR_mid)) &
           is.finite(LONG$ratio_pct), na.rm = TRUE)

wb_extra <- list(PH = PHTAB)
save_table(TAB, "ST_DeathTiming",
  title = "Supplementary Table | Sensitivity of the 5-year landmark analyses to the assignment of death times, with proportional-hazards diagnostics",
  footnotes = c(
    "All event times in the primary analyses lie on the biennial assessment grid; deaths are assigned to the visit at which death was recorded, the latest defensible time. The alternative assigns each death to the midpoint between the last visit at which the participant was known to be alive and the visit at which death was recorded, the earliest defensible time under interval observation. The two assignments bracket the possible dating of deaths.",
    "Hazard ratios are from Cox models with a single global reference (robust, highest intrinsic-capacity tertile), adjusted for age, sex, baseline comorbidity count, education, household income and area of residence, with follow-up censored at 5 years. Recovery is analysed among wave-1 disability carriers with death as a competing event (cause-specific models shown).",
    "Change is the difference in log hazard ratios expressed as a percentage of the primary log hazard ratio and is not shown where the primary estimate is negligible (HR between 0.95 and 1.05).",
    sprintf("Across all cells the median absolute change was %.1f%% and the maximum %.1f%%; %d estimate(s) changed direction.", med, mx, nfl),
    "The PH sheet reports scaled-Schoenfeld-residual tests (cox.zph) for the primary models: the exposure (cell) term and the global test.",
    "CI, confidence interval; HR, hazard ratio; PH, proportional hazards."))

## PH 시트를 같은 파일에 추가
if (requireNamespace("openxlsx", quietly = TRUE)) {
  fx <- file.path(FIG_DIR, "ST_DeathTiming.xlsx")
  if (file.exists(fx)) {
    wb <- openxlsx::loadWorkbook(fx)
    if (!("PH" %in% names(wb))) openxlsx::addWorksheet(wb, "PH")
    openxlsx::writeData(wb, "PH", PHTAB)
    openxlsx::saveWorkbook(wb, fx, overwrite = TRUE)
  }
}

cat("\n", strrep("=", 70), "\n요약\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  사망시점 배정: |변화| 중앙 %.1f%% · 최대 %.1f%% · 방향 역전 %d\n", med, mx, nfl))
cat("  Schoenfeld PH 검정:\n"); print(PHTAB, row.names = FALSE)
cat("=== 완료 ===  표:", file.path(FIG_DIR, "ST_DeathTiming.xlsx"), "\n")
