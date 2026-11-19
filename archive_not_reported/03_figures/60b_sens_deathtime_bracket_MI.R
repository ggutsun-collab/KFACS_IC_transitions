###############################################################################
## 60b_sens_deathtime_bracket_MI.R   (Table 3 · ED Table 3 — 사망시점 구간 괄호 민감도)
##                                                                        v261012
## 목적 : 60_ 은 "주분석 = 방문 격자(가장 늦은 시점)" 를 전제로 midpoint 를 대안으로
##        둡니다. 그러나 Table 3 은 실제로 보고 월/구간 중간점(followup_years)으로
##        날짜를 잡은 것으로 보입니다. 이 스크립트는
##          [0] 주분석 사망시각이 무엇인지 직접 판정하고
##          [1] 구간으로만 알려진 사망(death_time_method == "midpoint")을
##              구간 시작(마지막 생존 방문)과 구간 끝(사망 보고 방문)으로 옮겨
##              5년 사망·복합 악화를 m=20 MI, within-sex 척도로 다시 추정합니다.
##        보고 월이 있는 사망(reported_month)은 세 시나리오 모두 그대로 둡니다.
## 산출 : FIG_DIR/ST_DeathTimingBracket.xlsx, ST_DeathTimingBracket_values.csv
## 위치 : KFACS_IC_transitions/03_figures/  (33_figure3_strata_MI.R 과 같은 폴더)
## 실행 : setwd("C:/Users/user/OneDrive/바탕 화면/Nat_Aging/R/KFACS_IC_transitions")
##        source("03_figures/60b_sens_deathtime_bracket_MI.R", encoding = "UTF-8")
###############################################################################
message("\n=== 60b_sens_deathtime_bracket_MI v261012 ===")

## ── 헬퍼 적재 ──────────────────────────────────────────────────────────
## 스크립트 위치 -> 작업디렉터리 -> 프로젝트 고정 경로 순서로 _bootstrap.R 을 찾습니다.
## 찾지 못하면 어디를 찾았는지 출력하고 멈춥니다 (조용히 넘어가지 않음).
NAT_AGING_ROOT <- "C:/Users/user/OneDrive/바탕 화면/Nat_Aging/R/KFACS_IC_transitions"   # 프로젝트 루트 (R/_bootstrap.R 의 상위)
local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of)) { d <- dirname(of); break }
  }
  if (is.null(d) && requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && length(p) && nzchar(p)) d <- dirname(p)
  }
  rel  <- c("R", "../R", "../../R", ".", "..")
  home <- normalizePath("~/..", winslash = "/", mustWork = FALSE)          # C:/Users/user
  roots <- c(NAT_AGING_ROOT,
             file.path(home, c("OneDrive/바탕 화면", "OneDrive/Desktop", "Desktop", "바탕 화면"), "Nat_Aging"))
  cand <- unique(normalizePath(c(file.path(NAT_AGING_ROOT, "R"),
                                 if (!is.null(d)) file.path(d, rel), rel,
                                 file.path(roots, "R"), roots), winslash = "/", mustWork = FALSE))
  hit <- cand[file.exists(file.path(cand, "_bootstrap.R"))]
  if (!length(hit)) {                                                   # 마지막 수단: 33_ 이 있는 폴더를 직접 탐색
    message("[60b] 고정 경로에 없어 OneDrive·Desktop 아래를 탐색합니다 (1-2분 걸릴 수 있음) ...")
    base <- unique(c(file.path(home, "OneDrive"), file.path(home, "Desktop"), file.path(home, "Documents")))
    base <- base[dir.exists(base)]
    f <- unlist(lapply(base, function(b) list.files(b, pattern = "^_bootstrap\\.R$", recursive = TRUE,
                                                    full.names = TRUE, all.files = TRUE)))
    f <- f[!grepl("_retired|submission|\\.Rproj\\.user", f)]
    if (length(f) > 1) {                                                # 33_ 와 같은 폴더에 있는 것을 우선
      w <- file.exists(file.path(dirname(f), "33_figure3_strata_MI.R")); if (any(w)) f <- f[w] }
    hit <- dirname(f)
  }
  if (!length(hit))
    stop("_bootstrap.R 을 찾지 못했습니다. 확인한 위치:\n  ", paste(cand, collapse = "\n  "),
         "\n-> 33_figure3_strata_MI.R 이 있는 폴더에 이 파일을 두고 33_ 과 같은 방법으로 실행하십시오.")
  if (length(hit) > 1) message("[60b] _bootstrap.R 후보 ", length(hit), "개 중 첫 번째 사용:\n  ",
                               paste(hit, collapse = "\n  "))
  source(file.path(hit[1], "_bootstrap.R"), encoding = "UTF-8")
  message("[60b] _bootstrap.R 적재: ", hit[1])
})
if (!exists("kf_init", mode = "function"))
  stop("_bootstrap.R 은 읽었지만 kf_init() 이 정의되지 않았습니다. _bootstrap.R 판본을 확인하십시오.")
if (!exists("PROJECT_ROOT")) setwd(NAT_AGING_ROOT)
kf_init()
.need("survival"); suppressPackageStartupMessages(library(survival))

TAU <- 5; MIN_EV_HR <- 5; MIN_N_HR <- 30
MI_STEM <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")
MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))); SETS <- names(MI); rm(MI); invisible(gc())
SCEN <- c(primary = "As analysed (reported month / interval midpoint)",
          start   = "Interval-dated deaths moved to interval start",
          end     = "Interval-dated deaths moved to interval end")

## ── 33_ 과 동일한 within-sex 척도 ──────────────────────────────────────
apply_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]; z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE)) }
  d$gLIC_z <- z; d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS); d }
ADJ_WS <- setdiff(CFG$ADJ, "sex_f")

## ── 사망 구간 경계 (원자료: wave별 age_365 로 W1 이후 경과시간) ─────────
##  방문 시각(age_365)과 사망 날짜 정보는 대체되지 않는 관측값이므로 MI 세트와 무관합니다.
##  MI csv(imp_*.csv)에는 이 변수들이 없으므로 4-set 종단 RDS 의 MAIN 에서 한 번만 읽습니다.
DB_NEED <- c("id", "wave", "age_365", "death_event", "death_wave", "followup_years", "death_time_method")
find_bounds_source <- function() {
  c1 <- if (exists("IMP_DIR") && exists("IMP_STEM")) file.path(IMP_DIR, paste0(IMP_STEM, ".rds")) else character(0)
  roots <- unique(c(get0("PROJECT_ROOT", ifnotfound = NULL), NAT_AGING_ROOT,
                    dirname(dirname(NAT_AGING_ROOT))))
  c2 <- unlist(lapply(roots[dir.exists(roots)], function(r)
    list.files(r, pattern = "^KFACS_imputed_4sets_LONGITUDINAL.*\\.rds$", recursive = TRUE, full.names = TRUE)))
  c2 <- c2[!grepl("_retired|submission", c2)]
  if (length(c2)) c2 <- c2[order(file.info(c2)$mtime, decreasing = TRUE)]
  for (f in unique(c(c1, c2))) {
    if (!file.exists(f)) next
    x <- try(readRDS(f), silent = TRUE); if (inherits(x, "try-error")) next
    if (!is.data.frame(x)) x <- if ("MAIN" %in% names(x)) x$MAIN else x[[1]]
    if (is.data.frame(x) && all(DB_NEED %in% names(x))) { message("[60b] 사망 구간 정보: ", f); return(as.data.frame(x)) }
  }
  stop("age_365 · death_time_method 가 있는 4-set 종단 RDS 를 찾지 못했습니다.\n",
       "-> 콘솔에서 DB_FILE <- \"<KFACS_imputed_4sets_LONGITUDINAL_260730.rds 경로>\" 를 지정한 뒤 다시 source() 하십시오.")
}
death_bounds <- function(raw) {
  need <- c("id", "wave", "age_365", "death_event", "death_wave", "followup_years", "death_time_method")
  miss <- setdiff(need, names(raw))
  if (length(miss)) stop("원자료에 없는 변수: ", paste(miss, collapse = ", "))
  raw$id <- as.character(raw$id); raw$wave <- as.integer(as_num(raw$wave))
  a1 <- raw[raw$wave == 1, c("id", "age_365")]; names(a1)[2] <- "a1"
  r  <- merge(raw[, c("id", "wave", "age_365")], a1, by = "id")
  r$tw <- as_num(r$age_365) - as_num(r$a1)
  med <- tapply(r$tw, r$wave, stats::median, na.rm = TRUE)          # 결측 방문은 wave 중앙값
  tw  <- function(i, k) { v <- r$tw[r$id == i & r$wave == k]
    if (length(v) && is.finite(v[1])) v[1] else unname(med[as.character(k)]) }
  b1 <- raw[raw$wave == 1, ]
  b1 <- b1[as_num(b1$death_event) == 1 & is.finite(as_num(b1$death_wave)), ]
  data.frame(id = b1$id,
             fu  = as_num(b1$followup_years),
             mth = as.character(b1$death_time_method),
             lo  = mapply(tw, b1$id, as_num(b1$death_wave) - 1),
             hi  = mapply(tw, b1$id, as_num(b1$death_wave)),
             stringsAsFactors = FALSE)
}

## ── 사건 (사망 제외 악화 시각을 따로 둡니다) ────────────────────────────
mk_events <- function(dat) {
  dat <- dat[order(dat$id, dat$time), ]; sp <- split(seq_len(nrow(dat)), dat$id)
  do.call(rbind, lapply(sp, function(ix) {
    s <- dat$state[ix]; tm <- dat$time[ix]; ok <- is.finite(tm); s <- s[ok]; tm <- tm[ok]
    n <- length(tm); if (n < 1 || !is.finite(s[1])) return(NULL)
    s0 <- s[1]
    fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tw  <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] > s0 & s[-1] < 4) else NA_real_
    kd  <- which(is.finite(s) & s == 4); td <- if (length(kd)) tm[kd[1]] else NA_real_
    ka  <- which(is.finite(s) & s < 4);  ta <- if (length(ka)) max(tm[ka]) else NA_real_
    data.frame(id = as.character(dat$id[ix[1]]), s0 = s0, tmax = max(tm),
               t_wor = tw, t_dth = td, t_la = ta) }))
}
shift_death <- function(E, DB, scen) {
  E$t_new <- E$t_dth
  if (scen != "primary") {
    k <- match(E$id, DB$id); im <- !is.na(k) & DB$mth[k] %in% "midpoint" & is.finite(E$t_dth)
    tgt <- if (scen == "start") DB$lo[k] else DB$hi[k]
    E$t_new[im] <- pmax(tgt[im], E$t_la[im], na.rm = TRUE)       # 마지막 생존 관측 이전으로는 당기지 않음
  }
  E$tmax <- pmax(E$tmax, E$t_new, na.rm = TRUE); E
}
mk_surv <- function(E, out) { inf <- function(x) ifelse(is.finite(x), x, Inf)
  te <- if (out == "wor") pmin(inf(E$t_wor), inf(E$t_new)) else inf(E$t_new)
  ft <- pmin(te, E$tmax, TAU); ev <- as.integer(is.finite(te) & te <= ft)
  data.frame(id = E$id, s0 = E$s0, ft = pmax(ft, 1e-6), ev = ev) }
risk5 <- function(ft, ev) { if (!any(ev > 0)) return(c(p = 0, se = NA))
  s <- summary(survfit(Surv(ft, ev) ~ 1), times = TAU, extend = TRUE)
  p <- min(max(1 - s$surv[1], 1e-6), 1 - 1e-6); c(p = p, se = s$std.err[1] * s$surv[1] / (p * (1 - p))) }

###############################################################################
## 1. 세트별 추정
###############################################################################
DB_RAW <- if (exists("DB_FILE") && file.exists(DB_FILE)) {
  x <- readRDS(DB_FILE); if (!is.data.frame(x)) x <- if ("MAIN" %in% names(x)) x$MAIN else x[[1]]; as.data.frame(x)
} else find_bounds_source()
DB <- death_bounds(DB_RAW); rm(DB_RAW)
msg(sprintf("[60b] 사망 %d건 (보고 월 %d · 구간 중간점 %d) | 구간 끝-시작 중앙값 %.2f년",
            nrow(DB), sum(DB$mth %in% "reported_month"), sum(DB$mth %in% "midpoint"),
            stats::median(DB$hi - DB$lo, na.rm = TRUE)))
OUT <- list(); CHK <- list()
for (i in seq_along(SETS)) { set <- SETS[i]; msg(sprintf("[%s] %d/%d", set, i, length(SETS)))
  load_imputed(set, stem = MI_STEM); raw <- load_kfacs()
  d <- apply_scale(prep_long(raw)); d$id <- as.character(d$id); b <- get_baseline(d); b$id <- as.character(b$id)
  b$tert <- d$gLIC_tert_W1cut[match(b$id, d$id)]; b$z <- d$gLIC_z[match(b$id, d$id)]
  b$fr <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))
  b <- b[!is.na(b$tert) & !is.na(b$fr), ]; ADJ <- intersect(ADJ_WS, names(b))
  E0 <- mk_events(d)

  ## [0] 주분석 사망시각 판정: prep_long 의 사망행 시각 vs followup_years vs 방문 격자
  k <- match(E0$id, DB$id); ok <- !is.na(k) & is.finite(E0$t_dth)
  if (i == 1) {
    nd <- sum(is.finite(E0$t_dth)); nm <- sum(ok)
    msg(sprintf("[60b] 사망 id 매칭 %d / %d", nm, nd))
    if (nd > 0 && nm < 0.9 * nd) stop("MI 자료와 사망 구간 정보의 id 가 맞지 않습니다. 예: MI ",
      paste(head(E0$id[is.finite(E0$t_dth)], 3), collapse = ","), " / RDS ", paste(head(DB$id, 3), collapse = ","))
  }
  CHK[[set]] <- data.frame(set = set, n_death = sum(ok),
    med_abs_vs_fu   = stats::median(abs(E0$t_dth[ok] - DB$fu[k[ok]])),
    med_abs_vs_grid = stats::median(abs(E0$t_dth[ok] - DB$hi[k[ok]])),
    n_midpoint = sum(DB$mth[k[ok]] %in% "midpoint"))

  for (sc in names(SCEN)) { E <- shift_death(E0, DB, sc)
    for (out in c("dth", "wor")) {
      S <- merge(b[, c("id", "fr", "tert", "z", ADJ)], mk_surv(E, out), by = "id"); S <- S[S$s0 < 4, ]
      S$cell <- droplevels(factor(paste(S$fr, S$tert, sep = "|")))
      S$cell <- relevel(S$cell, ref = "Robust|T3")
      rhs <- paste(ADJ, collapse = " + ")
      cx <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ cell +", rhs)), data = S), silent = TRUE)
      for (cl in levels(S$cell)) { zc <- S[S$cell == cl, ]; ne <- sum(zc$ev); rk <- risk5(zc$ft, zc$ev)
        nm <- paste0("cell", cl); bh <- NA; sh <- NA
        if (cl != "Robust|T3" && !inherits(cx, "try-error") && nm %in% names(coef(cx)) &&
            ne >= MIN_EV_HR && nrow(zc) >= MIN_N_HR) { bh <- coef(cx)[nm]; sh <- sqrt(vcov(cx)[nm, nm]) }
        OUT[[length(OUT) + 1]] <- data.frame(set = set, scen = sc, outcome = out, term = cl,
          n = nrow(zc), events = ne, p = rk[["p"]], se_logit = rk[["se"]], b = bh, se = sh) }
      for (fs in c("Robust", "Pre-frail")) {                          # ED Table 3: 층 내 1 s.d. 당 HR
        R <- S[S$fr == fs, ]
        cz <- try(coxph(as.formula(paste("Surv(ft, ev == 1L) ~ z +", rhs)), data = R), silent = TRUE)
        bz <- if (inherits(cz, "try-error")) NA else coef(cz)["z"]
        sz <- if (inherits(cz, "try-error")) NA else sqrt(vcov(cz)["z", "z"])
        OUT[[length(OUT) + 1]] <- data.frame(set = set, scen = sc, outcome = out, term = paste0(fs, "|perSD"),
          n = nrow(R), events = sum(R$ev), p = NA, se_logit = NA, b = bz, se = sz) } } } }
L <- do.call(rbind, OUT); rownames(L) <- NULL
CK <- do.call(rbind, CHK)

###############################################################################
## 2. 풀링 (Rubin; 위험은 logit, HR 은 log 척도)
###############################################################################
M <- length(SETS)
pool <- function(x, se) { ok <- is.finite(x) & is.finite(se)
  if (sum(ok) < max(2, ceiling(0.8 * M))) return(c(NA, NA, NA, NA))
  r <- rubin_pool(x[ok], se[ok]); c(r["est"], r["lcl"], r["ucl"], r["p"]) }
P <- do.call(rbind, lapply(split(L, paste(L$scen, L$outcome, L$term)), function(z) {
  rr <- pool(stats::qlogis(pmin(pmax(z$p, 1e-6), 1 - 1e-6)), z$se_logit); rh <- pool(z$b, z$se)
  data.frame(scen = z$scen[1], outcome = z$outcome[1], term = z$term[1],
             n = round(stats::median(z$n)), events = round(stats::median(z$events)),
             risk = 100 * stats::plogis(rr[1]), lo = 100 * stats::plogis(rr[2]), hi = 100 * stats::plogis(rr[3]),
             hr = exp(rh[1]), hr_lo = exp(rh[2]), hr_hi = exp(rh[3]), hr_p = rh[4], stringsAsFactors = FALSE) }))
rownames(P) <- NULL
save_vals(P, "ST_DeathTimingBracket_values.csv", FIG_DIR)

###############################################################################
## 3. 표
###############################################################################
fh <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, h), "-")
fr <- function(e, l, h) ifelse(is.finite(e), sprintf("%.1f (%.1f-%.1f)", e, l, h), "-")
keys <- unique(P[P$scen == "primary", c("outcome", "term")])
keys <- keys[order(match(keys$outcome, c("dth", "wor")), grepl("perSD", keys$term), keys$term), ]
TAB <- data.frame(Outcome = c(dth = "Death", wor = "Composite worsening")[keys$outcome],
                  Term = sub("\\|perSD$", ", per +1 s.d. of IC", gsub("\\|", " / ", keys$term)),
                  check.names = FALSE, stringsAsFactors = FALSE)
for (sc in names(SCEN)) { z <- P[P$scen == sc, ]; k <- match(paste(keys$outcome, keys$term), paste(z$outcome, z$term))
  TAB[[paste0("Events, ", sc)]] <- z$events[k]
  TAB[[paste0("5-y risk %, ", sc)]] <- ifelse(grepl("perSD", keys$term), "", fr(z$risk[k], z$lo[k], z$hi[k]))
  TAB[[paste0("HR, ", sc)]] <- fh(z$hr[k], z$hr_lo[k], z$hr_hi[k]) }
TAB$Outcome[duplicated(TAB$Outcome)] <- ""

nmid <- round(stats::median(CK$n_midpoint))
save_table(TAB, "ST_DeathTimingBracket",
  title = "Supplementary Table | Sensitivity of the 5-year analyses to the dating of deaths known only to an interval",
  footnotes = c(
    sprintf(paste("Deaths were dated to the reported month or, when only the interval between visits was known (%d deaths),",
                  "to the midpoint of that interval. The alternatives move these interval-dated deaths to the start of the",
                  "interval (the last visit at which the participant was assessed) or to its end (the visit at which death",
                  "was recorded); deaths with a reported month are unchanged. The two alternatives bracket the possible dating."), nmid),
    paste("Cells are defined at wave 1 by the frailty phenotype and the within-sex tertile of intrinsic capacity; follow-up",
          "is censored at 5 years. Risks are one minus the Kaplan-Meier estimate at 5 years; hazard ratios are from Cox models",
          "with a single reference cell (robust, highest tertile) or per +1 s.d. of intrinsic capacity within the stratum,",
          "adjusted for age, baseline comorbidity count, education, household income and area of residence."),
    sprintf("Estimates are combined across %d multiply imputed datasets by Rubin's rules (logit scale for risks, log scale for hazard ratios).", M),
    "CI, confidence interval; HR, hazard ratio; IC, intrinsic capacity."))

###############################################################################
## 4. 판정
###############################################################################
cat("\n", strrep("=", 70), "\n[0] 주분석 사망시각은 무엇인가 (세트 중앙값)\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  |사망행 시각 - followup_years| 중앙값 : %.3f년\n", stats::median(CK$med_abs_vs_fu)))
cat(sprintf("  |사망행 시각 - 사망 보고 방문|  중앙값 : %.3f년\n", stats::median(CK$med_abs_vs_grid)))
cat(if (stats::median(CK$med_abs_vs_fu) < stats::median(CK$med_abs_vs_grid))
      "  -> 주분석은 보고 월/구간 중간점 날짜입니다. 60_ 각주의 '방문 격자' 서술은 사실과 다릅니다.\n"
    else "  -> 주분석은 방문 격자입니다. Fig. 3 legend 의 날짜 서술을 고쳐야 합니다.\n")
cat("\n[1] 핵심 추정치\n")
show <- P[P$term %in% c("Robust|perSD", "Robust|T1", "Pre-frail|T2") , ]
show <- show[order(show$outcome, show$term, match(show$scen, names(SCEN))),
             c("outcome", "term", "scen", "events", "risk", "hr", "hr_lo", "hr_hi")]
show[, 5:8] <- round(show[, 5:8], 2); print(show, row.names = FALSE)
r <- P[P$outcome == "dth" & P$term == "Robust|perSD", ]
cat(sprintf("\n  Robust 사망 per s.d. HR 상한이 1 미만인 시나리오: %d / %d\n", sum(r$hr_hi < 1, na.rm = TRUE), nrow(r)))
cat("=== 완료 ===  표:", file.path(FIG_DIR, "ST_DeathTimingBracket.xlsx"), "\n")
