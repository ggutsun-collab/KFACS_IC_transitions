###############################################################################
## 20_table1_baseline_MI.R   (Phase 4 — Table 1: within-sex IC 삼분위, m=20 MI 평균)
##                                                                        v260904
## 원본 20_table1_baseline.R 과 행 구성은 동일. 바뀐 것 두 가지:
##   ① 열 = within-sex 삼분위 (14_ 의 apply_scale 과 동일 정의: wave-1 성별 평균·SD,
##      성별 1/3 절단점). pooled 삼분위 판은 Extended Data 로 (원본 스크립트 그대로 사용).
##   ② 값 = 20개 완성 데이터의 평균 (평균·SD·빈도 모두 세트 평균; P 는 세트 중앙값).
##      Table 1 은 기술통계이므로 Rubin 풀링 대신 세트 평균을 씁니다 (표준 관행).
##
## 산출 (FIG_DIR): Table1_Baseline_MI.csv/.xlsx
## 실행: setwd("C:/Users/user/OneDrive/바탕 화면/Nat_Aging"); source("R/20_table1_baseline_MI.R", encoding="UTF-8")
###############################################################################
message("\n=== 20_table1_baseline_MI v260904 ===")

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
if (!exists("attach_extras")) stop("02_theme_and_output.R 이 구버전입니다 (attach_extras 없음).")

###############################################################################
## 1. 설정 · MI 적재
###############################################################################
MI_STEM <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
SETS <- names(MI); M <- length(SETS)
msg(sprintf("MI 세트 %d개", M))

## 14_ 와 동일한 within-sex 척도
apply_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]
  z <- d$gLIC; tt <- rep(NA_character_, nrow(d)); cuts <- list()
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE); cuts[[s]] <- ct
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE))
  }
  d$gLIC_z <- z
  d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS)
  attr(d, "ws_cuts") <- cuts
  d
}

n_switch <- function(v) { v <- v[!is.na(v)]; if (length(v) < 2) 0L else sum(v[-1] != v[-length(v)]) }

###############################################################################
## 2. 세트별 기저 데이터 (b_list)
###############################################################################
B <- vector("list", M); EXTRA <- NULL; CUTS <- list()
for (i in seq_along(SETS)) {
  set <- SETS[i]
  load_imputed(set, stem = MI_STEM)
  d <- apply_scale(prep_long(load_kfacs()))
  CUTS[[set]] <- attr(d, "ws_cuts")
  b <- get_baseline(d)
  b$gLIC_z <- d$gLIC_z[match(b$id, d$id)]                       # within-sex z (W1)
  b$gLIC_tert_W1cut <- d$gLIC_tert_W1cut[match(b$id, d$id)]
  if (is.null(EXTRA)) {                                          # 음주·사회활동: 관측값 1회만
    be <- try(attach_extras(b, verbose = TRUE), silent = TRUE)
    EXTRA <- if (!inherits(be, "try-error"))
      be[, intersect(c("id", "alcohol_heavy", "social_activity"), names(be)), drop = FALSE]
      else data.frame(id = b$id)
  }
  b <- merge(b, EXTRA, by = "id", all.x = TRUE, suffixes = c("", ".x"))
  b$educ_high   <- as.integer(as_num(b$educ_bl)     >= 5)
  b$income_high <- as.integer(as_num(b$income_3cat) == 3)
  b$married     <- as.integer(as_num(b$marri_status) == 1)
  b$female      <- as.integer(as_num(b$sex_unified) == 2)
  tr <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p)
    data.frame(id = p$id[1], n_tr_adl = n_switch(as_num(p$state)),
               n_tr_frailty = n_switch(as_num(p$frailty_3cat)))))
  b <- merge(b, tr, by = "id", all.x = TRUE)
  B[[i]] <- b
  msg(sprintf("[%s] %d/%d  N=%d  삼분위 %s", set, i, M, nrow(b),
              paste(table(b$gLIC_tert_W1cut), collapse = "/")))
}
LV <- CFG$TERT_LABELS

###############################################################################
## 3. 세트 평균 행 생성기
###############################################################################
.g <- function(b) droplevels(factor(b$gLIC_tert_W1cut, levels = LV))
.blank <- function(nm) { r <- data.frame(Characteristic = nm, Overall = "", stringsAsFactors = FALSE)
  for (l in LV) r[[l]] <- ""; r$P <- ""; r }
hdr <- function(nm) .blank(nm)
.medp <- function(ps) { ps <- ps[is.finite(ps)]; if (length(ps)) stats::median(ps) else NA }

## 연속형: 세트별 mean, sd 를 평균; P = Kruskal-Wallis 세트 중앙값
row_cont <- function(var, nm, dg = 1) {
  st <- lapply(B, function(b) { x <- as_num(b[[var]]); g <- .g(b)
    o <- c(mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
    tt <- unlist(lapply(LV, function(l) c(mean(x[g == l], na.rm = TRUE), stats::sd(x[g == l], na.rm = TRUE))))
    ok <- is.finite(x) & !is.na(g)
    p <- if (sum(ok) > 10) try(stats::kruskal.test(x[ok] ~ droplevels(g[ok]))$p.value, silent = TRUE) else NA
    c(o, tt, if (inherits(p, "try-error")) NA else p) })
  A <- colMeans(do.call(rbind, st)[, -length(st[[1]]), drop = FALSE])
  P <- .medp(vapply(st, function(z) z[length(z)], 0))
  r <- data.frame(Characteristic = nm, Overall = sprintf("%.*f (%.*f)", dg, A[1], dg, A[2]), stringsAsFactors = FALSE)
  for (k in seq_along(LV)) r[[LV[k]]] <- sprintf("%.*f (%.*f)", dg, A[1 + 2*k], dg, A[2 + 2*k])
  r$P <- fmt_p(P); r
}
## 이분형: 세트별 n, N 평균 (n 반올림); P = chi-square 세트 중앙값
row_bin <- function(var_or_x, nm, show_p = TRUE) {
  st <- lapply(B, function(b) {
    x <- if (is.character(var_or_x)) as.integer(as_num(b[[var_or_x]])) else as.integer(var_or_x(b))
    g <- .g(b)
    o <- c(sum(x == 1, na.rm = TRUE), sum(!is.na(x)))
    tt <- unlist(lapply(LV, function(l) { xi <- x[g == l]; c(sum(xi == 1, na.rm = TRUE), sum(!is.na(xi))) }))
    p <- NA
    if (show_p) { tb <- try(table(x, g), silent = TRUE)
      if (!inherits(tb, "try-error") && all(dim(tb) >= c(2, 2)))
        p <- try(suppressWarnings(stats::chisq.test(tb)$p.value), silent = TRUE) }
    c(o, tt, if (inherits(p, "try-error")) NA else p) })
  A <- colMeans(do.call(rbind, st)[, -length(st[[1]]), drop = FALSE])
  P <- .medp(vapply(st, function(z) z[length(z)], 0))
  f <- function(n, N) if (N > 0) sprintf("%d (%.1f)", round(n), 100 * n / N) else "-"
  r <- data.frame(Characteristic = nm, Overall = f(A[1], A[2]), stringsAsFactors = FALSE)
  for (k in seq_along(LV)) r[[LV[k]]] <- f(A[1 + 2*k], A[2 + 2*k])
  r$P <- if (show_p) fmt_p(P) else ""; r
}
row_multi <- function(var, nm, lev, lab = lev) {
  ps <- vapply(B, function(b) { x <- as.character(b[[var]]); x[!x %in% as.character(lev)] <- NA
    tb <- try(table(x, .g(b)), silent = TRUE)
    p <- if (!inherits(tb, "try-error") && all(dim(tb) >= c(2, 2)))
      try(suppressWarnings(stats::chisq.test(tb)$p.value), silent = TRUE) else NA
    if (inherits(p, "try-error")) NA_real_ else as.numeric(p) }, 0)
  h <- .blank(nm); h$P <- fmt_p(.medp(ps))
  rs <- lapply(seq_along(lev), function(k)
    row_bin(function(b) { x <- as.character(b[[var]]); as.integer(x == as.character(lev[k])) },
            paste0("    ", lab[k]), show_p = FALSE))
  do.call(rbind, c(list(h), rs))
}
has <- function(v) v %in% names(B[[1]]) && any(is.finite(as_num(B[[1]][[v]])))
omitted <- character(0)
opt_bin <- function(v, nm) if (has(v)) row_bin(v, nm) else { omitted <<- c(omitted, nm); NULL }

## 라벨 열(문자) — 첫 세트 기준으로 생성 (라벨은 세트 간 동일)
for (i in seq_along(B)) {
  B[[i]]$state_chr   <- as.character(B[[i]]$state_lab)
  B[[i]]$frailty_chr <- as.character(B[[i]]$frailty_lab)
  B[[i]]$area_chr    <- as.character(as_num(B[[i]]$area_3cat))
}

###############################################################################
## 4. 표 구성 (원본과 동일한 순서)
###############################################################################
R <- list(
  hdr("Demographics"),
  row_cont("age_unified", "Age at wave 1, years"),
  row_bin ("female",      "Female, n (%)"),
  row_cont("phys_bmi",    "Body mass index, kg/m^2"),
  hdr("Socioeconomic and social factors"),
  row_bin ("educ_high",   "Education >= middle school, n (%)"),
  row_bin ("married",     "Married / with partner, n (%)"),
  row_bin ("living_alone","Living alone, n (%)"),
  row_bin ("income_high", "Income, highest tertile, n (%)"),
  row_multi("area_chr", "Residential area", as.character(1:3),
            c("Metropolitan (city)", "Urban (dong)", "Rural (eup/myeon)")),
  hdr("Health and lifestyle"),
  row_cont("comorbid_count",  "Comorbidity count"),
  row_bin ("polypharmacy",    "Polypharmacy, n (%)"),
  row_bin ("smoking_current", "Current smoker, n (%)"),
  opt_bin ("alcohol_heavy",   "Alcohol >= 2-3 times/week, n (%)"),
  row_cont("pa_met_min_wk",   "Physical activity, MET-min/week", 0),
  row_bin ("strength_2plus",  "Resistance exercise >=2/week, n (%)"),
  opt_bin ("social_activity", "Social activity, n (%)"),
  hdr("Intrinsic capacity"),
  row_cont("gLIC_z", "IC level, wave-1 general-factor z-score (within sex)", 2),
  if (has("gLIC")) row_cont("gLIC", "IC level, unstandardized general-factor score", 2) else NULL,
  hdr("Health state at wave 1"),
  row_multi("state_chr",   "ADL disability state", CFG$STATE_LAB[1:3]),
  row_multi("frailty_chr", "Frailty phenotype",    CFG$FRAILTY_LAB[1:3]),
  hdr("Follow-up and observed transitions"),
  row_cont("followup_years", "Follow-up from wave 1, years"),
  row_cont("n_tr_adl",       "ADL state transitions, n"),
  row_cont("n_tr_frailty",   "Frailty state transitions, n"),
  row_bin ("death_event",    "Death during follow-up, n (%)"),
  local({
    py <- rowMeans(sapply(B, function(b) c(sum(as_num(b$followup_years), na.rm = TRUE),
      tapply(as_num(b$followup_years), .g(b), function(z) sum(z, na.rm = TRUE))[LV])))
    r <- data.frame(Characteristic = "Total person-years", Overall = format(round(py[1]), big.mark = ","),
                    stringsAsFactors = FALSE)
    for (k in seq_along(LV)) r[[LV[k]]] <- format(round(py[1 + k]), big.mark = ",")
    r$P <- ""; r }))
R <- Filter(Negate(is.null), R)
need <- c("Characteristic", "Overall", LV, "P")
TAB <- do.call(rbind, lapply(R, function(r) { for (n in setdiff(need, names(r))) r[[n]] <- ""; r[, need] }))

nn <- round(rowMeans(sapply(B, function(b) as.integer(table(.g(b))))))
TER <- c(T1 = "T1, low IC", T2 = "T2, mid", T3 = "T3, high IC")
names(TAB) <- c("Characteristic", sprintf("Overall (n=%d)", nrow(B[[1]])),
                sprintf("%s (n=%d)", TER[LV], nn), "P value")

###############################################################################
## 5. 저장
###############################################################################
cut_txt <- local({
  cm <- sapply(CUTS, function(cc) c(cc$Male, cc$Female)); cm <- rowMeans(cm)
  sprintf("men %.3f and %.3f; women %.3f and %.3f", cm[1], cm[2], cm[3], cm[4]) })
fn <- c(
  "Data are mean (s.d.) for continuous variables and n (%) for categorical variables, averaged over the 20 multiply imputed datasets (Methods); counts are rounded means.",
  "P values: Kruskal-Wallis test (continuous) or chi-square test (categorical) across intrinsic capacity (IC) tertiles; median across imputed datasets.",
  sprintf("IC tertiles are defined within sex from the wave-1 distribution of the general-factor score (mean cut-points across datasets: %s) and applied unchanged at every wave; each tertile therefore contains the same proportion of women by construction. Tertiles on the pooled (sex-combined) scale are shown in Extended Data.", cut_txt),
  "The within-sex z-score standardizes the general-factor score to the wave-1 mean and s.d. of each sex; the unstandardized score is shown for reference.",
  "Analytic sample = all participants with wave-1 data (n = 3,011).",
  "Person-years and state transitions are counted from wave 1 to the last observation or death; deaths are treated as an absorbing state.",
  "Education >= middle school corresponds to categories 5-8 of the original 8-level variable.",
  "Alcohol intake and social activity were not part of the imputation model (descriptive variables only) and are reported as observed at wave 1; percentages use the number of participants with a non-missing value as the denominator.",
  "MET, metabolic equivalent of task; s.d., standard deviation.")
if (length(omitted))
  fn <- append(fn, paste0("Not available in the analytic dataset and therefore omitted: ", paste(omitted, collapse = "; "), "."), after = 5)

save_table(TAB, "Table1_Baseline_MI",
  title = "Table 1 | Baseline characteristics of the longitudinal analytic sample, by within-sex intrinsic capacity tertile",
  footnotes = fn)
if (length(omitted)) message("[Table 1] 데이터에 없어 제외된 항목: ", paste(omitted, collapse = ", "))
cat("\n=== Table 1 (MI, within-sex) 완료 ===  행 수:", nrow(TAB), "\n")
print(TAB[TAB$Characteristic %in% c("Age at wave 1, years", "Female, n (%)",
        "IC level, wave-1 general-factor z-score (within sex)", "Death during follow-up, n (%)"), ],
      row.names = FALSE)
