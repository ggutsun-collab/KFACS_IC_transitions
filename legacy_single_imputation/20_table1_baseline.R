###############################################################################
## KF_Table1_Baseline_v260801.R
## Table 1 — Baseline characteristics by intrinsic-capacity tertile (wave 1)
## Nature Aging 서식: 편집가능 파일(csv/xlsx), 세로줄 없음, 각주 하단.
## 데이터: 종단 재대체 v10 · N = 3,011 (wave 1)
###############################################################################

###############################################################################
## 이 스크립트는 작업디렉터리와 무관하게 동작합니다.
## 자기 위치를 찾아 R/_bootstrap.R 을 불러오고, kf_init() 이 헬퍼를 적재합니다.
## (source() / RStudio Source 버튼 / Rscript 모두 지원)
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
if (!exists("attach_extras"))
  stop("02_theme_and_output.R 이 구버전입니다(음주/사회활동 병합 없음).\n",
       "  -> THEME_VERSION 2026-07-31a 이상으로 교체하십시오.")

load_imputed("MAIN")
d <- prep_long(load_kfacs()); b <- get_baseline(d)

###############################################################################
## 항목 구성은 업로드하신 tables_0727.xlsx 의 Table 1 을 그대로 따릅니다.
## 다만 이번 분석은 W1 앵커 · level-only 이므로 다음 3개 행은 뜻이 없어 제외했습니다.
##   · Age at Wave 2 (landmark)      -> W2 랜드마크를 쓰지 않음 (N 2,973 -> 3,011)
##   · IC, Wave-2 gLIC z-score       -> 노출은 rolling landmark 로 처리
##   · IC slope, per year (W1->W2)   -> slope 삭제(편집 방침)
## 데이터에 없는 항목(음주/사회활동)은 자동으로 건너뛰고 각주에 남깁니다.
###############################################################################

## ── 대체 데이터셋에서 빠진 기술변수 붙이기 ─────────────────────────────
## 음주 / 사회활동은 분석 공변량이 아니어서 종단 재대체 대상에서 빠졌습니다.
## 재대체는 필요 없고, wave-1 관측값을 원본에서 붙이면 됩니다.
##   · social_activity : KFACS_social_activity_*.xlsx (id, social_activity)
##   · alcohol         : KFACS_master_FINAL_preimput_*.xlsx 의 wave 1
## 파일은 data 폴더 이하에서 자동 탐색합니다. 못 찾으면 KF_theme 의
## SOCIAL_FILE / PREIMP_FILE 에 전체 경로를 직접 넣으십시오.
b <- attach_extras(b)

## ── 파생 변수 ────────────────────────────────────────────────────────────
b$educ_high   <- as.integer(as_num(b$educ_bl)     >= 5)   # 5-8 = 중학교 이상(확인 2026-07-29)
b$income_high <- as.integer(as_num(b$income_3cat) == 3)
b$married     <- as.integer(as_num(b$marri_status) == 1)
b$female      <- as.integer(as_num(b$sex_unified) == 2)

## 관찰된 상태전이 횟수 (사망 포함, 결측 구간은 건너뜀)
n_switch <- function(v) { v <- v[!is.na(v)]; if (length(v) < 2) 0L else sum(v[-1] != v[-length(v)]) }
tr <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p)
  data.frame(id = p$id[1],
             n_tr_adl     = n_switch(as_num(p$state)),
             n_tr_frailty = n_switch(as_num(p$frailty_3cat)))))
b <- merge(b, tr, by = "id", all.x = TRUE)

g <- droplevels(b$gLIC_tert_W1cut)
LV <- levels(g)

## ── 행 생성기 ────────────────────────────────────────────────────────────
.blank <- function(nm) { r <- data.frame(Characteristic = nm, Overall = "",
                                         stringsAsFactors = FALSE)
  for (l in LV) r[[l]] <- ""; r$P <- ""; r }
hdr <- function(nm) .blank(nm)

## 연속형: mean +- s.d. , P = Kruskal-Wallis (업로드 표의 각주와 동일)
row_cont <- function(x, nm, dg = 1) {
  x <- as_num(x)
  r <- data.frame(Characteristic = nm, Overall = fmt_m_sd(x, dg), stringsAsFactors = FALSE)
  for (l in LV) r[[l]] <- fmt_m_sd(x[g == l], dg)
  ok <- is.finite(x) & !is.na(g)
  p <- if (sum(ok) > 10 && length(unique(g[ok])) > 1)
         try(stats::kruskal.test(x[ok] ~ droplevels(g[ok]))$p.value, silent = TRUE) else NA
  r$P <- fmt_p(if (inherits(p, "try-error")) NA else p); r
}

## 이분형: n (%) , P = chi-square
row_bin <- function(x, nm, show_p = TRUE) {
  x <- as.integer(as_num(x))
  r <- data.frame(Characteristic = nm,
                  Overall = fmt_n_pct(sum(x == 1, na.rm = TRUE), sum(!is.na(x))),
                  stringsAsFactors = FALSE)
  for (l in LV) { xi <- x[g == l]; r[[l]] <- fmt_n_pct(sum(xi == 1, na.rm = TRUE), sum(!is.na(xi))) }
  p <- NA
  if (show_p) {
    tb <- try(table(x, g), silent = TRUE)
    if (!inherits(tb, "try-error") && all(dim(tb) >= c(2, 2)))
      p <- try(suppressWarnings(stats::chisq.test(tb)$p.value), silent = TRUE)
  }
  r$P <- if (show_p) fmt_p(if (inherits(p, "try-error")) NA else p) else ""
  r
}

## 다범주: 머리행에 전체 chi-square P, 아래에 들여쓴 범주행 (업로드 표와 동일한 모양)
row_multi <- function(x, nm, lev, lab = lev) {
  x <- as.character(x); x[!x %in% as.character(lev)] <- NA
  tb <- try(table(x, g), silent = TRUE)
  p  <- if (!inherits(tb, "try-error") && all(dim(tb) >= c(2, 2)))
          try(suppressWarnings(stats::chisq.test(tb)$p.value), silent = TRUE) else NA
  h <- .blank(nm); h$P <- fmt_p(if (inherits(p, "try-error")) NA else p)
  rs <- lapply(seq_along(lev), function(k)
    row_bin(as.integer(x == as.character(lev[k])), paste0("    ", lab[k]), show_p = FALSE))
  do.call(rbind, c(list(h), rs))
}

has <- function(v) v %in% names(b) && any(is.finite(as_num(b[[v]])))
omitted <- character(0)
opt_bin <- function(v, nm, f = function(z) z) {
  if (has(v)) row_bin(f(as_num(b[[v]])), nm) else { omitted <<- c(omitted, nm); NULL }
}

## ── 표 구성 (업로드 파일의 순서/그룹 그대로) ────────────────────────────
R <- list(
  hdr("Demographics"),
  row_cont(b$age_unified, "Age at wave 1, years"),
  row_bin (b$female,      "Female, n (%)"),
  row_cont(b$phys_bmi,    "Body mass index, kg/m^2"),

  hdr("Socioeconomic and social factors"),
  row_bin (b$educ_high,   "Education >= middle school, n (%)"),
  row_bin (b$married,     "Married / with partner, n (%)"),
  row_bin (b$living_alone,"Living alone, n (%)"),
  row_bin (b$income_high, "Income, highest tertile, n (%)"),
  row_multi(as_num(b$area_3cat), "Residential area", 1:3,
            c("Metropolitan (city)", "Urban (dong)", "Rural (eup/myeon)")),

  hdr("Health and lifestyle"),
  row_cont(b$comorbid_count,   "Comorbidity count"),
  row_bin (b$polypharmacy,     "Polypharmacy, n (%)"),
  row_bin (b$smoking_current,  "Current smoker, n (%)"),
  opt_bin ("alcohol_heavy",    "Alcohol >= 2-3 times/week, n (%)"),
  row_cont(b$pa_met_min_wk,    "Physical activity, MET-min/week", 0),
  row_bin (b$strength_2plus,   "Resistance exercise >=2/week, n (%)"),
  opt_bin ("social_activity",  "Social activity, n (%)"),

  hdr("Intrinsic capacity"),
  row_cont(b$gLIC_z, "IC level, wave-1 gLIC z-score", 1),

  hdr("Health state at wave 1"),
  row_multi(as.character(b$state_lab),   "ADL disability state", CFG$STATE_LAB[1:3]),
  row_multi(as.character(b$frailty_lab), "Frailty phenotype",    CFG$FRAILTY_LAB[1:3]),

  hdr("Follow-up and observed transitions"),
  row_cont(b$followup_years, "Follow-up from wave 1, years"),
  row_cont(b$n_tr_adl,       "ADL state transitions, n"),
  row_cont(b$n_tr_frailty,   "Frailty state transitions, n"),
  row_bin (b$death_event,    "Death during follow-up, n (%)"),
  local({
    py <- tapply(as_num(b$followup_years), g, function(z) sum(z, na.rm = TRUE))
    r <- data.frame(Characteristic = "Total person-years",
                    Overall = format(round(sum(as_num(b$followup_years), na.rm = TRUE)),
                                     big.mark = ","), stringsAsFactors = FALSE)
    for (l in LV) r[[l]] <- format(round(as.numeric(py[[l]])), big.mark = ",")
    r$P <- ""; r
  }))
R <- Filter(Negate(is.null), R)

need <- c("Characteristic", "Overall", LV, "P")
TAB <- do.call(rbind, lapply(R, function(r) { for (n in setdiff(need, names(r))) r[[n]] <- ""; r[, need] }))

nn  <- as.integer(table(g))
TER <- c(T1 = "T1, low IC", T2 = "T2, mid", T3 = "T3, high IC")
names(TAB) <- c("Characteristic", sprintf("Overall (n=%d)", nrow(b)),
                sprintf("%s (n=%d)", ifelse(LV %in% names(TER), TER[LV], LV), nn),
                "P value")

## ── 저장 ────────────────────────────────────────────────────────────────
sc <- attr(d, "scale_W1")
fn <- c(
  "Data are mean (s.d.) for continuous variables and n (%) for categorical variables.",
  "P values: Kruskal-Wallis test (continuous) or chi-square test (categorical) across intrinsic capacity (IC) tertiles.",
  sprintf("IC tertiles defined once from the wave-1 gLIC z-score distribution (cut-points %.3f and %.3f) and applied unchanged at every wave.",
          sc$cutpoints[1], sc$cutpoints[2]),
  "Analytic sample = all participants with wave-1 data (n = 3,011). Unlike the earlier wave-2-landmark table (n = 2,973), no landmark restriction is applied because IC slope is not used.",
  "Rows for age at the wave-2 landmark, wave-2 gLIC and IC slope are not shown: the present analysis is anchored at wave 1 and uses IC level only.",
  "Person-years and state transitions are counted from wave 1 to the last observation or death; deaths are treated as an absorbing state.",
  "Education >= middle school corresponds to categories 5-8 of the original 8-level variable.",
  "Alcohol intake and social activity were not part of the multiple-imputation model (they are descriptive variables, not analysis covariates) and are reported as observed at wave 1; percentages use the number of participants with a non-missing value as the denominator.",
  "gLIC, general latent intrinsic-capacity factor; MET, metabolic equivalent of task; s.d., standard deviation.")
if (length(omitted))
  fn <- append(fn, paste0("Not available in the imputed analytic dataset and therefore omitted: ",
                          paste(omitted, collapse = "; "), "."), after = 5)

save_table(TAB, "Table1_Baseline",
  title = "Table 1 | Baseline characteristics of the longitudinal analytic sample, by intrinsic capacity tertile",
  footnotes = fn)

if (length(omitted))
  message("[Table 1] 데이터에 없어 제외된 항목: ", paste(omitted, collapse = ", "))
cat("\n=== Table 1 완료 ===  행 수:", nrow(TAB), "\n")
