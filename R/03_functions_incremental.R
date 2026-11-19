###############################################################################
## 03_functions_incremental.R
## 증분가치(incremental value) 분석 공통 헬퍼
##   - Suppl Table 2 (Table 2 확장)  : KF_SupplTable2_TransitionsAdj_v260801.R
##   - Table 4 (감쇠표)              : KF_Table4_Attenuation_v260801.R
##
## 두 스크립트가 같은 계산을 쓰도록 여기에 한 번만 정의합니다.
## KF_common / KF_theme 를 먼저 source 한 뒤에 불러오십시오.
##
## ── 이 파일이 하는 일 ────────────────────────────────────────────────────
##  1) interval 시작 시점(t0)의 시간가변 공변량을 interval 표에 붙입니다.
##     (프레일티 표현형, CHS 연속점수 등 — 노출 IC 와 동일한 시점이어야
##      "같은 순간의 정보끼리" 비교하는 것이 됩니다)
##  2) 보정군을 바꿔 가며 전이별 IC 효과(IRR)를 재추정합니다.
##  3) 감쇠(attenuation)를 로그척도로 계산합니다.
###############################################################################

ADJ_VERSION <- "2026-07-31"

if (!exists("COMMON_VERSION") || COMMON_VERSION < "2026-07-30")
  stop("01_functions_common.R 을 먼저 source 하십시오.")

###############################################################################
## 1. interval 시작 시점 공변량 붙이기
##
##  build_intervals() 결과에는 t0(시작), t1(끝) 만 있고 시간가변 공변량은
##  gLIC_z / gLIC_tert 만 실려 옵니다. 프레일티·CHS 를 "같은 시점"에서
##  가져오려면 long 자료의 (id, time) 로 되찾아야 합니다.
##
##  ★ 시점 비교는 반드시 반올림해서 하십시오. t0 는 실수 연산의 결과라
##    d$time 과 비트 단위로 같지 않을 수 있습니다(조인이 통째로 실패).
###############################################################################
attach_start_cov <- function(iv, d, vars) {
  if (is.null(iv) || !nrow(iv)) return(iv)
  k_d  <- paste(d$id,  sprintf("%.3f", as_num(d$time)))
  keep <- !duplicated(k_d)
  k_iv <- paste(iv$id, sprintf("%.3f", as_num(iv$t0)))
  for (v in vars) {
    if (!v %in% names(d)) {
      iv[[v]] <- NA
      message("[attach_start_cov] '", v, "' 가 자료에 없습니다 -> NA")
      next
    }
    iv[[v]] <- unname(stats::setNames(d[[v]][keep], k_d[keep])[k_iv])
  }
  hit <- vapply(vars, function(v) mean(!is.na(iv[[v]])), numeric(1))
  message(sprintf("[attach_start_cov] 결합률: %s",
                  paste(sprintf("%s %.1f%%", vars, 100 * hit), collapse = ", ")))
  if (any(hit < 0.5))
    warning("시작시점 공변량 결합률이 50% 미만입니다 -> t0 / time 격자를 확인하십시오.")
  iv
}

###############################################################################
## 2. 보정군을 붙인 interval 표 만들기
##
##  frail_f  : interval 시작 시점의 프레일티 표현형 (Robust/Pre-frail/Frail)
##             ※ frailty 체계에서는 from-state 와 같아지므로 층 내 상수가 되고,
##               fit_transition_pois() 의 '상수 공변량 제거'가 자동으로 뺍니다.
##               (그래서 frailty 체계는 M1 = M0 이 됩니다. 정상 동작입니다.)
##  chs_z    : interval 시작 시점의 CHS 총점(0-5)을 표준화한 값
##             ※ "3범주가 거칠어서 IC 효과가 남는 것 아니냐"는 반론에 대한 답.
###############################################################################
prep_iv_adj <- function(d, state_var, exposure = "rolling") {
  iv <- build_intervals(d, state_var, exposure)
  if (is.null(iv)) return(NULL)
  iv <- attach_start_cov(iv, d, c("frail_state", "chs_total"))

  iv$frail_f <- factor(iv$frail_state, levels = 1:3,
                       labels = c("Robust", "Pre-frail", "Frail"))
  ## 사망행(4)은 interval 시작 시점에 올 수 없지만, 혹시 있으면 결측 처리
  iv$frail_f[!is.na(iv$frail_state) & iv$frail_state == 4] <- NA

  ch <- as_num(iv$chs_total)
  iv$chs_z <- if (all(is.na(ch))) NA_real_ else
    (ch - mean(ch, na.rm = TRUE)) / stats::sd(ch, na.rm = TRUE)
  iv
}

###############################################################################
## 3. 보정군별 전이 표
##  fit_transition_pois() 를 그대로 씁니다 -> 희소셀/분리/상수공변량 방어와
##  클러스터-로버스트 SE 가 모두 동일하게 적용됩니다.
###############################################################################
run_adj_set <- function(iv, model_lab, adj, exposure_terms = "gLIC_z",
                        transitions = NULL) {
  if (is.null(iv)) return(NULL)
  adj <- adj[vapply(adj, function(v) {
    vv <- all.vars(stats::as.formula(paste("~", v)))[1]
    vv %in% names(iv) && length(unique(stats::na.omit(iv[[vv]]))) > 1
  }, logical(1))]
  r <- run_transition_table(iv, transitions = transitions,
                            exposure_terms = exposure_terms, adj = adj)
  r$model <- model_lab
  r$adjset <- paste(adj, collapse = " + ")
  r
}

###############################################################################
## 4. 감쇠(attenuation)
##
##  ★ (IRR_adj - 1) / (IRR_base - 1) 로 계산하면 안 됩니다.
##    기저 IRR 이 1 근처인 전이(예: Severe -> Death, IRR 1.00)에서 분모가
##    0 에 가까워져 -1041%, +1948% 같은 값이 나옵니다(2026-07-31 확인).
##    효과크기는 로그척도가 자연척도이므로 log(IRR) 로 계산하고,
##    기저 효과 자체가 무시할 만하면(|log IRR| < LOGCUT) 보고하지 않습니다.
###############################################################################
ATT_LOGCUT <- 0.05          # |log IRR| < 0.05  (IRR 0.95-1.05) 이면 감쇠 미보고

## ★ 추가 규칙 (2026-07-31, 실행 결과 확인 후)
##   기저 추정치가 유의하지 않으면 감쇠를 보고하지 않습니다.
##   실행값에서 Severe -> Death 가 기저 0.88 (P=0.29, 즉 무효과)인데
##   "감쇠 +64%" 로 찍혔습니다. 무효과에서 무효과로 간 것을 '64% 감쇠'로
##   읽으면 오독입니다. 유의하지 않은 기저효과의 감쇠는 해석 대상이 아닙니다.
pct_attn <- function(irr_base, irr_adj, logcut = ATT_LOGCUT, p_base = NULL,
                     p_cut = 0.05) {
  b <- log(irr_base); a <- log(irr_adj)
  ok <- is.finite(b) & is.finite(a) & abs(b) >= logcut
  if (!is.null(p_base)) ok <- ok & is.finite(p_base) & p_base < p_cut
  ifelse(ok, 100 * (1 - a / b), NA_real_)
}

## ★ 보정이 실제로 적용된 행만 골라내기
##   프레일티 체계에서는 출발상태가 곧 프레일티라 보정항이 층 내 상수가 되어
##   자동으로 빠집니다 -> M1 == M0. 이런 행까지 넣어 중앙값을 내면
##   "감쇠 중앙값 0%" 가 되어 실제로 보정이 걸린 전이의 감쇠가 가려집니다.
attn_applied <- function(irr_base, irr_adj, tol = 1e-8) {
  is.finite(irr_base) & is.finite(irr_adj) &
    abs(log(irr_base) - log(irr_adj)) > tol
}
fmt_attn <- function(x) ifelse(is.na(x), "-", sprintf("%+.0f%%", x))

###############################################################################
## 5. 표시 헬퍼
###############################################################################
fmt_irr <- function(irr, lcl, ucl, note = NULL) {
  out <- ifelse(is.finite(irr) & is.finite(lcl) & is.finite(ucl),
                sprintf("%.2f (%.2f-%.2f)", irr, lcl, ucl), "NE")
  out
}

## 전이 행 순서: 악화 -> 회복 -> 사망, 그 안에서 중증도 순
SEV_RANK <- c(Normal = 1, Robust = 1, `No ADL disability` = 1,
              Mild = 2, `Pre-frail` = 2, `ADL disability` = 2,
              Severe = 3, Frail = 3, Death = 9)

kind_of <- function(from, to) {
  f <- SEV_RANK[as.character(from)]; t <- SEV_RANK[as.character(to)]
  ifelse(as.character(to) == "Death", "Death",
         ifelse(!is.na(f) & !is.na(t) & t > f, "Worsening", "Recovery"))
}

message("KF_adj_v260801 loaded  [ADJ_VERSION = ", ADJ_VERSION, "]")
