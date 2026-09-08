###############################################################################
## 01_functions_common.R  (KFACS intrinsic-capacity study - shared functions)
## KFACS Intrinsic Capacity (IC) longitudinal study - Study 2
## COMMON UTILITIES  (level-only; slope terms removed)
##
## 사용법: 다른 모든 스크립트 맨 위에서
##     source("01_functions_common.R")
##
## 설계 원칙
##  - 컬럼명은 인계장 부록 기준. 다만 실제 파일과 다를 수 있으므로 후보군을
##    두고 자동 탐지(pick_col). 없으면 명확한 에러 또는 NULL 반환(선택적 변수).
##  - 삼분위 절단점은 W1 gLIC에서 한 번 정의하여 모든 wave에 공통 적용.
##  - per +1 SD 의 SD 는 W1 gLIC 의 SD 로 고정(웨이브별 SD 표류 방지).
##    => rolling 분석에서도 "동일한 자(ruler)"로 비교됨. 논문 Methods에 명시할 것.
###############################################################################

## ---------------------------------------------------------------- 0. 옵션 ---
options(stringsAsFactors = FALSE, scipen = 999)

## 버전 태그 — 그림/분석 스크립트가 이 값을 확인합니다.
## 구버전 00_common 을 섞어 쓰면 default_transitions() 인자 오류나
## N=3,014 (정상은 3,011) 같은 증상이 납니다.
COMMON_VERSION <- "2026-07-30"

## VERSION 2 (2026-07-29) — 실제 KFACS_imputed_4sets_FINAL_260617.xlsx 구조 검증 반영
##   변경점: (1) 사망 이후 대체행 절단  (2) WAVE_GAP=2.0 로 시간축 생성
##           (3) frailty 체계에 Death 상태 주입  (4) W1 기준 z/삼분위
##           (5) 내장 gLIC_z·gLIC_tertile 과의 불일치 진단 출력

CFG <- list(
  DATA_FILE = "KFACS_imputed_4sets_FINAL_260617.xlsx",
  SHEET     = "MAIN",
  OUTDIR    = ".",
  SEED      = 20260728,
  ## KFACS 웨이브 간격(년). 중도절단자 wave5 시점 8.0 vs followup_years 중앙값 8.03
  ## 에서 역산한 값. 실제 방문일자 변수가 있으면 그것을 쓰는 편이 정확합니다.
  WAVE_GAP  = 2.0,
  ## 사망 이후(=death_wave 초과) 대체행 제거 여부. FALSE 로 두면 '부활' 행이 분석에 들어갑니다.
  DROP_POST_DEATH = TRUE,
  ## 범주형 노출(삼분위) 대비를 보고할 최소 셀 사건수. 어느 한 수준이라도
  ## 이보다 적으면 추정치를 보고하지 않고 NE 로 둡니다.
  MIN_CELL_EVENTS = 5,
  ## 삼분위 라벨: T1 = 최하위 IC, T3 = 최상위 IC
  TERT_LABELS = c("T1", "T2", "T3"),
  ## 기본 보정 공변량 (전이/생존 모형 공통)
  ##  age_c   : 기저연령 + 경과시간 (결정론적 -> 결측 0)
  ##  comorbid: 기저(W1) 값 고정 (W1 결측 0%. 시간가변 값은 후속 웨이브 결측이
  ##            사망과 상관되어 있어 complete-case 시 사망 사건의 29%가 소실됨)
  ADJ = c("age_c", "sex_f", "comorbid_bl", "edu_bl_f", "income_bl_f", "area_bl_f"),   ## primary (CORE) adjustment set
  ## 상태 라벨
  STATE_LAB   = c("Normal", "Mild", "Severe", "Death"),
  FRAILTY_LAB = c("Robust", "Pre-frail", "Frail", "Death")
)


## ------------------------------------------------------------ 1. 패키지 ---
.need <- function(pkgs, install = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    if (install) {
      message("[setup] installing: ", paste(miss, collapse = ", "))
      install.packages(miss, repos = "https://cloud.r-project.org")
    }
    miss <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) warning("[setup] still missing: ", paste(miss, collapse = ", "))
  }
  invisible(TRUE)
}

## 필수: readxl, openxlsx, survival  |  권장: sandwich, lmtest (robust SE)
.need(c("readxl", "openxlsx", "survival", "sandwich", "lmtest"))

suppressPackageStartupMessages(library(survival))
HAS_ROBUST  <- requireNamespace("sandwich", quietly = TRUE) &&
               requireNamespace("lmtest",  quietly = TRUE)
HAS_READXL  <- requireNamespace("readxl",   quietly = TRUE)
HAS_OPENXLSX<- requireNamespace("openxlsx", quietly = TRUE)
if (!HAS_ROBUST) message("[setup] sandwich/lmtest 없음 -> 모형기반 SE 로 대체(클러스터 보정 없음).")

set.seed(CFG$SEED)

## -------------------------------------------------------- 2. 소소한 유틸 ---
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

## 클러스터-로버스트 계수표 (sandwich/lmtest 없으면 모형기반 SE 로 자동 대체) ---
.robust_ct <- function(m, cluster) {
  vc <- NULL
  if (isTRUE(HAS_ROBUST)) vc <- try(sandwich::vcovCL(m, cluster = cluster, type = "HC0"), silent = TRUE)
  if (is.null(vc) || inherits(vc, "try-error")) vc <- stats::vcov(m)
  if (isTRUE(HAS_ROBUST)) {
    ct <- try(lmtest::coeftest(m, vcov. = vc), silent = TRUE)
    if (!inherits(ct, "try-error")) return(ct)
  }
  b  <- stats::coef(m); se <- sqrt(diag(vc))
  se <- se[names(b)]
  z  <- b / se
  out <- cbind(Estimate = b, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(out) <- names(b)
  out
}

## 컬럼 자동 탐지 --------------------------------------------------------------
pick_col <- function(df, candidates, required = TRUE, what = NULL) {
  what <- what %||% candidates[1]
  nm   <- names(df)
  ## 1) 정확 일치
  hit <- candidates[candidates %in% nm]
  if (length(hit)) return(hit[1])
  ## 2) 대소문자 무시
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit)) return(hit[1])
  ## 3) 부분 일치
  for (cand in candidates) {
    hit <- nm[grepl(paste0("^", cand), nm, ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  if (required) {
    stop(sprintf("[pick_col] 필수 컬럼을 찾지 못했습니다: %s\n  후보=%s\n  실제 컬럼=%s",
                 what, paste(candidates, collapse = "/"), paste(nm, collapse = ", ")))
  }
  message(sprintf("[pick_col] 선택적 컬럼 없음 -> 건너뜀: %s", what))
  NULL
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

## 안전한 반올림 포맷 ----------------------------------------------------------
fmt_est <- function(est, lo, hi, d = 2) {
  ifelse(is.na(est), "-",
         sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), est, lo, hi))
}
fmt_p <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_n_pct <- function(n, N, d = 1) sprintf("%d (%.*f)", n, d, 100 * n / N)
fmt_m_sd  <- function(x, d = 1) sprintf("%.*f (%.*f)", d, mean(x, na.rm = TRUE), d, sd(x, na.rm = TRUE))

## --------------------------------------------------------- 3. 데이터 적재 ---
load_kfacs <- function(path = CFG$DATA_FILE, sheet = CFG$SHEET) {
  if (!file.exists(path)) stop("데이터 파일을 찾을 수 없습니다: ", normalizePath(path, mustWork = FALSE))
  ## csv/rds 폴백 (엑셀 패키지가 없거나 사전 변환본을 쓰는 경우)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    d <- utils::read.csv(path, stringsAsFactors = FALSE)
    msg(sprintf("loaded csv: %s rows=%d cols=%d", basename(path), nrow(d), ncol(d)))
    return(d)
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.frame(readRDS(path)))
  if (!HAS_READXL) stop("readxl 패키지가 필요합니다 (또는 데이터를 csv 로 변환해 CFG$DATA_FILE 에 지정).")
  sh <- readxl::excel_sheets(path)
  if (!sheet %in% sh) {
    stop(sprintf("시트 '%s' 없음. 사용 가능한 시트: %s", sheet, paste(sh, collapse = ", ")))
  }
  d <- readxl::read_excel(path, sheet = sheet)
  d <- as.data.frame(d)
  msg(sprintf("loaded: %s [%s]  rows=%d  cols=%d", basename(path), sheet, nrow(d), ncol(d)))
  d
}

list_imputation_sheets <- function(path = CFG$DATA_FILE) {
  if (!HAS_READXL) stop("readxl 패키지가 필요합니다.")
  sh <- readxl::excel_sheets(path)
  imp <- sh[grepl("^(imp|IMP|set|SET|m)[ _-]?[0-9]+$", sh)]
  if (!length(imp)) imp <- setdiff(sh, "MAIN")
  imp
}

## ------------------------------------------- 4. 표준화 / 파생변수 만들기 ---
## 핵심 함수. 원자료(long) -> 분석용 long 으로 정규화.
prep_long <- function(raw) {

  v <- list(
    id       = pick_col(raw, c("id", "ID", "pid", "PID", "subject", "sid", "hhid_pid")),
    wave     = pick_col(raw, c("wave", "WAVE", "visit", "round", "w"), required = FALSE),
    time     = pick_col(raw, c("time", "years", "time_years", "t_years", "fu_time"), required = FALSE),
    state    = pick_col(raw, c("state", "ic_state", "state4")),
    frailty  = pick_col(raw, c("frailty_3cat", "frailty3", "frailty_cat")),
    gLIC     = pick_col(raw, c("gLIC", "glic", "gLIC_score")),
    gLIC_z   = pick_col(raw, c("gLIC_z", "glic_z"), required = FALSE),
    tert_in  = pick_col(raw, c("gLIC_tertile", "glic_tertile", "gLIC_tert"), required = FALSE),
    age      = pick_col(raw, c("age_unified", "age", "age_yr")),
    bage     = pick_col(raw, c("bage", "age_baseline", "age0", "baseline_age"), required = FALSE),
    sex      = pick_col(raw, c("sex_unified", "sex", "gender")),
    comorb   = pick_col(raw, c("comorbid_count", "comorbidity_count", "n_comorbid")),
    death    = pick_col(raw, c("death_event", "death", "died")),
    fu       = pick_col(raw, c("followup_years", "fu_years", "time_to_event", "survtime")),
    adl      = pick_col(raw, c("adl_state", "adl_cat", "adl", "kadl", "adl_disability"), required = FALSE)
  )

  d <- raw
  d$id             <- as.character(d[[v$id]])
  d$state          <- as_num(d[[v$state]])
  d$frailty_3cat   <- as_num(d[[v$frailty]])
  d$gLIC           <- as_num(d[[v$gLIC]])
  d$age_unified    <- as_num(d[[v$age]])
  d$sex_unified    <- as_num(d[[v$sex]])
  d$comorbid_count <- as_num(d[[v$comorb]])
  d$death_event    <- as_num(d[[v$death]])
  d$followup_years <- as_num(d[[v$fu]])
  d$adl_state      <- if (!is.null(v$adl)) as_num(d[[v$adl]]) else NA_real_

  ## wave -----------------------------------------------------------------
  if (!is.null(v$wave)) {
    d$wave <- as_num(d[[v$wave]])
  } else {
    d <- d[order(d$id), , drop = FALSE]
    d$wave <- ave(seq_len(nrow(d)), d$id, FUN = seq_along)
    message("[prep_long] wave 컬럼 없음 -> 행 순서로 wave 생성. 정렬 확인 필요.")
  }

  ## time(년) --------------------------------------------------------------
  if (!is.null(v$time)) {
    d$time <- as_num(d[[v$time]])
  } else {
    d$time <- (d$wave - min(d$wave, na.rm = TRUE)) * CFG$WAVE_GAP
    message(sprintf("[prep_long] time 컬럼 없음 -> (wave-1)*%.2f 년(격자)으로 생성.", CFG$WAVE_GAP))
  }
  d <- d[!is.na(d$followup_years), , drop = FALSE]   # 추적정보 없는 3명 제외

  ## =========================================================================
  ## ★ Death 를 흡수상태(absorbing)로 재구성
  ## -------------------------------------------------------------------------
  ## 구분해야 할 두 가지:
  ##  (1) 상태과정 : Death 는 흡수상태 — 한 번 들어가면 나올 수 없다.
  ##      현 파일은 balanced panel 이라 사망 이후 웨이브에 Normal/Mild/Severe
  ##      대체값이 들어 있어 '부활' 전이가 생긴다.
  ##  (2) 관측과정 : 사망 이후에는 측정 자체가 존재하지 않는다.
  ##      따라서 사망시점 이후의 gLIC/LIC_* 는 값이 있어서는 안 된다.
  ##
  ## 전이강도 모형에서는 Death 가 origin 이 되는 전이가 없으므로 (1)만 지키면
  ## 사후 행을 남기든 지우든 우도 기여가 같다. 그러나 (2)를 어기면 JM 의 lme
  ## 서브모형이 가짜 종단 측정치에 적합된다. 그래서 상태는 흡수시키고
  ## 측정값은 결측 처리한다.
  ##
  ## 또한 death 시점은 격자((wave-1)*GAP)가 아니라 followup_years 가 정확하다.
  ## 둘이 어긋나 격자시점 >= 사망시점이 되는 사례가 실제로 존재하며(KFACS 14명),
  ## 이를 처리하지 않으면 duration<=0 필터에서 사망 사건이 통째로 사라진다.
  ## =========================================================================
  meas_cols <- intersect(c("gLIC", "gLIC_z", "KICOPE", "KICOPE_z", "CIC",
                           grep("^(LIC_|cic_|kicope_)", names(d), value = TRUE)),
                         names(d))
  n0 <- nrow(d)
  d <- d[order(d$id, d$time), , drop = FALSE]
  sp <- split(seq_len(nrow(d)), d$id)

  built <- lapply(sp, function(ix) {
    p  <- d[ix, , drop = FALSE]
    fy <- p$followup_years[1]
    died <- isTRUE(p$death_event[1] == 1)
    obs <- p[is.na(p$state) | p$state != 4, , drop = FALSE]   # 원본 Death 행 제외
    if (!died) return(obs)                                    # 중도절단: 관측 그대로
    obs <- obs[obs$time < fy, , drop = FALSE]                 # 사망 이후 '관측'은 불가능
    if (!nrow(obs)) return(NULL)
    dr <- obs[nrow(obs), , drop = FALSE]                      # 실제 사망시점에 Death 행 생성
    dr$time  <- fy
    dr$state <- 4
    dr$wave  <- NA_real_
    for (cc in meas_cols) dr[[cc]] <- NA_real_                # 사망시점에 측정값 없음
    rbind(obs, dr)
  })
  d <- do.call(rbind, built[!vapply(built, is.null, logical(1))])
  rownames(d) <- NULL
  msg(sprintf("Death 흡수 재구성: %d -> %d 행 | 사망 사건 %d건 (전원 회수되어야 함)",
              n0, nrow(d), sum(d$state == 4, na.rm = TRUE)))

  ## 기저연령 --------------------------------------------------------------
  if (!is.null(v$bage)) {
    d$bage <- as_num(d[[v$bage]])
  } else {
    b <- d[order(d$id, d$time), c("id", "age_unified")]
    b <- b[!duplicated(b$id), ]
    names(b)[2] <- "bage"
    d <- merge(d, b, by = "id", all.x = TRUE)
  }
  d$age0 <- d$bage                     # Cox 공변량명 (JM 스크립트 요구)

  ## 정렬 ------------------------------------------------------------------
  d <- d[order(d$id, d$time), , drop = FALSE]
  rownames(d) <- NULL

  ## ---- W1 기준 스케일 고정 (핵심) ---------------------------------------
  w1 <- d[!duplicated(d$id), ]
  mu1 <- mean(w1$gLIC, na.rm = TRUE); sd1 <- sd(w1$gLIC, na.rm = TRUE)
  if (!is.finite(sd1) || sd1 == 0) stop("[prep_long] W1 gLIC 의 SD 를 계산할 수 없습니다.")
  d$gLIC_zW1 <- (d$gLIC - mu1) / sd1

  ## 원본 gLIC_z 가 있으면 보존하되, 분석에는 gLIC_zW1 사용
  d$gLIC_z_orig <- if (!is.null(v$gLIC_z)) as_num(d[[v$gLIC_z]]) else NA_real_
  d$gLIC_z <- d$gLIC_zW1

  cut1 <- unname(quantile(w1$gLIC, probs = c(1/3, 2/3), na.rm = TRUE))
  d$gLIC_tert_W1cut <- cut(d$gLIC, breaks = c(-Inf, cut1, Inf),
                           labels = CFG$TERT_LABELS, right = TRUE)
  d$gLIC_tert_W1cut <- factor(d$gLIC_tert_W1cut, levels = CFG$TERT_LABELS)

  ## 공변량 정규화 ---------------------------------------------------------
  ## ★ 결측이 결과(사망)와 상관된 시간가변 공변량을 그대로 쓰면 complete-case
  ##   분석에서 사망 사건이 선택적으로 소실된다. KFACS 에서는 시간가변
  ##   comorbid_count 를 쓸 때 사망 500건 중 146건(29%)이 사라진다.
  ##   -> 연령은 결정론적으로 재구성하고, 동반질환은 기저값으로 고정한다.
  d$age_visit <- d$bage + d$time                       # 결측 없음
  d$age_c     <- d$age_visit - mean(w1$age_unified, na.rm = TRUE)
  d$sex_f     <- factor(ifelse(d$sex_unified == 2, "Female", "Male"),
                        levels = c("Male", "Female"))
  cmb <- w1$comorbid_count[match(d$id, w1$id)]
  d$comorbid_bl <- cmb
  d$comorbid_tv <- d$comorbid_count                    # 시간가변 원본 보존(민감도용)

  ## ===== KF_ADJSET_PATCH (v260812) : 기저 사회경제 보정변수 ================
  ## 시간가변 값이 아니라 wave 1 값을 고정합니다 (comorbid_bl 과 같은 방침).
  .ADJ_ABSORB   <- FALSE   ## 결측을 기준수준에 흡수할 것인가
  .ADJ_ABSORB_N <- 30   ## 흡수할 때의 인원 기준
  .bl <- function(nm) {
    if (nm %in% names(w1)) as_num(w1[[nm]])[match(d$id, w1$id)]
    else rep(NA_real_, nrow(d))
  }
  ## 결측은 Unknown 이라는 별도 수준으로 둡니다. 기준수준에 흡수하면
  ## "결측자는 저학력이다" 라고 주장하는 셈입니다. 자료에 없는 것을
  ## 있는 것처럼 채우지 않습니다. (.ADJ_ABSORB 가 TRUE 일 때만 흡수)
  .fac <- function(v, lab_pos, lab_neg, pos_code) {
    x <- ifelse(is.na(v), "Unknown", ifelse(v %in% pos_code, lab_pos, lab_neg))
    if (isTRUE(.ADJ_ABSORB) && sum(x == "Unknown") > 0 &&
        sum(x == "Unknown") < .ADJ_ABSORB_N) {
      return(factor(ifelse(x == "Unknown", lab_neg, x), levels = c(lab_neg, lab_pos)))
    }
    droplevels(factor(x, levels = c(lab_neg, lab_pos, "Unknown")))
  }
  d$edu_bl_f    <- .fac(.bl("edu_high_bl"),  "MidSchoolPlus", "PrimaryOrLess", 1)
  d$income_bl_f <- .fac(.bl("income_high"),  "TopTertile",    "LowerTwo",      1)
  d$alone_bl_f  <- .fac(.bl("living_alone"), "LivingAlone",   "NotAlone",      2)
  ## ★ v260812: 거주지 결측을 조용히 Urban 으로 밀어 넣던 것을 고쳤습니다.
  ##   다른 변수와 같은 규칙 — 결측은 Unknown 으로 그대로 둡니다.
  .ar  <- .bl("area_3cat")
  .arx <- ifelse(is.na(.ar), "Unknown",
                 c("Metro", "Urban", "Rural")[pmax(1, pmin(3, round(.ar)))])
  d$area_bl_f <- {
    if (isTRUE(.ADJ_ABSORB) && sum(.arx == "Unknown") > 0 &&
        sum(.arx == "Unknown") < .ADJ_ABSORB_N) {
      factor(ifelse(.arx == "Unknown", "Metro", .arx),
             levels = c("Metro", "Urban", "Rural"))
    } else {
      droplevels(factor(.arx, levels = c("Metro", "Urban", "Rural", "Unknown")))
    }
  }
  ## ★ v260812: 지수에서 독거를 뺐습니다 (충돌 편향). 교육 + 소득, 0-2 점.
  d$ses_index   <- as.numeric(d$edu_bl_f    == "MidSchoolPlus") +
                   as.numeric(d$income_bl_f == "TopTertile")
  d$ses_index[is.na(d$ses_index)] <- 0
  msg(sprintf("[ADJSET v260812] 중졸이상 %.1f%% · 소득최상위 %.1f%% · 지수평균 %.2f",
              100 * mean(d$edu_bl_f == "MidSchoolPlus", na.rm = TRUE),
              100 * mean(d$income_bl_f == "TopTertile", na.rm = TRUE),
              mean(d$ses_index, na.rm = TRUE)))
  .unk <- function(f) {
    if (!("Unknown" %in% levels(f))) return("없음/흡수됨")
    sprintf("%d명", length(unique(d$id[f == "Unknown"])))
  }
  msg(sprintf("[ADJSET v260812] Unknown 수준 — 교육 %s · 소득 %s · 거주지 %s (흡수=%s)",
              .unk(d$edu_bl_f), .unk(d$income_bl_f), .unk(d$area_bl_f),
              as.character(.ADJ_ABSORB)))
  ## ===== /KF_ADJSET_PATCH ==================================================


  ## ---- ★ frailty 체계에 Death 상태 주입 ---------------------------------
  ## frailty_3cat 는 0/1/2 만 있고 사망 코드가 없어서, 사망 행에서도 직전
  ## frailty 값이 그대로 남습니다. 그대로 두면 Robust->Death 전이가 0건이 됩니다.
  ## (주의: 이 시점의 d 는 재정렬된 상태이므로 위쪽 is_death_row 인덱스를
  ##  재사용하면 안 됩니다. state 로 다시 판정합니다.)
  d$frail_state <- d$frailty_3cat + 1
  d$frail_state[!is.na(d$state) & d$state == 4] <- 4

  ## ---- ADL 체계 (선택) ---------------------------------------------------
  adl_sc <- pick_col(raw, c("adl_score"), required = FALSE)
  if (!is.null(adl_sc)) {
    x <- as_num(d[[adl_sc]])
    d$adl_state <- ifelse(x >= 1, 2, 1)
    d$adl_state[!is.na(d$state) & d$state == 4] <- 4
  }

  ## 라벨 ------------------------------------------------------------------
  d$state_lab   <- factor(CFG$STATE_LAB[d$state],        levels = CFG$STATE_LAB)
  d$frailty_lab <- factor(CFG$FRAILTY_LAB[d$frail_state], levels = CFG$FRAILTY_LAB)

  attr(d, "scale_W1") <- list(mean = mu1, sd = sd1, cutpoints = cut1,
                              n_W1 = nrow(w1))
  msg(sprintf("prep_long: N=%d, person-waves=%d | W1 gLIC mean=%.3f sd=%.3f | cut=%.3f / %.3f",
              length(unique(d$id)), nrow(d), mu1, sd1, cut1[1], cut1[2]))

  ## ---- 내장 변수와의 불일치 진단 ----------------------------------------
  if (!is.null(v$tert_in)) {
    tin <- as.character(d[[v$tert_in]])
    agree <- mean(tin == as.character(d$gLIC_tert_W1cut), na.rm = TRUE)
    msg(sprintf("내장 gLIC_tertile 과 W1기준 삼분위 일치율 = %.1f%%", 100 * agree))
    if (agree < 0.95)
      message("  ※ 내장 gLIC_tertile 은 전체 person-wave 풀 기준으로 만들어졌습니다.\n",
              "    (W1 에서 1/3씩 나뉘지 않음). 인계장 원칙(W1 절단점)과 다르므로\n",
              "    분석에는 gLIC_tert_W1cut 을 사용합니다.")
  }
  if (any(!is.na(d$gLIC_z_orig))) {
    r <- suppressWarnings(stats::cor(d$gLIC_z_orig, d$gLIC_z, use = "complete.obs"))
    sd_o <- stats::sd(d$gLIC, na.rm = TRUE)
    msg(sprintf("내장 gLIC_z 는 전체풀 SD(%.4f) 기준, 본 분석은 W1 SD(%.4f) 기준 (상관 %.3f)",
                sd_o, sd1, r))
  }
  msg(sprintf("보정 공변량 결측률: age_c %.1f%%, comorbid_bl %.1f%% (시간가변 comorbid 원본 %.1f%%)",
              100 * mean(is.na(d$age_c)), 100 * mean(is.na(d$comorbid_bl)),
              100 * mean(is.na(d$comorbid_tv))))
  d
}

## W1 기저 데이터셋 -----------------------------------------------------------
get_baseline <- function(d) {
  b <- d[!duplicated(d$id), , drop = FALSE]
  rownames(b) <- NULL
  b
}

## ------------------------------------------ 5. 전이 person-interval 구성 ---
## exposure = "rolling" : 각 interval 시작 wave(Wk)의 IC
## exposure = "W1"      : 모든 interval 에 W1 IC 를 적용
build_intervals <- function(d, state_var = c("state", "frailty_3cat", "adl_state"),
                            exposure = c("rolling", "W1"),
                            death_level = NULL) {

  state_var <- match.arg(state_var)
  exposure  <- match.arg(exposure)

  ## frailty_3cat 요청 시 Death 가 주입된 frail_state 를 사용
  src <- if (state_var == "frailty_3cat" && "frail_state" %in% names(d)) "frail_state" else state_var
  if (!src %in% names(d) || all(is.na(d[[src]]))) {
    message("[build_intervals] '", state_var, "' 이 없거나 모두 NA -> NULL 반환")
    return(NULL)
  }

  ## 상태를 1..K 정수로 재코딩
  s <- d[[src]]
  if (src == "frailty_3cat") s <- s + 1   # frail_state 가 없을 때만 해당
  d$.s <- s
  if (is.null(death_level)) death_level <- max(s, na.rm = TRUE)
  ## frailty 에는 Death 코드가 별도로 없을 수 있음 -> 4로 강제
  if (state_var == "frailty_3cat") death_level <- 4
  if (state_var == "state")        death_level <- 4

  d <- d[order(d$id, d$time), ]
  sp <- split(d, d$id)

  out <- lapply(sp, function(p) {
    n <- nrow(p)
    rows <- list()

    if (n >= 2) {
      for (k in seq_len(n - 1)) {
        dur <- p$time[k + 1] - p$time[k]
        if (!is.finite(dur) || dur <= 0) next
        rows[[length(rows) + 1]] <- data.frame(
          id = p$id[k], k = k,
          t0 = p$time[k], t1 = p$time[k + 1], dur = dur,
          from = p$.s[k], to = p$.s[k + 1],
          type = "obs",
          stringsAsFactors = FALSE)
      }
    }

    ## 사망 interval: 마지막 관측 wave 이후 사망이 발생했고,
    ## 이미 Death 상태 행으로 기록돼 있지 않은 경우에만 추가
    last_state <- p$.s[n]
    already_death <- isTRUE(last_state == death_level)
    if (!already_death && isTRUE(p$death_event[n] == 1)) {
      dur_d <- p$followup_years[n] - p$time[n]
      if (is.finite(dur_d) && dur_d > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          id = p$id[n], k = n,
          t0 = p$time[n], t1 = p$followup_years[n], dur = dur_d,
          from = last_state, to = death_level,
          type = "death", stringsAsFactors = FALSE)
      }
    }
    if (!length(rows)) return(NULL)
    res <- do.call(rbind, rows)
    ## exposure & 공변량 부착 (interval 시작 시점 wave = k행)
    idx <- if (exposure == "rolling") res$k else rep(1L, nrow(res))
    res$gLIC_z   <- p$gLIC_z[idx]
    res$gLIC_tert<- p$gLIC_tert_W1cut[idx]
    res$age_c    <- p$age_c[idx]
    res$sex_f    <- p$sex_f[idx]
    res$comorbid_count <- p$comorbid_count[idx]
    res$comorbid_bl    <- p$comorbid_bl[1]
    ## ===== KF_ADJSET_PATCH (v260812) =====
    res$edu_bl_f    <- p$edu_bl_f[1]
    res$income_bl_f <- p$income_bl_f[1]
    res$alone_bl_f  <- p$alone_bl_f[1]
    res$area_bl_f   <- p$area_bl_f[1]
    res$ses_index   <- p$ses_index[1]
    ## ===== /KF_ADJSET_PATCH =====
    res$comorbid_tv    <- p$comorbid_tv[idx]
    res$bage     <- p$bage[1]
    res
  })

  iv <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(iv) || !nrow(iv)) stop("[build_intervals] interval 이 하나도 생성되지 않았습니다.")
  rownames(iv) <- NULL
  iv <- iv[!is.na(iv$from) & !is.na(iv$to) & is.finite(iv$dur) & iv$dur > 0, ]
  iv$gLIC_tert <- factor(iv$gLIC_tert, levels = CFG$TERT_LABELS)
  iv$from_lab  <- .lab_for(state_var)[iv$from]
  iv$to_lab    <- .lab_for(state_var)[iv$to]
  attr(iv, "state_var") <- state_var
  attr(iv, "exposure")  <- exposure
  msg(sprintf("intervals[%s / %s]: n=%d, persons=%d, person-years=%.0f",
              state_var, exposure, nrow(iv), length(unique(iv$id)), sum(iv$dur)))
  iv
}

.lab_for <- function(state_var) {
  switch(state_var,
         state          = CFG$STATE_LAB,
         frailty_3cat   = CFG$FRAILTY_LAB,
         adl_state      = c("No ADL disability", "ADL disability", "Severe", "Death"))
}

## -------------------------------- 6. 전이별 Poisson (robust SE, cluster id)---
## origin state 에서 시작한 모든 interval 이 risk-time 을 기여,
## 결과는 특정 destination 으로의 전이 여부. offset = log(person-time).
fit_transition_pois <- function(iv, from_lab, to_lab,
                                exposure_terms = "gLIC_z",
                                adj = CFG$ADJ,
                                min_events = 5) {

  dsub <- iv[iv$from_lab == from_lab, , drop = FALSE]
  if (!nrow(dsub)) return(.empty_trans_row(from_lab, to_lab, exposure_terms))
  dsub$y <- as.integer(dsub$to_lab == to_lab)
  n_ev <- sum(dsub$y, na.rm = TRUE)
  py   <- sum(dsub$dur, na.rm = TRUE)

  if (n_ev < min_events) {
    r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
    r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
    r$note <- sprintf("events < %d (불안정) -> 추정 생략", min_events)
    return(r)
  }

  ## ---- 희소성/분리(separation) 사전 점검 --------------------------------
  ## 범주형 노출인데 어떤 수준에서 사건이 0이면 계수가 ±무한대로 발산한다.
  sparse_note <- ""
  for (x in exposure_terms) {
    v <- all.vars(stats::as.formula(paste("~", x)))[1]
    if (!v %in% names(dsub)) next
    if (is.factor(dsub[[v]]) || is.character(dsub[[v]])) {
      cnt <- tapply(dsub$y, droplevels(factor(dsub[[v]])), sum, na.rm = TRUE)
      if (any(cnt == 0, na.rm = TRUE)) {
        r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
        r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
        r$note <- sprintf("%d events overall, but %s stratum contributes 0 (%s) — tertile contrast not estimable",
                          n_ev, paste(names(cnt)[cnt == 0], collapse = "/"),
                          paste(sprintf("%s=%d", names(cnt), as.integer(cnt)), collapse = ", "))
        return(r)
      }
      ## ★ 사건 0 이 아니어도, 한 수준의 사건이 극소수면 점추정치가 폭주한다.
      ##    (예: Severe->Normal 의 T3 는 사건 1건 -> IRR 9.42, 95% CI 1.04–85.69)
      ##    이런 값은 표에 실리면 그 자체로 오독을 만들므로 보고하지 않는다.
      if (any(cnt < CFG$MIN_CELL_EVENTS, na.rm = TRUE)) {
        r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
        r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
        r$note <- sprintf("%d events overall, but %s stratum has < %d (%s) — tertile contrast not estimable",
                          n_ev, paste(names(cnt)[cnt < CFG$MIN_CELL_EVENTS], collapse = "/"),
                          CFG$MIN_CELL_EVENTS,
                          paste(sprintf("%s=%d", names(cnt), as.integer(cnt)), collapse = ", "))
        return(r)
      }
    }
  }

  rhs <- c(exposure_terms, adj)
  rhs <- rhs[vapply(rhs, function(x) {
    v <- all.vars(stats::as.formula(paste("~", x)))
    all(v %in% names(dsub)) && any(!is.na(dsub[[v[1]]]))
  }, logical(1))]
  ## 상수 공변량 제거 (예: 특정 origin 에 남성만 남는 경우)
  rhs <- rhs[vapply(rhs, function(x) {
    v <- all.vars(stats::as.formula(paste("~", x)))[1]
    length(unique(stats::na.omit(dsub[[v]]))) > 1
  }, logical(1))]

  ## ★ 노출항이 하나도 살아남지 못하면 공변량을 노출로 오보고하지 말 것
  if (!any(exposure_terms %in% rhs)) {
    r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
    r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
    r$note <- "노출 변수가 이 origin 에서 상수/결측 -> 추정 불가"
    return(r)
  }

  f <- stats::as.formula(paste0("y ~ ", paste(rhs, collapse = " + "), " + offset(log(dur))"))
  m <- try(stats::glm(f, family = stats::poisson(), data = dsub), silent = TRUE)
  if (inherits(m, "try-error") || !m$converged) {
    r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
    r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
    r$note <- "모형 수렴 실패"
    return(r)
  }

  ct <- .robust_ct(m, dsub$id)

  keep <- setdiff(rownames(ct), "(Intercept)")
  ## exposure_terms 에 해당하는 항만 추출 (공변량으로 대체 금지)
  want <- unlist(lapply(exposure_terms, function(x) {
    grep(paste0("^", gsub("([.\\[\\]()])", "\\\\\\1", x)), keep, value = TRUE)
  }))
  want <- unique(want)
  if (!length(want)) {
    r <- .empty_trans_row(from_lab, to_lab, exposure_terms)
    r$n_events <- n_ev; r$n_intervals <- nrow(dsub); r$pyears <- py
    r$note <- "노출항이 모형에서 추정되지 않음"
    return(r)
  }

  b   <- ct[want, "Estimate"]
  se  <- ct[want, "Std. Error"]
  ## ★ 수치적 발산 감지: SE 가 비정상적으로 크면 추정 자체를 신뢰할 수 없음
  bad <- !is.finite(b) | !is.finite(se) | se > 3 | abs(b) > 10
  est <- exp(b); lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se)
  p   <- ct[want, 4]
  est[bad] <- NA; lo[bad] <- NA; hi[bad] <- NA; p[bad] <- NA

  nt <- rep(sparse_note, length(want))
  nt[bad] <- "SE 발산(희소 데이터) -> 추정 억제"

  data.frame(
    transition  = paste(from_lab, "->", to_lab),
    from = from_lab, to = to_lab,
    term = want, IRR = as.numeric(est),
    lcl = as.numeric(lo), ucl = as.numeric(hi), p = as.numeric(p),
    n_events = n_ev, n_intervals = nrow(dsub), pyears = py,
    note = nt, stringsAsFactors = FALSE)
}

.empty_trans_row <- function(from_lab, to_lab, terms) {
  data.frame(transition = paste(from_lab, "->", to_lab),
             from = from_lab, to = to_lab, term = terms[1],
             IRR = NA_real_, lcl = NA_real_, ucl = NA_real_, p = NA_real_,
             n_events = 0L, n_intervals = 0L, pyears = 0,
             note = "해당 origin interval 없음", stringsAsFactors = FALSE)
}

## 관심 전이 목록 -------------------------------------------------------------
## ★ 하드코딩하지 말 것. 실제 관측된 모든 비대각(off-diagonal) 전이를
##    상태 라벨 순서대로 열거한다. (예전 버전은 Normal->Severe, Severe->Normal,
##    Frail->Robust 를 목록에서 누락시켜 표에서 통째로 빠져 있었다.)
default_transitions <- function(state_var, iv = NULL, min_events = 1) {
  lab <- .lab_for(state_var)
  if (is.null(iv)) {
    out <- list()
    for (f in lab) for (t in lab) if (f != t) out[[length(out) + 1]] <- c(f, t)
    return(out)
  }
  tb <- table(factor(iv$from_lab, levels = lab), factor(iv$to_lab, levels = lab))
  out <- list()
  for (f in lab) for (t in lab) {
    if (f == t) next
    if (tb[f, t] >= min_events) out[[length(out) + 1]] <- c(f, t)
  }
  out
}

run_transition_table <- function(iv, transitions = NULL,
                                 exposure_terms = "gLIC_z", adj = CFG$ADJ) {
  sv <- attr(iv, "state_var")
  transitions <- transitions %||% default_transitions(sv, iv)
  res <- lapply(seq_along(transitions), function(i) {
    tr <- transitions[[i]]
    r <- fit_transition_pois(iv, tr[1], tr[2], exposure_terms = exposure_terms, adj = adj)
    r$ord <- i                      # 표 행 순서 보존용 (알파벳 정렬 금지)
    r
  })
  out <- do.call(rbind, res)
  out$system   <- sv
  out$exposure <- attr(iv, "exposure")
  out
}

## ------------------------------------------------ 7. 보정 전이율 (/1000py)---
## 삼분위별 보정 전이율: Poisson 모형에서 공변량을 표본 평균으로 고정한 예측값
###############################################################################
## ★ ref 인자 (2026-08-01 추가) — 패널 간 비교를 가능하게 하는 핵심 수정
##
##  이전 판은 보정율을 "그 출발상태 표본(dsub)의 공변량 평균" 에서 예측했습니다.
##  그러면 Normal 패널은 Normal 인 사람들의 평균연령에서, Severe 패널은
##  Severe 인 사람들의 (훨씬 높은) 평균연령에서 예측하게 됩니다.
##  같은 그림에 놓고 세로축까지 통일해 두면 독자는 이것을 직접 비교하는데,
##  실제로는 서로 다른 기준인구의 율을 비교하는 셈이 됩니다.
##  -> Severe / Frail 패널의 율이 실제보다 높게 보입니다.
##
##  ref 에 코호트 전체의 기준 프로파일(list(age_c=..., sex_f=..., comorbid_bl=...))
##  을 넘기면 모든 전이를 같은 기준에서 예측합니다(직접표준화).
##  ref = NULL 이면 이전과 동일하게 동작합니다(하위호환).
###############################################################################
make_ref <- function(iv, adj = CFG$ADJ) {
  r <- list()
  for (a in intersect(adj, names(iv))) {
    v <- iv[[a]]
    r[[a]] <- if (is.numeric(v)) mean(v, na.rm = TRUE)
              else names(sort(table(v), decreasing = TRUE))[1]
  }
  r
}

adjusted_rates_by_tertile <- function(iv, from_lab, to_lab, adj = CFG$ADJ,
                                      min_events = 5, ref = NULL) {
  dsub <- iv[iv$from_lab == from_lab, , drop = FALSE]
  dsub$y <- as.integer(dsub$to_lab == to_lab)
  n_ev <- sum(dsub$y, na.rm = TRUE)
  base <- data.frame(transition = paste(from_lab, "->", to_lab),
                     tertile = CFG$TERT_LABELS,
                     events = NA_integer_, pyears = NA_real_,
                     crude_rate = NA_real_, adj_rate = NA_real_,
                     RR = NA_real_, lcl = NA_real_, ucl = NA_real_,
                     stringsAsFactors = FALSE)
  ## crude
  for (i in seq_along(CFG$TERT_LABELS)) {
    tt <- CFG$TERT_LABELS[i]
    dd <- dsub[which(dsub$gLIC_tert == tt), ]
    base$events[i] <- sum(dd$y, na.rm = TRUE)
    base$pyears[i] <- sum(dd$dur, na.rm = TRUE)
    base$crude_rate[i] <- 1000 * base$events[i] / base$pyears[i]
  }
  if (n_ev < min_events || length(unique(stats::na.omit(dsub$gLIC_tert))) < 2) return(base)

  dsub$gLIC_tert <- stats::relevel(factor(dsub$gLIC_tert, levels = CFG$TERT_LABELS), ref = "T1")
  rhs <- c("gLIC_tert", adj)
  rhs <- rhs[vapply(rhs, function(x) {
    v <- all.vars(stats::as.formula(paste("~", x)))[1]
    v %in% names(dsub) && length(unique(stats::na.omit(dsub[[v]]))) > 1
  }, logical(1))]
  f <- stats::as.formula(paste0("y ~ ", paste(rhs, collapse = " + "), " + offset(log(dur))"))
  m <- try(stats::glm(f, family = stats::poisson(), data = dsub), silent = TRUE)
  if (inherits(m, "try-error") || !m$converged) return(base)

  ct <- .robust_ct(m, dsub$id)

  for (i in seq_along(CFG$TERT_LABELS)) {
    tt <- CFG$TERT_LABELS[i]
    nm <- paste0("gLIC_tert", tt)
    if (tt == "T1") { base$RR[i] <- 1; base$lcl[i] <- NA; base$ucl[i] <- NA }
    else if (nm %in% rownames(ct)) {
      b <- ct[nm, "Estimate"]; s <- ct[nm, "Std. Error"]
      if (is.finite(b) && is.finite(s) && s <= 3 && abs(b) <= 10) {
        base$RR[i] <- exp(b); base$lcl[i] <- exp(b - 1.96 * s); base$ucl[i] <- exp(b + 1.96 * s)
      }   # 발산 시 NA 유지
    }
  }
  ## 보정 전이율: 공변량 평균 프로파일에서의 예측 rate
  nd <- dsub[rep(1, 3), , drop = FALSE]
  nd$gLIC_tert <- factor(CFG$TERT_LABELS, levels = CFG$TERT_LABELS)
  for (a in intersect(adj, names(nd))) {
    if (!is.null(ref) && !is.null(ref[[a]])) {
      ## 코호트 공통 기준 프로파일 (패널 간 비교 가능)
      nd[[a]] <- if (is.numeric(dsub[[a]])) as.numeric(ref[[a]])
                 else factor(as.character(ref[[a]]), levels = levels(factor(dsub[[a]])))
    } else if (is.numeric(dsub[[a]])) nd[[a]] <- mean(dsub[[a]], na.rm = TRUE)
    else nd[[a]] <- names(sort(table(dsub[[a]]), decreasing = TRUE))[1]
  }
  nd$dur <- 1
  pr <- try(stats::predict(m, newdata = nd, type = "response"), silent = TRUE)
  if (!inherits(pr, "try-error")) base$adj_rate <- 1000 * as.numeric(pr)
  ## 발산한 수준은 보정율도 신뢰 불가 -> 함께 억제 (crude_rate 는 그대로 보고)
  bad <- which(base$tertile != "T1" & (is.na(base$RR) | base$events == 0))
  if (length(bad)) base$adj_rate[bad] <- NA_real_
  base
}

## --------------------------------------------------- 8. Rubin 풀링 헬퍼 ---
## est/se 는 로그스케일(로그 IRR/HR) 로 넣을 것.
rubin_pool <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]; se <- se[ok]
  m <- length(est)
  if (m == 0) return(c(est = NA, se = NA, lcl = NA, ucl = NA, p = NA, m = 0))
  if (m == 1) {
    return(c(est = est, se = se,
             lcl = est - 1.96 * se, ucl = est + 1.96 * se,
             p = 2 * stats::pnorm(-abs(est / se)), m = 1))
  }
  qbar <- mean(est)
  ubar <- mean(se^2)
  bvar <- stats::var(est)
  tvar <- ubar + (1 + 1/m) * bvar
  df   <- (m - 1) * (1 + ubar / ((1 + 1/m) * bvar))^2
  crit <- stats::qt(0.975, df = df)
  se_t <- sqrt(tvar)
  c(est = qbar, se = se_t,
    lcl = qbar - crit * se_t, ucl = qbar + crit * se_t,
    p = 2 * stats::pt(-abs(qbar / se_t), df = df), m = m)
}

## 여러 대체셋에 동일 분석 반복 -> 결과 리스트 반환
run_over_imputations <- function(FUN, sheets = NULL, path = CFG$DATA_FILE, ...) {
  sheets <- sheets %||% list_imputation_sheets(path)
  msg("imputation sheets: ", paste(sheets, collapse = ", "))
  res <- lapply(sheets, function(sh) {
    msg("  -> ", sh)
    d <- prep_long(load_kfacs(path, sh))
    FUN(d, ...)
  })
  names(res) <- sheets
  res
}

## --------------------------------------------------------- 9. 저장 유틸 ---
save_xlsx <- function(sheets, file, outdir = CFG$OUTDIR) {
  stopifnot(is.list(sheets), !is.null(names(sheets)))
  if (!HAS_OPENXLSX) {                       # 폴백: 시트별 CSV
    stub <- sub("\\.xlsx$", "", file)
    for (nm in names(sheets)) {
      fp <- file.path(outdir, sprintf("%s__%s.csv", stub, gsub("[^A-Za-z0-9_]", "_", nm)))
      ## 문자열이 이미 UTF-8 이므로 재인코딩하지 말 것 (-> · – 가 깨짐)
      utils::write.csv(sheets[[nm]], fp, row.names = FALSE)
      msg("saved (csv fallback): ", fp)
    }
    return(invisible(outdir))
  }
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    s <- substr(gsub("[\\[\\]:*?/\\\\]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, s)
    openxlsx::writeData(wb, s, sheets[[nm]], withFilter = TRUE)
    openxlsx::freezePane(wb, s, firstRow = TRUE)
    openxlsx::setColWidths(wb, s, cols = seq_len(max(1, ncol(sheets[[nm]]))), widths = "auto")
  }
  fp <- file.path(outdir, file)
  openxlsx::saveWorkbook(wb, fp, overwrite = TRUE)
  msg("saved: ", fp)
  invisible(fp)
}


## ---- 보정군의 출판용 영문 표기 (각주 자동 생성: 코드가 곧 각주의 원천) ----
adj_phrase <- function(adj = CFG$ADJ) {
  lab <- c(age_c       = "age",
           sex_f       = "sex",
           comorbid_bl = "baseline comorbidity count",
           edu_bl_f    = "education",
           income_bl_f = "household income",
           area_bl_f   = "area of residence",
           alone_bl_f  = "living arrangement",
           frail_f     = "the concurrent frailty phenotype",
           chs_z       = "the continuous CHS frailty score",
           ses_index   = "the socioeconomic index")
  x <- unname(ifelse(adj %in% names(lab), lab[adj], adj))
  n <- length(x)
  if (n == 0) return("no covariates")
  if (n == 1) return(x)
  paste0(paste(x[-n], collapse = ", "), " and ", x[n])
}

msg(sprintf("01_functions_common.R loaded  [COMMON_VERSION = %s]", COMMON_VERSION))
