###############################################################################
## F2_main_JMdynpred_260730.R
## Figure 2 (MAIN) — Joint-model dynamic predictions  (VALUE-ONLY, slope 삭제)
##
##   행 a/b/c = IC 프로파일 3개 (하위 10 / 중앙 / 상위 90 백분위)
##   열       = 랜드마크 (자료의 추적기간에 맞춰 자동 결정, 기본 0 / 2 / 4 년)
##   좌축 = IC 궤적, 우축 = 예측 생존확률
##
##   ★ 기본값은 '가상 프로파일'입니다. 실제 참가자 3명을 쓰면 연령·동반질환
##     차이가 IC 효과를 덮어 Low IC 의 생존이 High IC 보다 높게 나옵니다
##     (2026-07-31 확인). PROFILE_MODE 참조.
##
## ── 지난 세션 버그 (재발 방지: 코드에 고정) ─────────────────────────────
##   predict(., process = "event") 는 newdata 에 Cox 적합에 쓰인 생존변수가
##   있어야 합니다. Cox 를 Surv(t, ev) 로 적합했으므로
##       nd$t <- max(nd$time);  nd$ev <- 0
##   또한 Cox 공변량명이 age0 이므로 long$age0 필요 (prep_long 이 생성).
##
## ── 랜드마크 선정 원칙 ──────────────────────────────────────────────────
##   고정 랜드마크(0/4/8년)를 쓰고, 각 랜드마크에서 event-free 이며 이후
##   예측 구간이 남아 있는 참가자만 고릅니다. '최종 측정' 을 랜드마크로 쓰면
##   예측할 구간이 남지 않고, 사건이 랜드마크보다 먼저 온 참가자가 섞입니다.
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
## ── 버전 확인 (구버전 파일이 섞이면 여기서 멈춥니다) ────────────────────
if (!exists("COMMON_VERSION") || COMMON_VERSION < "2026-07-30")
  stop("01_functions_common.R 이 아닙니다. 구버전 00_common 이 로드되었습니다.")
if (!exists("THEME_VERSION")  || THEME_VERSION  < "2026-08-01")
  stop("02_theme_and_output.R 이 아닙니다. 구버전 90_theme 이 로드되었습니다.")
if (!exists("load_imputed"))
  stop("load_imputed() 가 없습니다 -> 02_theme_and_output.R 을 다시 받으십시오.")
.need(c("JMbayes2", "nlme"))
suppressPackageStartupMessages({ library(JMbayes2); library(nlme) })

## 경로는 90_theme 의 IMP_DIR / IMP_STEM / FIG_DIR
###############################################################################
## ★ 랜드마크 / 예측지평 (2026-07-31 수정)
##
##  1차 렌더에서 landmark 4y·8y 패널의 예측곡선이 8.5년쯤에서 끊겼습니다.
##  x축 문제가 아니라 **데이터 문제**입니다. 이 코호트의 최대 추적기간이
##  약 8.6년이라, 그 너머는 기저위험이 정의되지 않아 JM 이 값을 내지 못합니다
##  (predict() 가 NA 를 반환 -> 필터에서 잘림).
##
##  따라서 랜드마크와 지평을 실제 추적기간에 맞춥니다.
##    · 요청값을 두고, 자료에서 계산한 TMAX 로 자동 조정
##    · 각 랜드마크가 최소 MIN_WINDOW 년의 예측구간을 갖도록 걸러냄
##    · 조정 내역은 콘솔에 출력 + 그림 각주에 반영
###############################################################################
## ★ 랜드마크 = "이 시점까지 관측된 IC 로 예측을 갱신한다" 는 뜻이므로,
##   반드시 **측정 wave 시점**이어야 합니다(중간 시점은 새 정보가 없습니다).
##   wave 간격 2년 -> W1=0, W2=2, W3=4, W4=6, W5=8년.
##   마지막 랜드마크는 MIN_WINDOW 년 이상 예측할 구간이 남아야 하므로
##   W5(8년)는 쓸 수 없습니다(추적 상한 8.2년 -> 남는 구간 0.2년).
##   기본값 0/2/4/6 = W1 / W1-W2 / W1-W3 / W1-W4 (쓸 수 있는 랜드마크 전부).
##   3열로 줄이려면 c(0, 4, 6) 으로 바꾸십시오.
LANDMARKS_REQ <- c(0, 2, 4, 6)   # 요청 랜드마크 (자료에 맞춰 자동 조정됨)
HORIZON_REQ   <- Inf          # Inf = 각 랜드마크에서 추적 종료(TMAX)까지 예측.
                              # 숫자를 넣으면 그 년수로 상한을 겁니다.
MIN_WINDOW    <- 2            # 이보다 짧은 예측구간만 남는 랜드마크는 버림

###############################################################################
## ★★ 프로파일 방식 (2026-07-31 수정) ★★
##
##  이전 판은 실제 참가자 3명(IC 10/50/90 백분위)을 골라 그렸습니다.
##  그 결과 "Low IC" 참가자의 예측 생존이 "High IC" 참가자보다 **높게** 나왔습니다.
##  IC 로만 고르고 연령·성별·동반질환을 통제하지 않아, 개인차가 IC 효과를
##  덮어버린 것입니다. 그림의 주장("IC 가 생존을 예측한다")과 정반대로 읽힙니다.
##
##  "synthetic" : 공변량을 코호트 대표값으로 고정하고 IC 궤적만 10/50/90
##                백분위로 달리한 가상 참가자 3명. 차이가 오직 IC 에서만 오므로
##                해석이 명확하고 순서가 뒤집히지 않습니다. (기본값·권장)
##  "observed"  : 실제 참가자. 단, 연령·동반질환이 중앙값 근처인 사람만 후보로
##                제한해 교란을 줄입니다.
###############################################################################
###############################################################################
## ★ "matched" 모드 추가 (2026-08-01) — 참고문헌(Stolz JGA 2022 Fig 3)과 동형
##
##  "synthetic" 의 문제: IC 궤적을 z0 + (평균기울기)*t 로 만들기 때문에 점이
##  완벽한 직선 위에 놓입니다. 실제 반복측정에서는 있을 수 없는 모양이고,
##  캡션은 그 점을 "IC trajectory supplied to the model" 이라고 부르므로
##  독자는 관측값으로 읽습니다. 제시상 정직하지 않습니다.
##
##  Stolz 는 실제 참가자 3명(모두 여성, 기저 74/81/77세)의 **원측정값**을
##  점으로 찍었습니다("Points are raw observations of IC"). 궤적이 들쭉날쭉한
##  것이 정상이고, 그 불규칙성이 곧 동적예측이 갱신되는 이유입니다.
##
##  "matched" : 실제 참가자 3명을 고르되 성별을 하나로 고정하고 연령·동반질환
##              창을 좁혀 교란을 줄입니다(창은 후보가 3명 이상 나올 때까지
##              단계적으로 넓히고, 무엇을 넓혔는지 콘솔에 남깁니다).
##              선택된 3명의 나이·성별·동반질환·측정횟수를 각주에 명시합니다.
##  "synthetic" : 가상 3명. 대비는 깨끗하지만 점이 직선이 됩니다.
##  "observed"  : 예전 방식(느슨한 제한).
###############################################################################
###############################################################################
## ★ MAKE_FIGURE 와 PROFILE_MODE 는 반드시 짝을 맞춰야 합니다 (2026-08-01)
##
##  이 스크립트는 기본적으로 그림을 그리지 않고 예측값(CSV)만 만들며,
##  본문 Figure 2 는 KF_Fig2ALT 가 그 CSV 를 읽어 그립니다.
##  그런데 ALT 는 **IC 궤적을 그리지 않는** 4패널 중첩곡선입니다. 즉
##  "점이 직선이라 가짜 같다" 는 문제가 아예 없습니다. 반대로 실제 참가자를
##  넣으면 ALT 가 세 개인의 예측을 그리게 되어, 그림의 주장(IC 차이)이 아니라
##  세 사람의 우연한 차이를 보여주게 됩니다. 실제로 그렇게 나왔습니다.
##
##    MAKE_FIGURE = FALSE (기본) -> ALT 가 본문 그림 -> PROFILE_MODE = synthetic
##    MAKE_FIGURE = TRUE         -> 개인 궤적 그림   -> PROFILE_MODE = matched
##
##  아래에서 자동으로 맞춥니다. 수동으로 바꾸려면 PROFILE_MODE 를 직접 쓰십시오.
###############################################################################
MAKE_FIGURE  <- FALSE                                     # 개인 궤적 그림 생성 여부
PROFILE_MODE <- if (MAKE_FIGURE) "matched" else "synthetic"
MATCH_MIN_OBS <- 4            # 궤적이 보이려면 최소 이만큼의 IC 측정이 필요
MATCH_AGE_W   <- c(2, 3, 5, 8)      # 연령 창(세)을 이 순서로 넓힙니다
MATCH_CMB_W   <- c(1, 2)            # 동반질환 창
PROFILE_Q    <- c(0.10, 0.50, 0.90)
TMAX_Q       <- 0.98         # 예측 가능 상한을 추적기간 분포의 몇 분위로 볼지
set.seed(CFG$SEED)

load_imputed("MAIN")
d <- prep_long(load_kfacs())

## ── 종단 / 생존 데이터 ──────────────────────────────────────────────────
long <- d[, c("id","time","gLIC","gLIC_z","age0","sex_f","comorbid_bl")]
long <- long[is.finite(long$time) & is.finite(long$gLIC_z), ]
long <- long[order(long$id, long$time), ]

sv <- do.call(rbind, lapply(split(d[order(d$id, d$time), ], d$id), function(p) {
  n <- nrow(p)
  data.frame(id = p$id[1], t = p$followup_years[n], ev = as.integer(p$death_event[n] == 1),
             age0 = p$age0[1], sex_f = p$sex_f[1], comorbid_bl = p$comorbid_bl[1])
}))
sv <- sv[is.finite(sv$t) & sv$t > 0 & !is.na(sv$ev), ]
sv <- sv[sv$id %in% unique(long$id), ]
long <- merge(long, sv[, c("id","t")], by = "id")
long <- long[long$time <= long$t, ]
long <- long[order(long$id, long$time), ]
sv <- sv[sv$id %in% unique(long$id), ]; sv <- sv[order(sv$id), ]
msg(sprintf("JM: N=%d, deaths=%d, long obs=%d", nrow(sv), sum(sv$ev), nrow(long)))

## ── 예측 가능 상한(TMAX) 과 랜드마크 확정 ───────────────────────────────
TMAX_ABS <- max(sv$t, na.rm = TRUE)
TMAX     <- as.numeric(stats::quantile(sv$t, TMAX_Q, na.rm = TRUE))
TMAX     <- min(TMAX, max(sv$t[sv$ev == 1], na.rm = TRUE))   # 마지막 사망시점 이내
## 관측 wave 격자 (예측에 쓸 수 있는 후보 시점)
WAVE_T <- sort(unique(round(long$time, 3)))
CAND   <- WAVE_T[WAVE_T + MIN_WINDOW <= TMAX]
LANDMARKS <- LANDMARKS_REQ[LANDMARKS_REQ %in% WAVE_T & LANDMARKS_REQ + MIN_WINDOW <= TMAX]
if (!length(LANDMARKS)) LANDMARKS <- if (length(CAND)) CAND[c(1, length(CAND))] else 0
notwave <- setdiff(LANDMARKS_REQ, WAVE_T)
if (length(notwave))
  message("[Fig2] 측정 wave 가 아닌 랜드마크는 제외했습니다(새 정보 없음): ",
          paste(notwave, collapse = ", "), "\n  사용 가능한 wave 시점: ",
          paste(WAVE_T, collapse = ", "), "년")
HORIZON   <- min(HORIZON_REQ, TMAX - min(LANDMARKS))   # 가장 긴 예측구간
if (!is.finite(HORIZON) || HORIZON < MIN_WINDOW) HORIZON <- MIN_WINDOW

dropped <- setdiff(LANDMARKS_REQ, LANDMARKS)
msg(sprintf("추적기간: 최대 %.1f년, 예측 상한(TMAX, q%.2f) %.1f년", TMAX_ABS, TMAX_Q, TMAX))
msg(sprintf("wave 시점 = %s년 | 랜드마크로 쓸 수 있는 시점 = %s년",
            paste(WAVE_T, collapse = ", "), paste(CAND, collapse = ", ")))
msg(sprintf("선택된 랜드마크 = %s (%s) | 예측: 각 랜드마크 -> %.1f년%s",
            paste(LANDMARKS, collapse = ", "),
            paste(sprintf("W%d", round(LANDMARKS / (if (exists("CFG") && is.finite(CFG$WAVE_GAP)) CFG$WAVE_GAP else 2)) + 1), collapse = ", "),
            TMAX,
            if (length(dropped))
              sprintf("  [제외됨: %s -- 남는 구간이 %.0f년 미만]",
                      paste(dropped, collapse = ", "), MIN_WINDOW) else ""))

## ── 적합 ────────────────────────────────────────────────────────────────
lmeFit <- try(nlme::lme(gLIC_z ~ time + age0 + sex_f, random = ~ time | id, data = long,
                        control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                                   msMaxIter = 200, returnObject = TRUE)),
              silent = TRUE)
if (inherits(lmeFit, "try-error")) {
  message("[JM] 랜덤기울기 실패 -> 랜덤절편")
  lmeFit <- nlme::lme(gLIC_z ~ time + age0 + sex_f, random = ~ 1 | id, data = long,
                      control = nlme::lmeControl(opt = "optim", returnObject = TRUE))
}
coxFit <- survival::coxph(survival::Surv(t, ev) ~ age0 + sex_f + comorbid_bl,
                          data = sv, x = TRUE, model = TRUE)

## VALUE-ONLY 를 명시 고정 (slope 가 섞이지 않도록)
## v260817: 적합 객체를 캐시하고 (재실행 시 재적합 생략), Rhat 수렴진단을 기록
JM_CACHE <- file.path(FIG_DIR, "F2_jmFit.rds")
if (file.exists(JM_CACHE)) {
  message("[JM] 캐시 발견 -> 재적합 생략: ", basename(JM_CACHE),
          "  (재적합하려면 이 파일을 지우십시오)")
  jmFit <- readRDS(JM_CACHE)
} else {
  jmFit <- try(JMbayes2::jm(coxFit, list(lmeFit), time_var = "time",
                            functional_forms = list("gLIC_z" = ~ value(gLIC_z)),
                            n_iter = 12000L, n_burnin = 4000L, n_chains = 3L,
                            n_thin = 5L, seed = CFG$SEED), silent = TRUE)
  if (inherits(jmFit, "try-error")) {
    message("[JM] functional_forms 실패 -> 기본(value) 설정")
    jmFit <- JMbayes2::jm(coxFit, list(lmeFit), time_var = "time",
                          n_iter = 12000L, n_burnin = 4000L, n_chains = 3L,
                          n_thin = 5L, seed = CFG$SEED)
  }
  saveRDS(jmFit, JM_CACHE)
  message("[JM] 적합 객체 저장: ", basename(JM_CACHE))
}
sm <- summary(jmFit)

## ── 수렴진단 (potential scale reduction factor, Rhat) ────────────────────
.rhat_all <- tryCatch({
  rr <- c()
  for (nm in names(sm)) {
    tb <- sm[[nm]]
    if (is.matrix(tb) && "Rhat" %in% colnames(tb)) rr <- c(rr, tb[, "Rhat"])
    else if (is.data.frame(tb) && "Rhat" %in% names(tb)) rr <- c(rr, tb$Rhat)
  }
  rr[is.finite(rr)]
}, error = function(e) numeric(0))
if (length(.rhat_all)) {
  .rh_max <- max(.rhat_all); .rh_med <- stats::median(.rhat_all)
  message(sprintf("[JM 수렴진단] Rhat 최대 %.3f · 중앙 %.3f (모수 %d개) -> %s",
                  .rh_max, .rh_med, length(.rhat_all),
                  if (.rh_max < 1.1) "수렴 양호 (<1.1)"
                  else "★ 1.1 초과 -- n_iter 를 늘려 재적합하십시오"))
  writeLines(c(sprintf("JM convergence diagnostics (%s)", format(Sys.Date())),
               sprintf("max Rhat    = %.4f", .rh_max),
               sprintf("median Rhat = %.4f", .rh_med),
               sprintf("n parameters monitored = %d", length(.rhat_all)),
               if (.rh_max < 1.1)
                 "Verdict: convergence satisfactory (all Rhat < 1.1)"
               else "Verdict: NOT converged -- increase n_iter and refit"),
             file.path(FIG_DIR, "F2_JM_convergence.txt"))
} else message("[JM 수렴진단] Rhat 컬럼을 찾지 못했습니다 -- summary 구조를 확인하십시오.")
ASSOC <- data.frame(term = rownames(sm$Survival), sm$Survival, check.names = FALSE)
ASSOC$HR <- exp(ASSOC$Mean); ASSOC$lo <- exp(ASSOC[["2.5%"]]); ASSOC$hi <- exp(ASSOC[["97.5%"]])
save_vals(ASSOC, "F2_JM_association.csv", FIG_DIR)
print(ASSOC[, c("term","HR","lo","hi", intersect("Rhat", names(ASSOC)))])

## ── 프로파일 정의 ──────────────────────────────────────────────────────
.lev <- function(v) if (is.factor(v)) levels(v) else sort(unique(v))
.mode <- function(v) { t <- table(v); names(t)[which.max(t)] }
AGE0_REF    <- mean(sv$age0, na.rm = TRUE)
COMORB_REF  <- mean(sv$comorbid_bl, na.rm = TRUE)
SEX_REF     <- factor(.mode(sv$sex_f), levels = .lev(sv$sex_f))
B_TIME      <- tryCatch(unname(nlme::fixef(lmeFit)["time"]), error = function(e) 0)
if (!is.finite(B_TIME)) B_TIME <- 0
Z_Q         <- as.numeric(stats::quantile(long$gLIC_z[long$time == min(long$time)],
                                          PROFILE_Q, na.rm = TRUE))
###############################################################################
## ★ 행 라벨 (2026-08-01 수정)
##  실제 참가자를 쓰면 "Low IC / Median IC / High IC" 라는 이름을 쓸 수 없습니다.
##  그 이름은 **wave 1 시점의 순위**일 뿐인데, 실제 사람의 IC 는 서로 교차합니다.
##  실제로 이 자료에서 'Median' 참가자의 IC 가 W2 이후 급락해 'Low' 참가자보다
##  낮아졌고, 그 결과 예측생존도 W2 랜드마크부터 역전됩니다.
##  이름을 그대로 두면 독자는 그림이 고장난 줄 압니다.
##  참고문헌(Stolz)도 같은 이유로 Participant A / B / C 로만 부릅니다.
###############################################################################
PROF_LAB <- if (identical(PROFILE_MODE, "synthetic"))
  c("a  Low IC (10th)", "b  Median IC", "c  High IC (90th)") else
  c("a  Participant A", "b  Participant B", "c  Participant C")

.mk_id <- function(pid)
  if (is.factor(long$id)) factor(pid, levels = c(levels(long$id), pid)) else as.character(pid)

if (identical(PROFILE_MODE, "synthetic")) {
  PARTS <- data.frame(id = paste0("PROFILE_", c("LOW","MID","HIGH")),
                      profile = PROF_LAB, z0 = Z_Q, stringsAsFactors = FALSE)
  msg(sprintf("프로파일(가상): gLIC_z 시작값 %s | 공변량 고정 age0=%.1f, sex=%s, comorbid=%.1f | 기울기 %.3f/년",
              paste(sprintf("%.2f", Z_Q), collapse = ", "),
              AGE0_REF, as.character(SEX_REF), COMORB_REF, B_TIME))
} else if (identical(PROFILE_MODE, "matched")) {
  ###########################################################################
  ## 실제 참가자 3명 — 성별 고정 + 연령·동반질환 창을 단계적으로 넓혀 선택
  ###########################################################################
  nobs   <- table(long$id)
  medage <- stats::median(sv$age0, na.rm = TRUE)
  medcmb <- stats::median(sv$comorbid_bl, na.rm = TRUE)
  base_ok <- sv$id %in% names(nobs[nobs >= MATCH_MIN_OBS]) &
             sv$t >= max(LANDMARKS) + MIN_WINDOW &
             as.character(sv$sex_f) == as.character(SEX_REF)
  elig <- character(0); usedA <- NA_real_; usedC <- NA_real_
  for (cw in MATCH_CMB_W) for (aw in MATCH_AGE_W) {
    k <- base_ok & abs(sv$age0 - medage) <= aw & abs(sv$comorbid_bl - medcmb) <= cw
    if (sum(k, na.rm = TRUE) >= 10) { elig <- sv$id[which(k)]; usedA <- aw; usedC <- cw; break }
  }
  if (!length(elig)) {                       # 그래도 부족하면 성별 제한만 유지
    elig <- sv$id[which(base_ok)]; usedA <- Inf; usedC <- Inf
  }
  if (length(elig) < 3) elig <- sv$id[which(sv$t >= max(LANDMARKS) + MIN_WINDOW &
                                            sv$id %in% names(nobs[nobs >= 3]))]
  z1 <- vapply(elig, function(i) {
    v <- long$gLIC_z[long$id == i]; if (length(v)) v[1] else NA_real_ }, numeric(1))
  ok <- is.finite(z1); elig <- elig[ok]; z1 <- z1[ok]
  ## 후보 분포의 10/50/90 백분위에 '가장 가까운' 사람을 고릅니다
  tgt  <- stats::quantile(z1, PROFILE_Q, na.rm = TRUE)
  pick <- integer(0)
  for (q in tgt) {
    cand <- setdiff(seq_along(z1), pick)
    pick <- c(pick, cand[which.min(abs(z1[cand] - q))])
  }
  PARTS <- data.frame(id = elig[pick], profile = PROF_LAB[seq_along(pick)],
                      z0 = z1[pick], stringsAsFactors = FALSE)
  info <- sv[match(PARTS$id, sv$id), c("age0","sex_f","comorbid_bl","t","ev")]
  PARTS$age0 <- info$age0; PARTS$sex <- as.character(info$sex_f)
  PARTS$comorbid_bl <- info$comorbid_bl; PARTS$fu <- info$t; PARTS$died <- info$ev
  PARTS$nobs <- as.integer(nobs[as.character(PARTS$id)])
  msg(sprintf("프로파일(실제·매칭): 후보 %d명 (성별 %s 고정, 연령 ±%s세, 동반질환 ±%s)",
              length(elig), as.character(SEX_REF),
              ifelse(is.finite(usedA), usedA, "제한없음"),
              ifelse(is.finite(usedC), usedC, "제한없음")))
  print(PARTS[, c("profile","id","z0","age0","sex","comorbid_bl","nobs","fu","died")],
        row.names = FALSE, digits = 3)
} else {
  ## observed 모드: 연령·동반질환이 중앙값 근처인 사람만 후보로 (교란 축소)
  nobs <- table(long$id)
  near <- abs(sv$age0 - stats::median(sv$age0, na.rm = TRUE)) <= 2 &
          abs(sv$comorbid_bl - stats::median(sv$comorbid_bl, na.rm = TRUE)) <= 1
  elig <- sv$id[near & sv$t >= max(LANDMARKS) + MIN_WINDOW &
                sv$id %in% names(nobs[nobs >= 4])]
  if (length(elig) < 3)
    elig <- sv$id[sv$t >= max(LANDMARKS) + MIN_WINDOW & sv$id %in% names(nobs[nobs >= 4])]
  if (length(elig) < 3) elig <- sv$id[sv$t > max(LANDMARKS)]
  z <- vapply(elig, function(i) long$gLIC_z[long$id == i][1], numeric(1))
  o <- elig[order(z)]
  idx <- unique(pmax(1, round(PROFILE_Q * length(o))))
  PARTS <- data.frame(id = o[idx], profile = PROF_LAB[seq_along(idx)],
                      z0 = z[order(z)][idx], stringsAsFactors = FALSE)
  msg(sprintf("프로파일(실제 참가자): %s", paste(PARTS$id, collapse = ", ")))
}
PARTS$profile <- factor(PARTS$profile, levels = PARTS$profile)
save_vals(PARTS, "F2_participants.csv", FIG_DIR)

## ── 예측용 newdata 구성 ────────────────────────────────────────────────
## 가상 프로파일: 관측 격자(0, gap, 2*gap, ... <= landmark)에서 IC 를
##   z(t) = z0 + (혼합모형의 평균 시간기울기) * t  로 둡니다.
##   세 프로파일이 평행하므로 차이는 오직 '수준'에서만 옵니다
##   (value-only JM 이 재는 것이 바로 이 수준입니다).
GAP <- if (!is.null(CFG$WAVE_GAP) && is.finite(CFG$WAVE_GAP)) CFG$WAVE_GAP else 2
mk_newdata <- function(k, lm_t) {
  if (identical(PROFILE_MODE, "synthetic")) {
    tt <- seq(0, lm_t, by = GAP); if (!length(tt)) tt <- 0
    data.frame(id = .mk_id(PARTS$id[k]), time = tt,
               gLIC_z = PARTS$z0[k] + B_TIME * tt,
               age0 = AGE0_REF, sex_f = SEX_REF, comorbid_bl = COMORB_REF,
               stringsAsFactors = FALSE)
  } else {
    long[long$id == PARTS$id[k] & long$time <= lm_t,
         c("id","time","gLIC_z","age0","sex_f","comorbid_bl"), drop = FALSE]
  }
}

## ── 동적예측 (★ t / ev 는 Cox 적합에 쓰인 생존변수라 반드시 필요) ───────
dynpred <- function(k, lm_t) {
  nd <- mk_newdata(k, lm_t)
  if (is.null(nd) || !nrow(nd)) return(NULL)
  nd$t  <- max(nd$time, na.rm = TRUE)   # 마지막 event-free 시점
  nd$ev <- 0                            # 아직 사건 없음
  ## ★ 예측 상한(TMAX) 을 넘기면 JM 이 NA 를 반환합니다 -> 반드시 잘라서 요청
  t0 <- max(nd$t); t1 <- min(t0 + HORIZON, TMAX)
  if (t1 - t0 < 0.5) { message("[dynpred] 예측구간 부족 k=", k, " lm=", lm_t); return(NULL) }
  tseq <- seq(t0, t1, length.out = 26)
  pr <- try(predict(jmFit, newdata = nd, process = "event",
                    times = tseq, return_newdata = TRUE), silent = TRUE)
  if (inherits(pr, "try-error")) {
    message("[dynpred] 실패 k=", k, " lm=", lm_t, " :: ",
            conditionMessage(attr(pr, "condition"))); return(NULL)
  }
  o <- as.data.frame(pr); cn <- names(o)
  g <- function(a, b) if (a %in% cn) a else grep(b, cn, value = TRUE)[1]
  data.frame(id = PARTS$id[k], landmark = lm_t,
             time = as_num(o[[g("times","^time")]]),
             cif  = as_num(o[[g("pred_CIF","^pred")]]),
             lo   = as_num(o[[g("low_CIF","^low")]]),
             hi   = as_num(o[[g("upp_CIF","^upp")]]), stringsAsFactors = FALSE)
}
PRED <- do.call(rbind, unlist(lapply(seq_len(nrow(PARTS)), function(k)
  lapply(LANDMARKS, function(l) dynpred(k, l))), recursive = FALSE))

## ★ 가상 프로파일 예측이 통째로 실패하면 실제 참가자 방식으로 자동 후퇴
if ((is.null(PRED) || !nrow(PRED)) && identical(PROFILE_MODE, "synthetic"))
  stop("가상 프로파일 예측이 모두 실패했습니다.\n",
       "  -> 파일 상단 PROFILE_MODE 를 \"observed\" 로 바꾸고 다시 실행하십시오.\n",
       "     (위 [dynpred] 실패 메시지에 원인이 찍혀 있습니다)")

PRED <- PRED[is.finite(PRED$time) & is.finite(PRED$cif), ]

###############################################################################
## ★ 순서 점검 — 실제 참가자를 쓰면 순서가 뒤집힐 수 있습니다.
##  가상 프로파일과 달리 실제 사람은 IC 외의 것도 다릅니다. 뒤집혔다면
##  그것은 오류가 아니라 사실이지만, 그림의 메시지와 어긋나므로 반드시
##  알고 넘어가야 합니다. 멈추지 않고 경고만 냅니다.
###############################################################################
CROSS_NOTE <- ""
if (nrow(PRED)) {
  E <- do.call(rbind, lapply(split(PRED, list(PRED$id, PRED$landmark), drop = TRUE),
    function(p) data.frame(id = p$id[1], lm = p$landmark[1],
                           surv_end = 1 - p$cif[which.max(p$time)],
                           stringsAsFactors = FALSE)))
  E$profile <- as.character(PARTS$profile[match(E$id, PARTS$id)])
  E$ord     <- match(E$profile, levels(PARTS$profile))
  cat("\n[Fig2] 추적 종료시점 예측생존 (랜드마크별)\n")
  for (lm in sort(unique(E$lm))) {
    s <- E[E$lm == lm, ]; s <- s[order(s$ord), ]
    cat(sprintf("   랜드마크 %.0f년 : %s %s\n", lm,
        paste(sprintf("%s %.3f", sub("^[a-z]\\s+", "", s$profile), s$surv_end), collapse = " | "),
        if (all(diff(s$surv_end) > 0)) "" else "  <-- wave-1 순서와 역전"))
  }
  ## 교차가 생겼으면 각주 문장을 자동 생성합니다 (숨기지 않고 설명합니다)
  cr <- vapply(sort(unique(E$lm)), function(lm) {
    s <- E[E$lm == lm, ]; s <- s[order(s$ord), ]; !all(diff(s$surv_end) > 0) }, logical(1))
  if (any(cr) && !identical(PROFILE_MODE, "synthetic")) {
    lm1 <- sort(unique(E$lm))[which(cr)[1]]
    CROSS_NOTE <- paste0(
      "The ordering of the predicted curves changes from the ",
      sprintf("%.0f", lm1), "-year landmark onwards. This is the point of a dynamic ",
      "prediction rather than a defect: the participants are ranked by their capacity at ",
      "wave 1, but the model uses the capacity measured most recently, so a participant ",
      "whose capacity falls steeply overtakes one who started lower and remained stable. ")
    msg("[Fig2] 랜드마크 ", lm1, "년부터 순서 역전 -> 각주에 설명 문장을 자동 삽입합니다.")
  }
}
PRED$surv    <- 1 - PRED$cif
PRED$surv_lo <- 1 - PRED$hi
PRED$surv_hi <- 1 - PRED$lo
PRED <- merge(PRED, PARTS[, c("id","profile")], by = "id")
save_vals(PRED, "F2_dynamic_predictions.csv", FIG_DIR)

## ★ 곡선이 실제로 어디까지 그려지는지 확인 (끊김 재발 감지)
chk <- do.call(rbind, lapply(split(PRED, list(PRED$landmark, PRED$id), drop = TRUE),
  function(z) data.frame(landmark = z$landmark[1], id = z$id[1],
                         from = min(z$time), to = max(z$time))))
print(chk[order(chk$landmark, chk$id), ], row.names = FALSE)
short <- chk[chk$to < TMAX - 0.3, ]
if (nrow(short))
  warning("추적 상한(", sprintf("%.1f", TMAX), "년)까지 도달하지 못한 곡선이 ",
          nrow(short), "개 있습니다. [dynpred] 메시지를 확인하십시오.")

## IC 궤적 (좌축) — 예측에 실제로 투입한 그 값을 그립니다.
## ★ 좌축 범위: 전체 range 를 쓰면 극단치 때문에 궤적이 납작해집니다(1차 렌더 문제).
##   5-95 백분위수를 기본으로 하고, 표시할 궤적은 반드시 포함되게 넓힙니다.
OBS <- do.call(rbind, lapply(seq_len(nrow(PARTS)), function(k) {
  z <- mk_newdata(k, max(LANDMARKS))
  if (is.null(z) || !nrow(z)) return(NULL)
  data.frame(id = PARTS$id[k], time = z$time, gLIC_z = z$gLIC_z,
             stringsAsFactors = FALSE)
}))
qz   <- stats::quantile(long$gLIC_z, c(0.05, 0.95), na.rm = TRUE)
rngz <- range(c(qz, OBS$gLIC_z), na.rm = TRUE)
rngz <- rngz + c(-1, 1) * 0.08 * diff(rngz)          # 여백
OBS <- merge(OBS, PARTS[, c("id","profile")], by = "id")
OBS <- do.call(rbind, lapply(LANDMARKS, function(l) { x <- OBS; x$landmark <- l; x }))
OBS$shown <- OBS$time <= OBS$landmark
save_vals(OBS, "F2_observed_IC.csv", FIG_DIR)

###############################################################################
## ── 그림 ────────────────────────────────────────────────────────────────
## ★ MAKE_FIGURE = FALSE 가 기본입니다 (2026-08-01).
##
##  이 스크립트는 이제 **조인트모형 적합과 예측값 산출까지만** 합니다.
##  본문 Figure 2 는 KF_Fig2ALT_JMdynpred_v260801.R 이 이 결과(CSV)를 읽어
##  그립니다. 개인 궤적 격자그림은 우리 자료로는 성립하지 않기 때문입니다:
##
##    참고문헌(Stolz)  1인당 IC 측정 최대 13회 / 18개월 간격 / 20년
##    본 연구          1인당 4-5회          / 2년 간격      / 8.2년
##
##  점 4-5개로는 개인 궤적이 사실상 측정오차입니다. 실제로 선택된 참가자
##  B 의 IC 가 6년간 1.65 SD 급락했다 반등했고, 그 결과 예측생존 순서가
##  W2 랜드마크부터 뒤집혔습니다(각주로 설명은 가능하나, 독자는 그림을
##  먼저 봅니다). 가상 프로파일로 바꾸면 점이 완벽한 직선이 되어 자료로
##  보이지 않습니다. 어느 쪽도 고칠 수 있는 문제가 아닙니다.
##
##  개인 궤적 그림이 다시 필요하면 아래를 TRUE 로 두십시오. 그러면
##  Figure2_indiv_JM_dynamic_predictions 라는 별도 이름으로 저장되어
##  ALT 가 만든 Figure 2 를 덮어쓰지 않습니다.
###############################################################################
if (!MAKE_FIGURE) {
  cat("\n=== 조인트모형 적합·예측 완료 (그림은 그리지 않음) ===\n")
  cat("  예측값 : ", file.path(FIG_DIR, "F2_dynamic_predictions.csv"), "\n", sep = "")
  cat("  관측 IC: ", file.path(FIG_DIR, "F2_observed_IC.csv"), "\n", sep = "")
  cat("  -> 본문 Figure 2 는 KF_Fig2ALT_JMdynpred_v260801.R 이 그립니다.\n")
  cat(">> F2_JM_association.csv 의 value(gLIC_z) HR 과 R-hat 을 반드시 확인하십시오.\n")
} else {

PRED$profile  <- factor(PRED$profile, levels = PARTS$profile)
OBS$profile   <- factor(OBS$profile,  levels = PARTS$profile)
## 열 라벨: 열이 4개면 짧게 (겹침 방지). wave 번호를 함께 적어 독자가 바로 압니다.
lab_lm <- if (length(LANDMARKS) >= 4)
  function(x) sprintf("W%d (%g y)", round(x / GAP) + 1, x) else
  function(x) sprintf("Landmark: wave %d (%g y)", round(x / GAP) + 1, x)
PRED$lmf <- factor(lab_lm(PRED$landmark), levels = lab_lm(LANDMARKS))
OBS$lmf  <- factor(lab_lm(OBS$landmark),  levels = lab_lm(LANDMARKS))

## ★ 색은 KF_theme 규칙을 따릅니다.
##   이 그림의 예측선은 '사망' 결과이므로 주홍(PAL_OUTCOME["Mortality"]) 입니다.
##   Suppl Fig 1 의 Death 행, Suppl Fig 2 의 Mortality 패널과 같은 색입니다.
COL_MORT <- PAL_OUTCOME[["Mortality"]]

## ★ x축: 모든 패널이 0 ~ TMAX 를 공유합니다(랜드마크 간 비교를 위해).
XMAX <- TMAX
## 열이 4개면 패널 폭이 45 mm 로 좁아집니다 -> 눈금을 성기게 해 라벨 겹침을 막습니다
.step <- if (XMAX > 12) 4 else if (length(LANDMARKS) >= 4) 4 else 2
XBRK  <- seq(0, floor(XMAX), by = .step)

## ★ 우축(생존확률): 예측 생존이 0.85-1.00 구간에 몰려 0-1 축에서는 차이가 안 보입니다.
##   관측된 하한에 맞춰 축을 잘라 확대하고, 각주에 잘린 사실을 명시합니다.
ZOOM_SURV <- TRUE      # FALSE 로 두면 생존축을 0-1 전체로 되돌립니다
S_LO <- if (ZOOM_SURV) min(PRED$surv_lo, na.rm = TRUE) else 0
S_LO <- max(0, min(0.9, floor(S_LO * 20) / 20))      # 0.05 단위로 내림
S_BRK <- pretty(c(S_LO, 1), n = 4); S_BRK <- S_BRK[S_BRK >= S_LO & S_BRK <= 1]

###############################################################################
## ★ 패널 아래가 비던 문제 — 진짜 원인 (2026-07-31 재수정)
##
##  1차 시도: 패널을 IC 띠 / 생존 띠 로 나눔  -> 실패.
##            축 눈금이 두 뭉치로 뭉쳐 오히려 그림이 망가졌습니다.
##
##  실제 원인: 예측지평이 4년으로 짧아 생존이 1.00 -> 0.88 밖에 안 내려가고,
##            곡선이 패널 위쪽에만 붙어 있었던 것입니다.
##            참고문헌(Stolz, J Gerontol 2022, Fig 3)은 각 랜드마크에서
##            **추적 종료 시점까지** 예측해 곡선이 패널을 가로지릅니다.
##
##  수정: 지평을 고정하지 않고 랜드마크 -> TMAX(추적 상한) 까지 예측합니다.
##        랜드마크가 뒤로 갈수록 예측구간이 짧아지는데, 이것이 바로
##        "정보가 쌓일수록 예측이 갱신된다"는 이 그림의 메시지입니다.
###############################################################################
to01 <- function(z) pmin(1, pmax(0, (z - rngz[1]) / diff(rngz)))
to_s <- function(s) pmin(1, pmax(0, (s - S_LO) / (1 - S_LO)))

Z_BRK <- pretty(rngz, n = 4); Z_BRK <- Z_BRK[Z_BRK >= rngz[1] & Z_BRK <= rngz[2]]

msg(sprintf("Fig2 축: x 0-%.1f년 | 생존축 %.2f-1.00 | IC축 %.2f~%.2f",
            XMAX, S_LO, rngz[1], rngz[2]))

fig <- ggplot() +
  geom_ribbon(data = PRED, aes(time, ymin = to_s(surv_lo), ymax = to_s(surv_hi)),
              fill = COL_MORT, alpha = 0.18) +
  geom_line(data = PRED, aes(time, to_s(surv)), colour = COL_MORT, linewidth = 0.5) +
  geom_point(data = subset(OBS, shown), aes(time, to01(gLIC_z)),
             size = 1.1, colour = "black") +
  geom_line(data = subset(OBS, shown), aes(time, to01(gLIC_z)),
            linewidth = LINE_PT, colour = "black") +
  geom_vline(aes(xintercept = landmark), data = unique(PRED[, c("landmark","lmf","profile")]),
             linetype = "22", linewidth = LINE_PT, colour = SEV_COL[["aux"]]) +
  facet_grid(profile ~ lmf, switch = "y") +
  scale_y_continuous(name = "Intrinsic capacity (gLIC z-score)",
                     limits = c(0, 1),
                     breaks = to01(Z_BRK), labels = sprintf("%.1f", Z_BRK),
                     sec.axis = sec_axis(~ ., name = "Predicted survival probability",
                                         breaks = to_s(S_BRK),
                                         labels = sprintf("%.2f", S_BRK))) +
  scale_x_continuous(limits = c(0, XMAX), breaks = XBRK, expand = expansion(mult = 0.02)) +
  labs(x = "Years from baseline") +
  theme_na() +
  theme(strip.placement   = "outside",
        strip.text.y.left  = element_text(angle = 0, hjust = 0, face = "bold"),
        panel.spacing.x    = unit(1.6, "mm"),
        panel.spacing.y    = unit(1.6, "mm"),
        axis.title.y.right = element_text(colour = COL_MORT),
        axis.text.y.right  = element_text(colour = COL_MORT),
        axis.ticks.y.right = element_line(colour = COL_MORT, linewidth = LINE_PT))

CAP_PROFILE <- if (identical(PROFILE_MODE, "synthetic")) {
  paste0("Rows are three hypothetical participants who differ only in intrinsic capacity ",
         "(gLIC z-score at the ", paste(sprintf("%.0fth", PROFILE_Q * 100), collapse = ", "),
         " centile of the wave-1 distribution, thereafter following the cohort-average ",
         "trajectory); age, sex and baseline comorbidity are held at cohort reference values ",
         "so that any difference between rows is attributable to IC alone. Black points and ",
         "lines: the IC trajectory supplied to the model up to the landmark ",
         "(dashed vertical line). ")
} else if (identical(PROFILE_MODE, "matched")) {
  paste0("Participants A, B and C are three real participants of the same sex (",
         tolower(as.character(SEX_REF)), ") and of similar age and baseline comorbidity, ",
         "selected to lie near the ",
         paste(sprintf("%.0fth", PROFILE_Q * 100), collapse = ", "),
         " centile of the wave-1 distribution of intrinsic capacity respectively. They were ",
         "aged ", paste(sprintf("%.0f", PARTS$age0), collapse = ", "),
         " years at baseline, with ", paste(sprintf("%.0f", PARTS$comorbid_bl), collapse = ", "),
         " chronic conditions and ", paste(PARTS$nobs, collapse = ", "),
         " repeated measurements. Black points are raw observations of intrinsic capacity, ",
         "joined by a line, up to the landmark (dashed vertical line); their irregularity is ",
         "the information the joint model uses to revise each prediction. ", CROSS_NOTE)
} else {
  paste0("Rows are three participants at the ",
         paste(sprintf("%.0fth", PROFILE_Q * 100), collapse = ", "),
         " centile of baseline IC, restricted to participants near the cohort median for age ",
         "and comorbidity. Black points and lines: observed IC trajectory up to the landmark ",
         "(dashed vertical line). ")
}

fig <- fig + annot_na(
  title = "Dynamic predictions of survival from the joint model (current IC level only)",
  caption = paste0(CAP_PROFILE,
                   "Vermillion line and band: predicted survival probability ",
                   "(95% credible interval) from the landmark to the end of follow-up (",
                   sprintf("%.1f", TMAX), " years), so the prediction window shortens as the ",
                   "landmark moves later and more IC measurements accrue. Waves were collected ",
                   "every ", sprintf("%.0f", GAP), " years, so the landmarks shown correspond to ",
                   paste(sprintf("wave %d", round(LANDMARKS / GAP) + 1), collapse = ", "),
                   ". The joint model uses the current value of IC only; no slope term is ",
                   "included and it cannot extrapolate beyond the observed follow-up. ",
                   "Left-hand axis, IC (gLIC z-score, 5th-95th percentile of the observed ",
                   "distribution); right-hand axis, survival probability, truncated at ",
                   sprintf("%.2f", S_LO), " to resolve differences between profiles. N = ",
                   nrow(sv), ", deaths = ", sum(sv$ev), "."))

save_na(fig, "Figure2_indiv_JM_dynamic_predictions", width_mm = NA_WIDTH$double,
        height_mm = 130, dir = FIG_DIR)
cat("\n=== 개인 궤적 그림 완료 (보조) ===\n")
cat(">> F2_JM_association.csv 의 value(gLIC_z) HR 과 R-hat 을 반드시 확인하십시오.\n")
}   # end if (MAKE_FIGURE)
