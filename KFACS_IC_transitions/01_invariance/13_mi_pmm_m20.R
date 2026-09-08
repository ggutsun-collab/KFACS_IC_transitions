###############################################################################
## 13_mi_pmm_m20.R   (Phase 2 — 다중대체 m=20 완성 데이터 생성)   v260903
##
## ── 왜 새로 만드는가 ─────────────────────────────────────────────────────
##  MAIN 세트는 VIM::kNN (k=5, Gower) 단일대체입니다. kNN 은 결정론적이라
##  seed 를 바꿔도 같은 값이 나오므로 "반복" 으로 MI 를 만들 수 없습니다.
##  Seq 세트의 VIM::hotdeck 도 ord_var 를 주면 순차(결정론적)입니다.
##  -> 표준적 해법: mice 의 predictive mean matching, m = 20.
##
## ── 설계 ─────────────────────────────────────────────────────────────────
##  · wide 형식(변수_w1..w5)으로 대체 -> 같은 변수의 다른 wave 값이 예측자가 됩니다.
##    (long 형식 pmm 은 개인 내 상관을 버립니다)
##  · 사망 이후 wave 의 셀도 mice 가 채우지만, long 으로 되돌릴 때 버립니다
##    (prep_long 이 어차피 사망 이후 행을 제거). 대체값이 예측자로만 잠깐 쓰입니다.
##  · 교육/소득/거주지 결측은 원고대로 Unknown 범주 유지 -> 대체하지 않고 예측자로만.
##  · frailty_3cat 은 대체하지 않고 chs_total(0-5) 을 대체한 뒤 파생(0 / 1-2 / 3-5).
##  · 완성 데이터마다 bifactor(특정요인 4개, sensory 는 g 만) 를 재적합해 gLIC 산출
##    — 00a 파이프라인의 postprocess() 와 동일.
##  · state / death 변수는 00a 와 동일한 규칙으로 복구.
##
## 출력: IMP_DIR/KFACS_mi_pmm_m20.rds  — list(MI01, ..., MI20) ; 각 원소는
##       MAIN 세트와 같은 열 구조(load_imputed(set) 로 그대로 읽힘)
##       IMP_DIR/KFACS_mi_pmm_m20_diagnostics.xlsx — 결측률·수렴 진단
##
## 실행: source("13_mi_pmm_m20.R", encoding = "UTF-8")   (30-90분)
###############################################################################

MI_VERSION <- "v260903"
message("\n=== 13_mi_pmm_m20 ", MI_VERSION, " ===")

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
.need(c("mice", "lavaan", "readxl", "openxlsx"))
suppressPackageStartupMessages({ library(mice); library(lavaan) })

###############################################################################
## 1. 설정
###############################################################################
M        <- 20
MAXIT    <- 10
SEED     <- 20260903
## 결측 보존 원본 (00b 가 state/FUP 를 붙인 판). 00_setup 의 PROJECT_ROOT 아래 data/.
PRE_FILE <- get0("PRE_FILE", ifnotfound = file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260618_state.xlsx"))  # 00_setup.R 에서 정함
PRE_SHEET <- "전체_long"
OUT_RDS  <- file.path(IMP_DIR, "KFACS_mi_pmm_m20.rds")

IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
TV_WAVE  <- c(IND17, "adl_score", "iadl_score", "chs_total", "comorbid_count",
              "phys_bmi", "pa_met_min_wk", "polypharmacy")          # wave 별 대체 대상
BL_NUM   <- c("age_unified")                                         # 기저 수치 예측자
BL_CAT   <- c("sex_unified", "educ_bl", "area_3cat", "income_3cat", "living_alone",
              "marri_status", "smoking_current", "strength_2plus", "pa_freq_3cat")
KEEP_ID  <- c("death_event", "death_wave", "followup_years", "death_date_new", "death_time_method")
BIFACTOR <- '
  g  =~ Balance+GS_ms+rev_CST+loss_of_Bwt+Appetite+EXH+HGS+rev_logMAR+rev_PTA+
        Orientation+Memory+Attention_Calculation+Language+Visuospatial+
        Negative_affect+Positive_affect+Motivation
  loco =~ Balance+GS_ms+rev_CST
  vita =~ loss_of_Bwt+Appetite+EXH+HGS
  cogn =~ Orientation+Memory+Attention_Calculation+Language+Visuospatial
  psyc =~ Negative_affect+Positive_affect+Motivation
  g ~~ 0*loco + 0*vita + 0*cogn + 0*psyc '

###############################################################################
## 2. 원본 로드 (한글 경로 회피)
###############################################################################
if (!file.exists(PRE_FILE)) stop("결측 보존 원본이 없습니다: ", PRE_FILE)
tmp <- file.path(tempdir(), "pre.xlsx"); file.copy(PRE_FILE, tmp, overwrite = TRUE)
pre <- as.data.frame(readxl::read_excel(tmp, sheet = PRE_SHEET))
pre$id <- as.character(pre$id); pre$wave <- as.integer(as_num(pre$wave))
msg(sprintf("pre-imp: %d행 x %d열, id %d, wave %s", nrow(pre), ncol(pre),
            length(unique(pre$id)), paste(sort(unique(pre$wave)), collapse = ",")))
need_cols <- c("id", "wave", "FUP", "state", TV_WAVE, BL_NUM, BL_CAT, KEEP_ID)
miss <- setdiff(need_cols, names(pre))
if (length(miss)) stop("원본에 없는 열: ", paste(miss, collapse = ", "))

## 동의철회 3명 제외 (분석 코호트 3,011)
pre <- pre[!pre$id %in% c("kf161224", "kf170602", "kf171295"), ]
ids <- sort(unique(pre$id)); waves <- sort(unique(pre$wave))
msg(sprintf("분석 코호트 id = %d", length(ids)))

## 사망 이후 / 사망 wave 의 측정은 존재하지 않음 -> 대체 대상에서 제외하기 위해 표시
dw <- tapply(pre$death_wave, pre$id, function(x) suppressWarnings(as.numeric(x[1])))
pre$post_death <- !is.na(dw[pre$id]) & pre$wave >= dw[pre$id]

###############################################################################
## 3. wide 구성
###############################################################################
bl <- pre[pre$wave == 1, c("id", BL_NUM, BL_CAT), drop = FALSE]
bl <- bl[!duplicated(bl$id), ]
for (v in BL_CAT) bl[[v]] <- factor(ifelse(is.na(bl[[v]]), "Unknown", as.character(bl[[v]])))
## Unknown 이 실제로 없는 변수는 수준 정리
for (v in BL_CAT) bl[[v]] <- droplevels(bl[[v]])

W <- bl
for (w in waves) for (v in TV_WAVE) {
  x <- pre[pre$wave == w, c("id", v)]; x <- x[!duplicated(x$id), ]
  W[[paste0(v, "_w", w)]] <- as_num(x[[v]])[match(W$id, x$id)]
}
rownames(W) <- W$id
X <- W[, setdiff(names(W), "id")]
## polypharmacy 는 0/1 -> pmm 으로 두어도 무방(정수 유지)
miss_rate <- sort(colMeans(is.na(X)), decreasing = TRUE)
cat("\n결측률 상위 15:\n"); print(round(head(miss_rate, 15), 3))

###############################################################################
## 4. mice — pmm, wide, m=20
###############################################################################
meth <- make.method(X)
meth[BL_CAT] <- ""                                  # Unknown 범주 유지, 대체 안 함
meth[names(meth) %in% names(X)[grepl("_w[0-9]$", names(X))]] <- "pmm"
meth["age_unified"] <- if (any(is.na(X$age_unified))) "pmm" else ""
pred <- quickpred(X, mincor = 0.15, minpuc = 0.2, method = "pearson")
## 기저 공변량은 항상 예측자로 포함
pred[, c(BL_NUM, BL_CAT)] <- 1
diag(pred) <- 0
cat(sprintf("\n예측자 행렬: 변수 %d개, 변수당 평균 예측자 %.1f개\n", ncol(X), mean(rowSums(pred))))

t0 <- Sys.time()
imp <- mice(X, m = M, method = meth, predictorMatrix = pred, maxit = MAXIT,
            seed = SEED, printFlag = TRUE)
msg(sprintf("mice 완료: %.1f분", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## 수렴 진단 (체인 평균의 wave 간 분산 비 — 간이 Rhat)
rh <- tryCatch({
  cm <- imp$chainMean            # var x iter x m
  vapply(dimnames(cm)[[1]], function(v) {
    z <- cm[v, , ]; if (all(is.na(z))) return(NA_real_)
    B <- var(colMeans(z, na.rm = TRUE)); Wv <- mean(apply(z, 2, var, na.rm = TRUE))
    if (!is.finite(Wv) || Wv == 0) return(NA_real_); sqrt((Wv + B) / Wv) }, numeric(1))
}, error = function(e) NULL)
if (!is.null(rh)) cat(sprintf("\n간이 Rhat: 중앙 %.3f · 최대 %.3f (%s)\n",
                              median(rh, na.rm = TRUE), max(rh, na.rm = TRUE), names(which.max(rh))))

###############################################################################
## 5. long 복원 + gLIC 재적합 + state/death 복구 (00a postprocess/restore 동형)
###############################################################################
death_vars <- pre[!duplicated(pre$id), c("id", KEEP_ID)]
fup <- pre[, c("id", "wave", "FUP")]
make_state <- function(d) {
  d$state <- ifelse(d$FUP %in% 4, 4L,
             ifelse(d$FUP %in% c(5, 6), 3L,
             ifelse(!is.na(d$adl_score) & d$adl_score >= 3, 3L,
             ifelse(!is.na(d$adl_score) & d$adl_score >= 1, 2L,
             ifelse(!is.na(d$adl_score) & d$adl_score == 0, 1L, NA_integer_)))))
  d$state_label <- factor(d$state, levels = 1:4, labels = c("Normal","Mild","Severe","Death"))
  d
}
to_long <- function(Xc) {
  Xc$id <- rownames(X)
  out <- do.call(rbind, lapply(waves, function(w) {
    cols <- paste0(TV_WAVE, "_w", w)
    d <- Xc[, c("id", BL_NUM, BL_CAT, cols)]
    names(d)[match(cols, names(d))] <- TV_WAVE
    d$wave <- w; d
  }))
  out <- merge(out, fup, by = c("id", "wave"), all.x = TRUE)
  out <- merge(out, death_vars, by = "id", all.x = TRUE)
  ## 사망 이후 행 제거 (사망 wave 행은 state=4 를 위해 남기되 측정값은 NA)
  dwv <- suppressWarnings(as.numeric(out$death_wave))
  out <- out[is.na(dwv) | out$wave <= dwv, ]
  dwv <- suppressWarnings(as.numeric(out$death_wave))
  is_dw <- !is.na(dwv) & out$wave == dwv
  out[is_dw, c(IND17, "adl_score", "iadl_score", "chs_total")] <- NA
  ## 파생
  out$chs_total    <- pmin(pmax(round(out$chs_total), 0), 5)
  out$frailty_3cat <- ifelse(is.na(out$chs_total), NA, ifelse(out$chs_total == 0, 0,
                             ifelse(out$chs_total <= 2, 1, 2)))
  out$adl_score  <- pmax(round(out$adl_score), 0)
  out$iadl_score <- pmax(round(out$iadl_score), 0)
  out$iadl_2cat  <- as.integer(out$iadl_score >= 1)
  out <- make_state(out)
  ## 기저 범주 예측자를 원 코딩으로 되돌림 (Unknown -> NA; prep_long 이 다시 Unknown 처리)
  for (v in BL_CAT) { z <- as.character(out[[v]]); z[z == "Unknown"] <- NA; out[[v]] <- as_num(z) }
  ## 01_functions_common 이 찾는 이름
  out$edu_high_bl <- out$educ_bl; out$income_high <- as.integer(out$income_3cat == 3)
  out[order(out$id, out$wave), ]
}
score_glic <- function(d, label) {
  ok <- stats::complete.cases(d[, IND17])
  Z <- d[ok, IND17]; Z[] <- lapply(Z, function(x) (x - mean(x)) / stats::sd(x))
  fit <- tryCatch(cfa(BIFACTOR, data = Z, std.lv = TRUE, estimator = "MLR"), error = function(e) NULL)
  if (is.null(fit) || !lavInspect(fit, "converged")) { warning(label, ": bifactor 실패"); d$gLIC <- NA; return(d) }
  g <- as.numeric(lavPredict(fit, newdata = Z)[, "g"])
  if (stats::cor(g, rowMeans(Z), use = "complete.obs") < 0) g <- -g       # 높을수록 좋음
  d$gLIC <- NA_real_; d$gLIC[ok] <- g
  d$gLIC_z <- as.numeric(scale(d$gLIC))
  d$gLIC_tertile <- cut(d$gLIC, stats::quantile(d$gLIC, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                        include.lowest = TRUE, labels = c("T1", "T2", "T3"))
  cat(sprintf("  [%s] bifactor CFI=%.3f RMSEA=%.3f | gLIC N=%d\n", label,
              fitMeasures(fit, "cfi.robust"), fitMeasures(fit, "rmsea.robust"), sum(ok)))
  d
}

MI <- vector("list", M); names(MI) <- sprintf("MI%02d", seq_len(M))
for (i in seq_len(M)) {
  Xc <- complete(imp, i)
  d  <- to_long(Xc)
  MI[[i]] <- score_glic(d, names(MI)[i])
}
## 검증 — MAIN 세트와 열·행 구조 비교
main_rds <- file.path(IMP_DIR, paste0(IMP_STEM, ".rds"))
if (file.exists(main_rds)) {
  mn <- readRDS(main_rds); mn <- if (is.data.frame(mn)) mn else mn$MAIN
  cat(sprintf("\nMAIN: %d행 | MI01: %d행 (사망 이후 행 제거로 MI 가 적습니다 — prep_long 결과는 같아야 함)\n",
              nrow(mn), nrow(MI[[1]])))
  cm <- merge(mn[, c("id", "wave", "gLIC")], MI[[1]][, c("id", "wave", "gLIC")], by = c("id", "wave"))
  cat(sprintf("gLIC MAIN vs MI01 상관 r = %.4f (관측 지표가 대부분이라 0.95 이상이어야 정상)\n",
              stats::cor(cm$gLIC.x, cm$gLIC.y, use = "complete.obs")))
}
saveRDS(MI, OUT_RDS)
msg("저장: ", OUT_RDS)

DIAG <- data.frame(variable = names(miss_rate), missing_pct = round(100 * miss_rate, 2),
                   rhat = if (!is.null(rh)) round(rh[names(miss_rate)], 3) else NA)
openxlsx::write.xlsx(list(missingness = DIAG,
                          settings = data.frame(key = c("m", "maxit", "seed", "method", "mincor", "source"),
                                                value = c(M, MAXIT, SEED, "pmm (wide format)", 0.15, basename(PRE_FILE)))),
                     sub("\\.rds$", "_diagnostics.xlsx", OUT_RDS), overwrite = TRUE)
cat("\n=== 완료 ===  다음: 14_mi_pool_tables.R\n")
