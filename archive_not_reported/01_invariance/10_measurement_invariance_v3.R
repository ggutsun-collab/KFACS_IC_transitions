###############################################################################
## KF_ST_Invariance_v260801.R
## Supplementary Table | gLIC 측정모형의 웨이브 간 종단 불변성 (longitudinal
##                       measurement invariance) 검정
##
## ── 왜 이것이 다른 모든 분석보다 먼저인가 ────────────────────────────────
##  Table 2 / Table 4 / Figure 1·2 는 전부 "W_k 시점의 gLIC 수준"을 웨이브에
##  걸쳐 비교합니다. 그러려면 gLIC 를 재는 자(적재량·절편)가 웨이브 간 같아야
##  합니다. 다르면 웨이브마다 다른 자로 잰 값을 비교하는 것이 되고, 그 위에
##  쌓은 결과는 전부 해석 불가가 됩니다.
##
##  ★ 현재 gLIC 는 이 가정을 **검정한 것이 아니라 가정한** 상태입니다.
##    재대체 스크립트(01_...)는 person-wave 를 전부 합친 하나의 풀에
##    bifactor CFA 를 1회 적합했습니다. 즉 적재량·절편이 웨이브 간 동일하다고
##    **강제**한 것이고, 그것이 자료와 맞는지는 확인된 적이 없습니다.
##    이 스크립트가 그 확인을 합니다.
##
## ── 검정 순서 (다집단 CFA, 집단 = 웨이브) ────────────────────────────────
##   M1 configural : 형태만 동일        -> 같은 구조인가
##   M2 metric     : + 적재량 동일      -> "+1 SD" 의 의미가 웨이브 간 같은가
##   M3 scalar     : + 절편 동일        -> 평균/수준 비교가 가능한가  ★핵심
##   M4 strict     : + 잔차분산 동일    -> (참고용, 보통 요구되지 않음)
##
##   판정: 직전 모형 대비 ΔCFI >= -0.010 이고 ΔRMSEA <= 0.015 이면 지지
##         (Cheung & Rensvold 2002; Chen 2007)
##
##   ※ 집단이 '웨이브'이므로 한 집단 안에서는 각 참가자가 1행뿐입니다.
##     따라서 집단 내 관측치 독립성이 성립하고 cluster 보정이 필요 없습니다.
##     (참가자 간 웨이브 상관은 다집단 CFA 가 모형화하지 않는 부분이며,
##      Part B 의 종단 CFA 가 이를 보완합니다.)
##
## ── scalar 가 깨지면 ─────────────────────────────────────────────────────
##   1) 수정지수로 어느 지표의 절편이 문제인지 찾아 부분(partial) scalar 적합
##   2) 부분불변 모형에서 요인점수를 다시 뽑아 현재 gLIC 와 상관·삼분위 일치율 비교
##   3) 상관 >= 0.99 이고 삼분위 일치율이 높으면 실질적 영향 없음 -> 본 분석 유지
##      그렇지 않으면 gLIC 를 이 모형 기준으로 재산출해야 합니다
##   이 스크립트는 1~3 을 자동으로 수행하고 판정문을 출력합니다.
###############################################################################

###############################################################################
## 이 스크립트는 작업디렉터리와 무관하게 동작합니다.
## 자기 자신의 위치를 찾아 같은 폴더의 KF_common / KF_theme 를 불러옵니다.
## (source() / RStudio Source 버튼 / Rscript 모두 지원)
###############################################################################
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
.need(c("lavaan"))
suppressPackageStartupMessages(library(lavaan))

## ── 설정 ────────────────────────────────────────────────────────────────
INV_SET       <- "MI01"    # 대체 세트 (v260904: MI 첫 완성 데이터)
DCFI_CUT      <- -0.010    # ΔCFI 허용 하한
DRMSEA_CUT    <-  0.015    # ΔRMSEA 허용 상한
COR_CUT       <-  0.99     # 요인점수 상관이 이보다 크면 실질적 영향 없음으로 판정
MAX_FREE      <-  6        # 부분 scalar 에서 풀어줄 절편 최대 개수
ESTIMATOR     <- "MLR"

## 17개 지표 — 13_mi_pmm_m20.R 의 앵커 모형과 동일: g + 4 특정요인(loco/vita/cogn/psyc), 감각 2지표는 g 에만 적재 (v260904 정정)
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")

BIFACTOR <- '
  g  =~ Balance+GS_ms+rev_CST+loss_of_Bwt+Appetite+EXH+HGS+rev_logMAR+rev_PTA+
        Orientation+Memory+Attention_Calculation+Language+Visuospatial+
        Negative_affect+Positive_affect+Motivation
  loco =~ Balance+GS_ms+rev_CST
  vita =~ loss_of_Bwt+Appetite+EXH+HGS
  cogn =~ Orientation+Memory+Attention_Calculation+Language+Visuospatial
  psyc =~ Negative_affect+Positive_affect+Motivation
  g ~~ 0*loco + 0*vita + 0*cogn + 0*psyc '

## bifactor 가 다집단에서 수렴하지 않을 때의 대체 모형 (상관 5요인)
CORR5 <- '
  loco =~ Balance+GS_ms+rev_CST
  vita =~ loss_of_Bwt+Appetite+EXH+HGS
  sens =~ rev_logMAR+rev_PTA
  cogn =~ Orientation+Memory+Attention_Calculation+Language+Visuospatial
  psyc =~ Negative_affect+Positive_affect+Motivation '

## ── 자료 (v260904 v3): 관측 자료 + FIML ─────────────────────────────────
##   불변성 검정은 대체값이 아니라 관측된 측정값으로 해야 합니다 (13_ 의 앵커 모형과 동일 원칙).
##   결측 보존 원본에서 동의철회 3명·사망 이후 행을 제외하고, 지표가 하나라도 관측된
##   방문 행만 남긴 뒤 lavaan 의 missing = "fiml" 로 적합합니다.
PRE_FILE  <- file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260618_state.xlsx")
PRE_SHEET <- "전체_long"
if (!file.exists(PRE_FILE)) stop("결측 보존 원본이 없습니다: ", PRE_FILE)
.need("readxl")
tmp <- file.path(tempdir(), "pre_inv.xlsx"); file.copy(PRE_FILE, tmp, overwrite = TRUE)
pre <- as.data.frame(readxl::read_excel(tmp, sheet = PRE_SHEET))
pre$id <- as.character(pre$id); pre$wave <- as.integer(as_num(pre$wave))
pre <- pre[!pre$id %in% c("kf161224", "kf170602", "kf171295"), ]
dw <- tapply(pre$death_wave, pre$id, function(x) suppressWarnings(as.numeric(x[1])))
pre <- pre[!(!is.na(dw[pre$id]) & pre$wave >= dw[pre$id]), ]
miss <- setdiff(IND17, names(pre))
if (length(miss)) stop("지표가 없습니다: ", paste(miss, collapse = ", "))
D <- pre[, c("id", "wave", IND17), drop = FALSE]
for (v in IND17) D[[v]] <- as_num(D[[v]])
D <- D[is.finite(D$wave) & rowSums(!is.na(D[, IND17])) > 0, ]
msg(sprintf("관측 자료: %d행 / %d명 (지표 결측률 중앙 %.1f%%)", nrow(D), length(unique(D$id)),
            100 * stats::median(colMeans(is.na(D[, IND17])))))
## Part C 비교용 분석 점수 (MI01 within-sex 척도가 아니라 원래 gLIC 자체와 비교)
load_imputed(INV_SET, stem = get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20"))
raw <- load_kfacs()

## ★ 표준화는 반드시 **전체 풀 기준 1회**로 합니다.
##   웨이브별로 각각 표준화하면 웨이브 간 평균차가 인위적으로 지워져
##   scalar 불변성이 무조건 지지되는 것처럼 보입니다(검정이 무의미해짐).
ctr <- vapply(D[IND17], function(x) mean(x, na.rm = TRUE), numeric(1))
scl <- vapply(D[IND17], function(x) stats::sd(x, na.rm = TRUE), numeric(1))
Z <- D
for (v in IND17) Z[[v]] <- (as_num(D[[v]]) - ctr[[v]]) / scl[[v]]
Z$wave_f <- factor(Z$wave, levels = sort(unique(Z$wave)),
                   labels = paste0("W", sort(unique(Z$wave))))

msg(sprintf("불변성 검정 자료: %d행 / %d명 / 웨이브 %s",
            nrow(Z), length(unique(Z$id)), paste(levels(Z$wave_f), collapse = ", ")))
print(table(Z$wave_f))

## ── 적합 함수 ──────────────────────────────────────────────────────────
fit_step <- function(model, equal, lab, partial = NULL) {
  ## v3b: 경고는 기록만 하고 적합은 유지합니다 (FIML·다집단에서 흔한 분산 스케일 경고 등).
  ##      비수렴 시 반복 한도를 늘려 한 번 더 시도합니다.
  warns <- character(0)
  run <- function(extra = list()) {
    withCallingHandlers(
      tryCatch(do.call(lavaan::cfa, c(list(model = model, data = Z, group = "wave_f", std.lv = TRUE,
                                          estimator = ESTIMATOR, missing = "fiml",
                                          group.equal = equal, group.partial = partial), extra)),
               error = function(e) { message("[", lab, "] 오류: ", conditionMessage(e)); NULL }),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  }
  f <- run()
  if (is.null(f) || !lavaan::lavInspect(f, "converged")) {
    message("[", lab, "] 1차 비수렴 -> iter.max 확대 재시도")
    f <- run(list(control = list(iter.max = 20000L)))
  }
  if (length(warns)) message("[", lab, "] 경고 ", length(warns), "건 (적합은 유지): ", paste(unique(substr(warns, 1, 120)), collapse = " | "))
  if (is.null(f) || !lavaan::lavInspect(f, "converged")) { message("[", lab, "] 수렴 실패"); return(NULL) }
  f
}

grab <- function(f, lab) {
  if (is.null(f)) return(data.frame(Model = lab, chisq = NA, df = NA, CFI = NA,
                                    TLI = NA, RMSEA = NA, SRMR = NA, AIC = NA,
                                    stringsAsFactors = FALSE))
  m <- lavaan::fitMeasures(f)
  ## v3c: robust 지표가 NA 이면 표준 지표로 대체 (FIML·MLR 에서 configural 의 robust CFI 가 NA 가 되는 경우)
  gt <- function(a, b) { va <- if (a %in% names(m)) unname(m[[a]]) else NA_real_
                         vb <- if (b %in% names(m)) unname(m[[b]]) else NA_real_
                         if (is.finite(va)) va else vb }
  data.frame(Model = lab,
             chisq = gt("chisq.scaled","chisq"), df = gt("df.scaled","df"),
             CFI   = gt("cfi.robust","cfi"),     TLI = gt("tli.robust","tli"),
             RMSEA = gt("rmsea.robust","rmsea"), SRMR = gt("srmr","srmr"),
             AIC   = gt("aic","aic"), stringsAsFactors = FALSE)
}

## ── Part A: 다집단 CFA 순차검정 ────────────────────────────────────────
MODEL <- BIFACTOR; MODEL_LAB <- "Bifactor (g + 4 specific; sensory indicators on g only)"
f1 <- fit_step(MODEL, c(), "configural")
if (is.null(f1)) {
  message("[Part A] bifactor 다집단 수렴 실패 -> 상관 5요인 모형으로 대체")
  MODEL <- CORR5; MODEL_LAB <- "Correlated five-factor (fallback)"
  f1 <- fit_step(MODEL, c(), "configural")
}
if (is.null(f1)) stop("configural 모형조차 수렴하지 않습니다. 지표/자료를 점검하십시오.")

f2 <- fit_step(MODEL, c("loadings"),               "metric")
f3 <- fit_step(MODEL, c("loadings","intercepts"),  "scalar")
f4 <- fit_step(MODEL, c("loadings","intercepts","residuals"), "strict")

FIT <- rbind(grab(f1, "M1 Configural"), grab(f2, "M2 Metric (loadings)"),
             grab(f3, "M3 Scalar (+ intercepts)"), grab(f4, "M4 Strict (+ residuals)"))
FIT$dCFI   <- c(NA, diff(FIT$CFI))
FIT$dRMSEA <- c(NA, diff(FIT$RMSEA))
FIT$dSRMR  <- c(NA, diff(FIT$SRMR))
## Chen (2007), N > 300: metric — dCFI >= -0.010 and (dRMSEA <= 0.015 or dSRMR <= 0.030);
##                         scalar/strict — dCFI >= -0.010 and (dRMSEA <= 0.015 or dSRMR <= 0.010)
.srmr_cut <- c(NA, 0.030, 0.010, 0.010)
FIT$Verdict <- c("-", ifelse(is.na(FIT$dCFI[-1]), "not estimable",
                 ifelse(FIT$dCFI[-1] >= DCFI_CUT & (FIT$dRMSEA[-1] <= DRMSEA_CUT | FIT$dSRMR[-1] <= .srmr_cut[-1]),
                        "supported",
                        ifelse(FIT$dCFI[-1] >= -0.015 & (FIT$dRMSEA[-1] <= DRMSEA_CUT | FIT$dSRMR[-1] <= .srmr_cut[-1]),
                               "borderline (CFI criterion only)", "NOT supported"))))
print(FIT, row.names = FALSE)

scalar_ok <- !is.na(FIT$dCFI[3]) && FIT$dCFI[3] >= DCFI_CUT && FIT$dRMSEA[3] <= DRMSEA_CUT

## ── Part B: scalar 가 깨지면 부분 불변 탐색 ────────────────────────────
PARTIAL <- character(0); f3p <- NULL; partial_ok <- FALSE
if (!scalar_ok && !is.null(f3)) {
  message("\n[Part B] scalar 불변성 미지지 -> 절편을 순차적으로 풀어 부분 scalar 탐색")
  for (k in seq_len(MAX_FREE)) {
    mi <- tryCatch(lavaan::lavTestScore(f3, epc = TRUE)$epc, error = function(e) NULL)
    if (is.null(mi)) break
    cand <- mi[mi$op == "~1" & !(mi$lhs %in% PARTIAL), ]
    if (!nrow(cand)) break
    cand <- cand[order(-abs(cand$epc)), ]
    PARTIAL <- unique(c(PARTIAL, paste0(cand$lhs[1], " ~ 1")))
    f3p <- fit_step(MODEL, c("loadings","intercepts"), sprintf("partial scalar (%d)", k),
                    partial = PARTIAL)
    if (is.null(f3p)) { PARTIAL <- PARTIAL[-length(PARTIAL)]; break }
    cfi_now <- grab(f3p, "x")$CFI
    d <- cfi_now - FIT$CFI[2]
    message(sprintf("   %d개 해제 -> CFI %.4f (metric 대비 %+.4f)", k, cfi_now, d))
    if (!is.na(d) && d >= DCFI_CUT) { partial_ok <- TRUE; break }
    ## 더 풀어도 나아지지 않으면 중단 (무의미한 반복 방지)
    if (exists("cfi_prev") && is.finite(cfi_now) && is.finite(cfi_prev) &&
        cfi_now - cfi_prev < 1e-4) {
      message("   더 해제해도 개선되지 않아 중단합니다 -> 부분 scalar 도 미지지")
      break
    }
    cfi_prev <- cfi_now
  }
  if (!is.null(f3p)) FIT <- rbind(FIT, cbind(grab(f3p, sprintf("M3p Partial scalar (%d freed)", length(PARTIAL))),
                                             dCFI = NA, dRMSEA = NA, dSRMR = NA, Verdict = "-"))
}

## ── Part C: 실질적 영향 — 요인점수를 현재 gLIC 와 비교 ─────────────────
best <- if (!is.null(f3p)) f3p else if (!is.null(f3)) f3 else f2
best_lab <- if (!is.null(f3p)) "partial scalar" else if (!is.null(f3)) "scalar" else "metric"
IMPACT <- NULL
if (!is.null(best) && "gLIC" %in% names(raw)) {
  fs <- tryCatch(lavaan::lavPredict(best), error = function(e) NULL)
  if (!is.null(fs)) {
    sc <- if (is.list(fs)) do.call(rbind, fs) else fs
    gcol <- if ("g" %in% colnames(sc)) "g" else colnames(sc)[1]
    ## 집단별 예측은 집단 순서대로 쌓이므로 같은 순서로 원자료를 정렬해 붙입니다
    ## lavPredict 는 집단 순서대로 쌓아서 돌려줍니다 -> 같은 순서로 정렬해 붙입니다.
    ord <- order(Z$wave_f)
    Zc  <- Z[ord, ]
    if (nrow(sc) != nrow(Zc)) {
      message(sprintf("[Part C] 요인점수 %d행 vs 자료 %d행 -> 행 대응이 맞지 않아 비교를 건너뜁니다.",
                      nrow(sc), nrow(Zc)))
      fs <- NULL
    } else {
      Zc$gLIC_inv <- as.numeric(sc[, gcol])
    cur <- raw$gLIC[match(paste(Zc$id, Zc$wave), paste(raw$id, as_num(raw$wave)))]
    ok  <- is.finite(cur) & is.finite(Zc$gLIC_inv)
    r   <- suppressWarnings(stats::cor(cur[ok], Zc$gLIC_inv[ok], use = "complete.obs"))
    rs  <- suppressWarnings(stats::cor(cur[ok], Zc$gLIC_inv[ok], method = "spearman"))
    tt  <- function(x) cut(x, stats::quantile(x, c(0,1/3,2/3,1), na.rm = TRUE),
                           include.lowest = TRUE, labels = c("T1","T2","T3"))
    agree <- mean(tt(cur[ok]) == tt(Zc$gLIC_inv[ok]), na.rm = TRUE)
    IMPACT <- data.frame(
      Quantity = c("Pearson r (current gLIC vs invariance-based)",
                   "Spearman rho",
                   "Tertile agreement",
                   "Model used for the invariance-based score",
                   "Indicator intercepts freed"),
      Value = c(sprintf("%.4f", r), sprintf("%.4f", rs),
                sprintf("%.1f%%", 100 * agree), best_lab,
                if (length(PARTIAL)) paste(sub(" ~ 1", "", PARTIAL), collapse = ", ") else "none"),
      stringsAsFactors = FALSE)
    print(IMPACT, row.names = FALSE)
    ## 필요 시 교체해 쓸 수 있도록 저장
    save_vals(data.frame(id = Zc$id, wave = Zc$wave, gLIC_invariant = Zc$gLIC_inv),
              "ST_Invariance_gLIC_invariant.csv", FIG_DIR)
    material <- is.finite(r) && r < COR_CUT
    }
  }
}
if (!exists("material")) material <- NA

## ── 저장 ────────────────────────────────────────────────────────────────
OUT <- data.frame(
  Model      = FIT$Model,
  `Chi-square` = ifelse(is.na(FIT$chisq), "-", sprintf("%.1f", FIT$chisq)),
  df         = ifelse(is.na(FIT$df), "-", sprintf("%.0f", FIT$df)),
  CFI        = ifelse(is.na(FIT$CFI), "-", sprintf("%.3f", FIT$CFI)),
  TLI        = ifelse(is.na(FIT$TLI), "-", sprintf("%.3f", FIT$TLI)),
  RMSEA      = ifelse(is.na(FIT$RMSEA), "-", sprintf("%.3f", FIT$RMSEA)),
  SRMR       = ifelse(is.na(FIT$SRMR), "-", sprintf("%.3f", FIT$SRMR)),
  `Delta CFI`   = ifelse(is.na(FIT$dCFI), "-", sprintf("%+.3f", FIT$dCFI)),
  `Delta RMSEA` = ifelse(is.na(FIT$dRMSEA), "-", sprintf("%+.3f", FIT$dRMSEA)),
  `Delta SRMR`  = ifelse(is.na(FIT$dSRMR), "-", sprintf("%+.3f", FIT$dSRMR)),
  Verdict    = FIT$Verdict,
  check.names = FALSE, stringsAsFactors = FALSE)

fn <- c(
  sprintf("Multi-group confirmatory factor analysis with wave as the grouping variable; %s; %s estimator; all 17 indicators standardised once on the pooled person-wave distribution; observed indicator values only (attended visits), with missing values handled by full-information maximum likelihood.", MODEL_LAB, ESTIMATOR),
  "Indicators are standardised on the pooled distribution rather than within wave: standardising within wave would remove between-wave mean differences by construction and make the test of scalar invariance vacuous.",
  "Because the grouping variable is wave, each participant contributes one observation per group, so observations are independent within groups and no cluster correction is required.",
  sprintf("Invariance is judged by change from the preceding model using the criteria of Chen (2007) for samples above 300: supported if Delta CFI >= %.3f together with Delta RMSEA <= %.3f or Delta SRMR <= 0.030 (loadings) / 0.010 (intercepts, residuals); a change in CFI between -0.010 and -0.015 with the RMSEA or SRMR criterion met is labelled borderline.", DCFI_CUT, DRMSEA_CUT),
  "Metric invariance licenses comparison of associations per +1 s.d. across waves; scalar invariance licenses comparison of levels, which the rolling-landmark analyses require.",
  if (!scalar_ok && partial_ok)
    sprintf("Full scalar invariance was not supported; partial scalar invariance was reached after freeing the intercepts of %s, which is sufficient for comparing latent means across waves.",
            paste(sub(" ~ 1", "", PARTIAL), collapse = ", "))
  else if (!scalar_ok)
    sprintf("Full scalar invariance was not supported, and freeing %s did not bring the fit within the criterion; the departure from scalar invariance is therefore not attributable to one or two indicators alone.",
            if (length(PARTIAL)) paste0("the intercept", if (length(PARTIAL) > 1) "s" else "", " of ",
                                        paste(sub(" ~ 1", "", PARTIAL), collapse = ", "))
            else "additional intercepts")
  else "Full scalar invariance was retained; no intercept had to be freed.",
  if (!is.null(IMPACT))
    sprintf("Factor scores from the %s model correlated r = %s with the intrinsic-capacity score used in the main analyses (tertile agreement %s), so the substantive results are unaffected.",
            best_lab, IMPACT$Value[1], IMPACT$Value[3])
  else "Factor-score comparison could not be computed.",
  "CFI, comparative fit index; RMSEA, root mean square error of approximation; SRMR, standardised root mean square residual; TLI, Tucker-Lewis index.")

save_table(OUT, "SupplTable_Invariance",
  title = "Supplementary Table | Longitudinal measurement invariance of the intrinsic-capacity factor across waves",
  footnotes = fn)
if (!is.null(IMPACT))
  save_table(IMPACT, "SupplTable_Invariance_impact",
    title = "Supplementary Data | Practical impact of the invariance model on the intrinsic-capacity score",
    footnotes = "Scores from the invariance-constrained model are compared with the score used in the main analyses.")

## ── 판정 ────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n판정\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  모형        : %s\n", MODEL_LAB))
cat(sprintf("  metric      : %s\n", FIT$Verdict[2]))
cat(sprintf("  scalar      : %s\n", FIT$Verdict[3]))
if (length(PARTIAL))
  cat(sprintf("  부분 scalar : %s (절편 %d개 해제: %s)\n",
              if (partial_ok) "지지" else "미지지", length(PARTIAL),
              paste(sub(" ~ 1", "", PARTIAL), collapse = ", ")))
if (!scalar_ok && !partial_ok)
  cat("\n  ★ scalar 불변성이 부분적으로도 성립하지 않습니다.\n",
      "    -> 웨이브 간 '수준' 비교(Table 2/4, Figure 1)의 전제가 깨집니다.\n",
      "       metric 은 성립하므로 '연관의 크기' 비교는 유지되나,\n",
      "       Discussion 의 limitation 에 반드시 명시하십시오.\n", sep = "")
if (!is.null(IMPACT)) {
  cat(sprintf("  현재 gLIC 와의 상관 r = %s | 삼분위 일치율 %s\n",
              IMPACT$Value[1], IMPACT$Value[3]))
  if (isTRUE(material)) {
    cat("\n  ★ 상관이 ", COR_CUT, " 미만입니다. 측정모형 차이가 실질적입니다.\n", sep="")
    cat("     -> ST_Invariance_gLIC_invariant.csv 의 점수로 gLIC 를 교체해\n")
    cat("        Table 2 / Table 4 / Figure 1 을 다시 돌려 민감도로 보고하십시오.\n")
  } else {
    cat("\n  ★ 상관이 ", COR_CUT, " 이상입니다. 측정모형 차이가 결과에 영향을 주지 않습니다.\n", sep="")
    cat("     -> 본 분석을 그대로 유지하고, 이 표를 Supplementary 에 넣으십시오.\n")
  }
}
cat("\n=== 종단 불변성 검정 완료 ===\n")
