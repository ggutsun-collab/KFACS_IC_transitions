###############################################################################
## KFACS v6 — 통합 imputation 파이프라인 (처음부터 끝까지 1회 실행)
##   입력 : KFACS_master_FINAL_preimput_260616.xlsx (전체_long, 결측 보존 원본)
##          KFACS_fup_w1w5.xlsx (FUP wave별 시트)
##   핵심 : ADL을 IADL 포함 거리로 대체 (전화조사 ADL 결측 → 같은시점 IADL 활용)
##   출력 : 4 imputed 세트 (MAIN/Seq/NoFrailty/MNAR) + FUP + 4-state + 복구변수
##   주의 : 한글경로 → setwd 후 상대 파일명 (인코딩 회피)
###############################################################################
suppressMessages({
  library(dplyr); library(tidyr); library(readxl); library(writexl)
  library(VIM); library(lavaan)
})

## ── 0. 설정 ──────────────────────────────────────────────────────────────
## 경로는 환경변수로 받습니다. 저장소에 개인 폴더 경로를 두지 않습니다.
##   Sys.setenv(KFACS_RAW_DIR = "D:/KF_DATA/build")
DIR_DATA <- Sys.getenv("KFACS_RAW_DIR", unset = ".")
OUT_DIR  <- Sys.getenv("KFACS_IMP_OUT", unset = file.path(DIR_DATA, "imputation_out"))
if (!dir.exists(DIR_DATA))
  stop("원자료 폴더가 없습니다: ", DIR_DATA,
       "\n  -> Sys.setenv(KFACS_RAW_DIR = \"...\") 로 지정하십시오.")
KNN_K      <- 5
MNAR_DELTA <- -0.5            # worse-shift (지표는 높을수록 좋음 → 음수가 악화)
setwd(DIR_DATA)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ── 변수 집합 정의 (실제 환경 확정값) ────────────────────────────────────
ind17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
lic_domain <- c("LIC_Locomotor","LIC_Vitality","LIC_Sensory","LIC_Cognition",
                "LIC_Psychological","LIC_Vision","LIC_Hearing")
disab_vars  <- c("adl_score","iadl_score")
frail_vars  <- c("chs_total","frailty_3cat")
comp_kicope <- c("kicope_locomotor","kicope_vitality","kicope_sensory",
                 "kicope_cognition","kicope_psychological")
comp_cic    <- c("cic_sppb20","cic_mna10","cic_hgs10","cic_va10","cic_ha10",
                 "cic_mmse20","cic_gds20")
match_full  <- c("age_unified","sex_unified","educ_bl","area_3cat","comorbid_count",
                 "phys_bmi","income_3cat","living_alone","marri_status","polypharmacy",
                 "pa_met_min_wk","pa_freq_3cat","smoking_current","strength_2plus")
# ★ ADL 대체 시 IADL을 거리변수로 추가 (전화조사 ADL 결측 핵심 보조변수)
match_adl   <- c(match_full, "iadl_score")

## ── 1. 원본 master 로드 ─────────────────────────────────────────────────
cat("===== 1. master 로드 =====\n")
master <- read_excel("KFACS_master_FINAL_preimput_260616.xlsx", sheet = "전체_long")
cat("master:", nrow(master), "x", ncol(master), "\n")

## ── 2. FUP 머지 (wave별 시트 → long → id+wave 결합) ─────────────────────
cat("\n===== 2. FUP 머지 =====\n")
fup_sheets <- c(`2`="W2(2018-2019)", `3`="W3(2020-2021)",
                `4`="W4(2022-2023)", `5`="W5(2024-2025)")
fup_long <- bind_rows(lapply(names(fup_sheets), function(w){
  d <- read_excel("KFACS_fup_w1w5.xlsx", sheet = fup_sheets[[w]])
  col <- grep("FU_reason", names(d), value = TRUE)[1]
  tibble(id = d[["No."]], wave = as.integer(w), FUP = as.numeric(d[[col]]))
}))
bl <- read_excel("KFACS_fup_w1w5.xlsx", sheet = "BL")
fup_long <- bind_rows(fup_long, tibble(id = bl[["No."]], wave = 1L, FUP = 1))

master <- master %>% mutate(wave = as.integer(wave)) %>%
  left_join(fup_long, by = c("id","wave"))
cat("FUP=4(사망):", sum(master$FUP == 4, na.rm = TRUE),
    "| FUP 결측:", sum(is.na(master$FUP)), "\n")

## ── 3. base 구성 (imputation 입력: 필요한 변수만, .rowid 부여) ───────────
cat("\n===== 3. base 구성 =====\n")
base_vars <- unique(c("id","wave","FUP", match_full, disab_vars, frail_vars,
                      ind17, lic_domain, comp_kicope, comp_cic))
base_vars <- base_vars[base_vars %in% names(master)]
base <- master %>% select(all_of(base_vars)) %>%
  mutate(.rowid = row_number())
# 범주형으로 다뤄야 매칭이 맞는 변수 (factor화)
cat_keys <- intersect(c("sex_unified","area_3cat","income_3cat","living_alone",
                        "marri_status","pa_freq_3cat","smoking_current",
                        "strength_2plus","polypharmacy","frailty_3cat"),
                      names(base))
base[cat_keys] <- lapply(base[cat_keys], factor)
cat("base:", nrow(base), "x", ncol(base), "\n")

## ── 공통 함수: kNN / sequential hot-deck / postprocess ──────────────────
run_knn <- function(df_in, target_vars, match_keys) {
  tv <- target_vars[target_vars %in% names(df_in)]
  tv <- tv[colSums(is.na(df_in[tv, drop = FALSE])) > 0]
  if (length(tv) == 0) return(df_in)
  dist_v <- intersect(unique(c("wave", match_keys)), names(df_in))
  VIM::kNN(df_in, variable = tv, dist_var = dist_v, k = KNN_K, imp_var = FALSE)
}
hd_block <- function(chunk, target_vars, match_keys) {
  tv <- target_vars[target_vars %in% names(chunk)]
  tv <- tv[colSums(is.na(chunk[tv, drop = FALSE])) > 0]
  if (length(tv) == 0) return(chunk)
  mk  <- match_keys[match_keys %in% names(chunk)]
  ord <- mk[sapply(chunk[mk], is.numeric)]
  VIM::hotdeck(chunk, variable = tv, domain_var = "sex_unified",
               ord_var = if (length(ord)) ord else NULL, imp_var = FALSE)
}
run_hotdeck <- function(df_in, target_vars, match_keys) {
  parts <- split(df_in, df_in$wave)
  out   <- lapply(parts, hd_block, target_vars = target_vars, match_keys = match_keys)
  res   <- do.call(rbind, out)
  res[order(res$.rowid), ]
}
postprocess <- function(d, label) {
  d$KICOPE <- rowSums(d[, comp_kicope], na.rm = FALSE)
  d$CIC    <- rowSums(d[, comp_cic],    na.rm = FALSE)
  # ★ bifactor: Vision/Hearing은 g에만, 특수요인은 loco/vita/cogn/psyc 4개 (sens 없음)
  bifactor <- '
    g  =~ Balance+GS_ms+rev_CST+loss_of_Bwt+Appetite+EXH+HGS+rev_logMAR+rev_PTA+
          Orientation+Memory+Attention_Calculation+Language+Visuospatial+
          Negative_affect+Positive_affect+Motivation
    loco =~ Balance+GS_ms+rev_CST
    vita =~ loss_of_Bwt+Appetite+EXH+HGS
    cogn =~ Orientation+Memory+Attention_Calculation+Language+Visuospatial
    psyc =~ Negative_affect+Positive_affect+Motivation
    g ~~ 0*loco + 0*vita + 0*cogn + 0*psyc '
  fit <- tryCatch(cfa(bifactor, data = d, std.lv = TRUE,
                      estimator = "MLR", missing = "ml"),
                  error = function(e) NULL)
  if (!is.null(fit) && lavInspect(fit, "converged")) {
    d$gLIC <- as.numeric(lavPredict(fit, newdata = d)[, "g"])
    cat(sprintf("  [%s] bifactor CFI=%.3f RMSEA=%.3f\n", label,
        fitMeasures(fit, "cfi.robust"), fitMeasures(fit, "rmsea.robust")))
  } else { warning(sprintf("[%s] bifactor 재적합 실패", label)); d$gLIC <- NA_real_ }
  d$gLIC_z   <- as.numeric(scale(d$gLIC))
  d$KICOPE_z <- as.numeric(scale(d$KICOPE))
  tert <- function(x) cut(x, quantile(x, c(0,1/3,2/3,1), na.rm = TRUE),
                          include.lowest = TRUE, labels = c("T1","T2","T3"))
  d$gLIC_tertile   <- tert(d$gLIC)
  d$KICOPE_tertile <- tert(d$KICOPE)
  d$CIC_tertile    <- tert(d$CIC)
  d
}

## ── 4. SET 1 — MAIN : kNN Gower, IADL→ADL 순서, frailty↔LIC 분리 ────────
cat("\n===== 4. SET1 MAIN =====\n")
m <- base
m <- run_knn(m, target_vars = "iadl_score",            match_keys = match_full)  # ① IADL 먼저
m <- run_knn(m, target_vars = c("adl_score", frail_vars), match_keys = match_adl) # ② ADL+frailty (IADL 포함 거리)
m <- run_knn(m, target_vars = c(ind17, lic_domain),    match_keys = c(match_full, disab_vars)) # ③ LIC
m <- run_knn(m, target_vars = c(comp_kicope, comp_cic),match_keys = c(match_full, disab_vars)) # ④ component
main_df <- postprocess(m, "MAIN")

## ── 5. SET 2 — Seq : sequential hot-deck (도너방식 민감도) ──────────────
cat("\n===== 5. SET2 Seq =====\n")
s <- base
s <- run_hotdeck(s, target_vars = "iadl_score",             match_keys = match_full)
s <- run_hotdeck(s, target_vars = c("adl_score", frail_vars), match_keys = match_adl)
s <- run_hotdeck(s, target_vars = c(ind17, lic_domain),     match_keys = c(match_full, disab_vars))
s <- run_hotdeck(s, target_vars = c(comp_kicope, comp_cic), match_keys = c(match_full, disab_vars))
seq_df <- postprocess(s, "Seq")

## ── 6. SET 3 — NoFrailty : LIC 매칭에 frailty/disab 배제 (민감도) ───────
cat("\n===== 6. SET3 NoFrailty =====\n")
f <- base
f <- run_knn(f, target_vars = "iadl_score",             match_keys = match_full)
f <- run_knn(f, target_vars = c("adl_score", frail_vars), match_keys = match_adl)
f <- run_knn(f, target_vars = c(ind17, lic_domain),     match_keys = match_full)  # ★ frailty/disab 배제
f <- run_knn(f, target_vars = c(comp_kicope, comp_cic), match_keys = match_full)
nofrail_df <- postprocess(f, "NoFrailty")

## ── 7. SET 4 — MNAR : MAIN + δ worse-shift (원결측 셀만, 연속 LIC지표) ──
cat("\n===== 7. SET4 MNAR =====\n")
mask <- as.data.frame(is.na(base))          # ★ matrix→df 고정 (원결측 위치)
mn   <- main_df
cont_ind <- ind17[sapply(ind17, function(v)
              v %in% names(base) &&
              length(unique(na.omit(base[[v]]))) > 10)]
for (v in cont_ind) {
  if (!v %in% names(mask) || !v %in% names(mn)) next
  sdv <- sd(base[[v]], na.rm = TRUE)
  idx <- mask[[v]]; idx[is.na(idx)] <- FALSE
  mn[[v]][idx] <- mn[[v]][idx] + MNAR_DELTA * sdv
}
mnar_df <- postprocess(mn, "MNAR")

## ── 8. 4-state 생성 (FUP=4 Death / FUP 5·6 or ADL≥3 Severe / 1-2 Mild / 0 Normal)
cat("\n===== 8. 4-state 생성 =====\n")
make_state <- function(d){
  d$state <- dplyr::case_when(
    d$FUP == 4                                ~ 4L,
    d$FUP %in% c(5,6)                         ~ 3L,
    !is.na(d$adl_score) & d$adl_score >= 3    ~ 3L,
    !is.na(d$adl_score) & d$adl_score >= 1    ~ 2L,
    !is.na(d$adl_score) & d$adl_score == 0    ~ 1L,
    TRUE                                      ~ NA_integer_
  )
  d$state_label <- factor(d$state, levels = 1:4,
                          labels = c("Normal","Mild","Severe","Death"))
  d
}

## ── 9. 빠진 변수 복구 (death류=id고정, iadl범주=imputed에서 재생성) ─────
death_vars <- master %>%
  select(id, death_event, death_wave, followup_years,
         death_date_new, death_time_method) %>%
  distinct(id, .keep_all = TRUE)

restore <- function(d){
  d <- make_state(d)
  d <- d %>% left_join(death_vars, by = "id") %>%
    mutate(iadl_2cat = as.integer(iadl_score >= 1),
           iadl_3cat = cut(iadl_score, c(-Inf,0,2,Inf), labels = c(0,1,2)))
  d
}
main_df    <- restore(main_df)
seq_df     <- restore(seq_df)
nofrail_df <- restore(nofrail_df)
mnar_df    <- restore(mnar_df)

## ── 10. 검증 ────────────────────────────────────────────────────────────
cat("\n===== 10. 검증 =====\n")
for (nm in c("main_df","seq_df","nofrail_df","mnar_df")) {
  d <- get(nm)
  cat(sprintf("[%s] nrow=%d | adl결측=%d | gLIC결측=%d | state결측=%d | cor(gLIC,chs)=%.3f\n",
      nm, nrow(d), sum(is.na(d$adl_score)), sum(is.na(d$gLIC)),
      sum(is.na(d$state)),
      cor(d$gLIC, as.numeric(d$chs_total), use = "complete.obs")))
}
cat("\nMAIN wave×state:\n")
print(addmargins(table(wave = main_df$wave, state = addNA(main_df$state_label))))
cat("\nMAIN 사망(death_event=1 @wave1):",
    sum(main_df$death_event == 1 & main_df$wave == 1, na.rm = TRUE), "\n")

## ── 11. 저장 ────────────────────────────────────────────────────────────
setwd(OUT_DIR)
write_xlsx(list(MAIN = main_df, Seq = seq_df,
                NoFrailty = nofrail_df, MNAR = mnar_df),
           "KFACS_imputed_4sets_FINAL_260617.xlsx")
saveRDS(list(MAIN = main_df, Seq = seq_df, NoFrailty = nofrail_df, MNAR = mnar_df),
        "KFACS_imputed_4sets_FINAL_260617.rds")
cat("\n저장 완료:", file.path(OUT_DIR, "KFACS_imputed_4sets_FINAL_260617.xlsx"), "\n")
cat("===== 파이프라인 종료 =====\n")
