###############################################################################
## gIC TRAJECTORY -- STEP 0.5: STATE 변수 정비 (pre-imp에 추가)
##
## 목적: pre-imp(관측치) 파일에 state 변수 추가. imputed는 이미 완비 -> 손대지 않음.
##
## ★ FUP는 원본 1차출처 KFACS_fup_w1w5.xlsx 사용 (imputed FUP와 100% 일치 검증).
##   코드정의(추적조사 상세 시트): 1센터본인 2가정방문 3전화응답 4사망
##   5장기요양시설입소 6병원입소 7연락안됨 8단순거절 9이사 10미참여
##   11센터대리 12전화대리 13동의철회 14가구대리. ★입소=5,6.
##
## 추가 컬럼 (5개):
##   FUP, state, state_label    <- 원본FUP + adl_score로 산출 (imputed와 정의 일치)
##   state_obs, state_obs_label <- 관측 기반 (사망후/미관측 NA, 민감도용)
##
## 검증(260618, Python 사전검증):
##   - 원본 FUP vs imputed FUP 일치율 100.00% (n=14,607, 불일치 0)
##   - id+wave 키 15,070 완전일치, 관측 adl_score 100% 동일
##   - state_obs vs imputed state 일치율 100.00% (관측구간, 불일치 0)
##
## ★ state 정의:
##   Death(4)  = w == death_wave (그 wave에만; 사망 후 행은 obs에서 NA)
##   Severe(3) = adl_score >= 3  OR  FUP in {5,6}
##   Mild(2)   = adl_score 1-2 ;  Normal(1) = adl_score 0
###############################################################################

## ---- PACKAGES -------------------------------------------------------------
if (!requireNamespace("openxlsx",quietly=TRUE)) install.packages("openxlsx")
if (!requireNamespace("dplyr",quietly=TRUE))    install.packages("dplyr")
suppressPackageStartupMessages({ library(openxlsx); library(dplyr) })

## ---- PATHS (Nat_Aging) ----------------------------------------------------
## 경로는 환경변수로 받습니다. 저장소에 개인 폴더 경로를 두지 않습니다.
##   Sys.setenv(KFACS_ROOT = "D:/path/to/Nat_Aging")   # data/ 가 있는 프로젝트 루트
DATA_DIR <- file.path(Sys.getenv("KFACS_ROOT", unset = "."), "data")
if (!dir.exists(DATA_DIR)) stop("data/ 폴더가 없습니다: ", DATA_DIR, "\n  -> Sys.setenv(KFACS_ROOT = \"...\") 로 지정하십시오.")
PATH_PRE <- file.path(DATA_DIR, "KFACS_master_FINAL_preimput_260616.xlsx")
PATH_FUP <- file.path(DATA_DIR, "KFACS_fup_w1w5.xlsx")   # ★ 원본 FUP
OUT_PRE  <- file.path(DATA_DIR, "KFACS_master_FINAL_preimput_260618_state.xlsx")

to_tmp <- function(src, tag){
  stopifnot(file.exists(src)); tmp <- file.path(tempdir(), paste0(tag,".xlsx"))
  file.copy(src, tmp, overwrite=TRUE); tmp
}

## ---- 1. 원본 FUP -> long 변환 ----------------------------------------------
fp <- to_tmp(PATH_FUP, "fup")
## BL(wave1): FU_reason 없음 -> 전원 1(센터본인)
bl <- read.xlsx(fp, sheet="BL")[, "No.", drop=FALSE]
names(bl) <- "id"; bl$wave <- 1; bl$FUP <- 1
mk_w <- function(sheet, col, wv){
  d <- read.xlsx(fp, sheet=sheet)[, c("No.", col)]
  names(d) <- c("id","FUP"); d$wave <- wv; d[, c("id","wave","FUP")]
}
fup_long <- rbind(
  bl[, c("id","wave","FUP")],
  mk_w("W2(2018-2019)","W2_sample_FU_reason",2),
  mk_w("W3(2020-2021)","W3_sample_FU_reason",3),
  mk_w("W4(2022-2023)","W4_sample_FU_reason",4),
  mk_w("W5(2024-2025)","W5_sample_FU_reason",5))
fup_long$wave <- as.numeric(fup_long$wave)
cat("FUP long:", nrow(fup_long), "rows |", dplyr::n_distinct(fup_long$id), "ids\n")

## ---- 2. pre-imp 로드 + FUP merge ------------------------------------------
pre <- read.xlsx(to_tmp(PATH_PRE,"pre"), sheet="전체_long")
pre$wave <- as.numeric(pre$wave)
pre2 <- pre %>% left_join(fup_long, by=c("id","wave"))
stopifnot(nrow(pre2)==nrow(pre))

## ---- 3. state (main) + state_obs (관측) -----------------------------------
mk_state <- function(adl, fup, dw, w, observed_only){
  out <- rep(NA_real_, length(adl))
  out[!is.na(adl) & adl==0] <- 1
  out[!is.na(adl) & adl>=1 & adl<=2] <- 2
  out[!is.na(adl) & adl>=3] <- 3
  out[!is.na(fup) & fup %in% c(5,6)] <- 3          # 입소
  out[!is.na(dw) & w==dw] <- 4                     # 사망 wave
  if (observed_only) out[!is.na(dw) & w>dw] <- NA_real_  # 사망후 미관측
  out
}
## state: imputed 정의와 동일하게 만들되, pre-imp 관측 adl 기반.
## (사망후 행은 pre-imp adl이 NA라 자동 NA; imputed처럼 채우지 않음)
pre2$state <- mk_state(pre2$adl_score, pre2$FUP, pre2$death_wave, pre2$wave, observed_only=FALSE)
pre2$state_obs <- mk_state(pre2$adl_score, pre2$FUP, pre2$death_wave, pre2$wave, observed_only=TRUE)
lab <- c("Normal","Mild","Severe","Death")
pre2$state_label     <- factor(pre2$state,     levels=1:4, labels=lab)
pre2$state_obs_label <- factor(pre2$state_obs, levels=1:4, labels=lab)

## ---- 4. VERIFY ------------------------------------------------------------
cat("\nstate_label 분포:\n");     print(table(pre2$state_label,     useNA="ifany"))
cat("state_obs_label 분포:\n"); print(table(pre2$state_obs_label, useNA="ifany"))
cat("입소(FUP 5,6) by wave:\n")
for (w in 1:5){
  s <- pre2[pre2$wave==w,]
  cat(sprintf("  wave %d: 입소=%d\n", w, sum(s$FUP %in% c(5,6), na.rm=TRUE)))
}

## ---- 5. SAVE --------------------------------------------------------------
write.xlsx(list("전체_long"=pre2), OUT_PRE, overwrite=TRUE)
cat("\n저장 완료:\n  ", OUT_PRE, "\n")
cat("추가 컬럼: FUP, state, state_label, state_obs, state_obs_label\n")
cat("※ imputed 파일은 state 완비 -> 수정 불필요.\n")
cat("※ 원본 FUP(KFACS_fup_w1w5.xlsx) 사용; imputed FUP와 100% 일치 검증됨.\n")
