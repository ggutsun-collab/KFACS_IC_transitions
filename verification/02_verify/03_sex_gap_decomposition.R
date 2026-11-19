# =============================================================================
# 11 — IC 점수의 성별 격차가 어디서 오는가
#
# 선생님 지적: 악력·보행속도·하지근력은 성별 절단점을 다르게 넣었다.
# 그렇다면 gLIC_z 가 왜 이렇게 성별로 갈리는지 다른 곳에서 설명되어야 합니다.
#
# 관찰된 사실 (이건 해석이 아니라 데이터입니다)
#   robust 층 성별 구성 : 남 715 / 여 634   <- 거의 균형. frailty 기준은 성별 보정됨.
#   그 안의 IC tertile  : T1 남10/여96,  T3 남586/여187   <- 극단적 분리
#   즉 frailty 는 성별 중립인데 IC 점수만 성별로 갈립니다.
#
# 이 스크립트는 그 격차가 어느 지표에서 오는지 지표 단위로 찾습니다.
# 제 이전 설명(수행기반 지표 미표준화)이 맞는지 틀리는지 여기서 결정됩니다.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(dplyr) })
OUT <- "sexgap_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "sexgap_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

D  <- d
W1 <- D %>% filter(!is.na(wave)) %>% group_by(id) %>% slice_min(wave, n=1, with_ties=FALSE) %>% ungroup()

## 성별 벡터 확보
sexv <- NULL
for (v in c("sex_unified","sex_f",".female")) if (v %in% names(W1)) { sexv <- W1[[v]]; sexnm <- v; break }
say("성별 변수: %s", sexnm); print(table(sexv, useNA="ifany"))
W1$.SEX <- as.character(sexv)

## ===========================================================================
rule("0. d 의 전체 열 이름 — 17개 IC 지표를 특정하기 위해")
print(matrix(c(names(D), rep("", (5 - length(names(D)) %% 5) %% 5)), ncol=5, byrow=TRUE), quote=FALSE)

## ===========================================================================
rule("1. gLIC 계열 변수의 성별 격차")
gl <- names(W1)[grepl("glic", tolower(names(W1)))]
say("gLIC 계열 열: %s", paste(gl, collapse=", "))
smd <- function(x, g) {
  x <- as.numeric(x); a <- x[g==unique(g)[1]]; b <- x[g==unique(g)[2]]
  (mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE)) /
    sqrt((var(a,na.rm=TRUE)+var(b,na.rm=TRUE))/2)
}
for (v in gl) {
  x <- W1[[v]]
  if (!is.numeric(x)) { say("\n%s (범주형):", v); print(table(x, W1$.SEX, useNA="ifany")); next }
  m <- tapply(x, W1$.SEX, function(z) c(n=sum(!is.na(z)), mean=mean(z,na.rm=TRUE), sd=sd(z,na.rm=TRUE)))
  say("\n%s :", v); print(round(do.call(rbind, m), 3))
  say("  표준화 평균차(SMD) = %.3f   |0.2| 이상이면 실질적 격차", smd(x, W1$.SEX))
}
say("\n>>> gLIC_z 의 SMD 가 크면, 성별 절단점을 쓴 지표가 있더라도")
say(">>> 최종 인자점수는 성별 중립이 아닙니다. 이것이 결정적 수치입니다.")

## ===========================================================================
rule("2. 어느 지표가 격차를 만드는가 — 모든 수치형 열의 성별 SMD")
num <- names(W1)[sapply(W1, is.numeric)]
num <- setdiff(num, c("wave","time","death_event","death_wave","death_wave_id","FUP"))
tab <- data.frame(var=num, stringsAsFactors=FALSE)
tab$n_ok <- sapply(num, function(v) sum(!is.na(W1[[v]])))
tab$SMD  <- sapply(num, function(v) suppressWarnings(smd(W1[[v]], W1$.SEX)))
tab <- tab %>% filter(n_ok > 1000, is.finite(SMD)) %>% arrange(desc(abs(SMD)))
say("성별 격차가 큰 순서 (상위 30개):")
print(head(as.data.frame(tab), 30), digits=3, row.names=FALSE)
write.csv(tab, file.path(OUT,"sex_SMD_all_vars.csv"), row.names=FALSE)
say("\n>>> 상위에 악력/보행속도/근력이 있으면 -> 수행지표가 여전히 성별 의존적입니다.")
say(">>> 상위에 인지/교육/우울이 있으면  -> 격차의 출처는 다른 영역입니다.")
say(">>> 어느 쪽인지에 따라 원고 Methods 와 Limitations 의 문장이 완전히 달라집니다.")

## ===========================================================================
rule("3. 교육으로 설명되는가 (한국 1930-40년대생 여성의 교육연수 격차)")
edu <- intersect(c("edu_bl","educ_bl","edu_high_bl","edu_bl_f"), names(W1))
if (length(edu)) {
  say("교육 변수: %s", edu[1])
  print(table(W1[[edu[1]]], W1$.SEX, useNA="ifany"))
  if ("gLIC_z" %in% names(W1)) {
    f0 <- lm(gLIC_z ~ .SEX, data=W1)
    f1 <- lm(as.formula(paste("gLIC_z ~ .SEX +", edu[1])), data=W1)
    say("\ngLIC_z 의 성별 계수:")
    say("  교육 미보정 : %.3f", coef(f0)[2])
    say("  교육 보정   : %.3f", coef(f1)[2])
    say("  감쇠 = %.1f%%", 100*(1-abs(coef(f1)[2])/abs(coef(f0)[2])))
    say(">>> 감쇠가 크면 성별 격차의 상당 부분이 교육을 통한 인지영역 차이입니다.")
  }
}

## ===========================================================================
rule("4. robust 층 안에서의 확인")
if (all(c("frailty_3cat","gLIC_z") %in% names(W1))) {
  RB <- W1 %>% filter(as.character(frailty_3cat)=="0")
  say("robust n = %d | 남 %d / 여 %d", nrow(RB),
      sum(RB$.SEX==unique(RB$.SEX)[1]), sum(RB$.SEX==unique(RB$.SEX)[2]))
  m <- tapply(RB$gLIC_z, RB$.SEX, function(z) c(n=length(z), mean=mean(z,na.rm=TRUE), sd=sd(z,na.rm=TRUE)))
  print(round(do.call(rbind,m),3))
  say("robust 층 내 gLIC_z 성별 SMD = %.3f", smd(RB$gLIC_z, RB$.SEX))
  cut_coh <- quantile(W1$gLIC_z, c(1/3,2/3), na.rm=TRUE)
  RB$tert <- cut(RB$gLIC_z, c(-Inf,cut_coh,Inf), labels=c("T1","T2","T3"))
  say("\ntertile x 성별:"); print(table(RB$tert, RB$.SEX))
  say("\n>>> frailty 기준은 성별 절단점을 써서 robust 층이 남634/여715 로 균형인데,")
  say(">>> 그 안에서 IC tertile 만 극단적으로 갈린다면 격차의 출처는 IC 점수 구성입니다.")
}
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
