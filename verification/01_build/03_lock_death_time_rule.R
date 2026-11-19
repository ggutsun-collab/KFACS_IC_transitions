# =============================================================================
# 09 — 사망시점 규칙 확정 (마지막 조각)
#
# 지금까지 확정된 것
#   n 3,011 / robust 1,349 / tertile 106-470-773 / 사망 500  <- 전부 원고와 일치
#   d$time 은 년 단위이고 최대 8.172 로 원고의 '최대 8.2년' 과 일치
#   d$wave 가 정확히 500 행에서 결측  -> 사망기록 행이 별도로 붙어 있을 가능성
#   연령차 기준 person-years 22,786 vs 원고 22,785  <- 사실상 일치
#
# 아직 안 맞는 것
#   원고 Table 3 의 5년 사망 2 / 23 / 26 (합 51)
#   midpoint 규칙으로는     5 / 35 / 44 (합 84)
#   -> 사망시점 배정 규칙이 다릅니다. 이 스크립트가 세 규칙을 전부 시험합니다.
#
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT <- "rebuild_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "death_time_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

D <- d
B <- readRDS(file.path(OUT, "wave1_frame_SOURCE.rds"))

## ===========================================================================
rule("1. wave 가 결측인 500 행의 정체")
NAW <- D %>% filter(is.na(wave))
say("wave 결측 행 = %d, 고유 id = %d", nrow(NAW), n_distinct(NAW$id))
say("그 행들의 death_event 합 = %d", sum(NAW$death_event, na.rm=TRUE))
say("그 행들의 time 요약:"); print(summary(NAW$time))
say("그 행들의 death_wave 분포:"); print(table(NAW$death_wave, useNA="ifany"))
say("그 행들의 death_time_method:"); print(table(NAW$death_time_method, useNA="ifany"))
say("\n>> 500 행 / 500 id / death_event 500 이면 이 행들이 사망기록 행이고,")
say(">> 그 행의 time 이 원고가 쓴 실제 사망시점입니다.")
say("\ndeath_wave 별 time 분포 (사망기록 행):")
if (nrow(NAW)) print(round(do.call(rbind, tapply(NAW$time, NAW$death_wave,
    function(z) c(n=length(z), min=min(z), med=median(z), max=max(z)))), 3))

## ===========================================================================
rule("2. 세 가지 사망시점 규칙 구성")
VIS <- D %>% filter(!is.na(wave))
S <- D %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event > 0, na.rm=TRUE)),
            dwave = suppressWarnings(max(death_wave, na.rm=TRUE)),
            mth = paste(sort(unique(na.omit(as.character(death_time_method)))), collapse="+"),
            t_rec = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)),      # 사망기록 행의 time
            t_lastvisit = suppressWarnings(max(time[!is.na(wave)], na.rm=TRUE)),
            .groups="drop")
for (v in c("dwave","t_rec","t_lastvisit")) S[[v]][!is.finite(S[[v]])] <- NA

S$R1 <- ifelse(S$ev==1, S$t_rec,               S$t_lastvisit)  # 규칙1: 사망기록 행의 time
S$R2 <- ifelse(S$ev==1, (S$dwave-1)*2 - 1,     S$t_lastvisit)  # 규칙2: 구간 midpoint
S$R3 <- ifelse(S$ev==1, (S$dwave-1)*2,         S$t_lastvisit)  # 규칙3: 사망이 기록된 방문시점

for (r in c("R1","R2","R3")) {
  v <- S[[r]]
  say("\n[%s] 결측 %d | 전체 중앙값 %.2f, 최대 %.2f | person-years %.0f  [원고 22,785]",
      r, sum(is.na(v)), median(v,na.rm=TRUE), max(v,na.rm=TRUE), sum(v,na.rm=TRUE))
  say("   사망자 시점: "); print(summary(v[S$ev==1]))
}

## ===========================================================================
rule("3. 어느 규칙이 원고 Table 3 (5년 사망 2 / 23 / 26) 을 재현하는가")
K0 <- B %>% select(id, frailty_3cat, gLIC_z,
                   any_of(c("age_c","bage","sex_f","comorbid_count","comorbid_bl",
                            "edu_bl_f","income_bl_f","area_bl_f"))) %>%
  left_join(S %>% select(id, ev, R1, R2, R3), by="id")
cuts <- quantile(K0$gLIC_z, c(1/3,2/3), na.rm=TRUE)
K0$tert <- cut(K0$gLIC_z, c(-Inf,cuts,Inf), labels=c("T1","T2","T3"))
RB0 <- K0 %>% filter(as.character(frailty_3cat)=="0")

for (r in c("R1","R2","R3")) {
  v <- RB0[[r]]
  tb <- RB0 %>% mutate(tt=v) %>% group_by(tert) %>%
    summarise(n=n(), d5=sum(ev[tt<=5], na.rm=TRUE), d8=sum(ev[tt<=8], na.rm=TRUE), .groups="drop")
  say("\n[%s]  5년 사망: %s   (원고 2 / 23 / 26, 합 51)", r,
      paste(tb$d5, collapse=" / "))
  say("      8년 사망: %s", paste(tb$d8, collapse=" / "))
  if (identical(as.integer(tb$d5), c(2L,23L,26L))) say("      >>>>>> 원고 Table 3 완전 재현")
}

## ===========================================================================
rule("4. 확정 규칙으로 최종 분석")
pick <- if (identical(as.integer((RB0 %>% mutate(tt=R1) %>% group_by(tert) %>%
              summarise(d5=sum(ev[tt<=5],na.rm=TRUE),.groups="drop"))$d5), c(2L,23L,26L))) "R1" else
        if (identical(as.integer((RB0 %>% mutate(tt=R3) %>% group_by(tert) %>%
              summarise(d5=sum(ev[tt<=5],na.rm=TRUE),.groups="drop"))$d5), c(2L,23L,26L))) "R3" else "R1"
say("사용 규칙: %s   (재현 규칙이 없으면 R1 = 사망기록 행의 time 을 씁니다)", pick)

RB <- RB0 %>% mutate(t_fu = .data[[pick]]) %>% filter(!is.na(t_fu), t_fu > 0, !is.na(gLIC_z))
say("robust n = %d, 사망 = %d", nrow(RB), sum(RB$ev))

## 성별 구성 — 지난번 계산이 틀렸으므로 원자료를 그대로 봅니다
say("\nsex_f 의 실제 값:"); print(table(RB$sex_f, useNA="ifany"))
say("tertile x sex_f 교차표:"); print(table(RB$tert, RB$sex_f, useNA="ifany"))

km <- function(dd,h){ if(!nrow(dd)) return(NA_real_)
  s <- survfit(Surv(pmin(t_fu,h), ifelse(t_fu>h,0L,ev))~1, data=dd)
  100*(1-summary(s,times=h,extend=TRUE)$surv[1]) }
tb <- RB %>% group_by(tert) %>%
  summarise(n=n(), d5=sum(ev[t_fu<=5]), d8=sum(ev[t_fu<=8]), .groups="drop")
tb$km5 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],5))
tb$km8 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],8))
print(as.data.frame(tb), digits=3, row.names=FALSE)
write.csv(tb, file.path(OUT,"robust_tertiles_LOCKED.csv"), row.names=FALSE)

RB$.z  <- (RB$gLIC_z - mean(K0$gLIC_z,na.rm=TRUE))/sd(K0$gLIC_z,na.rm=TRUE)
RB <- RB %>% group_by(sex_f) %>% mutate(.zs = as.numeric(scale(gLIC_z))) %>% ungroup()
cvs <- intersect(c("age_c","sex_f","comorbid_count","edu_bl_f","income_bl_f","area_bl_f"), names(RB))

fit <- function(z, adj, h) {
  x <- RB; x$.t <- pmin(x$t_fu,h); x$.e <- ifelse(x$t_fu>h,0L,x$ev); x$Z <- x[[z]]
  f <- as.formula(paste("Surv(.t,.e) ~ Z", if(adj && length(cvs)) paste("+",paste(cvs,collapse="+")) else ""))
  m <- coxph(f, data=x); s <- summary(m)
  sprintf("HR %.3f (%.3f-%.3f), 사망 %3d, P = %.3g", s$coef["Z","exp(coef)"],
          s$conf.int["Z","lower .95"], s$conf.int["Z","upper .95"], m$nevent, s$coef["Z","Pr(>|z|)"])
}
say("")
for (h in c(5,8)) {
  say("%d년 per +1 s.d. (코호트 척도), 미보정 : %s", h, fit(".z", FALSE, h))
  say("%d년 per +1 s.d. (코호트 척도), 보정   : %s", h, fit(".z", TRUE , h))
  say("%d년 per +1 s.d. (성별 내 척도), 보정  : %s", h, fit(".zs", TRUE, h))
}

## 단계별 보정 — 어느 변수가 연관성을 만드는지
say("\n보정 변수를 하나씩 추가할 때 8년 HR 변화:")
x <- RB; x$.t <- pmin(x$t_fu,8); x$.e <- ifelse(x$t_fu>8,0L,x$ev)
acc <- character(0)
for (v in c("(없음)", cvs)) {
  if (v != "(없음)") acc <- c(acc, v)
  f <- as.formula(paste("Surv(.t,.e) ~ .z", if(length(acc)) paste("+",paste(acc,collapse="+")) else ""))
  s <- summary(coxph(f, data=x))
  say("  + %-16s HR %.3f (%.3f-%.3f)", v, s$coef[".z","exp(coef)"],
      s$conf.int[".z","lower .95"], s$conf.int[".z","upper .95"])
}
say(">>> HR 을 크게 움직이는 변수가 곧 이 연관성의 열쇠입니다.")

saveRDS(RB, file.path(OUT,"robust_frame_LOCKED.rds"))
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
