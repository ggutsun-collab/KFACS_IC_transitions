# =============================================================================
# 08 — 검증 통과. 이제 남은 것은 '사망 시점' 하나입니다.
#
# 통과한 것 (d + iv 가 원고의 정본임이 확정):
#   n 3,011 / robust 1,349 / pre-frail 1,416 / frail 246 / 사망 500
#   tertile n = 106 / 470 / 773  <- 원고 Table 3 과 정확히 일치
#
# 깨진 것:
#   death_event 가 사람 단위 상수라서 모든 person-wave 에 1 로 복제되어 있습니다.
#   (14,582 행 중 2,027 행이 1, 그러나 고유 id 는 500)
#   그래서 min(time[death_event>0]) 이 전부 wave-1 의 0 으로 잡혔습니다.
#   d$time 은 명목 방문시점(0/2/4/6/8)이지 추적시간이 아닙니다.
#
# 이 스크립트는 FUP / death_date_new / death_wave / age 로 추적시간을 재구성하고,
# 원고의 '중앙값 8.0년, 최대 8.2년, 22,785 person-years' 로 검증한 뒤 분석합니다.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT <- "rebuild_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "time_fix_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

D <- d
B <- readRDS(file.path(OUT, "wave1_frame_SOURCE.rds"))   # 검증 통과 프레임

## ===========================================================================
rule("1. 시간 후보 열 진단")
peek <- function(x, nm) {
  say("%-16s class=%-10s NA=%d", nm, paste(class(x), collapse="/"), sum(is.na(x)))
  if (is.numeric(x)) print(summary(x)) else print(head(sort(unique(as.character(x))), 8))
}
for (v in c("FUP","time","death_date_new","death_wave","death_time_method",
            "age_unified","age_visit","age0","age_365","wave"))
  if (v %in% names(D)) peek(D[[v]], v)

if ("FUP" %in% names(D)) {
  f <- D %>% group_by(id) %>% summarise(n_u = n_distinct(FUP), mx = max(FUP, na.rm=TRUE), .groups="drop")
  say("\nFUP 가 사람 단위 상수입니까? 고유값 1개인 id = %d / %d", sum(f$n_u==1), nrow(f))
  say("FUP 최대값 요약:"); print(summary(f$mx))
  say(">> 값이 ~8 이면 년, ~96 이면 개월, ~2900 이면 일 단위입니다.")
}

## ===========================================================================
rule("2. 추적시간 재구성 후보")
S <- D %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event > 0, na.rm=TRUE)),
            fup_max = suppressWarnings(max(FUP, na.rm=TRUE)),
            dwave   = suppressWarnings(max(death_wave, na.rm=TRUE)),
            mth     = paste(sort(unique(na.omit(as.character(death_time_method)))), collapse="+"),
            last_visit = suppressWarnings(max(time, na.rm=TRUE)),
            age_first = suppressWarnings(min(age_unified, na.rm=TRUE)),
            age_last  = suppressWarnings(max(age_unified, na.rm=TRUE)),
            .groups="drop")
for (v in c("fup_max","dwave","last_visit","age_first","age_last"))
  S[[v]][!is.finite(S[[v]])] <- NA

## 단위 자동 판별
u <- median(S$fup_max, na.rm=TRUE)
scale_f <- if (is.na(u)) NA else if (u > 1000) 365.25 else if (u > 50) 12 else 1
say("FUP 중앙값 = %.2f  -> 나눌 값 %s", u, ifelse(is.na(scale_f),"?",scale_f))
S$tA <- S$fup_max / scale_f                                   # 후보 A: FUP
S$tB <- ifelse(S$ev==1, (S$dwave-1)*2 - 1, S$last_visit)      # 후보 B: death_wave midpoint
S$tC <- S$age_last - S$age_first                              # 후보 C: 연령차

cand <- list(A="FUP", B="death_wave midpoint", C="연령차")
for (k in names(cand)) {
  v <- S[[paste0("t", k)]]
  say("\n[%s] %s", k, cand[[k]])
  say("  전체 중앙값 %.2f, 최대 %.2f, 음수 %d, 결측 %d   [원고: 중앙값 8.0, 최대 8.2]",
      median(v,na.rm=TRUE), max(v,na.rm=TRUE), sum(v<0,na.rm=TRUE), sum(is.na(v)))
  say("  총 person-years = %.0f   [원고 22,785]", sum(v,na.rm=TRUE))
  say("  사망자 시점 요약:"); print(summary(v[S$ev==1]))
}

## death_wave 와의 정합 — 가장 강한 검증
say("\ndeath_wave 별 재구성 시점 (사망자만). wave2~0-2y, 3~2-4y, 4~4-6y, 5~6-8y 여야 정상")
for (k in names(cand)) {
  say("\n[%s]", k)
  print(round(tapply(S[[paste0("t",k)]][S$ev==1], S$dwave[S$ev==1],
                     function(z) c(n=length(z), med=median(z,na.rm=TRUE),
                                   min=min(z,na.rm=TRUE), max=max(z,na.rm=TRUE))) %>%
          do.call(rbind, .), 2))
}
say("\n원고 death_wave 분포: 2=38, 3=90, 4=164, 5=208 (합 500)")

## ===========================================================================
rule("3. 최적 후보로 분석 — 조/보정 대비를 나란히")
pick <- if (!is.na(scale_f) && abs(median(S$tA,na.rm=TRUE) - 8) < 1.5) "tA" else
        if (abs(median(S$tB,na.rm=TRUE) - 8) < 1.5) "tB" else "tC"
say("선택된 시간 변수: %s   (다르게 하시려면 pick 을 바꾸십시오)", pick)
S$t_fu <- S[[pick]]

K <- B %>% select(id, frailty_3cat, gLIC_z, any_of(c("age_c","bage","sex_f","comorbid_count",
        "comorbid_bl","edu_bl_f","income_bl_f","area_bl_f"))) %>%
  left_join(S %>% select(id, ev, t_fu), by="id") %>%
  filter(!is.na(t_fu), t_fu > 0, !is.na(gLIC_z))
say("분석 프레임 n = %d, 사망 = %d", nrow(K), sum(K$ev))

RB <- K %>% filter(as.character(frailty_3cat)=="0")
say("robust n = %d, 사망 = %d   [원고 1,349]", nrow(RB), sum(RB$ev))

cuts <- quantile(K$gLIC_z, c(1/3,2/3), na.rm=TRUE)
RB$tert <- cut(RB$gLIC_z, c(-Inf,cuts,Inf), labels=c("T1","T2","T3"))
RB$.z <- (RB$gLIC_z - mean(K$gLIC_z,na.rm=TRUE))/sd(K$gLIC_z,na.rm=TRUE)
RB <- RB %>% group_by(sex_f) %>% mutate(.zs = as.numeric(scale(gLIC_z))) %>% ungroup()

km <- function(dd,h){ if(!nrow(dd)) return(NA_real_)
  s <- survfit(Surv(pmin(t_fu,h), ifelse(t_fu>h,0L,ev))~1, data=dd)
  100*(1-summary(s,times=h,extend=TRUE)$surv[1]) }
tb <- RB %>% group_by(tert) %>%
  summarise(n=n(), pct_female=100*mean(as.numeric(as.character(sex_f))==1 |
              tolower(as.character(sex_f)) %in% c("f","female","여"), na.rm=TRUE),
            d5=sum(ev[t_fu<=5]), d8=sum(ev[t_fu<=8]), .groups="drop")
tb$km5 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],5))
tb$km8 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],8))
print(as.data.frame(tb), digits=3, row.names=FALSE)
write.csv(tb, file.path(OUT,"robust_tertiles_FINAL.csv"), row.names=FALSE)
say("원고 Table 3 : T1 n=106 d=2 (1.9%%) | T2 n=470 d=23 (4.9%%) | T3 n=773 d=26 (3.4%%)")

cvs <- intersect(c("age_c","sex_f","comorbid_count","edu_bl_f","income_bl_f","area_bl_f"), names(RB))
fit <- function(z, adj, h) {
  dta <- RB; dta$.t <- pmin(dta$t_fu,h); dta$.e <- ifelse(dta$t_fu>h,0L,dta$ev); dta$Z <- dta[[z]]
  f <- as.formula(paste("Surv(.t,.e) ~ Z", if(adj && length(cvs)) paste("+",paste(cvs,collapse="+")) else ""))
  s <- summary(coxph(f, data=dta))
  sprintf("HR %.3f (%.3f-%.3f), P = %.3g", s$coef["Z","exp(coef)"],
          s$conf.int["Z","lower .95"], s$conf.int["Z","upper .95"], s$coef["Z","Pr(>|z|)"])
}
say("\n%-46s %s", "5년 per +1 s.d., 미보정", fit(".z", FALSE, 5))
say("%-46s %s", "5년 per +1 s.d., 보정",   fit(".z", TRUE , 5))
say("%-46s %s", "8년 per +1 s.d., 미보정", fit(".z", FALSE, 8))
say("%-46s %s", "8년 per +1 s.d., 보정",   fit(".z", TRUE , 8))
say("%-46s %s", "8년 per +1 성별내 s.d., 보정", fit(".zs", TRUE, 8))
say("\n>>> 미보정과 보정의 차이가 크면, 이 연관성은 성별·연령 보정에 의존합니다.")
say(">>> tertile 여성비(위 표 pct_female)가 극단적이면 심사자가 반드시 지적합니다.")

saveRDS(RB, file.path(OUT,"robust_frame_FINAL.rds"))
saveRDS(K,  file.path(OUT,"wave1_frame_FINAL.rds"))
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
