# =============================================================================
# 05 — 찾았습니다. 원고의 1,349 = d$frailty_3cat == 0
#      이제 어느 정의가 옳은지 판정하고, 사망 500 규칙을 마저 찾습니다.
#
#   wave 1        robust   pre-frail   frail
#   frailty_3cat   1349      1416       246   (8.2% frail)   <- 원고
#   frailty_w1     1082      1469       460   (15.3% frail)  <- 재분석
#
# frail 유병률이 거의 2배 차이납니다. 둘 다 맞을 수는 없습니다.
#
# 그대로 source. 편집할 곳 없습니다.
# ※ 앞 스크립트가 에러로 죽으면서 sink 가 열려 있을 수 있습니다.
#    콘솔에 아무것도 안 보이면 sink() 를 두어 번 실행해 닫은 뒤 이 파일을 돌리십시오.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT <- "adjudicate_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "adjudicate_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")
hit  <- function(lab, n, target = 500)
  say("  %-52s = %4d %s", lab, n, if (isTRUE(n == target)) sprintf("  <<<<<< %d 일치", target) else "")

W <- dat_w1; W$id <- trimws(as.character(W$id))
D <- d;      D$id <- trimws(as.character(D$id))
w1 <- D %>% group_by(id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>% ungroup()

## ===========================================================================
rule("A. 두 frailty 정의의 정면 대조")

K <- W %>% select(id, d_death, t_death, gLIC_z, adl_w1, frailty_w1,
                  any_of(c("age","sex","comorb_n","educ","income","area"))) %>%
  left_join(w1 %>% select(id, frailty_3cat, any_of("chs_total")), by = "id")

say("교차표  (행: frailty_w1 / 열: frailty_3cat  0=robust 1=pre-frail 2=frail)")
print(table(K$frailty_w1, K$frailty_3cat, useNA = "ifany"))
say("\n완전일치율 = %.1f%%",
    100*mean(as.integer(factor(K$frailty_w1, levels=c("Robust","Pre-frail","Frail"))) - 1 ==
             as.numeric(as.character(K$frailty_3cat)), na.rm = TRUE))

## chs_total 로 직접 재구성 — 이것이 심판입니다
if ("chs_total" %in% names(K)) {
  say("\nwave-1 chs_total (Fried 기준 충족 개수) 분포:")
  print(table(K$chs_total, useNA = "ifany"))
  K$chs3 <- cut(as.numeric(K$chs_total), c(-0.5, 0.5, 2.5, 5.5),
                labels = c("Robust","Pre-frail","Frail"))
  say("\nchs_total 로 재구성한 분류 (0 / 1-2 / 3-5):")
  print(table(K$chs3, useNA = "ifany"))
  say("\nchs3 vs frailty_3cat:"); print(table(K$chs3, K$frailty_3cat, useNA = "ifany"))
  say("\nchs3 vs frailty_w1:");   print(table(K$chs3, K$frailty_w1,   useNA = "ifany"))
  a1 <- mean(as.integer(K$chs3) - 1 == as.numeric(as.character(K$frailty_3cat)), na.rm = TRUE)
  a2 <- mean(as.character(K$chs3) == as.character(K$frailty_w1), na.rm = TRUE)
  say("\nchs_total 재구성과의 일치율 : frailty_3cat %.1f%%  |  frailty_w1 %.1f%%", 100*a1, 100*a2)
  say(">>> 일치율이 높은 쪽이 Fried 기준을 그대로 구현한 변수입니다.")
  say(">>> 낮은 쪽은 결측 처리나 절단점이 다르게 들어간 파생변수일 가능성이 큽니다.")
}

say("\n참고: KFACS 및 지역사회 거주 70-84세 한국인 코호트의 기존 보고에서")
say("      Fried frailty 유병률은 대체로 7-10%% 대입니다.")
say("      frailty_3cat = 8.2%% / frailty_w1 = 15.3%%.")
say(">>> 외부 문헌과의 정합성은 frailty_3cat 쪽이 뚜렷이 높습니다.")

## 어느 쪽이 결측을 frail 로 밀어넣었는지
disagree <- K %>% filter(!is.na(frailty_3cat),
                         as.integer(factor(frailty_w1, levels=c("Robust","Pre-frail","Frail")))-1 !=
                         as.numeric(as.character(frailty_3cat)))
say("\n두 정의가 다른 인원 = %d 명", nrow(disagree))
if (nrow(disagree)) {
  print(table(disagree$frailty_w1, disagree$frailty_3cat))
  write.csv(disagree, file.path(OUT, "frailty_disagreement.csv"), row.names = FALSE)
}

## ===========================================================================
rule("B. 사망 500 규칙 재탐색  [tb 인덱싱 오류 수정]")

say("dat_w1 사망 = %d", sum(W$d_death))
tt <- sort(W$t_death[W$d_death == 1])
say("500번째 사망시점 = %.4f 년 / 501번째 = %.4f 년", tt[500], tt[min(501, length(tt))])
for (h in c(5, 6, 6.5, 7, 7.2, 7.5, 7.8, 8, 8.1, 8.2))
  hit(sprintf("t_death <= %.1f 년", h), sum(W$d_death[W$t_death <= h]))

say("\n사망시점 십분위별 건수:")
print(table(cut(W$t_death[W$d_death==1],
                breaks = quantile(W$t_death[W$d_death==1], 0:10/10), include.lowest = TRUE)))

P <- D %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event > 0, na.rm = TRUE)),
            mth = paste(sort(unique(na.omit(as.character(death_time_method)))), collapse="+"),
            dwav = suppressWarnings(max(as.numeric(death_wave), na.rm = TRUE)),
            has_date = any(!is.na(death_date_new)), .groups = "drop")
P$dwav[!is.finite(P$dwav)] <- NA

hit("d$death_event (id 단위)",           sum(P$ev))
hit("  + death_date_new 존재",           sum(P$ev & P$has_date))
hit("  + method = reported_month 만",    sum(P$ev & P$mth == "reported_month"))
hit("  + method = midpoint 만",          sum(P$ev & P$mth == "midpoint"))
hit("  + method 있음",                   sum(P$ev & nzchar(P$mth)))
say("\n사망자 method 조합:"); print(table(P$mth[P$ev==1], useNA="ifany"))
say("사망자 death_wave 분포:"); print(table(P$dwav[P$ev==1], useNA="ifany"))
dw <- sort(unique(na.omit(P$dwav)))
for (k in dw) hit(sprintf("death_wave <= %s", k), sum(P$ev & P$dwav <= k, na.rm=TRUE))

M <- W %>% select(id, d_death, t_death, frailty_w1) %>%
  left_join(P, by="id") %>% mutate(ev = ifelse(is.na(ev), 0L, ev))
extra <- M %>% filter(d_death == 1 & ev == 0)
say("\ndat_w1 에서만 사망인 id = %d 명", nrow(extra))
if (nrow(extra)) { print(summary(extra$t_death)); print(table(extra$frailty_w1, useNA="ifany"))
                   write.csv(extra, file.path(OUT,"extra_deaths.csv"), row.names=FALSE) }

## ===========================================================================
rule("C. 미리보기 — 원고 정의(frailty_3cat==0, n=1349)로 다시 계산하면")

RB <- K %>% filter(frailty_3cat == 0 | frailty_3cat == "0", !is.na(gLIC_z))
say("robust n = %d, 사망 = %d", nrow(RB), sum(RB$d_death))
cuts <- quantile(W$gLIC_z, c(1/3, 2/3), na.rm = TRUE)
RB$tert <- cut(RB$gLIC_z, c(-Inf, cuts, Inf), labels = c("T1","T2","T3"))
km <- function(dat,h){ s<-survfit(Surv(pmin(t_death,h), ifelse(t_death>h,0L,d_death))~1, data=dat)
                       100*(1-summary(s,times=h,extend=TRUE)$surv[1]) }
tb <- RB %>% group_by(tert) %>%
  summarise(n=n(), d5=sum(d_death[t_death<=5]), d8=sum(d_death[t_death<=8]), .groups="drop") %>%
  rowwise() %>% mutate(km5=km(RB[RB$tert==tert,],5), km8=km(RB[RB$tert==tert,],8)) %>% ungroup()
print(as.data.frame(tb), digits=3, row.names=FALSE)
say("원고 Table 3 : T1 n=106 d=2 (1.9%%) | T2 n=470 d=23 (4.9%%) | T3 n=773 d=26 (3.4%%)")
say(">>> n 은 맞아떨어집니까? 사망은? tertile 절단점까지 원고와 같은지 확인하십시오.")

cvs <- intersect(c("age","sex","comorb_n","educ","income","area"), names(RB))
RB$.z <- (RB$gLIC_z - mean(W$gLIC_z, na.rm=TRUE)) / sd(W$gLIC_z, na.rm=TRUE)
for (h in c(5,8)) {
  RB$.t <- pmin(RB$t_death, h); RB$.e <- ifelse(RB$t_death > h, 0L, RB$d_death)
  m <- coxph(as.formula(paste("Surv(.t,.e) ~ .z +", paste(cvs, collapse="+"))), data = RB)
  s <- summary(m)
  say("%d년 per +1 s.d. : HR %.3f (%.3f-%.3f), 사망 %d / %d, P = %.3g", h,
      s$coef[".z","exp(coef)"], s$conf.int[".z","lower .95"], s$conf.int[".z","upper .95"],
      m$nevent, m$n, s$coef[".z","Pr(>|z|)"])
}
say("\n※ 사망 500 규칙이 아직 미확정이므로 위 숫자도 잠정입니다.")
say("   다만 원고 정의로 바꾸었을 때 방향과 크기가 유지되는지는 여기서 이미 보입니다.")

rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
