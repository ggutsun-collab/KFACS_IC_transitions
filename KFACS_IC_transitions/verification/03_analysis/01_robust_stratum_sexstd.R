# =============================================================================
# 10 — 데이터 문제 종결. 최종 분석.
#
# 확정 (원고 완전 재현):
#   n 3,011 / robust 1,349 / tertile 106-470-773 / 사망 500
#   사망시점 = wave 결측 행(사망기록 행)의 time  -> Table 3 의 2 / 23 / 26 완전 재현
#
# 확정된 문제 (이번 분석의 핵심):
#   robust 층 내 IC-사망 연관성은 성별 보정에 전적으로 의존합니다.
#     미보정        HR 0.830 (0.650-1.059)  P = 0.13
#     + 연령        HR 0.919 (0.718-1.175)  P > 0.05
#     + 성별        HR 0.558 (0.425-0.733)  <-- 성별 하나가 전부
#   T1 남성 10 / 여성 96,  T3 남성 586 / 여성 187.
#   -> 성별 내 표준화 IC 를 primary 로 올려야 합니다.
#
# 이 스크립트가 원고 수정에 필요한 수치를 전부 생성합니다.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({
  for (p in c("survival","splines","dplyr")) library(p, character.only = TRUE) })
OUT <- "final_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "FINAL_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

## ===========================================================================
rule("0. 확정 프레임 재구성")
D <- d
W1 <- D %>% filter(!is.na(wave)) %>% group_by(id) %>% slice_min(wave, n=1, with_ties=FALSE) %>% ungroup()
S <- D %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event > 0, na.rm=TRUE)),
            t_death = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)),
            t_last  = suppressWarnings(max(time[!is.na(wave)], na.rm=TRUE)), .groups="drop") %>%
  mutate(t_death = ifelse(is.finite(t_death), t_death, NA),
         t_last  = ifelse(is.finite(t_last),  t_last,  NA),
         t_fu    = ifelse(ev==1, t_death, t_last))
IV1 <- iv %>% group_by(id) %>% slice_min(k, n=1, with_ties=FALSE) %>% ungroup() %>%
  select(id, gLIC_z, age_c, sex_f, comorbid_count, edu_bl_f, income_bl_f, area_bl_f)

## 복합악화: ADL 악화 또는 사망 중 먼저 오는 것
ADL <- D %>% filter(!is.na(wave)) %>% arrange(id, wave) %>% group_by(id) %>%
  mutate(base = first(adl_3cat)) %>%
  summarise(t_worse = { w <- which(adl_3cat > base & wave > min(wave))
                        if (length(w)) time[w[1]] else NA_real_ }, .groups="drop")

K <- W1 %>% select(id, frailty_3cat, adl_3cat) %>%
  left_join(S %>% select(id, ev, t_fu), by="id") %>%
  left_join(IV1, by="id") %>% left_join(ADL, by="id") %>%
  filter(!is.na(gLIC_z), !is.na(t_fu), t_fu > 0) %>%
  mutate(e_comp = as.integer(!is.na(t_worse) | ev==1),
         t_comp = pmin(ifelse(is.na(t_worse), Inf, t_worse), t_fu))
say("n = %d | robust = %d | 사망 = %d | 복합악화 = %d",
    nrow(K), sum(K$frailty_3cat==0), sum(K$ev), sum(K$e_comp))

## 표준화 두 가지
K$z_coh <- as.numeric(scale(K$gLIC_z))
K <- K %>% group_by(sex_f) %>% mutate(z_sex = as.numeric(scale(gLIC_z))) %>% ungroup()
cut_coh <- quantile(K$gLIC_z, c(1/3,2/3), na.rm=TRUE)
K$t_coh <- cut(K$gLIC_z, c(-Inf, cut_coh, Inf), labels=c("T1","T2","T3"))
K <- K %>% group_by(sex_f) %>%
  mutate(t_sex = cut(gLIC_z, c(-Inf, quantile(gLIC_z, c(1/3,2/3), na.rm=TRUE), Inf),
                     labels=c("T1","T2","T3"))) %>% ungroup()
CV <- c("age_c","sex_f","comorbid_count","edu_bl_f","income_bl_f","area_bl_f")
CVns <- setdiff(CV, "sex_f")   # 성별 내 척도에서는 성별을 층화로 다룸

## ===========================================================================
rule("1. robust 층 — 척도별 / 보정별 per-s.d. 추정")
RB <- K %>% filter(frailty_3cat == 0)
say("robust n = %d, 사망 = %d, 복합악화 = %d", nrow(RB), sum(RB$ev), sum(RB$e_comp))

est <- function(dat, z, covs, tv, evv, h) {
  x <- dat; x$.t <- pmin(x[[tv]], h); x$.e <- ifelse(x[[tv]] > h, 0L, x[[evv]]); x$Z <- x[[z]]
  f <- as.formula(paste("Surv(.t,.e) ~ Z", if(length(covs)) paste("+", paste(covs, collapse="+")) else ""))
  m <- coxph(f, data=x); s <- summary(m)
  data.frame(HR=s$coef["Z","exp(coef)"], lo=s$conf.int["Z","lower .95"],
             hi=s$conf.int["Z","upper .95"], ev=m$nevent, p=s$coef["Z","Pr(>|z|)"])
}
grid <- expand.grid(h=c(5,8), scale=c("z_coh","z_sex"), adj=c(FALSE,TRUE), stringsAsFactors=FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i,]; covs <- if (!g$adj) character(0) else if (g$scale=="z_sex") CVns else CV
  cbind(outcome="death", g, est(RB, g$scale, covs, "t_fu", "ev", g$h)) }))
res2 <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i,]; covs <- if (!g$adj) character(0) else if (g$scale=="z_sex") CVns else CV
  cbind(outcome="worsening", g, est(RB, g$scale, covs, "t_comp", "e_comp", g$h)) }))
ALL <- rbind(res, res2)
print(ALL, digits=3, row.names=FALSE)
write.csv(ALL, file.path(OUT,"perSD_estimates.csv"), row.names=FALSE)
say("\n>>> 성별 내 척도(z_sex)의 '미보정' 행을 보십시오.")
say(">>> 여기서도 유의하면, 임상 주장(진료실에서 보이는 관계)이 성립합니다.")
say(">>> 여기서도 비유의하면, 사망에 대한 임상 주장은 철회해야 합니다.")

## ===========================================================================
rule("2. Table 3 대체안 — 성별 내 tertile, 조위험과 보정 HR 병기")
km <- function(dd,tv,evv,h){ if(!nrow(dd)) return(NA_real_)
  x<-dd; x$.t<-pmin(x[[tv]],h); x$.e<-ifelse(x[[tv]]>h,0L,x[[evv]])
  s<-survfit(Surv(.t,.e)~1,data=x); 100*(1-summary(s,times=h,extend=TRUE)$surv[1]) }

for (sc in c("t_coh","t_sex")) {
  say("\n---- tertile 정의: %s ----", ifelse(sc=="t_coh","코호트 전체 절단점 (원고 방식)","성별 내 절단점 (제안)"))
  RB$G <- RB[[sc]]
  tb <- RB %>% group_by(G) %>%
    summarise(n=n(), male=sum(as.character(sex_f)=="Male"), female=sum(as.character(sex_f)=="Female"),
              d5=sum(ev[t_fu<=5]), d8=sum(ev[t_fu<=8]),
              w5=sum(e_comp[t_comp<=5]), .groups="drop")
  tb$km5  <- sapply(tb$G, function(g) km(RB[RB$G==g,],"t_fu","ev",5))
  tb$km8  <- sapply(tb$G, function(g) km(RB[RB$G==g,],"t_fu","ev",8))
  tb$kmw5 <- sapply(tb$G, function(g) km(RB[RB$G==g,],"t_comp","e_comp",5))
  print(as.data.frame(tb), digits=3, row.names=FALSE)
  RB$G <- relevel(factor(RB$G), ref="T3")
  for (h in c(5,8)) {
    x <- RB; x$.t <- pmin(x$t_fu,h); x$.e <- ifelse(x$t_fu>h,0L,x$ev)
    cu <- summary(coxph(Surv(.t,.e) ~ G, data=x))$conf.int
    covs <- if (sc=="t_sex") CVns else CV
    ad <- summary(coxph(as.formula(paste("Surv(.t,.e) ~ G +", paste(covs,collapse="+"))), data=x))$conf.int
    say("  %d년  T1 vs T3: 조 %.2f (%.2f-%.2f) | 보정 %.2f (%.2f-%.2f)", h,
        cu[1,1],cu[1,3],cu[1,4], ad[1,1],ad[1,3],ad[1,4])
    say("  %d년  T2 vs T3: 조 %.2f (%.2f-%.2f) | 보정 %.2f (%.2f-%.2f)", h,
        cu[2,1],cu[2,3],cu[2,4], ad[2,1],ad[2,3],ad[2,4])
  }
  write.csv(tb, file.path(OUT, paste0("table3_", sc, ".csv")), row.names=FALSE)
}

## ===========================================================================
rule("3. 성별 층화 — 남녀 각각에서 연관성이 있는가")
for (sx in unique(as.character(RB$sex_f))) {
  Z <- RB %>% filter(as.character(sex_f)==sx)
  e <- est(Z, "z_sex", setdiff(CVns,"sex_f"), "t_fu","ev",8)
  say("  %-8s n=%4d 사망=%3d  보정 HR %.3f (%.3f-%.3f) P=%.3g",
      sx, nrow(Z), sum(Z$ev), e$HR, e$lo, e$hi, e$p)
}
x <- RB; x$.t <- pmin(x$t_fu,8); x$.e <- ifelse(x$t_fu>8,0L,x$ev)
m0 <- coxph(as.formula(paste("Surv(.t,.e) ~ z_sex +", paste(CV,collapse="+"))), data=x)
m1 <- coxph(as.formula(paste("Surv(.t,.e) ~ z_sex*sex_f +", paste(CVns,collapse="+"))), data=x)
say("  IC x 성별 상호작용 P = %.3f", anova(m0,m1)$`Pr(>|Chi|)`[2])

## ===========================================================================
rule("4. 비선형성 (성별 내 척도, 8년, 보정)")
x <- RB; x$.t <- pmin(x$t_fu,8); x$.e <- ifelse(x$t_fu>8,0L,x$ev)
ml <- coxph(as.formula(paste("Surv(.t,.e) ~ z_sex +", paste(CVns,collapse="+"))), data=x)
mn <- coxph(as.formula(paste("Surv(.t,.e) ~ ns(z_sex,3) +", paste(CVns,collapse="+"))), data=x)
say("비선형성 LRT P = %.4f", anova(ml,mn)$`Pr(>|Chi|)`[2])
qs <- quantile(x$z_sex, c(.05,.10,.25,.50,.75,.90,.95), na.rm=TRUE)
nd <- data.frame(z_sex=as.numeric(qs))
for (v in CVns) { u<-x[[v]]; nd[[v]] <- if (is.numeric(u)) mean(u,na.rm=TRUE) else
  factor(names(sort(table(u),decreasing=TRUE))[1], levels=levels(factor(u))) }
lp <- predict(mn, newdata=nd, type="lp"); ref <- lp[4]
print(data.frame(pct=names(qs), HR_vs_median=round(exp(lp-ref),3)), row.names=FALSE)

## ===========================================================================
rule("5. PH 가정 · E-value · 절대위험")
print(cox.zph(ml))
ev8 <- est(RB, "z_sex", CVns, "t_fu","ev",8)
E <- function(r){ r<-ifelse(r<1,1/r,r); r+sqrt(r*(r-1)) }
say("E-value: 점추정 %.2f, CI 한계 %.2f", E(ev8$HR), E(ev8$hi))
sf <- survfit(ml, newdata=nd)
say("보정 8년 절대사망위험 (IC 백분위별):")
print(data.frame(pct=names(qs), risk=round(100*(1-summary(sf,times=8,extend=TRUE)$surv[1,]),1)), row.names=FALSE)

saveRDS(K,  file.path(OUT,"K_FINAL.rds")); saveRDS(RB, file.path(OUT,"RB_FINAL.rds"))
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
