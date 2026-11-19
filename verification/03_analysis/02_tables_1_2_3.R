# =============================================================================
# 12 — 성별 내 표준화 IC 로 전면 전환. Table 1/2/3 + 회복분석 일괄 재생성.
#
# 결정 사항 반영
#   primary 노출 = 성별 내 표준화 gLIC_z  (인자점수 단계에서 z 변환)
#   주 지평선    = 5년 (8년은 보조로 함께 출력)
#   측정모형     = 재적합하지 않음
#
# 이 스크립트 하나로 원고 수정에 필요한 모든 수치가 나옵니다. 마지막 실행입니다.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(survival); library(sandwich); library(lmtest); library(dplyr) })
OUT <- "tables_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "TABLES_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

## ===========================================================================
rule("0. 프레임 준비")
D <- d
W1 <- D %>% filter(!is.na(wave)) %>% group_by(id) %>% slice_min(wave,n=1,with_ties=FALSE) %>% ungroup()
S  <- D %>% group_by(id) %>%
  summarise(ev=as.integer(any(death_event>0,na.rm=TRUE)),
            td=suppressWarnings(max(time[is.na(wave)],na.rm=TRUE)),
            tl=suppressWarnings(max(time[!is.na(wave)],na.rm=TRUE)),.groups="drop") %>%
  mutate(td=ifelse(is.finite(td),td,NA), tl=ifelse(is.finite(tl),tl,NA),
         t_fu=ifelse(ev==1,td,tl))
ADL <- D %>% filter(!is.na(wave)) %>% arrange(id,wave) %>% group_by(id) %>%
  mutate(base=first(adl_3cat)) %>%
  summarise(t_worse={w<-which(adl_3cat>base & wave>min(wave)); if(length(w)) time[w[1]] else NA_real_},
            .groups="drop")
IV1 <- iv %>% group_by(id) %>% slice_min(k,n=1,with_ties=FALSE) %>% ungroup() %>%
  select(id,gLIC_z,age_c,sex_f,comorbid_count,edu_bl_f,income_bl_f,area_bl_f)

K <- W1 %>% select(id,frailty_3cat,adl_3cat) %>%
  left_join(S %>% select(id,ev,t_fu),by="id") %>% left_join(IV1,by="id") %>%
  left_join(ADL,by="id") %>% filter(!is.na(gLIC_z),!is.na(t_fu),t_fu>0) %>%
  mutate(e_comp=as.integer(!is.na(t_worse)|ev==1),
         t_comp=pmin(ifelse(is.na(t_worse),Inf,t_worse),t_fu),
         fr=factor(as.character(frailty_3cat),levels=c("0","1","2"),
                   labels=c("Robust","Pre-frail","Frail")))
K <- K %>% group_by(sex_f) %>%
  mutate(zS=as.numeric(scale(gLIC_z)),
         tS=cut(gLIC_z,c(-Inf,quantile(gLIC_z,c(1/3,2/3),na.rm=TRUE),Inf),labels=c("T1","T2","T3"))) %>%
  ungroup()
K$tS <- relevel(factor(K$tS), ref="T3")
say("n=%d | robust %d / pre-frail %d / frail %d | 사망 %d",
    nrow(K), sum(K$fr=="Robust"), sum(K$fr=="Pre-frail"), sum(K$fr=="Frail"), sum(K$ev))
say("성별 내 tertile 절단점 (성별로 다름) — 남/여 각각:")
print(K %>% group_by(sex_f) %>% summarise(c1=quantile(gLIC_z,1/3), c2=quantile(gLIC_z,2/3), .groups="drop"))

## ===========================================================================
rule("1. Table 1 대체 — 성별 내 tertile 의 구성")
t1 <- K %>% group_by(tS) %>% summarise(n=n(),
        male=sum(as.character(sex_f)=="Male"), female=sum(as.character(sex_f)=="Female"),
        pct_female=100*mean(as.character(sex_f)=="Female"),
        age=mean(age_c,na.rm=TRUE), comorb=mean(comorbid_count,na.rm=TRUE), .groups="drop")
print(as.data.frame(t1), digits=3, row.names=FALSE)
say(">>> 원고 Table 1 은 T1 여성 77.3%%, T3 20.9%% 였습니다. 위 pct_female 로 교체하십시오.")
write.csv(t1, file.path(OUT,"table1_tertile_composition.csv"), row.names=FALSE)

## ===========================================================================
rule("2. Table 3 전체 재생성 — 9셀 x 2결과, 5년 (8년 병기)")
km <- function(dd,tv,evv,h){ if(!nrow(dd)) return(c(NA,NA,NA))
  x<-dd; x$.t<-pmin(x[[tv]],h); x$.e<-ifelse(x[[tv]]>h,0L,x[[evv]])
  s<-summary(survfit(Surv(.t,.e)~1,data=x),times=h,extend=TRUE)
  100*c(1-s$surv[1], 1-s$upper[1], 1-s$lower[1]) }
CV <- c("age_c","comorbid_count","edu_bl_f","income_bl_f","area_bl_f")   # 성별 내 척도이므로 sex 제외

mk <- function(tv, evv, lab) {
  out <- K %>% group_by(fr, tS) %>%
    summarise(n=n(), ev5=sum(.data[[evv]][.data[[tv]]<=5]),
              ev8=sum(.data[[evv]][.data[[tv]]<=8]), .groups="drop")
  r5 <- t(mapply(function(f,g) km(K[K$fr==f & K$tS==g,],tv,evv,5), out$fr, out$tS))
  r8 <- t(mapply(function(f,g) km(K[K$fr==f & K$tS==g,],tv,evv,8), out$fr, out$tS))
  out$risk5 <- r5[,1]; out$lo5 <- r5[,2]; out$hi5 <- r5[,3]; out$risk8 <- r8[,1]
  # 전역 기준 = Robust/T3
  K$G <- interaction(K$fr, K$tS, sep="/"); K$G <- relevel(factor(K$G), ref="Robust/T3")
  x <- K; x$.t <- pmin(x[[tv]],5); x$.e <- ifelse(x[[tv]]>5,0L,x[[evv]])
  cu <- try(summary(coxph(Surv(.t,.e)~G, data=x))$conf.int, silent=TRUE)
  ad <- try(summary(coxph(as.formula(paste("Surv(.t,.e)~G +",paste(CV,collapse="+"))),data=x))$conf.int, silent=TRUE)
  gl <- paste0("G", levels(x$G)[-1])
  gethr <- function(m,nm) if (inherits(m,"try-error")||!(nm %in% rownames(m))) c(NA,NA,NA) else m[nm,c(1,3,4)]
  key <- paste(out$fr, out$tS, sep="/")
  H <- t(sapply(paste0("G",key), function(nm) gethr(cu,nm)))
  A <- t(sapply(paste0("G",key), function(nm) gethr(ad,nm)))
  out$HRcrude <- H[,1]; out$cl <- H[,2]; out$ch <- H[,3]
  out$HRadj   <- A[,1]; out$al <- A[,2]; out$ah <- A[,3]
  say("\n---- %s ----", lab)
  print(as.data.frame(out), digits=3, row.names=FALSE)
  write.csv(out, file.path(OUT, paste0("table3_", evv, ".csv")), row.names=FALSE)
  out
}
T3w <- mk("t_comp","e_comp","Composite worsening (5년 위험, Robust/T3 기준)")
T3d <- mk("t_fu",  "ev",    "Death, all cause (5년 위험, Robust/T3 기준)")

## ===========================================================================
rule("3. Table 2 재생성 — 성별 내 표준화 per-s.d. 로 18개 전이")
IV <- iv %>% left_join(K %>% select(id) %>% distinct(), by="id")
IV <- IV %>% group_by(sex_f) %>% mutate(zS=as.numeric(scale(gLIC_z))) %>% ungroup()
say("iv 행 %d | type: %s", nrow(IV), paste(unique(IV$type), collapse=", "))
say("from_lab -> to_lab 조합 수 = %d", nrow(distinct(IV, type, from_lab, to_lab)))

trs <- IV %>% distinct(type, from_lab, to_lab) %>% arrange(type, from_lab, to_lab)
res <- lapply(seq_len(nrow(trs)), function(i) {
  tr <- trs[i,]
  dat <- IV %>% filter(type==tr$type, from_lab==tr$from_lab)
  if (!nrow(dat)) return(NULL)
  dat$y <- as.integer(dat$to_lab == tr$to_lab)
  if (sum(dat$y) < 5) return(data.frame(type=tr$type, from=tr$from_lab, to=tr$to_lab,
      events=sum(dat$y), py=sum(dat$dur), IRR=NA, lo=NA, hi=NA, p=NA))
  f <- as.formula(paste("y ~ zS + age_c + comorbid_count + edu_bl_f + income_bl_f + area_bl_f + offset(log(dur))"))
  m <- try(glm(f, family=poisson, data=dat), silent=TRUE)
  if (inherits(m,"try-error")) return(NULL)
  ct <- coeftest(m, vcov=vcovCL(m, cluster=dat$id))
  b <- ct["zS","Estimate"]; se <- ct["zS","Std. Error"]
  data.frame(type=tr$type, from=tr$from_lab, to=tr$to_lab, events=sum(dat$y), py=sum(dat$dur),
             IRR=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=ct["zS","Pr(>|z|)"])
})
TAB2 <- bind_rows(res)
print(as.data.frame(TAB2), digits=3, row.names=FALSE)
write.csv(TAB2, file.path(OUT,"table2_perSD_sexstd.csv"), row.names=FALSE)
say(">>> 원고 Table 2 의 IRR 을 위 값으로 교체하십시오 (코호트 척도 -> 성별 내 척도).")

## ===========================================================================
rule("4. 회복 분석 — robust + ADL 장애 보유자")
RC <- K %>% filter(fr=="Robust", as.character(adl_3cat)!="0")
say("n = %d  [원고 103]", nrow(RC))
if (nrow(RC) > 20 && "t_worse" %in% names(RC)) {
  say("회복 분석은 별도 변수(t_recov/e_recov)가 필요합니다.")
  say("원고의 정의(wave-1 장애 보유자의 정상 ADL 복귀, 사망을 경쟁위험)를 쓰는 코드를")
  say("기존 스크립트에서 가져와 zS 로만 바꿔 다시 돌리십시오.")
  say("성별 구성: "); print(table(RC$sex_f))
}

## ===========================================================================
rule("5. 초록용 핵심 수치 요약")
RB <- K %>% filter(fr=="Robust")
e <- function(z,adj,tv,evv,h){ x<-RB; x$.t<-pmin(x[[tv]],h); x$.e<-ifelse(x[[tv]]>h,0L,x[[evv]]); x$Z<-x[[z]]
  f<-as.formula(paste("Surv(.t,.e)~Z",if(adj) paste("+",paste(CV,collapse="+")) else ""))
  s<-summary(coxph(f,data=x)); sprintf("%.2f (%.2f-%.2f), P=%.3g",
    s$coef["Z","exp(coef)"],s$conf.int["Z","lower .95"],s$conf.int["Z","upper .95"],s$coef["Z","Pr(>|z|)"]) }
say("robust, 5년 사망, per +1 성별내 s.d.  미보정 : HR %s", e("zS",FALSE,"t_fu","ev",5))
say("robust, 5년 사망, per +1 성별내 s.d.  보정   : HR %s", e("zS",TRUE ,"t_fu","ev",5))
say("robust, 5년 복합악화, per +1 s.d.     미보정 : HR %s", e("zS",FALSE,"t_comp","e_comp",5))
say("robust, 5년 복합악화, per +1 s.d.     보정   : HR %s", e("zS",TRUE ,"t_comp","e_comp",5))
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
