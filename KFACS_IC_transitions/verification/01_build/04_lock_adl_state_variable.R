# =============================================================================
# 14 — ADL 상태 변수 후보를 전부 원고 사건수에 대조
#      + 회복 패널의 사망 코딩 버그 수정
#
# 지난 실행 실패 원인: adl_3cat 에 결측 2,625건. 원고는 결측이 처리된 변수를 씁니다.
# 후보: iv$from/to (숫자), d$state, d$state_obs, d$adl_state, d$state_lab 등
# 목표: Normal->Mild 969 / Normal->Severe 141 / Normal->Death 361
#       Mild->Normal 651 / Mild->Severe 74 / Mild->Death 97
#       Severe->Normal 21 / Severe->Mild 30 / Severe->Death 42
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({
  for (p in c("survival","sandwich","lmtest","dplyr")) library(p, character.only=TRUE) })
OUT <- "adl_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "find_state_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

TGT <- c("Normal->Mild"=969,"Normal->Severe"=141,"Normal->Death"=361,
         "Mild->Normal"=651,"Mild->Severe"=74,"Mild->Death"=97,
         "Severe->Normal"=21,"Severe->Mild"=30,"Severe->Death"=42)
score <- function(fromv, tov, lab) {
  ok <- !is.na(fromv) & !is.na(tov)
  tb <- table(paste0(fromv[ok], "->", tov[ok]))
  hit <- sapply(names(TGT), function(k) { g <- if (k %in% names(tb)) as.integer(tb[[k]]) else 0L; g })
  n_ok <- sum(hit == TGT)
  say("  %-34s 일치 %d/9   %s", lab, n_ok,
      paste(sprintf("%s=%d", sub("->","→",names(TGT)), hit), collapse=" "))
  n_ok
}
norm <- function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% c("0","Normal","normal","None","none"), "Normal",
       ifelse(x %in% c("1","Mild","mild"), "Mild",
       ifelse(x %in% c("2","Severe","severe"), "Severe",
       ifelse(x %in% c("3","Death","death","Dead","dead"), "Death", x))))
  x
}

## ===========================================================================
rule("1. iv 의 from/to 숫자 열이 ADL 인가")
say("iv$from 값: %s", paste(names(table(iv$from, useNA='ifany')), collapse=", "))
say("iv$to   값: %s", paste(names(table(iv$to,   useNA='ifany')), collapse=", "))
say("iv$type 값: %s", paste(names(table(iv$type, useNA='ifany')), collapse=", "))
say("from x from_lab 교차표:"); print(table(iv$from, iv$from_lab, useNA="ifany"))
best <- list(n=-1)
n <- score(norm(iv$from), norm(iv$to), "iv$from / iv$to (그대로)")
if (n > best$n) best <- list(n=n, src="iv_raw")

## ===========================================================================
rule("2. d 의 상태 변수 후보를 구간에 얹어 대조")
V <- d %>% filter(!is.na(wave)) %>% arrange(id, wave) %>%
  group_by(id) %>% mutate(k_visit = row_number()) %>% ungroup()
cands <- intersect(c("state","state_obs","adl_state","state_lab","state_label",
                     "state_obs_label","adl_3cat","adl_2cat"), names(V))
say("후보 열: %s", paste(cands, collapse=", "))
for (v in cands) {
  say("\n%s 분포:", v); print(table(V[[v]], useNA="ifany"))
}
for (v in cands) {
  A <- iv %>%
    left_join(V %>% select(id, k_visit, .f = all_of(v)), by=c("id"="id","k"="k_visit")) %>%
    mutate(k2 = k + 1) %>%
    left_join(V %>% select(id, k_visit, .t = all_of(v)), by=c("id"="id","k2"="k_visit")) %>%
    mutate(.t = ifelse(as.character(to_lab)=="Death", "Death", as.character(.t)))
  n <- score(norm(A$.f), norm(A$.t), paste0("d$", v, " (구간에 얹음)"))
  if (n > best$n) best <- list(n=n, src=v, frame=A)
}

## in_risk 플래그를 걸어본 경우
if ("in_risk" %in% names(V)) {
  say("\nin_risk 분포:"); print(table(V$in_risk, useNA="ifany"))
}

rule("3. 판정")
say("최선 후보: %s  (일치 %d/9)", best$src, best$n)
if (best$n < 9) {
  say("!! 9/9 를 재현하는 변수가 없습니다.")
  say("   위 1절의 'from x from_lab 교차표' 와 각 후보의 분포를 보시고,")
  say("   원고 Table 2 를 만든 스크립트에서 상태 정의 부분을 알려주십시오.")
} else {
  say(">>> 재현 성공. 이 변수로 추정을 진행합니다.")
  A <- if (best$src == "iv_raw")
        iv %>% mutate(.f = norm(from), .t = norm(to)) else best$frame %>% mutate(.f = norm(.f), .t = norm(.t))
  A <- A %>% group_by(sex_f) %>% mutate(zS = as.numeric(scale(gLIC_z))) %>% ungroup()
  CV <- "age_c + comorbid_count + edu_bl_f + income_bl_f + area_bl_f"
  runtr <- function(f1, t1) {
    dat <- A %>% filter(.f == f1, !is.na(.t), !is.na(zS), dur > 0)
    dat$y <- as.integer(dat$.t == t1)
    if (sum(dat$y) < 5) return(data.frame(from=f1, to=t1, events=sum(dat$y), py=round(sum(dat$dur)),
                                          IRR=NA, lo=NA, hi=NA, p=NA))
    m <- glm(as.formula(paste("y ~ zS +", CV, "+ offset(log(dur))")), family=poisson, data=dat)
    ct <- coeftest(m, vcov=vcovCL(m, cluster=dat$id))
    b <- ct["zS","Estimate"]; se <- ct["zS","Std. Error"]
    data.frame(from=f1, to=t1, events=sum(dat$y), py=round(sum(dat$dur)),
               IRR=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=ct["zS","Pr(>|z|)"])
  }
  TAB <- bind_rows(lapply(strsplit(names(TGT),"->"), function(x) runtr(x[1], x[2])))
  print(as.data.frame(TAB), digits=3, row.names=FALSE)
  write.csv(TAB, file.path(OUT,"table2_ADL_sexstd_VERIFIED.csv"), row.names=FALSE)
}

## ===========================================================================
rule("4. 회복 패널 — 사망 코딩 수정본")
## 지난 실행 버그: t_end 를 t_lastvisit 으로 잘라서 사망(방문 이후 발생)이 0 이 되었습니다.
STV <- if (best$n == 9 && best$src != "iv_raw") best$src else "adl_3cat"
say("회복 정의에 사용할 상태 변수: %s", STV)
V$.st <- norm(V[[STV]])
S <- d %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event>0, na.rm=TRUE)),
            td = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)), .groups="drop") %>%
  mutate(td = ifelse(is.finite(td), td, NA_real_))
W1 <- V %>% group_by(id) %>% slice_min(k_visit, n=1, with_ties=FALSE) %>% ungroup() %>%
  select(id, st_w1 = .st, fr_w1 = frailty_3cat, gLIC_z) %>%
  left_join(iv %>% select(id, sex_f) %>% distinct(), by="id") %>%
  group_by(sex_f) %>%
  mutate(tS = cut(gLIC_z, c(-Inf, quantile(gLIC_z, c(1/3,2/3), na.rm=TRUE), Inf),
                  labels=c("T1","T2","T3"))) %>% ungroup()
REC <- V %>% arrange(id, k_visit) %>% group_by(id) %>%
  summarise(t_rec = { w <- which(.st == "Normal" & k_visit > 1); if (length(w)) time[w[1]] else NA_real_ },
            t_last = max(time), .groups="drop")
R <- W1 %>% left_join(REC, by="id") %>% left_join(S, by="id") %>%
  filter(st_w1 %in% c("Mild","Severe")) %>%
  mutate(t_end = pmin(ifelse(is.na(t_rec), Inf, t_rec),
                      ifelse(ev==1 & !is.na(td), td, Inf), na.rm=TRUE),
         t_end = ifelse(is.finite(t_end), t_end, t_last),
         status = ifelse(!is.na(t_rec) & t_rec <= t_end, 1L,
                  ifelse(ev==1 & !is.na(td) & td <= t_end + 1e-9, 2L, 0L))) %>%
  group_by(sex_f) %>% mutate(zS = as.numeric(scale(gLIC_z))) %>% ungroup()
say("wave-1 ADL 장애 보유자 = %d | robust = %d  [원고 103]",
    nrow(R), sum(as.character(R$fr_w1)=="0"))
print(table(R$fr_w1, R$status, dnn=c("frailty","0=cens 1=recov 2=death")))
say(">>> 사망 열이 0 이 아니어야 정상입니다.")

for (fr in c("0","1","2")) {
  Z <- R %>% filter(as.character(fr_w1)==fr); if (nrow(Z) < 15) next
  say("\n---- frailty_3cat = %s (n=%d, 회복 %d, 사망 %d) ----",
      fr, nrow(Z), sum(Z$status==1), sum(Z$status==2))
  Z$St <- factor(Z$status, 0:2, c("censor","recovery","death"))
  aj <- try(survfit(Surv(t_end, St) ~ tS, data=Z), silent=TRUE)
  if (!inherits(aj,"try-error")) {
    s5 <- summary(aj, times=5, extend=TRUE); i <- grep("recovery", colnames(s5$pstate))
    print(data.frame(group=s5$strata, CIF5_recovery=round(100*s5$pstate[,i],1)), row.names=FALSE)
  }
  cs <- try(summary(coxph(Surv(t_end, status==1) ~ zS + age_c + comorbid_count, data=Z)), silent=TRUE)
  if (!inherits(cs,"try-error")) say("  cause-specific HR = %.2f (%.2f-%.2f), P = %.3g",
      cs$coef["zS","exp(coef)"], cs$conf.int["zS","lower .95"],
      cs$conf.int["zS","upper .95"], cs$coef["zS","Pr(>|z|)"])
  fg <- try({ fd <- finegray(Surv(t_end, St) ~ ., data=Z, etype="recovery")
              summary(coxph(Surv(fgstart,fgstop,fgstatus) ~ zS, weight=fgwt, data=fd)) }, silent=TRUE)
  if (!inherits(fg,"try-error")) say("  Fine-Gray sHR  = %.2f (%.2f-%.2f), P = %.3g",
      fg$coef["zS","exp(coef)"], fg$conf.int["zS","lower .95"],
      fg$conf.int["zS","upper .95"], fg$coef["zS","Pr(>|z|)"])
  tb <- Z %>% group_by(tS) %>% summarise(n=n(), recov=sum(status==1), death=sum(status==2), .groups="drop")
  print(as.data.frame(tb), row.names=FALSE)
  write.csv(tb, file.path(OUT, paste0("recovery_fixed_fr", fr, ".csv")), row.names=FALSE)
}
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
