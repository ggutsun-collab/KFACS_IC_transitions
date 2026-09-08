# =============================================================================
# 13 — Table 2 의 ADL 전이 9개 + Table 3 의 회복 패널
#
# 새로 만들 것 없습니다. 원고 Fig 1 각주대로 ADL 패널과 frailty 패널은
# 같은 11,571 구간 / 22,730 person-years 를 씁니다. iv 의 구간에 ADL 라벨만 얹습니다.
#
# 검증: 원고 Table 2 의 ADL 사건수와 대조합니다.
#   Normal->Mild 969 | Normal->Severe 141 | Normal->Death 361
#   Mild->Normal 651 | Mild->Severe   74 | Mild->Death    97
#   Severe->Normal 21 | Severe->Mild   30 | Severe->Death  42
# 이 9개가 맞으면 프레임이 정확한 것이고, 그 위에서 성별 내 z 로 다시 추정합니다.
#
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({
  for (p in c("survival","sandwich","lmtest","dplyr")) library(p, character.only=TRUE) })
OUT <- "adl_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "adl_recovery_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")

D <- d
ADLLAB <- c("0"="Normal","1"="Mild","2"="Severe")

## ===========================================================================
rule("1. 방문 단위 ADL 상태 사다리 만들기")
V <- D %>% filter(!is.na(wave)) %>% arrange(id, wave) %>%
  group_by(id) %>% mutate(k_visit = row_number()) %>% ungroup() %>%
  select(id, wave, k_visit, time, adl_3cat, frailty_3cat, gLIC_z)
V$adl_lab <- ADLLAB[as.character(V$adl_3cat)]
say("방문 행 = %d, 고유 id = %d", nrow(V), n_distinct(V$id))
say("방문별 ADL 상태:"); print(table(V$adl_lab, useNA="ifany"))

## iv 의 구간에 시작/종료 방문을 붙입니다 (k = 구간 번호 = 시작 방문 번호)
IVa <- iv %>%
  left_join(V %>% select(id, k_visit, adl_from = adl_lab, gLIC_from = gLIC_z),
            by = c("id"="id", "k"="k_visit")) %>%
  mutate(k_end = k + 1) %>%
  left_join(V %>% select(id, k_visit, adl_end = adl_lab),
            by = c("id"="id", "k_end"="k_visit")) %>%
  mutate(adl_to = ifelse(as.character(to_lab) == "Death", "Death", adl_end))

say("\n구간 = %d  [원고 11,571] | person-years = %.0f  [원고 22,730]",
    nrow(IVa), sum(IVa$dur, na.rm=TRUE))
say("adl_from 결측 = %d | adl_to 결측 = %d", sum(is.na(IVa$adl_from)), sum(is.na(IVa$adl_to)))

## ===========================================================================
rule("2. 검증 — 원고 Table 2 의 ADL 사건수 재현")
TGT <- c("Normal->Mild"=969,"Normal->Severe"=141,"Normal->Death"=361,
         "Mild->Normal"=651,"Mild->Severe"=74,"Mild->Death"=97,
         "Severe->Normal"=21,"Severe->Mild"=30,"Severe->Death"=42)
obs <- IVa %>% filter(!is.na(adl_from), !is.na(adl_to)) %>%
  count(adl_from, adl_to) %>% mutate(key = paste0(adl_from,"->",adl_to))
for (k in names(TGT)) {
  got <- obs$n[match(k, obs$key)]; got <- ifelse(is.na(got), 0, got)
  say("  %-16s %5d   (원고 %4d) %s", k, got, TGT[[k]],
      if (isTRUE(got == TGT[[k]])) " OK" else " <<< 불일치")
}
say("\n>>> 9개가 모두 OK 면 프레임이 정확합니다. 불일치가 있으면 아래 추정은 쓰지 마십시오.")
say("전체 교차표:"); print(table(IVa$adl_from, IVa$adl_to, useNA="ifany"))

## ===========================================================================
rule("3. 성별 내 표준화 z 로 ADL 전이 9개 재추정")
IVa <- IVa %>% group_by(sex_f) %>% mutate(zS = as.numeric(scale(gLIC_z))) %>% ungroup()
CV <- "age_c + comorbid_count + edu_bl_f + income_bl_f + area_bl_f"

runtr <- function(fromst, tost) {
  dat <- IVa %>% filter(adl_from == fromst, !is.na(adl_to), !is.na(zS), dur > 0)
  if (!nrow(dat)) return(NULL)
  dat$y <- as.integer(dat$adl_to == tost)
  if (sum(dat$y) < 5) return(data.frame(from=fromst, to=tost, events=sum(dat$y),
        py=sum(dat$dur), IRR=NA, lo=NA, hi=NA, p=NA))
  m <- try(glm(as.formula(paste("y ~ zS +", CV, "+ offset(log(dur))")),
               family=poisson, data=dat), silent=TRUE)
  if (inherits(m,"try-error")) return(NULL)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = dat$id))
  b <- ct["zS","Estimate"]; se <- ct["zS","Std. Error"]
  data.frame(from=fromst, to=tost, events=sum(dat$y), py=round(sum(dat$dur)),
             IRR=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=ct["zS","Pr(>|z|)"])
}
TAB <- bind_rows(lapply(strsplit(names(TGT),"->"), function(x) runtr(x[1], x[2])))
print(as.data.frame(TAB), digits=3, row.names=FALSE)
write.csv(TAB, file.path(OUT,"table2_ADL_sexstd.csv"), row.names=FALSE)
say(">>> 이 값을 Manuscript 의 Table 2 상단 9행에 넣으십시오.")

## tertile 대비 (성별 내 tertile)
say("\n성별 내 tertile 대비 (T2 vs T1, T3 vs T1):")
W1 <- V %>% group_by(id) %>% slice_min(k_visit, n=1, with_ties=FALSE) %>% ungroup()
W1 <- W1 %>% left_join(iv %>% select(id, sex_f) %>% distinct(), by="id") %>%
  group_by(sex_f) %>%
  mutate(tS = cut(gLIC_z, c(-Inf, quantile(gLIC_z, c(1/3,2/3), na.rm=TRUE), Inf),
                  labels=c("T1","T2","T3"))) %>% ungroup()
IVt <- IVa %>% left_join(W1 %>% select(id, tS), by="id")
for (x in strsplit(names(TGT),"->")) {
  dat <- IVt %>% filter(adl_from==x[1], !is.na(adl_to), !is.na(tS), dur>0)
  dat$y <- as.integer(dat$adl_to==x[2])
  if (sum(dat$y) < 15) { say("  %-16s 사건 %d — 추정 불가(NE)", paste0(x[1],"->",x[2]), sum(dat$y)); next }
  m <- try(glm(as.formula(paste("y ~ tS +", CV, "+ offset(log(dur))")), family=poisson, data=dat), silent=TRUE)
  if (inherits(m,"try-error")) next
  ct <- coeftest(m, vcov=vcovCL(m, cluster=dat$id))
  f <- function(nm) if (nm %in% rownames(ct))
    sprintf("%.2f (%.2f-%.2f)", exp(ct[nm,1]), exp(ct[nm,1]-1.96*ct[nm,2]), exp(ct[nm,1]+1.96*ct[nm,2])) else "NE"
  say("  %-16s T2 vs T1 %s | T3 vs T1 %s", paste0(x[1],"->",x[2]), f("tST2"), f("tST3"))
}

## ===========================================================================
rule("4. 회복 패널 — wave-1 ADL 장애 보유자")
REC <- V %>% arrange(id, k_visit) %>% group_by(id) %>%
  summarise(adl_w1 = first(adl_lab), fr_w1 = first(frailty_3cat),
            t_rec = { w <- which(adl_lab=="Normal" & k_visit > 1)
                      if (length(w)) time[w[1]] else NA_real_ },
            t_lastvisit = max(time), .groups="drop")
S <- D %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event>0, na.rm=TRUE)),
            td = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)), .groups="drop") %>%
  mutate(td = ifelse(is.finite(td), td, NA))
R <- REC %>% left_join(S, by="id") %>% left_join(W1 %>% select(id, gLIC_z, tS, sex_f), by="id") %>%
  filter(adl_w1 != "Normal") %>%
  mutate(t_end = pmin(ifelse(is.na(t_rec), Inf, t_rec),
                      ifelse(ev==1 & !is.na(td), td, Inf), t_lastvisit, na.rm=TRUE),
         status = case_when(!is.na(t_rec) & t_rec <= t_end ~ 1L,     # 1 = 회복
                            ev==1 & !is.na(td) & td <= t_end ~ 2L,   # 2 = 사망(경쟁)
                            TRUE ~ 0L))
R <- R %>% group_by(sex_f) %>% mutate(zS = as.numeric(scale(gLIC_z))) %>% ungroup()
say("wave-1 ADL 장애 보유자 = %d", nrow(R))
say("그중 robust = %d   [원고 103]", sum(as.character(R$fr_w1)=="0"))
print(table(R$fr_w1, R$status, dnn=c("frailty","0=censor 1=recov 2=death")))

for (fr in c("0","1","2")) {
  Z <- R %>% filter(as.character(fr_w1)==fr)
  if (nrow(Z) < 15) { say("\nfrailty=%s : n=%d — 생략", fr, nrow(Z)); next }
  say("\n---- frailty_3cat = %s (n=%d) ----", fr, nrow(Z))
  Z$St <- factor(Z$status, 0:2, c("censor","recovery","death"))
  aj <- survfit(Surv(t_end, St) ~ tS, data=Z)
  s5 <- try(summary(aj, times=5, extend=TRUE), silent=TRUE)
  if (!inherits(s5,"try-error")) {
    idx <- grep("recovery", colnames(s5$pstate))
    print(data.frame(group=s5$strata, n=s5$n.risk[,1],
                     CIF5 = round(100*s5$pstate[, idx], 1)))
  }
  cs <- try(coxph(Surv(t_end, status==1) ~ zS + age_c + comorbid_count, data=Z), silent=TRUE)
  if (!inherits(cs,"try-error")) { s <- summary(cs)
    say("  cause-specific HR per +1 성별내 s.d. = %.2f (%.2f-%.2f), P = %.3g",
        s$coef["zS","exp(coef)"], s$conf.int["zS","lower .95"],
        s$conf.int["zS","upper .95"], s$coef["zS","Pr(>|z|)"]) }
  fg <- try({ fgd <- finegray(Surv(t_end, factor(status,0:2,c("censor","recovery","death"))) ~ .,
                              data=Z, etype="recovery")
              coxph(Surv(fgstart, fgstop, fgstatus) ~ zS, weight=fgwt, data=fgd) }, silent=TRUE)
  if (!inherits(fg,"try-error")) { s <- summary(fg)
    say("  Fine-Gray sHR = %.2f (%.2f-%.2f), P = %.3g", s$coef["zS","exp(coef)"],
        s$conf.int["zS","lower .95"], s$conf.int["zS","upper .95"], s$coef["zS","Pr(>|z|)"]) }
  tb <- Z %>% group_by(tS) %>% summarise(n=n(), recov=sum(status==1), death=sum(status==2), .groups="drop")
  print(as.data.frame(tb), row.names=FALSE)
  write.csv(tb, file.path(OUT, paste0("recovery_fr", fr, ".csv")), row.names=FALSE)
}
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
