# =============================================================================
# 07 — dat_w1 은 원자료와 연결되지 않습니다. 원자료(d, iv)에서 직접 다시 만듭니다.
#
# 확인된 사실
#   dat_w1$id = 1, 2, 3, ... 3011   (단순 행번호)
#   d$id      = "kf160001", ...     (실제 KFACS 식별자)
#   -> 두 파일은 어떤 방식으로도 조인되지 않습니다. 지문 매칭도 0건.
#
#   그리고 원자료 d 는 원고와 정확히 일치합니다:
#     death_event (id 단위) = 500        <- 원고 초록의 500
#     death_wave <= 5       = 500
#     frailty_3cat == 0     = 1,349      <- 원고 Table 3 의 robust
#   반면 dat_w1 은 556 / 1,082 로 둘 다 어긋납니다.
#
# 결론: d(+iv) 가 정본이고, dat_w1 은 출처 불명의 파생 프레임입니다.
#       HR 0.63 을 포함한 재분석 결과 전부 폐기 대상입니다.
#
# 이 스크립트는 d 와 iv 만으로 wave-1 분석 프레임을 다시 만들고,
# 원고의 4개 기준값(3,011 / 1,349 / 500 / 103)으로 검증한 뒤에만 분석을 진행합니다.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ library(survival); library(dplyr) })
OUT <- "rebuild_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "rebuild_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")
chk  <- function(lab, got, want)
  say("  %-40s %8s  (원고 %s) %s", lab, format(got, big.mark=","), format(want, big.mark=","),
      if (isTRUE(got == want)) "  OK" else "  <<< 불일치")

## ===========================================================================
rule("0. 원자료 열 탐색 — 무엇이 있는지 먼저 확인")
show_cols <- function(df, nm, patt) {
  h <- names(df)[grepl(patt, tolower(names(df)))]
  say("%-4s | %-14s : %s", nm, sub("\\|.*","",patt), if (length(h)) paste(head(h,10), collapse=", ") else "(없음)")
  invisible(h)
}
for (p in c("wave|visit","^time$|fu|dur|t0|t1","death|dead|mort","frail|fried|chs",
            "adl|disab","glic|^ic|intrinsic|capacity","age","sex|female","edu","income",
            "area|urban|resid","comorb|cci|charlson","date")) {
  show_cols(d, "d", p)
}
say("")
say("iv 열: %s", paste(names(iv), collapse=", "))

## ===========================================================================
rule("1. wave-1 프레임 재구성")

D <- d
W1 <- D %>% group_by(id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>% ungroup()
say("wave-1 행 = %d, 고유 id = %d", nrow(W1), length(unique(D$id)))

## 사망 여부와 사망시점
SURV <- D %>% group_by(id) %>%
  summarise(
    ev    = as.integer(any(death_event > 0, na.rm = TRUE)),
    t_ev  = suppressWarnings(min(time[death_event > 0], na.rm = TRUE)),
    t_last= suppressWarnings(max(time, na.rm = TRUE)),
    .groups = "drop") %>%
  mutate(t_ev = ifelse(is.finite(t_ev), t_ev, NA_real_),
         t_last = ifelse(is.finite(t_last), t_last, NA_real_),
         t_death = ifelse(ev == 1, t_ev, t_last))
say("사망 = %d, t_death 결측 = %d", sum(SURV$ev), sum(is.na(SURV$t_death)))
say("t_death 요약:"); print(summary(SURV$t_death))
say("사망자의 t_death 요약:"); print(summary(SURV$t_death[SURV$ev==1]))
say(">> time 열이 '년' 단위가 아니면(예: 일/월) 위 수치가 이상하게 보입니다. 확인하십시오.")

## 노출과 공변량 — iv 의 첫 구간에서
IV1 <- iv %>% group_by(id) %>% slice_min(k, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(id, gLIC_z, any_of(c("gLIC_tert","age_c","bage","sex_f","comorbid_count","comorbid_bl",
                              "edu_bl_f","income_bl_f","area_bl_f","alone_bl_f","ses_index")))
say("\niv 첫 구간 행 = %d, gLIC_z 결측 = %d", nrow(IV1), sum(is.na(IV1$gLIC_z)))

B <- W1 %>% select(id, wave, frailty_3cat,
                   any_of(c("frailty_lab","chs_total","age_unified","sex_unified"))) %>%
  left_join(SURV %>% select(id, ev, t_death), by = "id") %>%
  left_join(IV1, by = "id")
say("결합 후 행 = %d, frailty_3cat 결측 = %d, gLIC_z 결측 = %d",
    nrow(B), sum(is.na(B$frailty_3cat)), sum(is.na(B$gLIC_z)))

## ADL 상태
adl_col <- names(W1)[grepl("adl", tolower(names(W1)))]
say("\nADL 후보 열: %s", paste(adl_col, collapse=", "))
for (a in head(adl_col, 5)) { say("  %s:", a); print(table(W1[[a]], useNA="ifany")) }

## ===========================================================================
rule("2. 원고 기준값 검증 — 여기서 다 맞아야 분석을 진행합니다")
chk("전체 n",                 nrow(B), 3011)
chk("robust (frailty_3cat==0)", sum(as.character(B$frailty_3cat)=="0", na.rm=TRUE), 1349)
chk("pre-frail",              sum(as.character(B$frailty_3cat)=="1", na.rm=TRUE), 1416)
chk("frail",                  sum(as.character(B$frailty_3cat)=="2", na.rm=TRUE),  246)
chk("총 사망",                sum(B$ev, na.rm=TRUE), 500)
say("\n추적 중앙값 = %.2f, 최대 = %.2f   (원고: 중앙값 8.0, 최대 8.2)",
    median(B$t_death, na.rm=TRUE), max(B$t_death, na.rm=TRUE))

OKGO <- (nrow(B)==3011) && (sum(as.character(B$frailty_3cat)=="0",na.rm=TRUE)==1349) &&
        (sum(B$ev,na.rm=TRUE)==500)

## ===========================================================================
rule("3. robust 층 재분석  (검증 통과 시에만)")
if (!OKGO) {
  say("!! 기준값이 맞지 않아 분석을 중단합니다.")
  say("   위 2절에서 '<<< 불일치' 가 붙은 항목을 알려주십시오.")
  say("   특히 t_death 단위(년/일/월)와 time 열의 정의를 확인해야 합니다.")
} else {
  RB <- B %>% filter(as.character(frailty_3cat)=="0", !is.na(gLIC_z), !is.na(t_death))
  say("robust 분석 n = %d, 사망 = %d", nrow(RB), sum(RB$ev))

  SDc <- sd(B$gLIC_z, na.rm=TRUE); MUc <- mean(B$gLIC_z, na.rm=TRUE)
  RB$.z <- (RB$gLIC_z - MUc)/SDc
  cuts <- quantile(B$gLIC_z, c(1/3,2/3), na.rm=TRUE)
  RB$tert <- cut(RB$gLIC_z, c(-Inf,cuts,Inf), labels=c("T1","T2","T3"))
  say("코호트 tertile 절단점 = %.4f, %.4f", cuts[1], cuts[2])

  cvs <- intersect(c("age_c","bage","sex_f","comorbid_count","comorbid_bl",
                     "edu_bl_f","income_bl_f","area_bl_f"), names(RB))
  cvs <- cvs[sapply(cvs, function(v) length(unique(na.omit(RB[[v]])))>1)]
  say("보정변수: %s", paste(cvs, collapse=", "))

  km <- function(dd,h){ if(!nrow(dd)) return(NA_real_)
    s<-survfit(Surv(pmin(t_death,h), ifelse(t_death>h,0L,ev))~1, data=dd)
    100*(1-summary(s,times=h,extend=TRUE)$surv[1]) }
  tb <- RB %>% group_by(tert) %>%
    summarise(n=n(), d5=sum(ev[t_death<=5]), d8=sum(ev[t_death<=8]), .groups="drop")
  tb$km5 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],5))
  tb$km8 <- sapply(tb$tert, function(g) km(RB[RB$tert==g,],8))
  print(as.data.frame(tb), digits=3, row.names=FALSE)
  say("원고 Table 3 : T1 n=106 d=2 (1.9%%) | T2 n=470 d=23 (4.9%%) | T3 n=773 d=26 (3.4%%)")
  write.csv(tb, file.path(OUT,"robust_tertiles_SOURCE.csv"), row.names=FALSE)

  for (h in c(5,8)) {
    RB$.t <- pmin(RB$t_death,h); RB$.e <- ifelse(RB$t_death>h,0L,RB$ev)
    f <- as.formula(paste("Surv(.t,.e) ~ .z", if(length(cvs)) paste("+",paste(cvs,collapse="+")) else ""))
    m <- coxph(f, data=RB); s <- summary(m)
    say("%d년 per +1 s.d. : HR %.3f (%.3f-%.3f), 사망 %d / %d, P = %.3g", h,
        s$coef[".z","exp(coef)"], s$conf.int[".z","lower .95"], s$conf.int[".z","upper .95"],
        m$nevent, m$n, s$coef[".z","Pr(>|z|)"])
  }
  saveRDS(B,  file.path(OUT, "wave1_frame_SOURCE.rds"))
  saveRDS(RB, file.path(OUT, "robust_frame_SOURCE.rds"))
  say("\n검증 통과 프레임 저장: rebuild_out/wave1_frame_SOURCE.rds, robust_frame_SOURCE.rds")
}

rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
