# =============================================================================
# 15 — Table 2 의 tertile 대비 + Table 3 의 회복 패널을 두 척도로 모두 산출
#
# 본문 표를 원래 구조(전 열 / 전 패널)로 되돌리기 위해 비어 있던 칸을 채웁니다.
# 코호트 척도와 성별 내 척도를 나란히 내므로, 어느 쪽을 본문에 두든
# 나머지 하나는 그대로 Extended Data Table 1 이 됩니다.
#
# 확정 규칙 사용: 사망시점 = wave 결측 행의 time, 장애상태 = state_lab
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({
  for (p in c("survival","sandwich","lmtest","dplyr")) library(p, character.only=TRUE) })
OUT <- "fill_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "fill_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")
norm <- function(x) { x <- as.character(x)
  ifelse(x %in% c("0","Normal"), "Normal", ifelse(x %in% c("1","Mild"), "Mild",
  ifelse(x %in% c("2","Severe"), "Severe", ifelse(x %in% c("3","Death"), "Death", x)))) }

## ===========================================================================
rule("0. 프레임")
V <- d %>% filter(!is.na(wave)) %>% arrange(id, wave) %>%
  group_by(id) %>% mutate(kv = row_number()) %>% ungroup() %>%
  mutate(.st = norm(state_lab))
S <- d %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event>0, na.rm=TRUE)),
            td = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)),
            tl = suppressWarnings(max(time[!is.na(wave)], na.rm=TRUE)), .groups="drop") %>%
  mutate(td = ifelse(is.finite(td), td, NA), tl = ifelse(is.finite(tl), tl, NA),
         t_fu = ifelse(ev==1, td, tl))

W1 <- V %>% group_by(id) %>% slice_min(kv, n=1, with_ties=FALSE) %>% ungroup() %>%
  select(id, st_w1 = .st, fr_w1 = frailty_3cat) %>%
  left_join(iv %>% group_by(id) %>% slice_min(k, n=1, with_ties=FALSE) %>% ungroup() %>%
              select(id, gLIC_z, sex_f, age_c, comorbid_count, edu_bl_f, income_bl_f, area_bl_f),
            by="id") %>%
  left_join(S %>% select(id, ev, t_fu), by="id") %>%
  filter(!is.na(gLIC_z))

# 두 척도의 tertile
cut_coh <- quantile(W1$gLIC_z, c(1/3,2/3), na.rm=TRUE)
W1$T_coh <- cut(W1$gLIC_z, c(-Inf, cut_coh, Inf), labels=c("T1","T2","T3"))
W1 <- W1 %>% group_by(sex_f) %>%
  mutate(z_sex = as.numeric(scale(gLIC_z)),
         T_sex = cut(gLIC_z, c(-Inf, quantile(gLIC_z, c(1/3,2/3), na.rm=TRUE), Inf),
                     labels=c("T1","T2","T3"))) %>% ungroup()
W1$z_coh <- as.numeric(scale(W1$gLIC_z))
say("n = %d | robust %d | 사망 %d", nrow(W1), sum(W1$fr_w1==0), sum(W1$ev))
say("tertile n  코호트 척도: %s", paste(table(W1$T_coh), collapse=" / "))
say("tertile n  성별 내 척도: %s", paste(table(W1$T_sex), collapse=" / "))

## ===========================================================================
rule("1. Table 2 — 18개 전이의 tertile 대비 (두 척도)")
IVa <- iv %>%
  left_join(V %>% select(id, kv, .f = .st), by=c("id"="id","k"="kv")) %>%
  mutate(k2 = k+1) %>%
  left_join(V %>% select(id, kv, .e2 = .st), by=c("id"="id","k2"="kv")) %>%
  mutate(.t = ifelse(as.character(to_lab)=="Death", "Death", .e2)) %>%
  left_join(W1 %>% select(id, T_coh, T_sex), by="id")

TRS <- rbind(
  data.frame(sys="Disability", from=c("Normal","Normal","Normal","Mild","Mild","Mild",
                                      "Severe","Severe","Severe"),
             to=c("Mild","Severe","Death","Normal","Severe","Death","Normal","Mild","Death")),
  data.frame(sys="Frailty", from=c("Robust","Robust","Robust","Pre-frail","Pre-frail","Pre-frail",
                                   "Frail","Frail","Frail"),
             to=c("Pre-frail","Frail","Death","Robust","Frail","Death","Robust","Pre-frail","Death")))
CV <- "age_c + comorbid_count + edu_bl_f + income_bl_f + area_bl_f"

tert_row <- function(i, tv) {
  tr <- TRS[i,]
  dat <- if (tr$sys=="Disability") IVa %>% filter(.f==tr$from) else
                                   IVa %>% filter(as.character(from_lab)==tr$from)
  dat$y <- if (tr$sys=="Disability") as.integer(dat$.t==tr$to) else
                                     as.integer(as.character(dat$to_lab)==tr$to)
  dat$G <- dat[[tv]]
  dat <- dat %>% filter(!is.na(G), dur > 0)
  ev_by <- tapply(dat$y, dat$G, sum)
  out <- data.frame(sys=tr$sys, from=tr$from, to=tr$to, scale=tv,
                    ev_T1=ev_by[["T1"]], ev_T2=ev_by[["T2"]], ev_T3=ev_by[["T3"]],
                    T2vT1="NE", T3vT1="NE", stringsAsFactors=FALSE)
  if (min(ev_by) >= 5) {
    m <- try(glm(as.formula(paste("y ~ G +", CV, "+ offset(log(dur))")),
                 family=poisson, data=dat), silent=TRUE)
    if (!inherits(m,"try-error")) {
      ct <- coeftest(m, vcov=vcovCL(m, cluster=dat$id))
      f <- function(nm) if (nm %in% rownames(ct))
        sprintf("%.2f (%.2f-%.2f)", exp(ct[nm,1]), exp(ct[nm,1]-1.96*ct[nm,2]),
                exp(ct[nm,1]+1.96*ct[nm,2])) else "NE"
      out$T2vT1 <- f("GT2"); out$T3vT1 <- f("GT3")
    }
  }
  out
}
TAB2 <- bind_rows(lapply(seq_len(nrow(TRS)), function(i)
          rbind(tert_row(i,"T_sex"), tert_row(i,"T_coh"))))
print(as.data.frame(TAB2), row.names=FALSE)
write.csv(TAB2, file.path(OUT,"table2_tertile_contrasts_BOTH.csv"), row.names=FALSE)
say(">>> scale = T_sex 는 본문(또는 ED)용, T_coh 는 대조용입니다.")
say(">>> 사건 5건 미만 tertile 이 있으면 NE 로 표기했습니다 (원고 규칙과 동일).")

## ===========================================================================
rule("2. Table 3 — 회복 패널 (두 척도)")
REC <- V %>% arrange(id, kv) %>% group_by(id) %>%
  summarise(t_rec = { w <- which(.st=="Normal" & kv > 1); if (length(w)) time[w[1]] else NA_real_ },
            t_last = max(time), .groups="drop")
R <- W1 %>% left_join(REC, by="id") %>% filter(st_w1 %in% c("Mild","Severe")) %>%
  mutate(t_end = pmin(ifelse(is.na(t_rec), Inf, t_rec),
                      ifelse(ev==1 & !is.na(t_fu), t_fu, Inf), na.rm=TRUE),
         t_end = ifelse(is.finite(t_end), t_end, t_last),
         status = ifelse(!is.na(t_rec) & t_rec <= t_end, 1L,
                  ifelse(ev==1 & !is.na(t_fu) & t_fu <= t_end + 1e-9, 2L, 0L)),
         St = factor(status, 0:2, c("censor","recovery","death")),
         fr = factor(as.character(fr_w1), c("0","1","2"), c("Robust","Pre-frail","Frail")))
say("wave-1 ADL 장애 보유자 = %d (robust %d)  [원고 337 / 103]",
    nrow(R), sum(R$fr=="Robust"))

for (tv in c("T_sex","T_coh")) {
  say("\n---- tertile 정의: %s ----", tv)
  R$G <- R[[tv]]
  R$cell <- interaction(R$fr, R$G, sep="/")
  R$cell <- relevel(factor(R$cell), ref="Robust/T3")
  tb <- R %>% group_by(fr, G) %>%
    summarise(n=n(), recov=sum(status==1), death=sum(status==2), .groups="drop")
  cif <- sapply(seq_len(nrow(tb)), function(i) {
    Z <- R %>% filter(fr==tb$fr[i], G==tb$G[i])
    if (nrow(Z) < 3) return(NA_real_)
    s <- try(summary(survfit(Surv(t_end, St) ~ 1, data=Z), times=5, extend=TRUE), silent=TRUE)
    if (inherits(s,"try-error")) return(NA_real_)
    100*s$pstate[1, grep("recovery", colnames(s$pstate))] })
  tb$CIF5 <- round(cif, 1)
  cs <- try(summary(coxph(Surv(t_end, status==1) ~ cell + age_c + comorbid_count, data=R)), silent=TRUE)
  fgd <- try(finegray(Surv(t_end, St) ~ ., data=R, etype="recovery"), silent=TRUE)
  fg <- if (!inherits(fgd,"try-error"))
          try(summary(coxph(Surv(fgstart,fgstop,fgstatus) ~ cell, weight=fgwt, data=fgd)), silent=TRUE) else NULL
  gethr <- function(o, nm) if (is.null(o) || inherits(o,"try-error") ||
                               !(nm %in% rownames(o$conf.int))) "NE" else
    sprintf("%.2f (%.2f-%.2f)", o$conf.int[nm,1], o$conf.int[nm,3], o$conf.int[nm,4])
  key <- paste0("cell", tb$fr, "/", tb$G)
  tb$HR_cs  <- sapply(key, function(k) gethr(cs, k))
  tb$sHR_fg <- sapply(key, function(k) gethr(fg, k))
  tb$HR_cs[tb$fr=="Robust" & tb$G=="T3"]  <- "1.00 (reference)"
  tb$sHR_fg[tb$fr=="Robust" & tb$G=="T3"] <- "1.00 (reference)"
  print(as.data.frame(tb), row.names=FALSE)
  write.csv(tb, file.path(OUT, paste0("table3_recovery_", tv, ".csv")), row.names=FALSE)
  # 층별 per-s.d.
  for (f in levels(R$fr)) {
    Z <- R %>% filter(fr==f); if (nrow(Z) < 20) next
    zz <- if (tv=="T_sex") "z_sex" else "z_coh"; Z$Z <- Z[[zz]]
    o <- try(summary(coxph(Surv(t_end, status==1) ~ Z + age_c + comorbid_count, data=Z)), silent=TRUE)
    if (!inherits(o,"try-error"))
      say("  %-10s per +1 s.d. cause-specific HR %.2f (%.2f-%.2f), P = %.3g", f,
          o$coef["Z","exp(coef)"], o$conf.int["Z","lower .95"],
          o$conf.int["Z","upper .95"], o$coef["Z","Pr(>|z|)"])
  }
}
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
