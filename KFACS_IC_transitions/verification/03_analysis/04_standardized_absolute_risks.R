# =============================================================================
# 17 — 보정 절대위험 (한 파일 완결판)
#
# 돌고 있는 15번은 중단하셔도 됩니다. 이 파일 하나면 끝납니다.
#   - 프레임 구성부터 다시 하므로 K 가 없어도 됩니다 (d, iv 만 있으면 됩니다)
#   - 15b 의 빠른 구현을 내장했습니다 (survfit 반복 호출 제거)
#   - 부트스트랩 1000회 x 4조합, 보통 1-2분
#   - 느린 구현이 이미 낸 값과 자동 대조하여 구현이 맞는지 먼저 확인합니다
#
# 산출: 조위험(KM) 과 보정 절대위험(g-computation)을
#       코호트 척도 / 성별 내 척도 x 사망 / 복합악화 의 4조합으로.
# 그대로 source.
# =============================================================================
while (sink.number() > 0) sink()
suppressPackageStartupMessages({ for (p in c("survival","dplyr")) library(p, character.only=TRUE) })
OUT <- "stdrisk_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "stdrisk_ALL_log.txt"), split = TRUE)
say  <- function(...) cat(sprintf(...), "\n", sep="")
rule <- function(t) cat("\n", strrep("=",76), "\n", t, "\n", strrep("=",76), "\n", sep="")
norm <- function(x){ x <- as.character(x)
  ifelse(x %in% c("0","Normal"),"Normal", ifelse(x %in% c("1","Mild"),"Mild",
  ifelse(x %in% c("2","Severe"),"Severe", ifelse(x %in% c("3","Death"),"Death", x)))) }

H  <- 5      # 지평선(년)
B  <- 1000   # 부트스트랩 반복수

## ===========================================================================
rule("1. 분석 프레임 구성  (확정 규칙: 사망시점 = wave 결측 행의 time, 상태 = state_lab)")
V <- d %>% filter(!is.na(wave)) %>% arrange(id, wave) %>%
  group_by(id) %>% mutate(kv = row_number()) %>% ungroup() %>% mutate(.st = norm(state_lab))
S <- d %>% group_by(id) %>%
  summarise(ev = as.integer(any(death_event>0, na.rm=TRUE)),
            td = suppressWarnings(max(time[is.na(wave)], na.rm=TRUE)),
            tl = suppressWarnings(max(time[!is.na(wave)], na.rm=TRUE)), .groups="drop") %>%
  mutate(td = ifelse(is.finite(td), td, NA), tl = ifelse(is.finite(tl), tl, NA),
         t_fu = ifelse(ev==1, td, tl))
ORD <- c("Normal","Mild","Severe")
ADL <- V %>% arrange(id, kv) %>% group_by(id) %>%
  mutate(base = first(.st)) %>%
  summarise(t_worse = { w <- which(kv > 1 & match(.st, ORD) > match(base, ORD))
                        if (length(w)) time[w[1]] else NA_real_ }, .groups="drop")

K <- V %>% group_by(id) %>% slice_min(kv, n=1, with_ties=FALSE) %>% ungroup() %>%
  select(id, fr_w1 = frailty_3cat) %>%
  left_join(iv %>% group_by(id) %>% slice_min(k, n=1, with_ties=FALSE) %>% ungroup() %>%
              select(id, gLIC_z, sex_f, age_c, comorbid_count, edu_bl_f, income_bl_f, area_bl_f),
            by="id") %>%
  left_join(S %>% select(id, ev, t_fu), by="id") %>% left_join(ADL, by="id") %>%
  filter(!is.na(gLIC_z), !is.na(t_fu), t_fu > 0) %>%
  mutate(e_comp = as.integer(!is.na(t_worse) | ev == 1),
         t_comp = pmin(ifelse(is.na(t_worse), Inf, t_worse), t_fu),
         fr = factor(as.character(fr_w1), c("0","1","2"), c("Robust","Pre-frail","Frail")))
cut_coh <- quantile(K$gLIC_z, c(1/3,2/3), na.rm=TRUE)
K$T_coh <- cut(K$gLIC_z, c(-Inf, cut_coh, Inf), labels=c("T1","T2","T3"))
K <- K %>% group_by(sex_f) %>%
  mutate(T_sex = cut(gLIC_z, c(-Inf, quantile(gLIC_z, c(1/3,2/3), na.rm=TRUE), Inf),
                     labels=c("T1","T2","T3"))) %>% ungroup()
say("n = %d | 사망 %d | 복합악화 %d   [원고 3,011 / 500]",
    nrow(K), sum(K$ev), sum(K$e_comp))
say("frailty: %s", paste(names(table(K$fr)), table(K$fr), sep="=", collapse="  "))
say("코호트 척도 tertile: %s | 성별 내 척도 tertile: %s",
    paste(table(K$T_coh), collapse="/"), paste(table(K$T_sex), collapse="/"))
assign("K", K, envir = .GlobalEnv)

## ===========================================================================
## 2. 보정 절대위험 (g-computation, 빠른 구현)
##    risk_g = 1 - mean_i exp( -H0(H) * exp( lp_i with cell set to g ) )
## ===========================================================================
stdrisk <- function(dat, tv, ev, tertv, covs) {
  dat$.t <- pmin(dat[[tv]], H); dat$.e <- ifelse(dat[[tv]] > H, 0L, dat[[ev]])
  dat$cell <- droplevels(interaction(dat$fr, dat[[tertv]], sep = "/"))
  m  <- coxph(as.formula(paste("Surv(.t,.e) ~ cell +", paste(covs, collapse="+"))), data = dat)
  bh <- basehaz(m, centered = TRUE)
  H0 <- if (any(bh$time <= H)) max(bh$hazard[bh$time <= H]) else 0
  lp0 <- predict(m, type = "lp"); b <- coef(m); lv <- levels(dat$cell)
  cf <- setNames(c(0, b[paste0("cell", lv[-1])]), lv); cf[is.na(cf)] <- 0
  own <- cf[as.character(dat$cell)]
  vapply(lv, function(g) 100*(1 - mean(exp(-H0 * exp(lp0 - own + cf[[g]])))), numeric(1))
}
boot_ci <- function(dat, tv, ev, tertv, covs, seed = 20260818) {
  set.seed(seed); pt <- stdrisk(dat, tv, ev, tertv, covs)
  bs <- matrix(NA_real_, B, length(pt), dimnames = list(NULL, names(pt)))
  for (i in seq_len(B)) {
    j <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    r <- try(stdrisk(dat[j, ], tv, ev, tertv, covs), silent = TRUE)
    if (!inherits(r, "try-error")) bs[i, names(r)] <- r }
  data.frame(cell = names(pt), risk = round(pt, 1),
             lo = round(apply(bs, 2, quantile, .025, na.rm=TRUE), 1),
             hi = round(apply(bs, 2, quantile, .975, na.rm=TRUE), 1), row.names = NULL)
}
crude <- function(dat, tv, ev, tertv) {
  dat$cell <- droplevels(interaction(dat$fr, dat[[tertv]], sep="/"))
  do.call(rbind, lapply(levels(dat$cell), function(g) {
    z <- dat[dat$cell == g, , drop=FALSE]; if (!nrow(z)) return(NULL)
    z$.t <- pmin(z[[tv]], H); z$.e <- ifelse(z[[tv]] > H, 0L, z[[ev]])
    s <- summary(survfit(Surv(.t,.e) ~ 1, data=z), times=H, extend=TRUE)
    data.frame(cell=g, n=nrow(z), events=sum(z$.e), crude=round(100*(1-s$surv[1]),1)) }))
}

CV_coh <- c("age_c","comorbid_count","edu_bl_f","income_bl_f","area_bl_f","sex_f")
CV_sex <- c("age_c","comorbid_count","edu_bl_f","income_bl_f","area_bl_f")

## ===========================================================================
rule("2. 구현 검증 — 느린 구현이 낸 값과 대조 (코호트 척도 · 사망)")
chk <- round(stdrisk(K, "t_fu", "ev", "T_coh", CV_coh), 1)
ref <- c("Robust/T1"=3.4,"Robust/T2"=6.9,"Robust/T3"=2.8,
         "Pre-frail/T1"=9.1,"Pre-frail/T2"=6.2,"Pre-frail/T3"=4.6,
         "Frail/T1"=14.8,"Frail/T2"=2.8,"Frail/T3"=33.9)
cmp <- data.frame(cell=names(chk), fast=as.numeric(chk), slow=ref[names(chk)], row.names=NULL)
cmp$diff <- round(cmp$fast - cmp$slow, 1)
print(cmp, row.names = FALSE)
ok <- all(abs(cmp$diff) <= 0.3, na.rm = TRUE)
say(ifelse(ok, ">>> 일치. 아래 결과를 그대로 쓰셔도 됩니다.",
               ">>> 불일치. 아래 결과를 쓰지 마시고 이 표를 알려주십시오."))

## ===========================================================================
res <- list()
for (tertv in c("T_coh","T_sex")) {
  covs <- if (tertv=="T_coh") CV_coh else CV_sex
  rule(sprintf("3. %s — 조위험(KM) 대 보정 절대위험 (%d년, 부트스트랩 %d회)",
       ifelse(tertv=="T_coh","코호트 척도 tertile (원고 Table 3)","성별 내 척도 tertile"), H, B))
  for (o in list(c("t_fu","ev","Death"), c("t_comp","e_comp","Composite worsening"))) {
    say("\n---- %s ----", o[3])
    z <- merge(crude(K,o[1],o[2],tertv), boot_ci(K,o[1],o[2],tertv,covs), by="cell", all.x=TRUE)
    z$shift <- round(z$risk - z$crude, 1)
    z <- z[order(z$cell), ]; print(z, row.names = FALSE)
    write.csv(z, file.path(OUT, sprintf("stdrisk_%s_%s.csv", tertv, o[2])), row.names=FALSE)
    res[[paste(tertv, o[3])]] <- z
  }
}

## ===========================================================================
rule("4. 판정")
for (k in names(res)) {
  z <- res[[k]]
  g <- function(pat) { x <- z[grepl(pat, z$cell), ]; x[order(x$cell), ] }
  rb <- g("^Robust"); pf <- g("^Pre-frail")
  say("[%s]", k)
  say("  robust    조 %-18s -> 보정 %-18s (%s)",
      paste(sprintf("%.1f", rb$crude), collapse="/"),
      paste(sprintf("%.1f", rb$risk),  collapse="/"),
      ifelse(all(diff(rb$risk) < 0), "단조", "비단조"))
  say("            사건수 %s", paste(rb$events, collapse=" / "))
  say("  pre-frail 조 %-18s -> 보정 %-18s (%s)",
      paste(sprintf("%.1f", pf$crude), collapse="/"),
      paste(sprintf("%.1f", pf$risk),  collapse="/"),
      ifelse(all(diff(pf$risk) < 0), "단조", "비단조"))
}
say("\n>>> 판단 기준")
say("   1) 보정으로 pre-frail 층에 gradient 가 생기는가  -> Table 3 위험 열 교체 근거")
say("   2) robust 층이 보정 후에도 비단조인가            -> 사건 수 문제이며 각주로 처리")
say("   3) 성별 내 척도에서 robust 사건 수가 몇 배인가   -> Extended Data 정당화")
rule("끝"); say("출력: %s", normalizePath(OUT)); sink()
