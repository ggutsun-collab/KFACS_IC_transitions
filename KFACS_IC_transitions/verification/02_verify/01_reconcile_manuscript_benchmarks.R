# =============================================================================
# 02 — 원고 숫자 vs 실제 데이터 대조
#
# 환경 스캔에서 드러난 것: 재분석은 새 숫자를 준 것이 아니라,
# 원고의 Table 1 / Table 3 / 초록 분모가 현재 데이터와 맞지 않는다는 것을 드러냈습니다.
#
#   원고                          dat_w1 (실제)
#   robust 44.8% (= 1,349)   ->   Robust 1,082 (35.9%)
#   deaths 500 (16.6%)       ->   d_death 합계 556 (18.5%)
#   robust 5년 사망 51건     ->   R8 기준 120건으로 보임
#   ADL 장애 robust 103명    ->   rc/csd n = 109
#
# 이 스크립트는 위 네 가지를 확정하고, 원인이 사망시점 배정 규칙인지 확인합니다.
# 그대로 source 하십시오. 편집할 곳 없습니다.
# =============================================================================

suppressPackageStartupMessages({ library(survival); library(dplyr) })

OUT <- "reconcile_out"; dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "reconcile_log.txt"), split = TRUE)

say  <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function(t) cat("\n", strrep("=", 76), "\n", t, "\n", strrep("=", 76), "\n", sep = "")
dump1 <- function(nm) {
  if (!exists(nm, envir = .GlobalEnv)) { cat("\n[", nm, "] 없음\n", sep = ""); return(invisible()) }
  cat("\n---- ", nm, " ----\n", sep = "")
  print(get(nm, envir = .GlobalEnv), digits = 4)
}

## ---------------------------------------------------------------------------
rule("0. 기존 결과 객체 전체 출력 (해석은 아래에서)")
for (o in c("cnt","tab3","R1","R3","R3b","R4","R5","R6","R8","nn","lrt","ref")) dump1(o)
if (exists("spline_curve")) {
  sc <- spline_curve
  say("\nspline_curve: %d 행, ic 범위 %.2f ~ %.2f", nrow(sc), min(sc$ic), max(sc$ic))
  write.csv(sc, file.path(OUT, "spline_curve_full.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------------------
rule("1. 코호트 기본 수치 — 원고 Table 1 / 초록과 대조")
W <- dat_w1
say("n = %d", nrow(W))
print(table(W$frailty_w1, useNA = "ifany"))
say("robust 비율 = %.1f%%   [원고: 44.8%%]", 100*mean(W$frailty_w1 == "Robust", na.rm = TRUE))
say("총 사망 = %d (%.1f%%)   [원고 초록: 500 (16.6%%)]", sum(W$d_death), 100*mean(W$d_death))
say("추적기간 중앙값 = %.2f년, 최대 = %.2f년", median(W$t_death), max(W$t_death))
say("ADL 장애 보유 & robust (원고 103명):")
if ("adl_w1" %in% names(W)) {
  tb <- W %>% filter(frailty_w1 == "Robust") %>% count(adl_w1)
  print(as.data.frame(tb))
  say("  -> robust 중 adl_w1 != 정상 인원 = %d   [원고: 103, rc/csd: 109]",
      sum(W$frailty_w1 == "Robust" & !(W$adl_w1 %in% c(0, "0", "none", "None", "Normal", "normal")), na.rm = TRUE))
}

## ---------------------------------------------------------------------------
rule("2. 5년 vs 8년 사망 — frailty 층별 (원고 Table 3 재현 시도)")
km_risk <- function(d, h) {
  s <- survfit(Surv(pmin(t_death, h), ifelse(t_death > h, 0L, d_death)) ~ 1, data = d)
  100 * (1 - summary(s, times = h, extend = TRUE)$surv[1])
}
tabf <- W %>% group_by(frailty_w1) %>%
  summarise(n = n(),
            deaths_5y = sum(d_death[t_death <= 5]),
            deaths_8y = sum(d_death[t_death <= 8]),
            deaths_all = sum(d_death), .groups = "drop") %>%
  rowwise() %>%
  mutate(km5 = km_risk(W[W$frailty_w1 == frailty_w1, ], 5),
         km8 = km_risk(W[W$frailty_w1 == frailty_w1, ], 8)) %>% ungroup()
print(as.data.frame(tabf), digits = 3)
write.csv(tabf, file.path(OUT, "deaths_by_frailty.csv"), row.names = FALSE)
say("\n>> 원고 Table 3 의 robust 층 5년 사망은 2+23+26 = 51 건입니다.")
say(">> 위 표의 robust deaths_5y 와 비교하십시오. 2배 이상 차이 나면 Table 3 은 재생성 대상입니다.")
say(">> 51/1349 = 3.8%% = 연 0.78%% 는 평균 76세 코호트에서 생물학적으로 불가능합니다.")

## ---------------------------------------------------------------------------
rule("3. robust 층 IC tertile — 원고 Table 3 직접 재현")
RB <- W %>% filter(frailty_w1 == "Robust")
cuts <- quantile(W$gLIC_z, c(1/3, 2/3), na.rm = TRUE)
RB$tert <- cut(RB$gLIC_z, c(-Inf, cuts, Inf), labels = c("T1","T2","T3"))
t3 <- RB %>% group_by(tert) %>%
  summarise(n = n(), d5 = sum(d_death[t_death <= 5]), d8 = sum(d_death[t_death <= 8]), .groups="drop") %>%
  rowwise() %>% mutate(km5 = km_risk(RB[RB$tert == tert, ], 5),
                       km8 = km_risk(RB[RB$tert == tert, ], 8)) %>% ungroup()
print(as.data.frame(t3), digits = 3)
write.csv(t3, file.path(OUT, "robust_tertiles_reproduced.csv"), row.names = FALSE)
say("\n원고 Table 3 : T1 n=106 d=2 (1.9%%) | T2 n=470 d=23 (4.9%%) | T3 n=773 d=26 (3.4%%)")
say(">> n 도 사망도 재현되지 않으면, Table 3 은 다른 데이터 빌드에서 나온 것입니다.")
say(">> 그 경우 T1 의 '2건 이상'은 통계적 발견이 아니라 빌드 불일치의 산물입니다.")
say("\n5년 단조성 : %s", paste(sprintf("%s=%.1f%%", t3$tert, t3$km5), collapse = "  "))
say("8년 단조성 : %s", paste(sprintf("%s=%.1f%%", t3$tert, t3$km8), collapse = "  "))

## ---------------------------------------------------------------------------
rule("4. 사망시점 배정 규칙이 원인인지 검정")
if (exists("d") && "death_time_method" %in% names(d)) {
  say("d$death_time_method 분포:")
  print(table(d$death_time_method, useNA = "ifany"))
}
if (exists("d") && all(c("death_wave","death_date_new") %in% names(d))) {
  say("\n>> 사망시점을 '사망이 기록된 방문'에 배정하면, wave 3-4 사이(4~6년차) 사망이")
  say(">> 6년차 방문에 배정되어 5년 창 밖으로 밀려납니다. 이것이 Table 3 의 5년 사망을")
  say(">> 체계적으로 과소계수한 원인일 가능성이 큽니다.")
}
say("\n같은 데이터에서 배정 규칙만 바꾼 5년 사망 건수 비교:")
if ("t_death" %in% names(W)) {
  for (shift in c(0, -1, 1)) {
    tt <- W$t_death + shift
    say("  t_death %+d년 이동 시 robust 5년 사망 = %d",
        shift, sum(W$d_death[W$frailty_w1 == "Robust" & tt <= 5]))
  }
  say("  (1년 이동만으로 건수가 크게 바뀌면 5년 창은 배정규칙에 극도로 민감합니다)")
}

## ---------------------------------------------------------------------------
rule("5. 스플라인 형태 — 저IC 구간이 실제로 더 위험한가")
if (exists("spline_curve")) {
  sc <- spline_curve
  qq <- quantile(rb$.ic_coh, c(.05,.10,.25,.50,.75,.90,.95), na.rm = TRUE)
  idx <- sapply(qq, function(v) which.min(abs(sc$ic - v)))
  pt <- data.frame(pct = names(qq), ic = round(sc$ic[idx],3),
                   HR = sc$HR[idx], lo = sc$lo[idx], hi = sc$hi[idx])
  print(pt, digits = 3, row.names = FALSE)
  write.csv(pt, file.path(OUT, "spline_percentiles.csv"), row.names = FALSE)
  hr5 <- sc$HR[idx[1]]; hr25 <- sc$HR[idx[3]]; hr50 <- sc$HR[idx[4]]
  say("")
  if (hr5 > hr25 && hr25 > hr50) {
    say(">>> 판정: 단조 증가. 저IC 로 갈수록 위험이 더 가파릅니다.")
    say(">>> 비선형성은 '선형 HR 0.63 이 오히려 보수적'이라는 뜻이며, 논문에 유리합니다.")
    say(">>> Results 에 '위험 증가가 낮은 용량 구간에서 가팔라졌다'로 서술하십시오.")
  } else if (hr5 < hr25) {
    say(">>> 판정: 저IC 끝에서 평탄화 또는 역전. 8년으로 늘려도 해결되지 않았습니다.")
    say(">>> HR 0.63 은 선형가정이 이 구간을 덮은 오설정 요약치입니다.")
    say(">>> 원인 후보: (a) 생존자 선택, (b) 저IC-robust 표현형의 희소성,")
    say(">>>            (c) 측정 바닥효과. 어느 쪽이든 Limitations 에 명시해야 합니다.")
  } else {
    say(">>> 판정 보류. 곡선 PDF 를 직접 보십시오.")
  }
} else say("spline_curve 객체가 없습니다.")

## ---------------------------------------------------------------------------
rule("6. rb 프레임 무결성 확인")
if (exists("rb")) {
  say("rb: n = %d, 사망 = %d", nrow(rb), sum(rb$d_death))
  say("rb 는 dat_w1 의 robust 전체와 일치합니까? %s",
      if (nrow(rb) == sum(W$frailty_w1 == "Robust")) "예 (제외 없음, 완전사례 편향 없음)" else "아니오 — 제외 발생")
  for (v in c(".ic_coh", ".ic_rob", ".ic_sex")) if (v %in% names(rb))
    say("  %s : mean %.3f, sd %.3f", v, mean(rb[[v]], na.rm=TRUE), sd(rb[[v]], na.rm=TRUE))
  say("gLIC_z 전체코호트 sd = %.4f, robust 층 sd = %.4f (비 %.3f)",
      sd(W$gLIC_z, na.rm=TRUE), sd(rb$gLIC_z, na.rm=TRUE),
      sd(rb$gLIC_z, na.rm=TRUE)/sd(W$gLIC_z, na.rm=TRUE))
}

rule("끝")
say("출력: %s", normalizePath(OUT))
sink()
