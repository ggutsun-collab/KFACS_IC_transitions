###############################################################################
## KF_Table2_Transitions_v260801.R
## Table 2 - Intrinsic-capacity level and state transitions
##   전이별 보정 IRR: per +1 SD  및 삼분위 대비(T2 vs T1, T3 vs T1)
##   slope 항 없음. rolling 노출(각 interval 시작 시점의 IC).
###############################################################################

###############################################################################
## 이 스크립트는 작업디렉터리와 무관하게 동작합니다.
## 자기 위치를 찾아 R/_bootstrap.R 을 불러오고, kf_init() 이 헬퍼를 적재합니다.
## (source() / RStudio Source 버튼 / Rscript 모두 지원)
###############################################################################
## 위치 탐색과 헬퍼 로드는 R/_bootstrap.R 한 곳에서만 정의합니다.
local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of)) { d <- dirname(of); break }
  }
  cand <- c(if (!is.null(d)) file.path(d, c("R", "../R", "../../R", ".", "..")),
            "R", "../R", "../../R", ".", "..")
  for (p in cand)
    if (file.exists(file.path(p, "_bootstrap.R"))) {
      source(file.path(p, "_bootstrap.R"), encoding = "UTF-8"); break }
})
kf_init()
load_imputed("MAIN")
d <- prep_long(load_kfacs())

###############################################################################
## ★ 체계 이름 정정 (2026-07-31)
##
##  이전 라벨: state -> "Intrinsic-capacity state"  ← 틀렸습니다.
##  state 는 IC 가 아니라 **ADL 장애 중증도**입니다. 데이터로 확인했습니다:
##      state 1 = adl_3cat 0 (Normal, 9,913행)
##      state 2 = adl_3cat 1 (Mild,   1,405행)
##      state 3 = adl_3cat 2 (Severe,   139행)
##  IC 는 결과가 아니라 **노출**(gLIC_z / gLIC_tert)입니다.
##
##  adl_state 는 같은 ADL 을 이분화한 것(adl_score >= 1)이라
##  "Any ADL disability" 로 이름을 명확히 했습니다.
##  -> 라벨만 바뀝니다. 모형·사건수·person-years·추정치는 그대로입니다.
###############################################################################
blocks <- list(c("state",        "ADL disability, 3 states (normal / mild / severe)"),
               c("frailty_3cat", "Frailty phenotype (robust / pre-frail / frail)"),
               c("adl_state",    "Any ADL disability (binary)"))
LONG <- do.call(rbind, lapply(blocks, function(bk) {
  iv <- build_intervals(d, bk[1], "rolling")
  if (is.null(iv)) return(NULL)
  r <- rbind(run_transition_table(iv, exposure_terms = "gLIC_z"),
             run_transition_table(iv, exposure_terms = "gLIC_tert"))
  r$sysname <- bk[2]; r
}))
LONG$est <- fmt_est(LONG$IRR, LONG$lcl, LONG$ucl)
LONG$P   <- fmt_p(LONG$p)
LONG$col <- ifelse(LONG$term == "gLIC_z", "per +1 SD",
             ifelse(LONG$term == "gLIC_tertT2", "T2 vs T1",
             ifelse(LONG$term == "gLIC_tertT3", "T3 vs T1", LONG$term)))

W <- unique(LONG[, c("sysname", "ord", "transition", "n_events", "n_intervals", "pyears")])
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$ord), ]
for (cc in c("per +1 SD", "T2 vs T1", "T3 vs T1")) {
  s <- LONG[LONG$col == cc, c("sysname", "transition", "est", "P")]
  s <- s[!duplicated(s[, 1:2]), ]; names(s)[3:4] <- c(paste0("e_", cc), paste0("p_", cc))
  W <- merge(W, s, by = c("sysname", "transition"), all.x = TRUE, sort = FALSE)
}
W <- W[order(match(W$sysname, vapply(blocks, `[`, character(1), 2)), W$ord), ]
for (cc in grep("^e_", names(W), value = TRUE)) {
  pc <- sub("^e_", "p_", cc)
  bad <- is.na(W[[cc]]) | W[[cc]] %in% c("", "-")
  W[[cc]][bad] <- "NE"; W[[pc]][bad] <- "-"
}
NEROW <- LONG[LONG$term == "gLIC_tert", c("sysname", "transition", "note")]
NEROW <- NEROW[!duplicated(NEROW[, 1:2]) & nzchar(NEROW$note), ]

TAB <- data.frame(
  System = W$sysname, Transition = W$transition,
  `No. of events` = W$n_events, `Person-years` = round(W$pyears),
  `IRR per +1 s.d. (95% CI)` = W[["e_per +1 SD"]], `P` = W[["p_per +1 SD"]],
  `T2 vs T1 (95% CI)` = W[["e_T2 vs T1"]], `P ` = W[["p_T2 vs T1"]],
  `T3 vs T1 (95% CI)` = W[["e_T3 vs T1"]], `P  ` = W[["p_T3 vs T1"]],
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

save_table(TAB, "Table2_Transitions",
  title = "Table 2 | Intrinsic-capacity level and state-transition rates",
  footnotes = c(
    "Incidence rate ratios from transition-specific Poisson models with a log(person-time) offset, adjusted for age, sex, baseline comorbidity count, education, household income and area of residence; standard errors are clustered by participant.",
    "Exposure is the intrinsic capacity measured at the start of each interval (rolling landmark), expressed per +1 s.d. of the wave-1 distribution and as tertiles of that distribution.",
    "All observed off-diagonal transitions are listed. Rows are ordered by state, not alphabetically.",
    sprintf("NE, not estimable for the tertile contrast only: a tertile contributed fewer than %d events to that transition, so the coefficient diverges rather than converging. The total number of events is given and the per-s.d. estimate remains valid. Affected transitions: %s.",
            CFG$MIN_CELL_EVENTS,
            if (nrow(NEROW)) paste(unique(NEROW$transition), collapse = "; ") else "none"),
    "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation."))
save_vals(LONG, "Table2_values_long.csv", FIG_DIR)
cat("\n=== Table 2 완료 ===\n")
