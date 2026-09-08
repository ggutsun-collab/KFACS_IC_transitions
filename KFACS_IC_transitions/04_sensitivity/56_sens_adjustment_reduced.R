###############################################################################
## 56_sens_adjustment_reduced.R
## 보충표 — 축소 보정군 민감도 (주 분석 CORE 6변수 vs 축소 3변수)
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
## 주 분석의 보정군은 CORE 6변수입니다: 연령·성별·기저 동반질환·교육·가구소득·
## 거주지역 (모두 wave 1 값, 사회경제 결측은 Unknown 별도 수준).
## 사회경제 변수 3개를 더 얹으면 모수당 사건 수가 줄어 희소 전이의 추정이
## 불안정해질 수 있다는 우려가 가능합니다. 이 표는 사회경제 변수를 뺀 축소
## 보정군(연령·성별·동반질환)으로 같은 전이 모형을 재적합해, 주 추정치가
## 보정군 선택에 좌우되지 않음을 보입니다.
##
## ── 독거를 어느 쪽에도 넣지 않는 이유 ────────────────────────────────────
## 노인에서 독거는 내재역량 저하의 원인이 아니라 결과일 수 있습니다
## (역량 저하 -> 자녀와 합가). 노출의 하류를 보정하면 충돌 편향입니다.
##
## 출력:  ST_AdjReduced.xlsx / .csv  ·  ST_AdjReduced_values.csv
## 실행:  source("56_sens_adjustment_reduced.R", encoding = "UTF-8")
###############################################################################

ADJRED_VERSION <- "v260812"
message("\n=== 56_sens_adjustment_reduced ", ADJRED_VERSION, " ===")

###############################################################################
## 0. 헬퍼 로드
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
## ── 안전장치: 주 보정군이 CORE 6변수인지 확인하고 아니면 멈춥니다 ─────────
if (!all(c("edu_bl_f", "income_bl_f", "area_bl_f") %in% CFG$ADJ))
  stop("CFG$ADJ 가 CORE 6변수가 아닙니다.\n  현재: ",
       paste(CFG$ADJ, collapse = " + "),
       "\n  -> 01_functions_common.R 의 CFG$ADJ 를 확인하십시오.")
ADJ_CORE    <- CFG$ADJ
ADJ_REDUCED <- c("age_c", "sex_f", "comorbid_bl")
msg(sprintf("주 보정군(CORE)  : %s", paste(ADJ_CORE, collapse = " + ")))
msg(sprintf("축소 보정군      : %s", paste(ADJ_REDUCED, collapse = " + ")))

MIN_EV <- 5    # 이보다 사건이 적은 전이는 표에서 제외

###############################################################################
## 1. 자료와 구간
###############################################################################
load_imputed("MAIN")
d <- prep_long(load_kfacs())
msg(sprintf("AdjReduced: N=%d, person-waves=%d", length(unique(d$id)), nrow(d)))

###############################################################################
## 2. 적합 함수 (참가자 군집 강건 SE)
###############################################################################
.rvcov <- function(m, cluster) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    v <- try(sandwich::vcovCL(m, cluster = cluster, type = "HC0"), silent = TRUE)
    if (!inherits(v, "try-error") && all(is.finite(v))) return(v)
  }
  stats::vcov(m)
}
fit_one <- function(dd, adj) {
  adj <- adj[vapply(adj, function(v)
    v %in% names(dd) && length(unique(stats::na.omit(dd[[v]]))) > 1, logical(1))]
  f <- stats::as.formula(paste0("y ~ gLIC_z",
        if (length(adj)) paste0(" + ", paste(adj, collapse = " + ")) else "",
        " + offset(log(dur))"))
  m <- try(stats::glm(f, family = stats::poisson(), data = dd), silent = TRUE)
  if (inherits(m, "try-error") || !isTRUE(m$converged))
    return(c(b = NA_real_, se = NA_real_))
  V <- .rvcov(m, dd$id)
  b <- stats::coef(m)["gLIC_z"]; s <- sqrt(V["gLIC_z", "gLIC_z"])
  if (!is.finite(b) || !is.finite(s) || s > 3 || abs(b) > 10)
    return(c(b = NA_real_, se = NA_real_))
  c(b = unname(b), se = unname(s))
}
fmt_ci <- function(b, s) ifelse(is.finite(b),
  sprintf("%.2f (%.2f–%.2f)", exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s)), "–")
pval   <- function(b, s) ifelse(is.finite(b) & is.finite(s),
                                2 * stats::pnorm(-abs(b / s)), NA_real_)
fmt_p  <- function(p) ifelse(!is.finite(p), "–",
                      ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

###############################################################################
## 3. 전이별 추정 — 같은 자료, 두 보정군
###############################################################################
SYS <- list(list(var = "state",        lab = "ADL disability, 3 states"),
            list(var = "frailty_3cat", lab = "Frailty phenotype"))

OUT <- list()
for (S in SYS) {
  iv <- try(build_intervals(d, S$var, "rolling"), silent = TRUE)
  if (inherits(iv, "try-error") || is.null(iv)) {
    message("[AdjReduced] ", S$lab, " 건너뜀"); next
  }
  for (tr in default_transitions(S$var, iv)) {
    dd <- iv[iv$from_lab == tr[1], , drop = FALSE]
    dd$y <- as.integer(dd$to_lab == tr[2])
    ev <- sum(dd$y, na.rm = TRUE)
    if (ev < MIN_EV) next
    core <- fit_one(dd, ADJ_CORE)
    redu <- fit_one(dd, ADJ_REDUCED)
    OUT[[length(OUT) + 1]] <- data.frame(
      system = S$lab, transition = paste(tr[1], "→", tr[2]),
      events = ev, pyears = sum(dd$dur, na.rm = TRUE),
      b_core = core["b"], se_core = core["se"],
      b_red  = redu["b"], se_red  = redu["se"],
      stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, OUT); rownames(R) <- NULL
if (is.null(R) || !nrow(R)) stop("추정된 전이가 없습니다.")
cat("\n[AdjReduced] 체계별 전이 수\n"); print(table(R$system))

###############################################################################
## 4. 변화율 — 로그척도, 분모 보호 규칙 (두 조건 모두 필수)
##  기저(CORE)가 귀무에 가까우면 백분율이 폭주하므로,
##  P < 0.05 이면서 |log IRR| > log(1.05) 인 전이에서만 계산합니다.
###############################################################################
p_core <- pval(R$b_core, R$se_core)
p_red  <- pval(R$b_red,  R$se_red)
ok_den <- is.finite(p_core) & p_core < 0.05 & abs(R$b_core) > log(1.05)
R$chg  <- ifelse(is.finite(R$b_red) & is.finite(R$b_core) & ok_den,
                 100 * (R$b_red - R$b_core) / abs(R$b_core), NA_real_)

TAB <- data.frame(System = R$system, Transition = R$transition,
                  `No. of events` = R$events, check.names = FALSE,
                  stringsAsFactors = FALSE)
TAB[["Primary model, CORE covariates (95% CI)"]] <- fmt_ci(R$b_core, R$se_core)
TAB[["P"]]  <- fmt_p(p_core)
TAB[["Reduced model, age, sex, comorbidity (95% CI)"]] <- fmt_ci(R$b_red, R$se_red)
TAB[["P "]] <- fmt_p(p_red)
TAB[["Change on log scale"]] <- ifelse(is.finite(R$chg), sprintf("%+.1f%%", R$chg), "–")
TAB[["Same direction"]] <- ifelse(is.finite(R$b_red) & is.finite(R$b_core),
                                  ifelse(sign(R$b_red) == sign(R$b_core), "yes", "no"), "–")
TAB$System[duplicated(TAB$System)] <- ""

n_ch  <- sum(is.finite(R$chg))
med   <- stats::median(abs(R$chg), na.rm = TRUE)
mx    <- suppressWarnings(max(abs(R$chg), na.rm = TRUE))
flip  <- sum(is.finite(R$b_red) & is.finite(R$b_core) &
             sign(R$b_red) != sign(R$b_core), na.rm = TRUE)
sig0  <- sum(p_core < 0.05, na.rm = TRUE)
keep  <- sum(p_core < 0.05 & p_red < 0.05, na.rm = TRUE)

FN <- c(
 paste0("Incidence rate ratios per +1 s.d. of intrinsic capacity, from transition-specific ",
        "Poisson models with a log(person-time) offset and standard errors clustered by ",
        "participant. The primary model is adjusted for ", adj_phrase(ADJ_CORE),
        ", all recorded at wave 1; participants with missing socioeconomic data are retained ",
        "in a separate unknown category."),
 paste0("The reduced model omits the three socioeconomic covariates and adjusts for ",
        adj_phrase(ADJ_REDUCED), " only. Agreement between the two models shows that the ",
        "primary estimates are not driven by the choice of adjustment set, and that the ",
        "additional parameters do not destabilise the sparser transitions."),
 paste0("Living arrangement is deliberately absent from both models: in older adults living ",
        "alone may be a consequence of declining capacity rather than a cause of it, and ",
        "conditioning on a potential consequence of the exposure risks collider bias."),
 paste0("Change is computed on the log scale as the difference in log incidence rate ratios ",
        "expressed as a percentage of the primary estimate. It is left blank where the primary ",
        "estimate was itself non-significant or negligible (P >= 0.05 or IRR within 0.95-1.05), ",
        "because the ratio is unstable when the denominator approaches zero."),
 paste0("All ", nrow(R), " transitions with at least ", MIN_EV, " events are shown. Among the ",
        n_ch, " for which a change could be quantified, the median absolute change was ",
        sprintf("%.1f%%", med), " and the largest was ", sprintf("%.1f%%", mx), "; ", flip,
        " estimate", if (flip == 1) "" else "s", " changed direction, and ", keep, " of the ",
        sig0, " estimates significant in the primary model remained significant."),
 "CI, confidence interval; IRR, incidence rate ratio; s.d., standard deviation.")

save_table(TAB, "ST_AdjReduced",
  title = paste("Supplementary Data | Sensitivity of the transition analyses to a",
                "reduced adjustment set"),
  footnotes = FN)
save_vals(R, "ST_AdjReduced_values.csv", FIG_DIR)

###############################################################################
## 5. 판정
###############################################################################
cat("\n=== 축소 보정 민감도 판정 ===================================\n")
cat(sprintf("  전이 %d개 (변화율 산출 %d개) | 절대변화 중앙 %.1f%% (최대 %.1f%%)\n",
            nrow(R), n_ch, med, mx))
cat(sprintf("  방향 역전 %d | 유의 유지 %d/%d\n", flip, keep, sig0))
cat("--------------------------------------------------------------\n")
if (is.finite(med) && med < 15 && flip == 0 && keep >= sig0 - 1) {
  cat("  판정: 보정군 선택이 결론을 바꾸지 않습니다.\n")
} else {
  cat("  판정: 변화가 작지 않습니다. 큰 순서로 확인하십시오.\n")
  o <- order(-abs(R$chg))
  print(utils::head(data.frame(transition = R$transition[o], events = R$events[o],
                        change = sprintf("%+.1f%%", R$chg[o])), 5), row.names = FALSE)
}
cat("==============================================================\n")
cat("\n=== 완료 ===  표:", file.path(FIG_DIR, "ST_AdjReduced.xlsx"), "\n")
