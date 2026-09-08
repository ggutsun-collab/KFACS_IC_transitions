###############################################################################
## KF_CC_Sensitivity_v260812.R
## 완전자료 민감도 — 무응답 처리 방식이 결과를 만들었는가
##
## ── 왜 다시 재는가 ───────────────────────────────────────────────────────
## 패치의 [6b] 가 "계수 변화 중앙값 -0.3%, 최대 -27.9%" 를 냈습니다.
## 그런데 그 -27.9% 는 자료가 아니라 나눗셈이 만든 숫자일 수 있습니다.
##
## 이 프로젝트에서 같은 함정을 이미 두 번 겪었습니다 (Table 4 감쇠율의
## -74.3%, ST_AdjSensitivity 의 동일 문제). 그때 정한 규칙은
##
##     기저 추정치가 귀무에 가까우면 변화율을 계산하지 않는다
##     ok_den <- p0 < 0.05 & abs(b0) > log(1.05)      ← 두 조건 모두
##
## 였는데, [6b] 코드에는 뒤쪽 절반만 들어갔습니다. Severe -> Death 는
## IRR 0.81 (log = -0.211) 이라 0.81 -> 0.86 이라는 작은 이동도 -28% 로
## 찍힙니다.
##
## 그러나 다른 가능성도 있습니다. 완전자료는 192명을 빼는데, 사건이
## 21건뿐인 전이에서 1-2건이 빠지면 계수가 실제로 움직입니다. 그건
## 편향이 아니라 소표본 잡음입니다.
##
## 어느 쪽인지는 추측이 아니라 전이 이름과 사건 수를 봐야 압니다.
## 이 스크립트는 그것만 합니다. 아무것도 바꾸지 않습니다.
##
## 전제: KF_PATCH_adjset_v260812 를 적용해 CFG$ADJ 가 CORE 인 상태
## 실행: source("KF_CC_Sensitivity_v260812.R", encoding = "UTF-8")
###############################################################################

message("\n=== 완전자료 민감도 v260812 ===")
MIN_EV   <- 5
OUT_STEM <- "ST_CompleteCase"

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
if (!any(c("edu_bl_f", "income_bl_f", "area_bl_f", "ses_index") %in% CFG$ADJ))
  stop("CFG$ADJ 에 사회경제 변수가 없습니다. 먼저 KF_PATCH_adjset_v260812 를 적용하십시오.\n",
       "  현재: ", paste(CFG$ADJ, collapse = " + "))
cat(sprintf("\n[0] 보정군 : %s\n", paste(CFG$ADJ, collapse = " + ")))

OUT_DIR <- FIG_DIR   ## 모든 산출물은 00_setup.R 의 단일 출력 폴더로

###############################################################################
## 1. 적합
###############################################################################
.rvcov <- function(m, cl) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    v <- try(sandwich::vcovCL(m, cluster = cl, type = "HC0"), silent = TRUE)
    if (!inherits(v, "try-error") && all(is.finite(v))) return(v)
  }
  stats::vcov(m)
}
fit1 <- function(dd) {
  adj <- CFG$ADJ[vapply(CFG$ADJ, function(v)
    v %in% names(dd) && length(unique(stats::na.omit(dd[[v]]))) > 1, logical(1))]
  f <- stats::as.formula(paste0("y ~ gLIC_z",
        if (length(adj)) paste0(" + ", paste(adj, collapse = " + ")) else "",
        " + offset(log(dur))"))
  m <- try(suppressWarnings(stats::glm(f, family = stats::poisson(), data = dd)),
           silent = TRUE)
  if (inherits(m, "try-error") || !isTRUE(m$converged)) return(c(b = NA, se = NA))
  V <- .rvcov(m, dd$id)
  b <- stats::coef(m)["gLIC_z"]
  s <- if ("gLIC_z" %in% rownames(V)) sqrt(V["gLIC_z", "gLIC_z"]) else NA_real_
  if (!is.finite(b) || !is.finite(s) || s > 3 || abs(b) > 10) return(c(b = NA, se = NA))
  c(b = unname(b), se = unname(s))
}
## Unknown 이 하나라도 있는 행을 제외 = 완전자료
keep_cc <- function(dd) {
  k <- rep(TRUE, nrow(dd))
  for (v in CFG$ADJ)
    if (is.factor(dd[[v]]) && "Unknown" %in% levels(dd[[v]])) k <- k & dd[[v]] != "Unknown"
  k
}

load_imputed("MAIN")
d <- prep_long(load_kfacs())

SYS <- list(list(var = "state", lab = "ADL disability"),
            list(var = "frailty_3cat", lab = "Frailty phenotype"))
ROWS <- list()
for (S in SYS) {
  iv <- try(build_intervals(d, S$var, "rolling"), silent = TRUE)
  if (inherits(iv, "try-error")) next
  for (tr in default_transitions(S$var, iv)) {
    dd <- iv[iv$from_lab == tr[1], , drop = FALSE]
    dd$y <- as.integer(dd$to_lab == tr[2])
    ev  <- sum(dd$y, na.rm = TRUE)
    if (ev < MIN_EV) next
    k   <- keep_cc(dd)
    a   <- fit1(dd)
    cc  <- fit1(dd[k, , drop = FALSE])
    ROWS[[length(ROWS) + 1]] <- data.frame(
      system = S$lab, transition = paste(tr[1], "\u2192", tr[2]),
      ev_all = ev, ev_cc = sum(dd$y[k], na.rm = TRUE),
      n_all  = nrow(dd), n_cc = sum(k),
      b_all = a["b"], se_all = a["se"], b_cc = cc["b"], se_cc = cc["se"],
      stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, ROWS); rownames(R) <- NULL
if (is.null(R) || !nrow(R)) stop("추정된 전이가 없습니다.")

###############################################################################
## 2. 변화율 — 프로젝트 표준 규칙을 그대로 적용합니다
##    기저가 귀무에 가까우면(P >= 0.05 이거나 IRR 이 0.95-1.05) 비웁니다.
##    분모가 0 에 가까울 때 나오는 큰 백분율은 자료가 아니라 나눗셈입니다.
###############################################################################
pval <- function(b, s) ifelse(is.finite(b) & is.finite(s), 2 * stats::pnorm(-abs(b / s)), NA)
R$p_all  <- pval(R$b_all, R$se_all)
R$ok_den <- is.finite(R$p_all) & R$p_all < 0.05 & abs(R$b_all) > log(1.05)
R$chg    <- ifelse(is.finite(R$b_cc) & is.finite(R$b_all) & R$ok_den,
                   100 * (R$b_cc - R$b_all) / abs(R$b_all), NA_real_)
R$irr_all <- exp(R$b_all); R$irr_cc <- exp(R$b_cc)
R$ev_lost <- R$ev_all - R$ev_cc
R$n_lost  <- R$n_all - R$n_cc

fci <- function(b, s) ifelse(is.finite(b),
  sprintf("%.2f (%.2f\u2013%.2f)", exp(b), exp(b - 1.96*s), exp(b + 1.96*s)), "\u2014")
TAB <- data.frame(
  System      = R$system,
  Transition  = R$transition,
  `Events (all / complete)` = sprintf("%d / %d", R$ev_all, R$ev_cc),
  `Events lost`             = R$ev_lost,
  `All records`             = fci(R$b_all, R$se_all),
  `Complete cases`          = fci(R$b_cc,  R$se_cc),
  `Change`  = ifelse(is.finite(R$chg), sprintf("%+.1f%%", R$chg), "\u2014 (기저 귀무)"),
  `Same direction` = ifelse(is.finite(R$b_cc) & is.finite(R$b_all),
                            ifelse(sign(R$b_cc) == sign(R$b_all), "yes", "no"), "\u2014"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

cat("\n[1] 전이별 비교\n")
print(TAB, row.names = FALSE)

###############################################################################
## 3. 판정
###############################################################################
cat("\n[2] 판정\n")
cat(sprintf("  무응답으로 제외된 구간 : %d개 (전체의 %.1f%%)\n",
            max(R$n_lost), 100 * max(R$n_lost) / max(R$n_all)))
ok <- is.finite(R$chg)
cat(sprintf("  변화율을 계산한 전이 : %d개 / %d개\n", sum(ok), nrow(R)))
cat(sprintf("  계산하지 않은 전이 %d개 — 기저가 귀무라 분모가 작습니다: %s\n",
            sum(!ok), paste(R$transition[!ok], collapse = ", ")))
if (any(ok)) {
  cat(sprintf("  변화율 중앙값 %+.1f%% · 최대 %+.1f%% (%s)\n",
              stats::median(R$chg[ok]),
              R$chg[ok][which.max(abs(R$chg[ok]))],
              R$transition[ok][which.max(abs(R$chg[ok]))]))
}
flip <- sum(is.finite(R$b_cc) & is.finite(R$b_all) & sign(R$b_cc) != sign(R$b_all))
sigc <- sum(is.finite(R$b_cc) & is.finite(R$b_all) &
            ((R$p_all < .05) != (pval(R$b_cc, R$se_cc) < .05)), na.rm = TRUE)
cat(sprintf("  방향 역전 %d개 · 유의성 변화 %d개\n", flip, sigc))

## 귀무 전이는 백분율 대신 IRR 의 절대 이동으로 봅니다
cat("\n[3] 기저가 귀무인 전이 — 백분율 대신 IRR 자체를 보십시오\n")
nn <- R[!R$ok_den & is.finite(R$b_all) & is.finite(R$b_cc), , drop = FALSE]
if (nrow(nn)) {
  for (i in seq_len(nrow(nn)))
    cat(sprintf("  %-22s IRR %.2f -> %.2f  (차이 %+.3f · 기저 P = %.3f)\n",
                nn$transition[i], nn$irr_all[i], nn$irr_cc[i],
                nn$irr_cc[i] - nn$irr_all[i], nn$p_all[i]))
  cat("  이 전이들은 애초에 귀무이므로 백분율 변화가 의미를 갖지 않습니다.\n")
} else cat("  없습니다.\n")

cat("\n=== 결론 ====================================================\n")
big <- any(ok) && max(abs(R$chg[ok])) >= 10
if (!big && flip == 0 && sigc == 0) {
  cat("  무응답 처리 방식이 결과를 만들지 않았습니다.\n")
  cat("  Discussion 에 쓸 문장:\n")
  cat("   \"Results were essentially unchanged in a complete-case analysis\n")
  cat("    restricted to participants with complete socioeconomic data\n")
  cat("    (n excluded = ", max(R$n_lost), " intervals).\"\n", sep = "")
} else {
  cat("  ★ 완전자료에서 결과가 움직입니다. 아래를 순서대로 확인하십시오.\n")
  cat("   1) [1] 표에서 어느 전이인지, 사건이 몇 건 빠졌는지\n")
  cat("   2) 사건이 적은 전이라면 소표본 잡음입니다 — 각주로 처리 가능\n")
  cat("   3) 사건이 많은 전이인데 움직인다면 무응답이 정보성입니다\n")
  cat("      -> 그때는 소득을 CORE 에서 빼거나 다중대체를 검토해야 합니다\n")
}
cat("=============================================================\n")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  ad <- function(nm, x) { openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, x)
                          openxlsx::setColWidths(wb, nm, 1:ncol(x), "auto") }
  ad("complete_case", TAB); ad("raw", R)
  openxlsx::saveWorkbook(wb, file.path(OUT_DIR, paste0(OUT_STEM, ".xlsx")), overwrite = TRUE)
  cat(sprintf("\n저장: %s\n", file.path(OUT_DIR, paste0(OUT_STEM, ".xlsx"))))
} else {
  utils::write.csv(TAB, file.path(OUT_DIR, paste0(OUT_STEM, ".csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n저장: CSV\n")
}
cat("\n=== 완료 ===\n")
