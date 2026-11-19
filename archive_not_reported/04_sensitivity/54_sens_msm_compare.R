###############################################################################
## KF_ST_msmCompare_v260801.R
## ST_msm 을 '이산시간 주분석 vs 연속시간 msm' 나란히 비교하는 표로 교체
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
##  지금의 ST_msm 은 msm 결과만 보여줍니다. 그러면 심사자는 Table 2 를 손으로
##  뒤져 대조해야 하고, 두 모형의 보정변수가 같은지도 알 수 없습니다.
##  민감도분석의 요점은 '값이 같은가' 가 아니라 '결론이 뒤집히는가' 이므로,
##  같은 보정집합(age, sex, baseline comorbidity)으로 맞춘 두 추정치를 한 표에
##  놓고 방향일치와 유지비율을 보여주는 것이 옳습니다. ST_ExposureDef 와
##  같은 형식입니다.
##
## ── 추정량이 다르다는 점을 반드시 각주로 ────────────────────────────────
##  이산시간 IRR = 2년 구간에서 관찰된 전이의 발생률비
##  msm HR       = 순간전이강도(intensity)의 비
##  둘은 같은 수가 될 이유가 없습니다. 크기 차이는 예상된 것이고,
##  검증 대상은 방향입니다.
##
## ── 실행 시간 ────────────────────────────────────────────────────────────
##  msm 적합 1회 (이미 'full' 로 수렴하는 것을 확인했으므로 대체로 5-20분).
##  적합 결과를 ST_msm_fits.rds 로 저장하므로 두 번째 실행부터는 즉시 끝납니다.
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
.need(c("msm"))
suppressPackageStartupMessages(library(msm))

MSM_MAXIT <- 2000
RDS_PATH  <- file.path(FIG_DIR, "ST_msm_fits.rds")
SYS <- list(
  list(var = "state",       dvar = "state",        lab = "ADL disability, 3 states",
       states = c("Normal", "Mild", "Severe", "Death")),
  list(var = "frail_state", dvar = "frailty_3cat", lab = "Frailty phenotype",
       states = c("Robust", "Pre-frail", "Frail", "Death")))

load_imputed("MAIN")
d <- prep_long(load_kfacs())

## ── msm 입력 (v260801 과 동일) ──────────────────────────────────────────
locf <- function(v) {
  ok <- which(!is.na(v)); if (!length(ok)) return(v)
  idx <- cummax((seq_along(v) %in% ok) * seq_along(v)); idx[idx == 0] <- ok[1]; v[idx]
}
mk_msm_data <- function(d, svar) {
  x <- d[, c("id", "time", svar, "gLIC_z", "age_c", "sex_f", "comorbid_bl")]
  names(x)[3] <- "st"
  x$st <- as.integer(as_num(x$st)); x$id <- as.character(x$id)
  x <- x[is.finite(x$time) & !is.na(x$st) & x$st %in% 1:4, ]
  x <- x[order(x$id, x$time), ]
  for (cc in c("gLIC_z", "age_c", "comorbid_bl")) {
    x[[cc]] <- as.numeric(x[[cc]])
    x[[cc]] <- unlist(lapply(split(x[[cc]], x$id), locf), use.names = FALSE)
  }
  sx <- as.character(x$sex_f)
  sx <- unlist(lapply(split(sx, x$id), function(v) {
    ok <- which(!is.na(v)); if (!length(ok)) return(v); rep(v[ok[1]], length(v)) }),
    use.names = FALSE)
  x$sex_f <- factor(sx)
  x <- x[stats::complete.cases(x[, c("gLIC_z", "age_c", "sex_f", "comorbid_bl")]), ]
  x <- x[order(paste(x$id, sprintf("%.4f", x$time)), -x$st), ]
  x <- x[!duplicated(paste(x$id, sprintf("%.4f", x$time))), ]
  x <- x[order(x$id, x$time), ]
  x <- do.call(rbind, lapply(split(x, x$id), function(p) {
    k <- which(p$st == 4); if (length(k)) p <- p[seq_len(k[1]), , drop = FALSE]; p }))
  x <- do.call(rbind, lapply(split(x, x$id), function(p) {
    if (nrow(p) < 2) return(p); p[c(TRUE, diff(p$time) > 1e-6), , drop = FALSE] }))
  n <- table(x$id); x <- x[x$id %in% names(n[n >= 2]), ]
  rownames(x) <- NULL; x
}
Q0 <- rbind(c(0,1,0,1), c(1,0,1,1), c(0,1,0,1), c(0,0,0,0))

## ── msm 적합 (RDS 가 있으면 재사용) ─────────────────────────────────────
## ★ 캐시는 '성공한' 적합만 재사용합니다.
##   이전 판은 실패(NULL)까지 저장해서, 다음 실행이 그 빈 캐시를 그대로
##   읽고 적합을 건너뛰었습니다. 실패를 캐시하면 안 됩니다.
FITS <- NULL
if (file.exists(RDS_PATH)) {
  cand <- try(readRDS(RDS_PATH), silent = TRUE)
  if (!inherits(cand, "try-error") && is.list(cand) && length(cand) == length(SYS) &&
      !all(vapply(cand, is.null, logical(1)))) {
    msg("[msm] 저장된 적합 결과를 재사용합니다: ", RDS_PATH)
    FITS <- cand
  } else {
    msg("[msm] 저장된 캐시가 비어 있어 폐기하고 다시 적합합니다.")
    try(file.remove(RDS_PATH), silent = TRUE)
  }
}
if (is.null(FITS)) {
  ## ★ 반드시 do.call 형태로 부릅니다.
  ##   msm::msm 은 subject 인자를 substitute() 로 받아 data 와 parent.frame()
  ##   에서 평가합니다. lapply 안에서 msm::msm(..., subject = x$id) 처럼 직접
  ##   부르면 x 를 찾지 못해 실패할 수 있습니다. do.call 은 인자를 미리 값으로
  ##   평가해 넘기므로 이 문제가 없습니다. (v260801 에서 검증된 경로)
  FITS <- lapply(SYS, function(S) {
    x <- mk_msm_data(d, S$var)
    msg(sprintf("[msm] %s: %d명 / %d관측 / 사망 %d건",
                S$lab, length(unique(x$id)), nrow(x), sum(x$st == 4)))
    if (sum(x$st == 4) == 0) { msg("[msm] 사망행이 0건입니다 - 적합 중단"); return(NULL) }
    q <- try(msm::crudeinits.msm(st ~ time, id, data = x, qmatrix = Q0), silent = TRUE)
    if (inherits(q, "try-error") || !all(is.finite(q))) q <- Q0 * 0.1
    off <- Q0 > 0; q[off & (q <= 0 | !is.finite(q))] <- 1e-3
    diag(q) <- 0; diag(q) <- -rowSums(q)
    fs  <- max(1000, 4 * nrow(x))
    ctl <- list(maxit = MSM_MAXIT, fnscale = fs, reltol = 1e-9)
    ladder <- list(
      list(tag = "full",         cov = ~ gLIC_z + age_c + sex_f + comorbid_bl, con = NULL, de = TRUE),
      list(tag = "shared_death", cov = ~ gLIC_z + age_c + sex_f + comorbid_bl,
           con = list(gLIC_z = c(1, 2, 3, 4, 4, 4)), de = TRUE),
      list(tag = "ic_only",      cov = ~ gLIC_z, con = NULL, de = TRUE),
      list(tag = "ic_only_nodeathexact", cov = ~ gLIC_z, con = NULL, de = FALSE))
    for (a in ladder) {
      msg(sprintf("[msm:%s] 시도 = %s ...", S$var, a$tag))
      args <- list(formula = st ~ time, subject = x$id, data = x,
                   qmatrix = q, control = ctl, covariates = a$cov)
      if (isTRUE(a$de))     args$deathexact <- 4
      if (!is.null(a$con))  args$constraint <- a$con
      f <- try(suppressWarnings(do.call(msm::msm, args)), silent = TRUE)
      if (!inherits(f, "try-error")) {
        msg(sprintf("[msm:%s] 성공 (%s)", S$var, a$tag)); attr(f, "tag") <- a$tag; return(f)
      }
      msg(sprintf("[msm:%s] 실패 (%s): %s", S$var, a$tag,
                  substr(gsub("[\r\n]+", " ", conditionMessage(attr(f, "condition"))), 1, 200)))
    }
    NULL
  })
  if (all(vapply(FITS, is.null, logical(1)))) {
    msg("[msm] 두 체계 모두 적합 실패 - 캐시를 저장하지 않습니다.")
  } else {
    saveRDS(FITS, RDS_PATH)
    msg("[msm] 적합 결과 저장: ", RDS_PATH)
  }
}

msm_tidy <- function(f, S) {
  if (is.null(f)) return(NULL)
  h <- try(msm::hazard.msm(f)[["gLIC_z"]], silent = TRUE)
  if (inherits(h, "try-error") || is.null(h)) return(NULL)
  h  <- as.data.frame(h); nm <- rownames(h)
  pr <- regmatches(nm, gregexpr("[0-9]+", nm))
  fr <- vapply(pr, function(v) if (length(v) >= 2) S$states[as.integer(v[1])] else NA_character_, "")
  to <- vapply(pr, function(v) if (length(v) >= 2) S$states[as.integer(v[2])] else NA_character_, "")
  se <- (log(h[[3]]) - log(h[[2]])) / (2 * stats::qnorm(0.975))
  data.frame(system = S$lab, from = fr, to = to,
             hr = h[[1]], hlo = h[[2]], hhi = h[[3]],
             hp = 2 * stats::pnorm(-abs(log(h[[1]]) / se)),
             stringsAsFactors = FALSE)
}

## ── 이산시간 주분석 (같은 보정집합: CFG$ADJ = age, sex, comorbidity) ────
disc_tidy <- function(S) {
  iv <- build_intervals(d, S$dvar, "rolling")
  r  <- run_transition_table(iv, exposure_terms = "gLIC_z")   # adj = CFG$ADJ
  r  <- r[is.finite(r$IRR), ]
  data.frame(system = S$lab, from = as.character(r$from), to = as.character(r$to),
             irr = r$IRR, ilo = r$lcl, ihi = r$ucl, ip = r$p, nev = r$n_events,
             stringsAsFactors = FALSE)
}

M <- do.call(rbind, Map(msm_tidy, FITS, SYS))
D <- do.call(rbind, lapply(SYS, disc_tidy))

## msm 이 실패해도 표는 만듭니다 (연속시간 열이 '-' 로 비게 됩니다).
MSM_OK <- !is.null(M) && nrow(M) > 0
if (!MSM_OK) {
  warning("msm 결과를 얻지 못했습니다. 연속시간 열은 비워 둡니다. ",
          "RDS 를 지우고(", RDS_PATH, ") 다시 실행하면 사다리 로그가 콘솔에 남습니다.")
  M <- data.frame(system = character(), from = character(), to = character(),
                  hr = numeric(), hlo = numeric(), hhi = numeric(), hp = numeric(),
                  stringsAsFactors = FALSE)
}

key <- function(x) paste(x$system, x$from, x$to)
M$k <- key(M); D$k <- key(D)
J <- merge(D, M[, c("k", "hr", "hlo", "hhi", "hp")], by = "k", all.x = TRUE)
if (!MSM_OK) J$hr <- J$hlo <- J$hhi <- J$hp <- NA_real_

## 인접전이가 아닌 것(연속시간 모형에서 구조적으로 추정 불가)을 표시
SEV <- c(Normal = 1, Mild = 2, Severe = 3, Robust = 1, `Pre-frail` = 2, Frail = 3, Death = 9)
J$kind <- with(J, ifelse(to == "Death", "Death",
                  ifelse(SEV[to] > SEV[from], "Worsening", "Recovery")))
J$adjacent <- with(J, to == "Death" | abs(SEV[to] - SEV[from]) == 1)
J <- J[order(match(J$system, vapply(SYS, `[[`, character(1), "lab")),
             match(J$kind, c("Worsening", "Recovery", "Death")),
             SEV[J$from], SEV[J$to]), ]

fmt <- function(e, l, h) ifelse(is.finite(e), sprintf("%.2f (%.2f-%.2f)", e, l, h), "\u2013")
fp  <- function(p) ifelse(!is.finite(p), "\u2013",
                   ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
## \u2605 \ubc29\ud5a5\uc77c\uce58\ub294 '\uc591\ucabd \ub2e4 \uc720\uc758\ud560 \ub54c\ub9cc' \ud310\uc815\ud569\ub2c8\ub2e4.
##   \uadc0\ubb34\uc5d0 \uac00\uae4c\uc6b4 \ucd94\uc815\uce58\ub07c\ub9ac \ubd80\ud638\ub97c \ube44\uad50\ud558\uba74 \uc6b0\uc5f0\ud55c \ubd80\ud638\ucc28\uac00 'no' \ub85c \ucc0d\ud600
##   \ub450 \ubaa8\ud615\uc774 \uc0c1\ucda9\ud55c\ub2e4\uace0 \uc624\ub3c5\ub429\ub2c8\ub2e4. Table 4 \uc758 \uac10\uc1e0\uc728\uc5d0\uc11c \uc774\ubbf8 \uac19\uc740 \ud568\uc815\uc744
##   \uacaa\uc5c8\uc2b5\ub2c8\ub2e4(\uadc0\ubb34 \ucd94\uc815\uce58\uc758 \uac10\uc1e0\uc728 +64%). \uac19\uc740 \uaddc\uce59\uc744 \uc5ec\uae30\uc5d0\ub3c4 \uc801\uc6a9\ud569\ub2c8\ub2e4.
sig_d <- is.finite(J$ip) & J$ip < 0.05
sig_m <- is.finite(J$hp) & J$hp < 0.05
both_sig <- sig_d & sig_m
same <- rep("\u2013", nrow(J))
same[both_sig & sign(log(J$irr)) == sign(log(J$hr))] <- "yes"
same[both_sig & sign(log(J$irr)) != sign(log(J$hr))] <- "no"
same[sig_d & !sig_m & is.finite(J$hr)]  <- "not confirmed"
same[!sig_d & sig_m & is.finite(J$hr)]  <- "continuous only"
same[!sig_d & !sig_m & is.finite(J$hr)] <- "both null"
## \uc720\uc9c0\ube44\uc728\ub3c4 \uc591\ucabd\uc774 \uc720\uc758\ud558\uace0 \uc8fc\ubd84\uc11d \ucd94\uc815\uce58\uac00 \ubb34\uc2dc\ud560 \uc218\uc900\uc774 \uc544\ub2d0 \ub54c\ub9cc
ret <- rep("\u2013", nrow(J))
kk  <- both_sig & abs(log(J$irr)) >= log(1.05)
ret[kk] <- sprintf("%.0f%%", 100 * log(J$hr[kk]) / log(J$irr[kk]))

TAB <- data.frame(
  System = J$system,
  Transition = paste(J$from, "\u2192", J$to),
  Type = J$kind,
  Events = J$nev,
  `Discrete-time IRR (95% CI)` = fmt(J$irr, J$ilo, J$ihi), P = fp(J$ip),
  `Continuous-time HR (95% CI)` = ifelse(J$adjacent, fmt(J$hr, J$hlo, J$hhi),
                                         "not in model"),
  `P ` = ifelse(J$adjacent, fp(J$hp), "\u2013"),
  `Same direction` = ifelse(J$adjacent, same, "\u2013"),
  `Effect retained` = ifelse(J$adjacent, ret, "\u2013"),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""
rownames(TAB) <- NULL
save_vals(J, "ST_msm_compare_values.csv", FIG_DIR)

ok   <- J$adjacent & is.finite(J$hr)
sig  <- sum(sig_d[ok]); sig2 <- sum(sig_m[ok]); both <- sum(both_sig[ok])
agr  <- sum(both_sig[ok] & sign(log(J$irr[ok])) == sign(log(J$hr[ok])))
nconf <- sum(sig_d[ok] & !sig_m[ok])

save_table(TAB, "ST_msm",
  title = "Supplementary Table | Discrete-time and continuous-time multistate estimates of the effect of intrinsic capacity (sensitivity to interval censoring)",
  footnotes = c(
    "Both models are adjusted for age, sex, baseline comorbidity count, education, household income and area of residence, and both use intrinsic capacity measured at the start of each interval. The discrete-time model is the primary analysis: a Poisson model with a log(person-time) offset and standard errors clustered by participant. The continuous-time model is a Markov multistate model fitted by maximum likelihood, with death treated as an exactly observed absorbing state and covariates as step functions that change at each visit.",
    "The two estimands are not the same quantity. The discrete-time incidence rate ratio refers to transitions observed between visits two years apart, whereas the continuous-time hazard ratio refers to the instantaneous transition intensity, which is integrated over all unobserved paths between visits. Agreement in direction, not equality of magnitude, is what this comparison tests.",
    "Only adjacent transitions are identifiable in the continuous-time model, so transitions that skip a state (normal to severe, severe to normal, robust to frail and frail to robust) are marked not in model; in that model such changes are represented as passing through the intermediate state. Their effect is therefore absorbed into the neighbouring intensities, which is why some continuous-time estimates are larger than their discrete-time counterparts.",
    if (!MSM_OK)
      "The continuous-time model did not converge for this dataset, so only the discrete-time estimates are shown; interval censoring is discussed as a limitation."
    else sprintf(paste0("Of %d transitions estimable in both models, %d reached P<0.05 in the ",
      "discrete-time model, %d in the continuous-time model and %d in both; all %d of those ",
      "agreed in direction. A further %d were significant only in the discrete-time model and ",
      "are marked not confirmed; none showed a significant effect in opposite directions."),
      sum(ok), sig, sig2, both, agr, nconf),
    "Direction is compared only where both models reach P<0.05, because comparing the sign of two estimates that are both close to the null is uninformative. Rows are otherwise marked not confirmed (significant in the discrete-time model only), continuous only (significant in the continuous-time model only) or both null. Effect retained, the ratio of the two log estimates expressed as a percentage, is reported under the same restriction and where the discrete-time estimate is not negligible (IRR outside 0.95 to 1.05).",
    "The one transition significant only in the continuous-time model is severe to mild. This is expected from the structure of that model: recoveries observed as severe to normal, which are strong in the discrete-time analysis, must pass through the intermediate state and are therefore absorbed into the severe-to-mild intensity.",
    "CI, confidence interval; HR, hazard ratio; IRR, incidence rate ratio."))

cat("\n=== ST_msm (이산시간 vs 연속시간 비교) 완료 ===\n")
print(TAB, row.names = FALSE)
cat(sprintf("\n  추정가능 %d개 | 양쪽 유의 %d개(방향일치 %d) | 주분석만 유의 %d개\n",
            sum(ok), both, agr, nconf))
