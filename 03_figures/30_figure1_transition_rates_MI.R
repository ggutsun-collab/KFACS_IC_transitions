###############################################################################
## KF_Fig1_TransitionRates_v260811.R
## Figure 1 (MI, within-sex) — Adjusted state-transition rates by within-sex IC tertile
## v260904: 30_figure1_transition_rates.R 의 데이터 블록만 m=20 MI × within-sex 로 교체
##
## ── v260811 에서 바뀐 것: "숨긴 칸이 하나도 없게" ───────────────────────
## 이전 판에는 그림에 나타나지 않는 정보가 세 종류 있었습니다.
##
##   ① 사건 0건인 칸       -> 아예 그려지지 않아 독자가 존재 자체를 모릅니다
##   ② 희소 칸의 신뢰구간  -> 축이 망가진다는 이유로 감췄습니다
##   ③ 보정율 추정 실패 칸 -> 역시 그려지지 않았습니다
##
## 이제 셋 다 그립니다.
##
##   ① 사건 0건 : '3의 법칙'(rule of three) 상한을 계산해 그 지점에 아래를
##      향한 빈 삼각형을 찍습니다. 삼각형이 아래를 가리키는 것이 "참값은 이
##      아래 어딘가" 라는 뜻이며, 0 을 로그축에 찍을 수 없는 문제를 우회합니다.
##      0/T 의 단측 97.5% 상한 = 3/T. 1,000 인년당으로 환산해 3000/pyears.
##      ★ 인년이 매우 적은 칸은 상한이 관측된 다른 율보다 높게 나옵니다.
##        그것이 정확한 표현입니다 — 그 칸에 대해 우리는 아무것도 모릅니다.
##
##   ② 희소 칸 : 신뢰구간을 점선으로 가늘게 그립니다. 감추지 않습니다.
##      구간이 크게 벌어진다는 사실 자체가 "믿지 말라"는 가장 정직한 신호입니다.
##      y축 범위는 이 구간들까지 포함해 다시 계산합니다.
##
##   ③ 보정율 실패 : 조율(crude rate)을 x 기호로 찍어 위치만 알려 줍니다.
##
## ── 자주 오해되는 점 ─────────────────────────────────────────────────────
## 속 빈 기호는 "통계적으로 유의하지 않다" 는 뜻이 아닙니다.
## 이 그림은 가설검정이 아니라 발생률과 그 정밀도를 보여 줍니다.
## 유의성 검정은 Table 2 의 IRR 이 합니다. 속 빈 기호는 오직
## "사건 5건 미만 또는 인년 50 미만" 이라는 자료량의 문제를 뜻합니다.
##
## 실행:  source("KF_Fig1_TransitionRates_v260811.R", encoding = "UTF-8")
###############################################################################

FIG1_SCRIPT_VERSION <- "v260904 MI within-sex (모든 칸 표시)"
message("\n=== KF_Fig1_TransitionRates ", FIG1_SCRIPT_VERSION, " ===")

## ── 헬퍼 적재 ──────────────────────────────────────────────────────────
## 위치 탐색과 헬퍼 로드는 R/_bootstrap.R 한 곳에서만 정의합니다.
## 모든 경로(PROJECT_ROOT / IMP_DIR / FIG_DIR)는 R/00_setup.R 에서만 정합니다.
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

###############################################################################
## 자가 교정 — KF_common 이 구버전이어도 공통 기준 프로파일을 씁니다
###############################################################################
.KF_PATCHED <- character(0)
if (!exists("make_ref", mode = "function")) {
  make_ref <- function(iv, adj = CFG$ADJ) {
    r <- list()
    for (a in intersect(adj, names(iv))) {
      v <- iv[[a]]
      r[[a]] <- if (is.numeric(v)) mean(v, na.rm = TRUE)
                else names(sort(table(v), decreasing = TRUE))[1]
    }
    r
  }
  .KF_PATCHED <- c(.KF_PATCHED, "make_ref() 정의")
}
if (!exists(".robust_ct", mode = "function")) {
  .robust_ct <- function(m, cluster) {
    has <- requireNamespace("sandwich", quietly = TRUE) &&
           requireNamespace("lmtest", quietly = TRUE)
    vc <- if (has) try(sandwich::vcovCL(m, cluster = cluster, type = "HC0"), silent = TRUE) else NULL
    if (is.null(vc) || inherits(vc, "try-error")) vc <- stats::vcov(m)
    if (has) { ct <- try(lmtest::coeftest(m, vcov. = vc), silent = TRUE)
               if (!inherits(ct, "try-error")) return(ct) }
    b <- stats::coef(m); se <- sqrt(diag(vc))[names(b)]; z <- b / se
    out <- cbind(Estimate = b, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    rownames(out) <- names(b); out
  }
  .KF_PATCHED <- c(.KF_PATCHED, ".robust_ct() 정의")
}
if (!("ref" %in% names(formals(adjusted_rates_by_tertile)))) {
  adjusted_rates_by_tertile <- function(iv, from_lab, to_lab, adj = CFG$ADJ,
                                        min_events = 5, ref = NULL) {
    dsub <- iv[iv$from_lab == from_lab, , drop = FALSE]
    dsub$y <- as.integer(dsub$to_lab == to_lab)
    n_ev <- sum(dsub$y, na.rm = TRUE)
    base <- data.frame(transition = paste(from_lab, "->", to_lab),
                       tertile = CFG$TERT_LABELS, events = NA_integer_,
                       pyears = NA_real_, crude_rate = NA_real_, adj_rate = NA_real_,
                       RR = NA_real_, lcl = NA_real_, ucl = NA_real_,
                       stringsAsFactors = FALSE)
    for (i in seq_along(CFG$TERT_LABELS)) {
      tt <- CFG$TERT_LABELS[i]; dd <- dsub[which(dsub$gLIC_tert == tt), ]
      base$events[i] <- sum(dd$y, na.rm = TRUE)
      base$pyears[i] <- sum(dd$dur, na.rm = TRUE)
      base$crude_rate[i] <- 1000 * base$events[i] / base$pyears[i]
    }
    if (n_ev < min_events || length(unique(stats::na.omit(dsub$gLIC_tert))) < 2) return(base)
    dsub$gLIC_tert <- stats::relevel(factor(dsub$gLIC_tert, levels = CFG$TERT_LABELS), ref = "T1")
    rhs <- c("gLIC_tert", adj)
    rhs <- rhs[vapply(rhs, function(x) {
      v <- all.vars(stats::as.formula(paste("~", x)))[1]
      v %in% names(dsub) && length(unique(stats::na.omit(dsub[[v]]))) > 1 }, logical(1))]
    f <- stats::as.formula(paste0("y ~ ", paste(rhs, collapse = " + "), " + offset(log(dur))"))
    m <- try(stats::glm(f, family = stats::poisson(), data = dsub), silent = TRUE)
    if (inherits(m, "try-error") || !m$converged) return(base)
    ct <- .robust_ct(m, dsub$id)
    for (i in seq_along(CFG$TERT_LABELS)) {
      tt <- CFG$TERT_LABELS[i]; nm <- paste0("gLIC_tert", tt)
      if (tt == "T1") base$RR[i] <- 1
      else if (nm %in% rownames(ct)) {
        b <- ct[nm, "Estimate"]; s <- ct[nm, "Std. Error"]
        if (is.finite(b) && is.finite(s) && s <= 3 && abs(b) <= 10) {
          base$RR[i] <- exp(b); base$lcl[i] <- exp(b - 1.96 * s); base$ucl[i] <- exp(b + 1.96 * s)
        }
      }
    }
    nd <- dsub[rep(1, 3), , drop = FALSE]
    nd$gLIC_tert <- factor(CFG$TERT_LABELS, levels = CFG$TERT_LABELS)
    for (a in intersect(adj, names(nd))) {
      if (!is.null(ref) && !is.null(ref[[a]]))
        nd[[a]] <- if (is.numeric(dsub[[a]])) as.numeric(ref[[a]])
                   else factor(as.character(ref[[a]]), levels = levels(factor(dsub[[a]])))
      else if (is.numeric(dsub[[a]])) nd[[a]] <- mean(dsub[[a]], na.rm = TRUE)
      else nd[[a]] <- names(sort(table(dsub[[a]]), decreasing = TRUE))[1]
    }
    nd$dur <- 1
    pr <- try(stats::predict(m, newdata = nd, type = "response"), silent = TRUE)
    if (!inherits(pr, "try-error")) base$adj_rate <- 1000 * as.numeric(pr)
    bad <- which(base$tertile != "T1" & (is.na(base$RR) | base$events == 0))
    if (length(bad)) base$adj_rate[bad] <- NA_real_
    base
  }
  .KF_PATCHED <- c(.KF_PATCHED, "adjusted_rates_by_tertile() 를 ref= 지원판으로 교체")
}
if (length(.KF_PATCHED)) {
  message("\n[Fig1] KF_common 자동 교정:")
  for (s in .KF_PATCHED) message("       - ", s)
} else message("[Fig1] KF_common 최신판 확인")

###############################################################################
## 설정
###############################################################################
MIN_EV  <- 5      # 사건이 이보다 적은 칸 -> 속 빈 기호
MIN_PY  <- 50     # 인년이 이보다 적은 칸 -> 사건수와 무관하게 속 빈 기호
FIG_H   <- 168    # 그림 높이(mm). 범례 2줄 기준
Y_SCALE <- "auto" # "auto" | "log" | "linear" | "free"

###############################################################################
## 데이터 — m=20 완성 데이터 × within-sex 척도, 칸별 Rubin 풀링   (v260904 MI 판)
##   · 세트마다 14_ 와 동일한 within-sex z·삼분위를 적용하고 세트별 보정율을 구한 뒤
##     log(rate) 척도에서 Rubin 규칙으로 합칩니다 (within = 1/events, between = 세트 분산).
##   · events = 세트 중앙값(반올림), person-years = 세트 평균.
##   · within-sex 척도에서는 성별이 표준화에 흡수되므로 보정군에서 sex_f 를 뺍니다.
###############################################################################
MI_STEM <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
N_SETS  <- NA          # NA = 전체 세트. 빠른 점검은 3.
MIN_M_OK <- 0.5      # v260904b: 세트의 절반 이상에서 추정되면 풀링 (x 표시 최소화)
ADJ_WS  <- setdiff(CFG$ADJ, "sex_f")
adj_phrase_orig <- adj_phrase
adj_phrase <- function(adj = ADJ_WS) adj_phrase_orig(adj)

apply_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]
  z <- d$gLIC; tt <- rep(NA_character_, nrow(d))
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    mu <- mean(w1$gLIC[k1], na.rm = TRUE); sdv <- stats::sd(w1$gLIC[k1], na.rm = TRUE)
    z[k] <- (d$gLIC[k] - mu) / sdv
    ct <- stats::quantile(w1$gLIC[k1], c(1/3, 2/3), na.rm = TRUE)
    tt[k] <- as.character(cut(d$gLIC[k], c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE))
  }
  d$gLIC_z <- z; d$gLIC_tert_W1cut <- factor(tt, levels = CFG$TERT_LABELS); d
}

MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds")))
SETS <- names(MI); if (is.finite(N_SETS)) SETS <- SETS[seq_len(min(N_SETS, length(SETS)))]
M <- length(SETS); rm(MI); invisible(gc())
msg(sprintf("Figure 1 (MI): 세트 %d개, within-sex 척도, 보정 = %s", M, paste(ADJ_WS, collapse = "+")))

collect_one <- function(d, state_var, sys_label) {
  iv <- build_intervals(d, state_var, "rolling")
  if (is.null(iv)) return(NULL)
  REF <- make_ref(iv, adj = ADJ_WS)
  tl <- default_transitions(state_var, iv)
  ## v260904e: 삼분위 칸별 적합 대신 per-s.d. 연속 모형(Table 2 의 주 모형)에서
  ##   IC = -1 / 0 / +1 s.d. 의 보정 전이율을 예측합니다. 칸이 비어도 모형 기반 추정치가
  ##   모든 위치에 나오며, 자료 지지(그 삼분위 칸의 사건·인년)는 기호로 표시합니다.
  Z_AT <- c(T1 = -1, T2 = 0, T3 = 1)
  rates_by_z <- function(iv, from_lab, to_lab) {
    dsub <- iv[iv$from_lab == from_lab, , drop = FALSE]
    dsub$y <- as.integer(dsub$to_lab == to_lab)
    base <- data.frame(transition = paste(from_lab, "->", to_lab), tertile = names(Z_AT), z = unname(Z_AT),
                       events = NA_integer_, pyears = NA_real_, crude_rate = NA_real_,
                       adj_rate = NA_real_, se_log = NA_real_, RR = NA_real_, lcl = NA_real_, ucl = NA_real_,
                       stringsAsFactors = FALSE)
    for (i in seq_along(Z_AT)) { dd <- dsub[which(dsub$gLIC_tert == names(Z_AT)[i]), ]
      base$events[i] <- sum(dd$y, na.rm = TRUE); base$pyears[i] <- sum(dd$dur, na.rm = TRUE)
      base$crude_rate[i] <- 1000 * base$events[i] / base$pyears[i] }
    n_ev <- sum(dsub$y, na.rm = TRUE)
    if (n_ev < 5) return(base)
    adj <- ADJ_WS[vapply(ADJ_WS, function(v) v %in% names(dsub) && length(unique(stats::na.omit(dsub[[v]]))) > 1, logical(1))]
    fit_one <- function(adj) {
      f <- stats::as.formula(paste0("y ~ gLIC_z", if (length(adj)) paste0(" + ", paste(adj, collapse = " + ")) else "", " + offset(log(dur))"))
      m <- try(stats::glm(f, family = stats::poisson(), data = dsub), silent = TRUE)
      if (inherits(m, "try-error") || !m$converged || any(is.na(stats::coef(m)))) return(NULL)
      V <- if (requireNamespace("sandwich", quietly = TRUE)) try(sandwich::vcovCL(m, cluster = dsub$id, type = "HC0"), silent = TRUE) else NULL
      if (is.null(V) || inherits(V, "try-error") || any(!is.finite(V))) V <- stats::vcov(m)
      nd <- dsub[rep(1, length(Z_AT)), , drop = FALSE]; nd$gLIC_z <- unname(Z_AT); nd$dur <- 1
      for (a in intersect(adj, names(nd))) {
        nd[[a]] <- if (!is.null(REF[[a]])) { if (is.numeric(dsub[[a]])) as.numeric(REF[[a]]) else factor(as.character(REF[[a]]), levels = levels(factor(dsub[[a]]))) }
                   else if (is.numeric(dsub[[a]])) mean(dsub[[a]], na.rm = TRUE) else names(sort(table(dsub[[a]]), decreasing = TRUE))[1] }
      X <- try(stats::model.matrix(stats::delete.response(stats::terms(m)), nd), silent = TRUE)
      if (inherits(X, "try-error") || ncol(X) != length(stats::coef(m))) return(NULL)
      eta <- as.numeric(X %*% stats::coef(m)); se <- sqrt(pmax(diag(X %*% V %*% t(X)), 0))
      b <- stats::coef(m)["gLIC_z"]; sb <- sqrt(V["gLIC_z", "gLIC_z"])
      if (!is.finite(b) || !is.finite(sb) || sb > 3) return(NULL)
      list(rate = 1000 * exp(eta), se = se, b = b, sb = sb)
    }
    r <- fit_one(adj); used <- "full"
    if (is.null(r)) { r <- fit_one(intersect(c("age_c", "comorbid_bl"), adj)); used <- "reduced (age, comorbidity)" }
    if (is.null(r)) { r <- fit_one(character(0)); used <- "capacity only" }
    if (is.null(r)) return(base)
    base$adj_rate <- r$rate; base$se_log <- r$se
    base$RR <- exp(r$b * base$z); base$lcl <- exp((r$b - 1.96 * r$sb) * base$z); base$ucl <- exp((r$b + 1.96 * r$sb) * base$z)
    base$adj_used <- used
    base
  }
  r <- do.call(rbind, lapply(tl, function(tr) rates_by_z(iv, tr[1], tr[2])))
  if (!"adj_used" %in% names(r)) r$adj_used <- NA_character_
  p <- strsplit(as.character(r$transition), " -> ", fixed = TRUE)
  r$from <- vapply(p, `[`, character(1), 1); r$to <- vapply(p, `[`, character(1), 2)
  r$system <- sys_label; r
}
PER <- vector("list", M)
for (i in seq_along(SETS)) {
  load_imputed(SETS[i], stem = MI_STEM)
  d <- apply_scale(prep_long(load_kfacs()))
  if (i == 1) msg(sprintf("Figure 1: N=%d (W1), person-waves=%d", length(unique(d$id)), nrow(d)))
  r <- rbind(collect_one(d, "state", "a  ADL disability"),
             collect_one(d, "frailty_3cat", "b  Frailty phenotype"))
  r$set <- SETS[i]; PER[[i]] <- r
  msg(sprintf("[%s] %d/%d 완료", SETS[i], i, M))
}
PER <- do.call(rbind, PER)
save_vals(PER, "F1_transition_rates_MI_perset.csv", FIG_DIR)

## 칸별 풀링
.pool_log <- function(x, se) {                # x, se: log 척도
  ok <- is.finite(x) & is.finite(se); if (sum(ok) < max(2, ceiling(MIN_M_OK * M))) return(c(NA, NA))
  x <- x[ok]; se <- se[ok]; m <- length(x)
  W <- mean(se^2); Bv <- if (m > 1) stats::var(x) else 0
  c(mean(x), sqrt(W + (1 + 1/m) * Bv))
}
key <- unique(PER[, c("system", "transition", "from", "to", "tertile")])
RATE <- do.call(rbind, lapply(seq_len(nrow(key)), function(j) {
  k <- key[j, ]
  s <- PER[PER$system == k$system & PER$transition == k$transition & PER$tertile == k$tertile, ]
  ev <- round(stats::median(s$events, na.rm = TRUE)); py <- mean(s$pyears, na.rm = TRUE)
  pr <- .pool_log(log(s$adj_rate), s$se_log)
  prr <- if (k$tertile == "T1") c(0, 0) else .pool_log(log(s$RR), (log(s$ucl) - log(s$lcl)) / 3.92)
  au <- names(sort(table(s$adj_used), decreasing = TRUE))[1]
  data.frame(transition = k$transition, tertile = k$tertile, events = as.integer(ev), pyears = py, adj_used = au,
             crude_rate = 1000 * ev / py, adj_rate = exp(pr[1]), se_log = pr[2],
             RR = exp(prr[1]), lcl = exp(prr[1] - 1.96 * prr[2]), ucl = exp(prr[1] + 1.96 * prr[2]),
             from = k$from, to = k$to, system = k$system, m_used = sum(is.finite(s$adj_rate)),
             stringsAsFactors = FALSE)
}))
RATE$lo <- RATE$adj_rate * exp(-1.96 * RATE$se_log)
RATE$hi <- RATE$adj_rate * exp(+1.96 * RATE$se_log)

###############################################################################
## ★ 칸 분류 — 네 가지. 어느 것도 그림에서 빠지지 않습니다.
###############################################################################
##  ok     : 사건 >= MIN_EV 이고 인년 >= MIN_PY. 채운 기호 + 실선 CI
##  sparse : 사건이 있지만 위 기준 미달. 속 빈 기호 + 점선 CI
##  zero   : 사건 0건. 3의 법칙 상한에 아래를 향한 빈 삼각형 + 바닥까지 선
##  nonest : 사건은 있으나 보정율 추정 실패. 조율 위치에 x 기호
###############################################################################
RATE$cls <- "ok"
RATE$cls[!is.na(RATE$adj_rate) & (RATE$events < MIN_EV | RATE$pyears < MIN_PY)] <- "sparse"
RATE$cls[!is.na(RATE$adj_rate) & (RATE$events == 0 | is.na(RATE$events))] <- "extrap"
RATE$cls[is.na(RATE$adj_rate) & RATE$events > 0] <- "nonest"
RATE$cls[is.na(RATE$adj_rate) & (RATE$events == 0 | is.na(RATE$events))] <- "zero"

## 사건 0건의 단측 97.5% 상한 (3의 법칙), 1,000 인년당
RATE$rule3 <- ifelse(RATE$cls == "zero" & is.finite(RATE$pyears) & RATE$pyears > 0,
                     3000 / RATE$pyears, NA_real_)
## 그림에 찍을 y 값
RATE$yval <- RATE$adj_rate
RATE$yval[RATE$cls == "zero"]   <- RATE$rule3[RATE$cls == "zero"]
RATE$yval[RATE$cls == "nonest"] <- RATE$crude_rate[RATE$cls == "nonest"]

## 축 순서
SEV <- c(Normal = 1, Robust = 1, Mild = 2, `Pre-frail` = 2, Severe = 3, Frail = 3, Death = 4)
RATE$from_rank <- SEV[as.character(RATE$from)]
RATE$to_rank   <- SEV[as.character(RATE$to)]
lev_all <- unique(c(CFG$STATE_LAB, CFG$FRAILTY_LAB))
RATE$from <- factor(RATE$from, levels = lev_all)
RATE$to   <- factor(RATE$to,   levels = lev_all)
LAB_TERTILE <- c(T1 = "Low IC (-1 s.d.)", T2 = "Average IC (0)", T3 = "High IC (+1 s.d.)")
PAL_TERTILE_LAB <- stats::setNames(unname(PAL_TERTILE[c("T1","T2","T3")]), unname(LAB_TERTILE))
RATE$tertile <- factor(RATE$tertile, levels = names(LAB_TERTILE), labels = unname(LAB_TERTILE))
fac_lev <- paste("Transition from", c("Normal","Mild","Severe","Robust","Pre-frail","Frail"))
RATE$facet <- factor(paste("Transition from", RATE$from), levels = fac_lev)

save_vals(RATE[, c("system","transition","from","to","tertile","adj_used","events","pyears",
                   "crude_rate","adj_rate","lo","hi","RR","cls","rule3","yval")],
          "F1_transition_rates_MI.csv", FIG_DIR)

###############################################################################
## y축 — 이제 희소 칸의 CI 와 3의 법칙 상한까지 포함해 계산합니다
###############################################################################
.pts <- RATE$yval[is.finite(RATE$yval) & RATE$yval > 0]
.ci  <- c(RATE$lo, RATE$hi)
.ci  <- .ci[is.finite(.ci) & .ci > 0]
ALLV <- c(.pts, .ci)
Y_SPREAD <- max(.pts) / min(.pts)
Y_LOG <- switch(Y_SCALE, log = TRUE, linear = FALSE, free = FALSE, Y_SPREAD > 10)
Y_LIM <- c(0, max(ALLV) * 1.05)
if (Y_LOG) Y_LIM <- c(min(ALLV) * 0.85, max(ALLV) * 1.15)
msg(sprintf("Fig1 y축: 공통 %s | 점 범위 %.1f~%.0f (배율 %.0f) | 축 %.2f~%.0f",
            if (Y_LOG) "로그" else "선형", min(.pts), max(.pts), Y_SPREAD, Y_LIM[1], Y_LIM[2]))

###############################################################################
## 그림
###############################################################################
HAS_AXES <- utils::packageVersion("ggplot2") >= "3.5.0"
fw <- function(sc) {
  a <- list(facets = ~ facet, nrow = 1, scales = sc)
  if (HAS_AXES) a$axes <- "all_y"
  do.call(ggplot2::facet_wrap, a)
}

CLS_LEV <- c("ok", "sparse", "extrap", "zero", "nonest")
## ★ 범례 라벨은 짧게. 긴 문구를 넣었더니 범례가 지면(180 mm)보다 넓어져
##   좌우로 넘치고 패널까지 밀려 잘렸습니다. 자세한 정의는 캡션이 합니다.
CLS_LAB <- c(ok     = "Adequate data",
             sparse = paste0("<", MIN_EV, " events or <", MIN_PY, " person-years"),
             extrap = "Model-based (no events in this tertile)",
             zero   = "No events (upper 97.5% limit)",
             nonest = "Crude rate (adjusted not estimable)")
CLS_SHP <- c(ok = 21, sparse = 21, extrap = 23, zero = 6, nonest = 4)
## ★ v260811b — 선 종류를 범례에 넣습니다.
##  이전에는 점선의 뜻이 캡션 본문에만 있었습니다. 범례에는 기호 세 개만
##  있어서, 그림만 보는 독자는 점선을 '유의하지 않음' 으로 읽습니다.
##  실제로 저자 본인이 두 번 물었습니다 — 심사자는 반드시 묻습니다.
##  linetype 을 cls 에 매핑하고, shape 척도와 제목·라벨을 완전히 같게 두어
##  ggplot 이 두 범례를 하나로 합치게 합니다. 그러면 키 하나에
##  '기호 + 그 기호가 쓰는 선' 이 같이 그려집니다.
CLS_LTY <- c(ok = "solid", sparse = "22", extrap = "12", zero = "blank", nonest = "solid")
LEG_LINE <- TRUE   # FALSE 로 두면 예전처럼 선 종류를 범례에서 뺍니다
## ★ 자료에 실제로 존재하는 분류만 범례에 넣습니다. 없는 기호를 범례에 남기면
##   독자가 "그런 칸이 있는데 안 보인다"고 오해합니다. 두 패널이 같은 범례를
##   쓰도록 패널별이 아니라 RATE 전체에서 계산합니다.
CLS_USE <- CLS_LEV[CLS_LEV %in% unique(RATE$cls)]
if (!length(CLS_USE)) CLS_USE <- "ok"
FILL_OV <- ifelse(CLS_USE == "ok", "black", NA)
LEG_ROW <- 1   # 세로로 쌓으므로 각 범례는 한 줄
msg(sprintf("[Fig1] 범례에 표시할 분류 %d개: %s", length(CLS_USE),
            paste(CLS_USE, collapse = ", ")))

## leg = TRUE 인 패널에서만 범례를 만듭니다.
## ★ 두 패널이 각각 범례를 만들면 patchwork 가 둘을 합치지 못하고 세 줄로
##   쌓는 일이 있습니다(실제로 'Adequate data' 줄이 두 번 나왔습니다).
##   한쪽에서만 만들면 확실합니다.
mk <- function(sys, leg = FALSE) {
  x <- RATE[RATE$system == sys & is.finite(RATE$yval), ]
  x$facet <- droplevels(x$facet)
  ## ★ x축 순서 (v260811 수정)
  ##  이전에는 tapply+sort 로 순위를 뽑았는데, 사건 0건 칸이 있는 패널에서
  ##  순서가 어긋났습니다(Severe 패널이 Death, Normal, Mild 로 나왔습니다).
  ##  SEV 에서 직접 순서를 만들어 factor levels 로 넣습니다. 스케일에
  ##  limits 를 주면 free_x 가 무력화되므로 levels 만 씁니다.
  lv <- intersect(names(sort(SEV)), unique(as.character(x$to)))
  x$to <- factor(as.character(x$to), levels = lv)
  x$cls <- factor(x$cls, levels = CLS_USE)

  ## CI: ok 는 실선, sparse 는 점선. 없는 칸은 접어서 감춥니다.
  x$lo_p <- ifelse(is.finite(x$lo), x$lo, x$yval)
  x$hi_p <- ifelse(is.finite(x$hi), x$hi, x$yval)
  x$show <- factor(is.finite(x$lo) & x$cls %in% c("ok", "sparse"),
                   levels = c("FALSE", "TRUE"))
  ## 0건 칸: 상한에서 축 바닥까지 내려가는 선
  ## ★ 그 패널에 0건 칸이 없으면 z 가 0행이 됩니다. 길이 1 을 0행에 대입하면
  ##   "replacement has 1 row, data has 0" 오류가 납니다. rep() 로 길이를 맞춥니다.
  ###########################################################################
  ## ★ 0건 칸의 '점선 꼬리' 레이어를 제거했습니다 (v260811 최종)
  ##  이 레이어 하나가 세 가지 문제를 동시에 만들고 있었습니다.
  ##
  ##  (1) x축 순서가 뒤바뀜
  ##      facet 의 scales="free_x" 에서 패널의 이산 축 순서는 factor 의
  ##      levels 가 아니라 '레이어를 훑으며 값이 처음 나타난 순서' 로 정해집니다.
  ##      꼬리 레이어가 첫 레이어였고 그 안에는 0건 칸의 도착 상태 하나만
  ##      들어 있었으므로, 그 상태가 축 맨 앞으로 나왔습니다.
  ##      (Mild 패널 -> Severe 가 앞, Severe 패널 -> Death 가 앞)
  ##
  ##  (2) 꼬리와 삼각형의 가로 위치가 어긋남
  ##      position_dodge 는 '그 레이어에 실제로 있는 group 수' 로 슬롯을
  ##      나눕니다. 꼬리 레이어에는 삼분위가 하나뿐이라 가운데로 왔고,
  ##      점 레이어(3군)의 삼각형은 오른쪽 슬롯에 놓였습니다.
  ##      (예전에 CI 와 점이 어긋났던 것과 정확히 같은 원인입니다)
  ##
  ##  (3) 이웃한 삼분위의 신뢰구간으로 오독됨
  ##      Mild -> Severe 에서는 T2 의 CI 아래에 붙은 점선처럼,
  ##      Severe -> Death 에서는 T2 위로 뻗은 파란 점선처럼 보였습니다.
  ##
  ##  삼각형만으로 상한을 표시하고, "그 아래 어딘가" 라는 뜻은 캡션이
  ##  정확하게 전달합니다. 정보 손실은 없습니다.
  ###########################################################################
  ## 채우기는 ok 에만
  x$fillv <- ifelse(x$cls == "ok", as.character(x$tertile), NA_character_)

  pd <- position_dodge(width = 0.65, preserve = "single")

  ## ★ 범례 키에 기호와 선을 같이 그리기 위한 공통 guide 객체.
  ##  shape 와 linetype 두 척도가 '제목 없음 + 같은 라벨' 이면 ggplot 이
  ##  둘을 한 범례로 합치고, 키 안에 선을 깔고 그 위에 기호를 얹습니다.
  ##  override.aes 의 colour="black" 은 키를 삼분위 색이 아닌 중립색으로
  ##  만들고, linewidth 는 키 안의 선 굵기입니다.
  gl_cls <- if (leg) guide_legend(order = 1, nrow = LEG_ROW,
                                  override.aes = list(colour = "black", fill = FILL_OV,
                                                      size = 1.5, linewidth = 0.40))
            else "none"

  p <- ggplot(x, aes(to, yval, colour = tertile, group = tertile)) +
    ## 신뢰구간 — 선 종류가 곧 칸의 분류입니다 (ok 실선 / sparse 점선 / zero 없음)
    geom_errorbar(aes(ymin = lo_p, ymax = hi_p, alpha = show, linetype = cls),
                  width = 0.38, position = pd, linewidth = 0.40, na.rm = TRUE,
                  show.legend = if (LEG_LINE)
                    c(linetype = TRUE, colour = FALSE, alpha = FALSE) else FALSE) +
    scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1), guide = "none", drop = FALSE) +
    scale_linetype_manual(values = unname(CLS_LTY[CLS_USE]), limits = CLS_USE,
                          labels = unname(CLS_LAB[CLS_USE]), name = NULL, drop = FALSE) +
    geom_point(aes(shape = cls, fill = fillv), position = pd,
               size = 1.15, stroke = 0.35, na.rm = TRUE) +
    scale_colour_manual(values = PAL_TERTILE_LAB, name = "Intrinsic capacity (within-sex s.d.)", drop = FALSE) +
    scale_fill_manual(values = PAL_TERTILE_LAB, guide = "none",
                      na.value = NA, drop = FALSE) +
    scale_shape_manual(values = CLS_SHP[CLS_USE], limits = CLS_USE,
                       labels = unname(CLS_LAB[CLS_USE]), name = NULL, drop = FALSE) +
    guides(shape    = gl_cls,
           linetype = if (LEG_LINE) gl_cls else "none",
           colour   = if (leg) guide_legend(order = 2,
                        override.aes = list(shape = 21, fill = unname(PAL_TERTILE_LAB)))
                      else "none") +
    ## ★ scale_x_discrete(limits=) 는 쓰지 않습니다.
    ##  limits 를 주면 facet 의 scales="free_x" 를 덮어써서 모든 패널에
    ##  체계 전체의 상태가 다 나타납니다(출발 상태 자신까지 빈 칸으로).
    ##  factor 의 levels 만 올바르면 free_x 가 패널마다 없는 수준을 떨어뜨리고
    ##  순서는 그대로 유지합니다.
    fw(if (identical(Y_SCALE, "free")) "free" else "free_x") +
    { if (Y_LOG) scale_y_continuous(trans = "log10", limits = Y_LIM,
                                    breaks = c(1, 3, 10, 30, 100, 300, 1000),
                                    labels = c("1","3","10","30","100","300","1,000"))
      else if (identical(Y_SCALE, "free")) scale_y_continuous()
      else scale_y_continuous(limits = Y_LIM) } +
    labs(x = "State at end of transition",
         y = "Transitions per 1,000 person-years",
         title = sys) +
    theme_na() +
    theme(panel.spacing.x = unit(3.2, "mm"),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "plain"))
  p
}

pa <- mk("a  ADL disability",   leg = TRUE)
pb <- mk("b  Frailty phenotype", leg = FALSE)

###############################################################################
## 캡션 — 표준이 아닌 칸을 하나도 빠짐없이 열거합니다
###############################################################################
en <- function(sel, f) {
  s <- RATE[sel, ]
  if (!nrow(s)) return("")
  paste(unique(vapply(seq_len(nrow(s)), function(i) f(s[i, ]), "")), collapse = "; ")
}
sp_txt <- en(RATE$cls == "sparse", function(r)
  sprintf("%s, %s: %d events in %.0f person-years", r$transition, r$tertile, r$events, r$pyears))
z_txt <- en(RATE$cls == "zero", function(r)
  sprintf("%s, %s: 0 events in %.0f person-years, upper limit %.0f", r$transition, r$tertile,
          r$pyears, r$rule3))
ne_txt <- en(RATE$cls == "nonest", function(r)
  sprintf("%s, %s: %d events in %.0f person-years", r$transition, r$tertile, r$events, r$pyears))

CAP <- paste0(
  "Poisson models with a log(person-time) offset, adjusted for ", adj_phrase(),
  "; standard errors clustered by participant. Capacity enters each model as a continuous ",
  "within-sex standardized score measured at the start of the interval (rolling landmark), and rates are ",
  "predicted at -1, 0 and +1 s.d.; estimates were obtained in each of 20 multiply imputed datasets and ",
  "combined by Rubin's rules on the log-rate scale. Symbols reflect the data supporting each position: ",
  "event counts (medians across datasets) and person-years (means) in the corresponding within-sex ",
  "tertile of the origin state. Rates are predicted at a single reference profile common to every transition, so ",
  "that panels can be compared with one another as well as within themselves; participants ",
  "occupying the more impaired states are older, and predicting each panel at its own mean ",
  "would otherwise make those panels appear higher for reasons of age rather than of state. ",
  "All panels share one vertical axis",
  if (Y_LOG) paste0("; the axis is logarithmic because the plotted values span a factor of about ",
    sprintf("%.0f", Y_SPREAD), ", and on a linear axis the confidence intervals of the least ",
    "frequent transitions would be shorter than the plotting symbol. ") else ". ",
  "The horizontal axis differs between panels because the destination states differ. ",
  "Vertical bars are 95% confidence intervals computed from the number of events. ",
  "Every tertile cell of every transition is shown, using four symbols. ",
  "Filled circles with solid intervals are positions whose tertile contributed at least ", MIN_EV, " events and at least ",
  MIN_PY, " person-years. ",
  if (any(RATE$cls == "extrap")) paste0("Open diamonds with dash-dot intervals are model-based predictions at positions whose tertile ",
    "contributed no events (", en(RATE$cls == "extrap", function(r) sprintf("%s, %s", r$transition, r$tertile)),
    "); these extrapolate the fitted capacity gradient and should be read as such. ") else "",
  if (nzchar(sp_txt)) paste0(
    "Open circles with dotted intervals are cells that fall below one or both of those ",
    "thresholds; both criteria are applied because a rate estimated from a handful of events ",
    "over a very short exposure time is unstable however many events it happens to contain. ",
    "The intervals of these cells are shown rather than omitted, because their width is the ",
    "most direct statement of how little the cell supports. They are not interpreted (",
    sp_txt, "). ") else "",
  if (nzchar(z_txt)) paste0(
    "Cells in which no transition was observed cannot be placed on a rate axis, and are shown ",
    "as a downward-pointing open triangle at the one-sided 97.5% upper limit obtained from the ",
    "rule of three (three events divided by the accrued person-time); the triangle points downwards ",
    "because the true rate lies anywhere below it. Where the accrued person-time ",
    "is very small this limit lies above the rates observed elsewhere, which is the correct ",
    "representation: the cell is uninformative rather than reassuring (", z_txt, "). ") else "",
  if (nzchar(ne_txt)) paste0(
    "Crosses mark cells in which events occurred but the adjusted rate could not be estimated; ",
    "the crude rate is plotted in its place (", ne_txt, "). ") else "",
  local({ z <- unique(RATE[RATE$adj_used != "full", c("transition", "adj_used")])
    if (nrow(z)) paste0("Where the fully adjusted model could not be estimated because a covariate level or a tertile ",
      "contributed no events, the adjustment was reduced for that transition only (",
      paste(sprintf("%s: %s", z$transition, z$adj_used), collapse = "; "), "). ") else "" }),
  "Symbol shape reflects the amount of data supporting a cell and carries no information about ",
  "statistical significance; hypothesis tests for the effect of intrinsic capacity are reported ",
  "as incidence rate ratios in Table 2.")

## ★ 범례를 가로로 한 줄에 몰면 폭이 넘칩니다. 세로로 쌓습니다.
fig <- (pa / pb) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.margin = margin(t = 1, b = 0, unit = "mm"),
        legend.spacing.y = unit(0.6, "mm"),
        legend.key.size = unit(3.4, "mm"),
        legend.text = element_text(size = 5.8),
        legend.title = element_text(size = 5.8))
.af <- names(formals(annot_na))
if (all(c("title", "caption") %in% .af)) {
  fig <- fig + annot_na(width_mm = NA_WIDTH$double,
    title = "Adjusted state-transition rates by intrinsic-capacity tertile", caption = CAP)
} else {
  fig <- fig + annot_na(width_mm = NA_WIDTH$double)
  writeLines(c("Figure 1 | Adjusted state-transition rates by intrinsic-capacity tertile", "", CAP),
             file.path(FIG_DIR, "Figure1_TransitionRates_MI_caption.txt"), useBytes = TRUE)
  message("[Fig1] 캡션을 직접 저장했습니다.")
}
save_na(fig, "Figure1_TransitionRates_MI", width_mm = NA_WIDTH$double,
        height_mm = FIG_H, dir = FIG_DIR)

###############################################################################
## 자체 검증
###############################################################################
tb <- table(factor(RATE$cls, levels = CLS_LEV))
cat("\n=== Figure 1 자체 검증 =====================================\n")
cat(sprintf("  전체 칸 %d개 = 표준 %d + 희소 %d + 사건0 %d + 추정실패 %d\n",
            nrow(RATE), tb["ok"], tb["sparse"], tb["zero"], tb["nonest"]))
cat(sprintf("  그림에 찍힌 칸 %d개  (누락 %d개)\n",
            sum(is.finite(RATE$yval)), sum(!is.finite(RATE$yval))))
if (any(!is.finite(RATE$yval))) {
  cat("  ★ 아래 칸이 여전히 그려지지 않습니다 — 인년이 0이라 3의 법칙 상한도 없습니다.\n")
  print(RATE[!is.finite(RATE$yval), c("transition","tertile","events","pyears")], row.names = FALSE)
}
cat(sprintf("  y축: 공통 %s, 배율 %.0f, 범위 %.2f~%.0f\n",
            if (Y_LOG) "로그" else "선형", Y_SPREAD, Y_LIM[1], Y_LIM[2]))
if (tb["zero"] > 0) {
  cat("  사건 0건 칸의 3의 법칙 상한:\n")
  zz <- RATE[RATE$cls == "zero", c("transition","tertile","pyears","rule3")]
  zz$rule3 <- round(zz$rule3, 1); zz$pyears <- round(zz$pyears, 1)
  print(zz, row.names = FALSE)
}
cat("------------------------------------------------------------\n")
cat(sprintf("  범례에 %d개 분류가 표시됩니다: %s\n", length(CLS_USE), paste(CLS_USE, collapse = ", ")))
cat("  자료에 없는 분류는 범례에서도 빼므로, 위 개수와 그림의 범례가 일치해야 합니다.\n")

## ★ 범례 병합 점검 (v260811b)
##  shape 와 linetype 이 한 범례로 합쳐졌는지 셉니다. 기대값은 2줄입니다
##  (분류 범례 1 + 삼분위 범례 1). 3이 나오면 병합이 실패해 예전처럼
##  범례가 한 줄 더 생긴 것이므로, 위쪽 LEG_LINE <- FALSE 로 되돌리십시오.
if (isTRUE(LEG_LINE)) {
  ng <- try({
    g  <- ggplot2::ggplotGrob(pa)
    bx <- g$grobs[grep("guide-box", g$layout$name)]
    n  <- 0L
    for (bb in bx) if (!is.null(bb$layout))
      n <- n + sum(grepl("^guides", bb$layout$name))
    n
  }, silent = TRUE)
  if (inherits(ng, "try-error") || !is.finite(ng) || ng == 0L) {
    cat("  선 종류 범례: 자동 점검 불가 (ggplot2 판이 달라 구조가 다릅니다).\n")
    cat("    PDF 를 열어 범례가 두 줄인지 눈으로 확인하십시오.\n")
  } else if (ng <= 2L) {
    cat(sprintf("  선 종류 범례: 병합 성공 (범례 %d줄). 키 안에 기호와 선이 함께 그려집니다.\n", ng))
  } else {
    cat(sprintf("  ★ 선 종류 범례: 병합 실패 (범례 %d줄). LEG_LINE <- FALSE 로 되돌리십시오.\n", ng))
  }
}
cat("  패널별 x축 순서 (출발 상태 자신은 없어야 합니다):\n")
for (.f in levels(droplevels(RATE$facet))) {
  .st <- unique(as.character(RATE$to[RATE$facet == .f & is.finite(RATE$yval)]))
  cat(sprintf("    %-26s %s\n", .f,
              paste(intersect(names(sort(SEV)), .st), collapse = " -> ")))
}
cat("============================================================\n")
cat("\n=== Figure 1 완료 ===  근거수치:", file.path(FIG_DIR, "F1_transition_rates_MI.csv"), "\n")
