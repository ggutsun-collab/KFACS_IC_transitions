###############################################################################
## KF_Table3_v260811.R
## Table 3 재작성 — 복합 악화 / 사망 / 복합 회복
##
## ── 이 스크립트가 바꾸는 것 ──────────────────────────────────────────────
## (1) worsening 을 "복합 악화" 로 명시적으로 정의합니다.
##     복합 악화 = ADL 상태가 처음으로 나빠지는 시점. 사망(state 4)은 상태
##     번호가 가장 크므로 자동으로 포함됩니다. 별도 항을 더하지 않습니다.
##     이전 각주는 "worse intrinsic-capacity state" 라고 적혀 있었으나 코드는
##     ADL 상태를 세고 있었습니다. 각주를 코드에 맞춥니다.
##
## (2) 복합 회복 열을 추가합니다.
##     복합 회복 = ADL 상태가 처음으로 좋아지는 시점. 사망은 회복을 영원히
##     막으므로 경쟁위험입니다. 카플란–마이어를 쓰면 회복 확률이 부풀려집니다.
##     Aalen–Johansen 누적발생함수로 5년 절대위험을 내고, 위험비는
##     원인별 Cox 로 냅니다. Fine–Gray 부분분포 위험비는 보조로 붙입니다.
##
## (3) 회복 열의 층을 줄입니다.
##     진단 결과 Frail × T2 = 2명, Frail × T3 = 1명이었습니다. 표본이
##     적어서가 아니라 "쇠약하면서 내재역량 최상위이면서 ADL 장애" 가 거의
##     성립하지 않는 조합이기 때문입니다. Frail 행은 삼분위를 합칩니다.
##
## (4) 쇠약 기준 회복을 보충표로 따로 계산합니다.
##     Pre-frail 층 안에서 T1 220 / T2 250 / T3 133 사건으로 검정력이 충분하고,
##     "내재역량이 쇠약표현형보다 호전을 잘 포착한다" 는 주장의 실제 근거가
##     여기 있습니다. 본표에 넣지 않는 이유는 기준군이 달라지기 때문입니다.
##
## 실행:  source("KF_Table3_v260811.R", encoding = "UTF-8")
###############################################################################

message("\n=== Table 3 재작성 v260811 ===")

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
TAU        <- 5      # 추적 절단 (년)

## ── 칸의 신뢰성 규칙 (v260811c 개정) ────────────────────────────────────────
## 이전 판은 "사건 10건 미만이면 위험비를 내지 않는다" 하나로 처리했습니다.
## 그 규칙은 106명·9사건인 Robust/T1 과 12명·9사건인 회복표의 Pre-frail/T3 를
## 똑같이 지웠습니다. 두 칸은 전혀 다른 상황입니다. 앞은 사건이 드문 것이고
## 뒤는 사람이 없는 것입니다. 사건 수와 인원을 따로 봅니다.
## ── 위험비 인쇄 정책 (v260813: '가능하면 인쇄, 주의 표시' 로 완화) ─────────
##  숫자를 최대한 보여주되, '추정치'라 부를 수 없는 것만 막습니다.
##   · 인쇄 + 주의표시(†): 사건 < FLAG_EV 또는 인원 < FLAG_N  -> 각주로 주의 안내
##   · 미인쇄(—)         : 사건 < MIN_EV_HR 또는 인원 < MIN_N_HR 또는 수치 발산
##     (사건 1-2건의 Cox HR 은 추정이 아니라 잡음입니다. 이 프로젝트에서 실제로
##      사건 2건짜리 칸이 HR 2.10 (1.13-3.93) 같은 성립 불가능한 값으로 인쇄된
##      전례가 있어, 그 선만은 지킵니다.)
MIN_EV_HR  <- 3      # 사건이 이보다 적으면 위험비를 내지 않습니다
MIN_N_HR   <- 10     # 칸 인원이 이보다 적으면 위험비를 내지 않습니다
FLAG_EV    <- 10     # 사건이 이보다 적으면 † (주의) 표시
FLAG_N     <- 30     # 인원이 이보다 적으면 † (주의) 표시
MIN_N_RISK <- 10     # 인원이 이보다 적으면 절대위험조차 인쇄하지 않습니다
                     #   (n = 2 짜리 칸이 "50.0%" 로 인쇄되는 것을 막습니다)
OUT_STEM   <- "Table3"

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
need <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop("패키지가 필요합니다: install.packages(\"", p, "\")")
need("survival")
has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
if (!has_cmprsk) message("  * cmprsk 없음 — Fine-Gray 열은 비웁니다 (필수 아님)")

## ★ 저장 위치를 고정합니다 (v260811c)
##  이전 판은 스크립트가 놓인 폴더를 따라갔고, 스크립트가 Downloads 에 있어서
##  결과가 Downloads 로 갔습니다. 결과 폴더를 직접 지정하고, 없을 때만
##  스크립트 위치로 물러납니다.
OUT_DIR <- FIG_DIR   ## 모든 산출물은 00_setup.R 의 단일 출력 폴더로
message("[Table3] 저장 폴더: ", OUT_DIR)

## ── 1. 자료 ─────────────────────────────────────────────────────────────────
load_imputed("MAIN")
d <- prep_long(load_kfacs())
b <- get_baseline(d)

b$tert <- factor(b$gLIC_tert_W1cut, levels = c("T1", "T2", "T3"))
b$fr   <- droplevels(factor(b$frailty_lab, levels = c("Robust", "Pre-frail", "Frail")))
b      <- b[!is.na(b$tert) & !is.na(b$fr), , drop = FALSE]

ADJ <- intersect(CFG$ADJ, names(b))
cat(sprintf("\n[0] 보정군 : %s\n", adj_phrase(ADJ)))
if (!setequal(ADJ, CFG$ADJ))
  warning("보정변수 일부가 baseline 에 없습니다: ",
          paste(setdiff(CFG$ADJ, ADJ), collapse = ", "))

## ── 2. 개인별 사건 시각 ─────────────────────────────────────────────────────
## 각 참여자에 대해 (기저상태, 마지막 관찰시각, 회복시각, 악화시각, 사망시각)
mk_events <- function(dat, sysvar) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  sp  <- split(seq_len(nrow(dat)), dat$id)
  res <- lapply(sp, function(ix) {
    s  <- dat[[sysvar]][ix]
    tm <- dat$time[ix]
    ok <- is.finite(tm)
    s  <- s[ok]; tm <- tm[ok]
    n  <- length(tm)
    if (n < 1) return(NULL)
    s0 <- s[1]
    if (!is.finite(s0)) return(NULL)
    tmax <- max(tm)
    fst <- function(cond) { k <- which(cond); if (length(k)) tm[-1][k[1]] else NA_real_ }
    tr <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] <  s0 & s[-1] >= 1) else NA_real_
    tw <- if (n >= 2) fst(is.finite(s[-1]) & s[-1] >  s0)              else NA_real_
    kd <- which(is.finite(s) & s == 4)
    td <- if (length(kd)) tm[kd[1]] else NA_real_

    ## 확정 악화 — 나빠진 상태가 다음 wave 에서도 유지되어야 사건입니다.
    ##   · 사망은 되돌아올 수 없으므로 확인 없이 곧바로 사건입니다.
    ##   · 일과성 악화(다음 wave 에 기저 이하로 복귀)는 건너뛰고 계속 추적합니다.
    ##   · 마지막 wave 에서 처음 나빠지면 확인할 wave 가 없으므로 사건이 아니라
    ##     그 시점의 중도절단입니다 (tmax 와 같은 시각이라 자동으로 처리됩니다).
    ts <- NA_real_
    if (n >= 2) {
      cand <- which(is.finite(s[-1]) & s[-1] > s0) + 1L
      for (jj in cand) {
        if (s[jj] == 4) { ts <- tm[jj]; break }
        if (jj < n && is.finite(s[jj + 1]) && s[jj + 1] > s0) { ts <- tm[jj]; break }
        if (jj == n) break
      }
    }
    data.frame(id = dat$id[ix[1]], s0 = s0, tmax = tmax,
               t_rec = tr, t_wor = tw, t_sus = ts, t_dth = td)
  })
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

## 절대위험/위험비를 위한 (시간, 사건) 만들기
mk_surv <- function(E, which_out) {
  inf <- function(x) ifelse(is.finite(x), x, Inf)
  if (which_out == "rec") {
    ft <- pmin(inf(E$t_rec), inf(E$t_dth), E$tmax, TAU)
    ev <- rep(0L, nrow(E))
    ev[is.finite(E$t_rec) & E$t_rec <= ft] <- 1L                       # 회복
    ev[ev == 0L & is.finite(E$t_dth) & E$t_dth <= ft] <- 2L            # 사망(경쟁)
  } else if (which_out == "wor") {
    ft <- pmin(inf(E$t_wor), E$tmax, TAU)                              # 사망 포함
    ev <- as.integer(is.finite(E$t_wor) & E$t_wor <= ft)
  } else if (which_out == "sus") {
    ft <- pmin(inf(E$t_sus), E$tmax, TAU)                              # 확정 악화
    ev <- as.integer(is.finite(E$t_sus) & E$t_sus <= ft)
  } else {
    ft <- pmin(inf(E$t_dth), E$tmax, TAU)
    ev <- as.integer(is.finite(E$t_dth) & E$t_dth <= ft)
  }
  ft <- pmax(ft, 1e-6)
  data.frame(id = E$id, s0 = E$s0, ft = ft, ev = ev)
}

## ── 3. 절대위험 ─────────────────────────────────────────────────────────────
## 회복 : Aalen-Johansen (사망 경쟁)   /  악화·사망 : 1 - KM
## ★ 신뢰구간은 survfit 이 계산한 변환 구간을 씁니다 (v260811c)
##  이전 판은 추정치 ± 1.96 × 표준오차 를 직접 계산했습니다. 위험이 0 이나
##  1 에 가까우면 정규근사가 경계를 뚫고 나가고, 저는 그것을 잘라내고
##  있었습니다. Robust/T3 회복이 "94.9 (87.9-100.0)" 으로 인쇄된 이유가
##  이것입니다. survfit 의 로그(또는 로그-로그) 변환 구간은 경계를 넘지
##  않으므로 그대로 씁니다.
risk_at <- function(ft, ev, competing, tt = TAU) {
  out <- c(est = NA_real_, lo = NA_real_, hi = NA_real_)
  ok <- is.finite(ft) & is.finite(ev)
  ft <- ft[ok]; ev <- ev[ok]
  if (!length(ft)) return(out)
  gcol <- function(s, nm, j) {
    m <- s[[nm]]; if (is.null(m)) return(NA_real_)
    if (is.matrix(m)) as.numeric(m[1, j]) else as.numeric(m[1])
  }
  if (competing) {
    if (!any(ev == 1L)) return(c(est = 0, lo = NA_real_, hi = NA_real_))
    ef <- factor(ev, levels = 0:2, labels = c("censor", "event", "compete"))
    f  <- try(survival::survfit(survival::Surv(ft, ef) ~ 1, conf.type = "log-log"),
              silent = TRUE)
    if (inherits(f, "try-error"))
      f <- try(survival::survfit(survival::Surv(ft, ef) ~ 1), silent = TRUE)
    if (inherits(f, "try-error")) return(out)
    s <- try(summary(f, times = tt, extend = TRUE), silent = TRUE)
    if (inherits(s, "try-error")) return(out)
    j <- match("event", f$states); if (is.na(j)) return(out)
    p  <- gcol(s, "pstate", j)
    lo <- gcol(s, "lower",  j)
    hi <- gcol(s, "upper",  j)
    if (!is.finite(lo) || !is.finite(hi)) {         # 하위 호환 — 정규근사로 물러남
      sd <- gcol(s, "std.err", j)
      lo <- p - 1.96 * sd; hi <- p + 1.96 * sd
    }
    out <- c(est = p, lo = max(0, lo), hi = min(1, hi))
  } else {
    if (!any(ev > 0)) return(c(est = 0, lo = NA_real_, hi = NA_real_))
    f <- try(survival::survfit(survival::Surv(ft, ev) ~ 1, conf.type = "log-log"),
             silent = TRUE)
    if (inherits(f, "try-error"))
      f <- try(survival::survfit(survival::Surv(ft, ev) ~ 1), silent = TRUE)
    if (inherits(f, "try-error")) return(out)
    s <- try(summary(f, times = tt, extend = TRUE), silent = TRUE)
    if (inherits(s, "try-error")) return(out)
    su <- s$surv[1]; lu <- s$lower[1]; uu <- s$upper[1]
    if (!is.finite(lu) || !is.finite(uu)) { lu <- su; uu <- su }
    out <- c(est = 1 - su, lo = max(0, 1 - min(1, uu)), hi = min(1, 1 - max(0, lu)))
  }
  out
}

## ── 4. 셀 정의 ──────────────────────────────────────────────────────────────
## 악화·사망 : 9칸 (Robust/Pre-frail/Frail × T1/T2/T3), 기준 Robust|T3
## 회복      : Robust·Pre-frail 은 삼분위, Frail 은 통합 (진단 결과 반영)
cell_full <- function(fr, tert) factor(paste(fr, tert, sep = "|"),
  levels = as.vector(t(outer(levels(fr), levels(tert), paste, sep = "|"))))

cell_rec <- function(fr, tert) {
  lab <- ifelse(as.character(fr) == "Frail", "Frail|all",
                paste(as.character(fr), as.character(tert), sep = "|"))
  lv <- c(paste("Robust",    c("T1","T2","T3"), sep = "|"),
          paste("Pre-frail", c("T1","T2","T3"), sep = "|"), "Frail|all")
  factor(lab, levels = lv)
}

## ── 5. 본표 (ADL 기준) ──────────────────────────────────────────────────────
Eadl <- mk_events(d, "state")
B    <- b[, c("id", "fr", "tert", ADJ), drop = FALSE]

fit_block <- function(E, B, which_out, eligible, cellf, ref) {
  S <- mk_surv(E, which_out)
  M <- merge(B, S, by = "id")
  M <- M[eligible(M$s0), , drop = FALSE]
  M$cell <- cellf(M$fr, M$tert)
  M <- M[!is.na(M$cell), , drop = FALSE]
  M$cell <- droplevels(M$cell)
  if (!ref %in% levels(M$cell)) stop("기준칸이 없습니다: ", ref)
  lv0 <- levels(M$cell)                    # 표에 찍을 자연 순서 (기준칸을 앞으로 빼기 전)
  M$cell <- stats::relevel(M$cell, ref = ref)

  comp <- which_out == "rec"
  lv   <- levels(M$cell)
  R <- do.call(rbind, lapply(lv, function(cl) {
    z  <- M[M$cell == cl, , drop = FALSE]
    ne <- sum(z$ev == 1L)
    rk <- risk_at(z$ft, z$ev, competing = comp)
    data.frame(cell = cl, n = nrow(z), events = ne,
               pyrs = sum(z$ft),
               risk = 100 * rk["est"], risk_lo = 100 * rk["lo"], risk_hi = 100 * rk["hi"],
               stringsAsFactors = FALSE)
  }))

  ## 원인별 Cox — 회복이면 사망은 중도절단
  fml <- stats::as.formula(paste("survival::Surv(ft, ev == 1L) ~ cell",
                                 if (length(ADJ)) paste("+", paste(ADJ, collapse = " + ")) else ""))
  cx <- try(survival::coxph(fml, data = M), silent = TRUE)
  R$HR <- NA_real_; R$HR_lo <- NA_real_; R$HR_hi <- NA_real_; R$P <- NA_real_
  R$se <- NA_real_
  if (!inherits(cx, "try-error")) {
    cf <- summary(cx)$coefficients
    ci <- summary(cx)$conf.int
    for (i in seq_len(nrow(R))) {
      nmi <- paste0("cell", R$cell[i])
      if (R$cell[i] == ref) { R$HR[i] <- 1; next }
      if (nmi %in% rownames(cf)) {
        R$HR[i]    <- ci[nmi, "exp(coef)"]
        R$HR_lo[i] <- ci[nmi, "lower .95"]
        R$HR_hi[i] <- ci[nmi, "upper .95"]
        R$P[i]     <- cf[nmi, ncol(cf)]
        R$se[i]    <- cf[nmi, "se(coef)"]
      }
    }
  }
  ## ★ 신뢰성 규칙 — 사건 수와 인원을 따로 봅니다 (v260811c)
  ##  drop : 사건이 너무 적거나(<5) 사람이 너무 적으면(<30) 위험비를 내지 않습니다
  ##  flag : 위험비는 내되 구간이 넓다고 표시합니다 (사건 <10)
  ##  hide : 인원이 한 자리면 절대위험도 인쇄하지 않습니다 (n = 2 짜리 칸 방지)
  unstable <- is.finite(R$se) & (R$se > 3 | abs(log(pmax(R$HR, 1e-12))) > 10)
  drop <- (R$events < MIN_EV_HR | R$n < MIN_N_HR | unstable) & R$cell != ref
  flag <- !drop & (R$events < FLAG_EV | R$n < FLAG_N) & R$cell != ref
  hide <- R$n < MIN_N_RISK
  R$HR[drop] <- NA; R$HR_lo[drop] <- NA; R$HR_hi[drop] <- NA; R$P[drop] <- NA
  R$risk[hide] <- NA; R$risk_lo[hide] <- NA; R$risk_hi[hide] <- NA
  R$HR[hide]   <- NA; R$HR_lo[hide]   <- NA; R$HR_hi[hide]   <- NA; R$P[hide] <- NA
  R$note <- ""
  R$note[flag] <- sprintf("사건 %d, n %d — 주의(†)", R$events[flag], R$n[flag])
  R$note[drop] <- sprintf("사건 %d, n %d — 추정 안 함", R$events[drop], R$n[drop])
  R$note[drop & unstable] <- "수치 발산 — 추정 억제"
  R$note[hide] <- sprintf("n = %d — 인원이 너무 적어 제시하지 않음", R$n[hide])
  R$flag <- flag; R$drop <- drop; R$hide <- hide

  ## Fine-Gray (회복만, 보조)
  R$sHR <- NA_real_
  if (comp && has_cmprsk) {
    X <- stats::model.matrix(stats::as.formula(paste("~ cell",
          if (length(ADJ)) paste("+", paste(ADJ, collapse = " + ")) else "")), data = M)[, -1, drop = FALSE]
    fg <- try(cmprsk::crr(M$ft, M$ev, X, failcode = 1, cencode = 0), silent = TRUE)
    if (!inherits(fg, "try-error")) {
      co <- exp(fg$coef); nmf <- names(fg$coef)
      for (i in seq_len(nrow(R))) {
        if (R$cell[i] == ref) { R$sHR[i] <- 1; next }
        k <- match(paste0("cell", R$cell[i]), nmf)
        if (!is.na(k) && !drop[i] && !hide[i]) R$sHR[i] <- co[k]
      }
    }
  }
  ## 기준칸을 맨 위로 올린 것은 모형 사정일 뿐입니다 — 표는 원래 순서로 되돌립니다
  R <- R[match(lv0, R$cell), , drop = FALSE]
  rownames(R) <- NULL
  attr(R, "n_model") <- nrow(M)
  R
}

###############################################################################
## 층내 추세검정 — 삼분위를 서수(1-3)로 넣은 층별 Cox 의 Wald P
##  희소 칸 하나에 기대지 않고 층 전체 자료를 쓰므로, 개별 칸이 신뢰성 규칙에
##  걸려 비어 있어도 층의 기울기는 검정력 있게 잡힙니다 (v260813 추가).
###############################################################################
trend_block <- function(E, B, which_out, eligible,
                        strata = c("Robust", "Pre-frail", "Frail")) {
  S <- mk_surv(E, which_out)
  M <- merge(B, S, by = "id")
  M <- M[eligible(M$s0), , drop = FALSE]
  out <- vapply(strata, function(fr) {
    z <- M[as.character(M$fr) == fr & !is.na(M$tert), , drop = FALSE]
    if (nrow(z) < MIN_N_HR || sum(z$ev == 1L) < MIN_EV_HR ||
        length(unique(z$tert)) < 2) return(NA_real_)
    z$tt <- as.numeric(z$tert)
    adj2 <- ADJ[vapply(ADJ, function(a)
      a %in% names(z) && length(unique(stats::na.omit(z[[a]]))) > 1, logical(1))]
    fml <- stats::as.formula(paste("survival::Surv(ft, ev == 1L) ~ tt",
             if (length(adj2)) paste("+", paste(adj2, collapse = " + ")) else ""))
    cx <- try(survival::coxph(fml, data = z), silent = TRUE)
    if (inherits(cx, "try-error")) return(NA_real_)
    cf <- try(summary(cx)$coefficients, silent = TRUE)
    if (inherits(cf, "try-error") || !"tt" %in% rownames(cf)) return(NA_real_)
    unname(cf["tt", ncol(cf)])
  }, numeric(1))
  stats::setNames(out, strata)
}

cat("\n[1] 복합 악화 (ADL 상태 악화 또는 사망)\n")
T_wor <- fit_block(Eadl, B, "wor", function(s0) is.finite(s0) & s0 < 4,
                   cell_full, "Robust|T3")
print(T_wor[, c("cell","n","events","risk","risk_lo","risk_hi","HR","HR_lo","HR_hi","P")],
      row.names = FALSE, digits = 3)

cat("\n[2] 사망\n")
T_dth <- fit_block(Eadl, B, "dth", function(s0) is.finite(s0) & s0 < 4,
                   cell_full, "Robust|T3")
print(T_dth[, c("cell","n","events","risk","risk_lo","risk_hi","HR","HR_lo","HR_hi","P")],
      row.names = FALSE, digits = 3)

cat("\n[3] 복합 회복 (ADL 상태 호전, 사망은 경쟁위험) — wave 1 Mild/Severe 만\n")
T_rec <- fit_block(Eadl, B, "rec", function(s0) is.finite(s0) & s0 %in% c(2, 3),
                   cell_full, "Robust|T3")   # v260813: Frail 도 T1/T2/T3 분리
                                             # (희소 칸은 신뢰성 규칙이 자동 처리)
print(T_rec[, c("cell","n","events","risk","risk_lo","risk_hi","HR","HR_lo","HR_hi","P","note")],
      row.names = FALSE, digits = 3)

## ── 6. 보충표 : 쇠약 기준 회복 ──────────────────────────────────────────────
T_rec_fr <- NULL
if ("frail_state" %in% names(d)) {
  Efr <- mk_events(d, "frail_state")
  cell_frrec <- function(fr, tert) {
    lab <- paste(as.character(fr), as.character(tert), sep = "|")
    lv  <- c(paste("Pre-frail", c("T1","T2","T3"), sep = "|"),
             paste("Frail",     c("T1","T2","T3"), sep = "|"))
    factor(ifelse(as.character(fr) == "Robust", NA, lab), levels = lv)
  }
  cat("\n[4] 보충 — 쇠약 기준 복합 회복 (wave 1 Pre-frail/Frail 만, 기준 Pre-frail|T3)\n")
  T_rec_fr <- fit_block(Efr, B, "rec", function(s0) is.finite(s0) & s0 %in% c(2, 3),
                        cell_frrec, "Pre-frail|T3")
  print(T_rec_fr[, c("cell","n","events","risk","risk_lo","risk_hi","HR","HR_lo","HR_hi","P","note")],
        row.names = FALSE, digits = 3)
} else {
  cat("\n[4] frail_state 가 없어 보충표를 건너뜁니다.\n")
}

## ── 6b. 보충표 : 확정 악화 (민감도) ─────────────────────────────────────────
## 이 코호트의 ADL 장애는 상당 부분 일과성입니다. 복합 악화는 "한 번이라도
## 넘어섰는가" 를 묻기 때문에 다음 wave 에 바로 복귀한 악화도 영구 악화와
## 똑같이 1건으로 셉니다. 일과성 악화가 내재역량과 약하게만 연관된다면
## 주 결과의 위험비는 영값 쪽으로 희석됩니다. 확정 악화는 그 희석을 걷어낸
## 값이므로, 주 결과가 유지되면 결론이 단단해지고 위험비가 커지면 그 자체가
## 본문에 쓸 숫자입니다.
cat("\n[5] 보충 — 확정 악화 (다음 wave 에서 유지된 악화, 사망은 확인 없이 사건)\n")
T_sus <- fit_block(Eadl, B, "sus", function(s0) is.finite(s0) & s0 < 4,
                   cell_full, "Robust|T3")
print(T_sus[, c("cell","n","events","risk","risk_lo","risk_hi","HR","HR_lo","HR_hi","P","note")],
      row.names = FALSE, digits = 3)

cat("\n  주 정의 대비 사건 수 변화\n")
CMP <- data.frame(
  cell        = T_wor$cell,
  ev_any      = T_wor$events,
  ev_sus      = T_sus$events[match(T_wor$cell, T_sus$cell)],
  HR_any      = T_wor$HR,
  HR_sus      = T_sus$HR[match(T_wor$cell, T_sus$cell)],
  stringsAsFactors = FALSE)
CMP$ev_drop_pct <- ifelse(CMP$ev_any > 0, 100 * (CMP$ev_any - CMP$ev_sus) / CMP$ev_any, NA_real_)
CMP$HR_ratio    <- ifelse(is.finite(CMP$HR_any) & CMP$HR_any > 0, CMP$HR_sus / CMP$HR_any, NA_real_)
print(CMP, row.names = FALSE, digits = 3)
cat(sprintf("  전체 사건 %d → %d (%.1f%% 감소)\n",
            sum(CMP$ev_any), sum(CMP$ev_sus, na.rm = TRUE),
            100 * (sum(CMP$ev_any) - sum(CMP$ev_sus, na.rm = TRUE)) / max(1, sum(CMP$ev_any))))
.n_unconf <- sum(is.finite(Eadl$t_wor) & !is.finite(Eadl$t_sus) & Eadl$t_wor <= TAU)
cat(sprintf("  확인되지 않아 사건에서 빠진 사람 %d명 (일과성 또는 마지막 wave 발생)\n", .n_unconf))

## ── 6c. 층내 추세검정 (P for trend) ────────────────────────────────────────
.el_all <- function(s0) is.finite(s0) & s0 < 4
.el_rec <- function(s0) is.finite(s0) & s0 %in% c(2, 3)
TRw <- trend_block(Eadl, B, "wor", .el_all)
TRd <- trend_block(Eadl, B, "dth", .el_all)
TRs <- trend_block(Eadl, B, "sus", .el_all)
TRr <- trend_block(Eadl, B, "rec", .el_rec, strata = c("Robust", "Pre-frail"))
TRrf <- if (exists("Efr") && !is.null(T_rec_fr))
  trend_block(Efr, B, "rec", .el_rec, strata = c("Pre-frail", "Frail")) else NULL
cat("\n[5b] 층내 추세검정 P (서수 삼분위, 층별 Cox)\n")
.pt <- function(nm, v) cat(sprintf("  %-12s %s\n", nm,
  paste(sprintf("%s %s", names(v),
                ifelse(is.finite(v), sprintf("%.4f", v), "-")), collapse = " | ")))
.pt("악화",      TRw); .pt("사망", TRd); .pt("확정악화", TRs)
.pt("회복(ADL)", TRr); if (!is.null(TRrf)) .pt("회복(쇠약)", TRrf)

## ── 7. 표 조립 ──────────────────────────────────────────────────────────────
fmt_r <- function(e, lo, hi) {
  o <- rep("\u2014", length(e))
  k <- is.finite(e) & is.finite(lo) & is.finite(hi)
  o[k] <- sprintf("%.1f (%.1f\u2013%.1f)", e[k], lo[k], hi[k])
  j <- is.finite(e) & !k
  o[j] <- sprintf("%.1f", e[j])
  o
}
fmt_h <- function(h, lo, hi)
  ifelse(!is.finite(h), "\u2014",
         ifelse(h == 1, "1.00 (reference)", sprintf("%.2f (%.2f\u2013%.2f)", h, lo, hi)))
fmt_p <- function(p) ifelse(!is.finite(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

key <- T_wor$cell
TAB <- data.frame(
  Stratum          = sub("\\|.*", "", key),
  Tertile          = sub(".*\\|", "", key),
  n                = T_wor$n,
  Worsen_events    = T_wor$events,
  Worsen_risk      = fmt_r(T_wor$risk, T_wor$risk_lo, T_wor$risk_hi),
  Worsen_HR        = paste0(fmt_h(T_wor$HR, T_wor$HR_lo, T_wor$HR_hi),
                            ifelse(T_wor$flag, "\u2020", "")),
  Worsen_P         = fmt_p(T_wor$P),
  Death_events     = T_dth$events,
  Death_risk       = fmt_r(T_dth$risk, T_dth$risk_lo, T_dth$risk_hi),
  Death_HR         = paste0(fmt_h(T_dth$HR, T_dth$HR_lo, T_dth$HR_hi),
                            ifelse(T_dth$flag, "\u2020", "")),
  Death_P          = fmt_p(T_dth$P),
  stringsAsFactors = FALSE)

rk <- T_rec$cell
RECTAB <- data.frame(
  Stratum        = sub("\\|.*", "", rk),
  Tertile        = sub(".*\\|", "", rk),
  n              = T_rec$n,
  Recov_events   = T_rec$events,
  Recov_risk     = fmt_r(T_rec$risk, T_rec$risk_lo, T_rec$risk_hi),
  Recov_HR       = paste0(fmt_h(T_rec$HR, T_rec$HR_lo, T_rec$HR_hi),
                          ifelse(T_rec$flag, "\u2020", "")),
  Recov_P        = fmt_p(T_rec$P),
  Recov_sHR      = ifelse(is.finite(T_rec$sHR), sprintf("%.2f", T_rec$sHR), ""),
  Note           = T_rec$note,
  stringsAsFactors = FALSE)

sk <- T_sus$cell
SUSTAB <- data.frame(
  Stratum      = sub("\\|.*", "", sk),
  Tertile      = sub(".*\\|", "", sk),
  n            = T_sus$n,
  Any_events   = T_wor$events[match(sk, T_wor$cell)],
  Any_HR       = fmt_h(T_wor$HR, T_wor$HR_lo, T_wor$HR_hi)[match(sk, T_wor$cell)],
  Sust_events  = T_sus$events,
  Sust_risk    = fmt_r(T_sus$risk, T_sus$risk_lo, T_sus$risk_hi),
  Sust_HR      = paste0(fmt_h(T_sus$HR, T_sus$HR_lo, T_sus$HR_hi),
                        ifelse(T_sus$flag, "\u2020", "")),
  Sust_P       = fmt_p(T_sus$P),
  Note         = T_sus$note,
  stringsAsFactors = FALSE)

RECFR <- NULL
if (!is.null(T_rec_fr)) {
  kf <- T_rec_fr$cell
  RECFR <- data.frame(
    Stratum      = sub("\\|.*", "", kf),
    Tertile      = sub(".*\\|", "", kf),
    n            = T_rec_fr$n,
    Recov_events = T_rec_fr$events,
    Recov_risk   = fmt_r(T_rec_fr$risk, T_rec_fr$risk_lo, T_rec_fr$risk_hi),
    Recov_HR     = paste0(fmt_h(T_rec_fr$HR, T_rec_fr$HR_lo, T_rec_fr$HR_hi),
                          ifelse(T_rec_fr$flag, "\u2020", "")),
    Recov_P      = fmt_p(T_rec_fr$P),
    Note         = T_rec_fr$note,
    stringsAsFactors = FALSE)
}

## ── 8. 각주 (코드와 일치하도록 새로 작성) ───────────────────────────────────
.eg <- b[is.finite(b$state) & b$state %in% c(2, 3), , drop = FALSE]
.tt <- table(factor(.eg$fr, levels = levels(b$fr)), factor(.eg$tert, levels = c("T1","T2","T3")))
.n_fr_t2 <- as.integer(.tt["Frail", "T2"])
.n_fr_t3 <- as.integer(.tt["Frail", "T3"])
.n_sev_bl <- sum(is.finite(b$state) & b$state == 3)

FOOT <- c(
 sprintf("Strata are defined at wave 1 by frailty phenotype and by tertile of the general intrinsic-capacity factor, tertile cut-points being taken from the wave-1 distribution. Follow-up is censored at %g years.", TAU),
 "Composite worsening is the first transition to a worse ADL disability state; death is the most severe state and is therefore included in this outcome by construction. Absolute risks are one minus the Kaplan-Meier estimate.",
 "Death is all-cause mortality.",
 "Composite recovery is the first transition to a better ADL disability state and is defined only for participants with mild or severe ADL disability at wave 1. Death precludes recovery and is treated as a competing event: absolute risks are Aalen-Johansen cumulative incidences, hazard ratios are from cause-specific Cox models, and subdistribution hazard ratios from Fine-Gray models are given for comparison.",
 sprintf("For recovery, only %d and %d participants were frail at wave 1 and in the middle or the highest capacity tertile while carrying ADL disability; those two cells are therefore reported as counts only.",
         .n_fr_t2, .n_fr_t3),
 sprintf("Hazard ratios marked \u2020 come from cells with fewer than %d events or fewer than %d participants and should be interpreted with caution: the point estimate is unstable and the interval wide. A dash is shown only where an estimate would be numerically meaningless (fewer than %d events, fewer than %d participants, or a divergent fit); cells with fewer than %d participants are reported as counts only.",
         FLAG_EV, FLAG_N, MIN_EV_HR, MIN_N_HR, MIN_N_RISK),
 sprintf("All models are adjusted for %s.", adj_phrase(ADJ)),
 "Reference cell is the robust stratum, highest capacity tertile, for every outcome in this table.",
 sprintf("Participants with severe ADL disability at wave 1 (n = %d) can worsen only by dying, so their worsening column is a mortality column; their number is too small to affect the estimates.", .n_sev_bl))

FOOT_SUS <- c(
 "Sustained worsening requires the worse ADL disability state to persist at the following wave. Death is irreversible and is counted without confirmation. A transient worsening that has returned to the wave-1 state or better by the next wave is not counted and the participant remains at risk.",
 "Participants whose first worsening occurs at their final observed wave cannot be confirmed and are censored at that wave rather than counted as events; this is the conventional treatment of confirmed-progression endpoints and it uses information from the confirming wave, which may fall beyond the 5-year analysis window.",
 "This analysis addresses the possibility that the primary definition, which counts any worsening at any wave, dilutes the association with intrinsic capacity by including short-lived ADL disability. In this cohort a large fraction of ADL disability resolves, so the concern is not hypothetical.",
 sprintf("All models are adjusted for %s. The reference cell is unchanged.", adj_phrase(ADJ)))

cat("\n[6] 각주\n"); for (i in seq_along(FOOT)) cat(sprintf("  %d. %s\n", i, FOOT[i]))
cat("\n  확정 악화 보충표 각주\n")
for (i in seq_along(FOOT_SUS)) cat(sprintf("  S%d. %s\n", i, FOOT_SUS[i]))

## ── 본표 병합: 결과별 세로 블록 (v260813 — 한 표로 통합) ────────────────
.sec_row <- function(lab) data.frame(
  S = lab, T = "", n = "", E = "", R = "", H = "", P = "", sH = "",
  stringsAsFactors = FALSE)
.blk_row <- function(S, T, n, E, R, H, P, sH = "") data.frame(
  S = S, T = T, n = as.character(n), E = as.character(E),
  R = blank(R), H = blank(H), P = blank(P), sH = blank(sH),
  stringsAsFactors = FALSE)
blank <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "\u2014", as.character(x))
.blank_keep <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "", as.character(x))
MAIN_PANEL <- rbind(
  .sec_row("Composite worsening (ADL state or death)"),
  .blk_row(TAB$Stratum, TAB$Tertile, TAB$n, TAB$Worsen_events,
           TAB$Worsen_risk, TAB$Worsen_HR, TAB$Worsen_P),
  .sec_row("Death (all cause)"),
  .blk_row(TAB$Stratum, TAB$Tertile, TAB$n, TAB$Death_events,
           TAB$Death_risk, TAB$Death_HR, TAB$Death_P),
  .sec_row("Recovery from ADL disability (death as competing risk)"),
  .blk_row(RECTAB$Stratum, RECTAB$Tertile, RECTAB$n, RECTAB$Recov_events,
           RECTAB$Recov_risk, RECTAB$Recov_HR, RECTAB$Recov_P, RECTAB$Recov_sHR))
## 구획행은 대시 대신 빈칸으로
.is_sec <- MAIN_PANEL$T == ""
MAIN_PANEL[.is_sec, -1] <- ""
names(MAIN_PANEL) <- c("Frailty phenotype", "IC tertile", "n", "Events",
                       "5-yr risk, % (95% CI)", "HR (95% CI)", "P",
                       "Subdistribution HR")

## ── 9. 저장 ─────────────────────────────────────────────────────────────────
outx <- file.path(OUT_DIR, paste0(OUT_STEM, ".xlsx"))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  add <- function(nm, x) { if (is.null(x)) return(invisible());
    openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, x)
    openxlsx::setColWidths(wb, nm, 1:ncol(x), "auto") }
  add("Table3_main", MAIN_PANEL)
  add("Table3_worsen_death", TAB)
  add("Table3_recovery_ADL", RECTAB)
  add("ST_sustained_worsen", SUSTAB)
  add("ST_recovery_frailty", RECFR)
  add("footnotes", data.frame(
        table = c(rep("Table 3", length(FOOT)), rep("ST sustained worsening", length(FOOT_SUS))),
        no    = c(seq_along(FOOT), seq_along(FOOT_SUS)),
        text  = c(FOOT, FOOT_SUS)))
  add("cmp_any_vs_sustained", CMP)
  add("raw_trend", data.frame(
        outcome = rep(c("worsening","death","sustained"), each = length(TRw)),
        stratum = rep(names(TRw), 3),
        p_trend = c(TRw, TRd, TRs)))
  add("raw_worsen",  T_wor); add("raw_death", T_dth); add("raw_recov", T_rec)
  add("raw_sustained", T_sus)
  if (!is.null(T_rec_fr)) add("raw_recov_frailty", T_rec_fr)
  openxlsx::saveWorkbook(wb, outx, overwrite = TRUE)
  cat(sprintf("\n저장: %s\n", outx))
} else {
  utils::write.csv(TAB,    file.path(OUT_DIR, paste0(OUT_STEM, "_worsen_death.csv")), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(RECTAB, file.path(OUT_DIR, paste0(OUT_STEM, "_recovery.csv")),     row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(SUSTAB, file.path(OUT_DIR, paste0(OUT_STEM, "_sustained.csv")),    row.names = FALSE, fileEncoding = "UTF-8")
  if (!is.null(RECFR))
    utils::write.csv(RECFR, file.path(OUT_DIR, paste0(OUT_STEM, "_recovery_frailty.csv")), row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n저장: CSV (openxlsx 없음)\n")
}

## ── 9b. Word 출력 (Nature Aging 서식) ───────────────────────────────────────
## 세로줄 없음 · 표 위와 머리글 아래와 표 아래에만 가로줄 · Arial 8pt ·
## 표 제목은 "Table 3 | ..." 형태 · 각주는 표 바로 아래 단락.
## 표가 넓어 가로 방향 절(section)에 넣습니다.
TIT3  <- "Table 3 | Five-year risk of worsening, death and recovery according to frailty phenotype and intrinsic-capacity tertile at wave 1"
TITS  <- "Supplementary Table | Sustained worsening: worsening required to persist at the following wave"
TITF  <- "Supplementary Table | Recovery defined on the frailty phenotype rather than on ADL disability"

blank <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "—", as.character(x))


A_PANEL <- data.frame(
  `Frailty phenotype` = TAB$Stratum,
  `IC tertile`        = TAB$Tertile,
  `n`                 = TAB$n,
  `Worsening: events`               = TAB$Worsen_events,
  `Worsening: 5-yr risk, % (95% CI)`= blank(TAB$Worsen_risk),
  `Worsening: HR (95% CI)`          = blank(TAB$Worsen_HR),
  `Worsening: P`                    = blank(TAB$Worsen_P),
  `Death: events`                   = TAB$Death_events,
  `Death: 5-yr risk, % (95% CI)`    = blank(TAB$Death_risk),
  `Death: HR (95% CI)`              = blank(TAB$Death_HR),
  `Death: P`                        = blank(TAB$Death_P),
  check.names = FALSE, stringsAsFactors = FALSE)

B_PANEL <- data.frame(
  `Frailty phenotype` = RECTAB$Stratum,
  `IC tertile`        = RECTAB$Tertile,
  `n`                 = RECTAB$n,
  `Events`            = RECTAB$Recov_events,
  `5-year cumulative incidence, % (95% CI)` = blank(RECTAB$Recov_risk),
  `Cause-specific HR (95% CI)` = blank(RECTAB$Recov_HR),
  `P`                 = blank(RECTAB$Recov_P),
  `Subdistribution HR` = blank(RECTAB$Recov_sHR),
  check.names = FALSE, stringsAsFactors = FALSE)

S_PANEL <- data.frame(
  `Frailty phenotype` = SUSTAB$Stratum,
  `IC tertile`        = SUSTAB$Tertile,
  `n`                 = SUSTAB$n,
  `Any worsening, events`       = SUSTAB$Any_events,
  `Any worsening, HR (95% CI)`  = blank(SUSTAB$Any_HR),
  `Sustained, events`           = SUSTAB$Sust_events,
  `Sustained, 5-year risk, % (95% CI)` = blank(SUSTAB$Sust_risk),
  `Sustained, HR (95% CI)`      = blank(SUSTAB$Sust_HR),
  `P`                           = blank(SUSTAB$Sust_P),
  check.names = FALSE, stringsAsFactors = FALSE)

F_PANEL <- if (!is.null(RECFR)) data.frame(
  `Frailty phenotype` = RECFR$Stratum,
  `IC tertile`        = RECFR$Tertile,
  `n`                 = RECFR$n,
  `Events`            = RECFR$Recov_events,
  `5-year cumulative incidence, % (95% CI)` = blank(RECFR$Recov_risk),
  `Cause-specific HR (95% CI)` = blank(RECFR$Recov_HR),
  `P`                 = blank(RECFR$Recov_P),
  check.names = FALSE, stringsAsFactors = FALSE) else NULL

outd <- file.path(OUT_DIR, paste0(OUT_STEM, ".docx"))
ok_docx <- FALSE

if (requireNamespace("officer", quietly = TRUE)) {
  ok_docx <- tryCatch({
    has_flx <- requireNamespace("flextable", quietly = TRUE)
    doc <- officer::read_docx()

    ttl <- function(d, t) officer::body_add_par(d, t, style = "Normal")
    .sub_par <- function(d, t) officer::body_add_par(d, t, style = "Normal")

    addtab <- function(d, df) {
      if (has_flx) {
        ft <- flextable::flextable(df)
        ft <- flextable::theme_booktabs(ft)          # 가로줄만, 세로줄 없음
        .si <- which(!nzchar(as.character(df[[2]])) & nzchar(as.character(df[[1]])))
        if (length(.si)) ft <- flextable::bold(ft, i = .si, part = "body")
        ft <- flextable::font(ft, fontname = "Arial", part = "all")
        ft <- flextable::fontsize(ft, size = 8, part = "all")
        ft <- flextable::padding(ft, padding = 2, part = "all")
        ft <- flextable::align(ft, j = 3:ncol(df), align = "center", part = "all")
        ft <- flextable::autofit(ft)
        flextable::body_add_flextable(d, ft)
      } else officer::body_add_table(d, df, header = TRUE)
    }
    foot <- function(d, v, pre = "") {
      for (i in seq_along(v))
        d <- officer::body_add_par(d, paste0(pre, i, ". ", v[i]), style = "Normal")
      d
    }

    doc <- ttl(doc, TIT3)
    doc <- addtab(doc, MAIN_PANEL)
    doc <- officer::body_add_par(doc, "")
    doc <- foot(doc, FOOT)
    doc <- officer::body_add_break(doc)

    doc <- ttl(doc, TITS)
    doc <- addtab(doc, S_PANEL)
    doc <- officer::body_add_par(doc, "")
    doc <- foot(doc, FOOT_SUS, pre = "S")

    if (!is.null(F_PANEL)) {
      doc <- officer::body_add_break(doc)
      doc <- ttl(doc, TITF)
      doc <- addtab(doc, F_PANEL)
      doc <- officer::body_add_par(doc, "")
      doc <- officer::body_add_par(doc,
        "1. Recovery is the first transition to a better frailty state and is defined only for participants who were pre-frail or frail at wave 1; the robust stratum is therefore absent. Death is treated as a competing event.", style = "Normal")
      doc <- officer::body_add_par(doc,
        "2. The reference cell is the pre-frail stratum, highest capacity tertile, and therefore differs from Table 3. This analysis is reported separately for that reason.", style = "Normal")
    }
    doc <- officer::body_end_section_landscape(doc)
    print(doc, target = outd)
    TRUE
  }, error = function(e) { message("  officer 실패: ", conditionMessage(e)); FALSE })
}

if (!ok_docx) {
  ## 대체 경로 — Word 가 그대로 여는 HTML 을 .doc 로 씁니다.
  esc <- function(x) { x <- as.character(x)
    x <- gsub("&","&amp;",x,fixed=TRUE); x <- gsub("<","&lt;",x,fixed=TRUE)
    gsub(">","&gt;",x,fixed=TRUE) }
  htab <- function(df) paste0(
    "<table cellspacing='0' cellpadding='3' style='border-collapse:collapse;font-family:Arial;font-size:8pt'>",
    "<tr>", paste0("<th style='border-top:1px solid #000;border-bottom:1px solid #000;text-align:center'>",
                   esc(names(df)), "</th>", collapse=""), "</tr>",
    paste0(apply(df, 1, function(r)
      paste0("<tr>", paste0("<td style='text-align:center'>", esc(r), "</td>", collapse=""), "</tr>")),
      collapse=""),
    "<tr><td colspan='", ncol(df), "' style='border-top:1px solid #000'></td></tr></table>")
  fl <- function(v, pre="") paste0("<p style='font-family:Arial;font-size:7.5pt;margin:2pt 0'>",
                                   pre, seq_along(v), ". ", esc(v), "</p>", collapse="")
  h <- c("<html xmlns:o='urn:schemas-microsoft-com:office:office'><head><meta charset='utf-8'>",
    "<style>@page{size:A4 landscape;margin:2cm} p{font-family:Arial}</style></head><body>",
    paste0("<p style='font-family:Arial;font-size:9pt'><b>", esc(TIT3), "</b></p>"),
    htab(MAIN_PANEL), fl(FOOT),
    "<br style='page-break-before:always'>",
    paste0("<p style='font-family:Arial;font-size:9pt'><b>", esc(TITS), "</b></p>"), htab(S_PANEL), fl(FOOT_SUS, "S"),
    if (!is.null(F_PANEL)) c("<br style='page-break-before:always'>",
      paste0("<p style='font-family:Arial;font-size:9pt'><b>", esc(TITF), "</b></p>"), htab(F_PANEL)) else "",
    "</body></html>")
  outd <- file.path(OUT_DIR, paste0(OUT_STEM, ".doc"))
  writeLines(h, outd, useBytes = TRUE)
  cat("\n  ※ officer 패키지가 없어 Word 가 여는 .doc(HTML) 로 저장했습니다.\n")
  cat("     진짜 .docx 를 원하시면: install.packages(c(\"officer\",\"flextable\"))\n")
}
cat(sprintf("저장: %s\n", outd))

## ── 10. 본문에 쓸 문장용 숫자 ───────────────────────────────────────────────
cat("\n[7] 본문 인용용\n")
pick <- function(R, cl) R[R$cell == cl, , drop = FALSE]
if (!is.null(T_rec_fr)) {
  for (cl in paste("Pre-frail", c("T1","T2","T3"), sep = "|")) {
    z <- pick(T_rec_fr, cl)
    if (nrow(z)) cat(sprintf("  쇠약기준 회복 %-14s : %d/%d, 5년 누적발생 %.1f%% (%.1f-%.1f)\n",
                             cl, z$events, z$n, z$risk, z$risk_lo, z$risk_hi))
  }
}
for (cl in c("Robust|T1","Robust|T3")) {
  z <- pick(T_wor, cl); y <- pick(T_sus, cl)
  if (nrow(z)) cat(sprintf("  복합악화     %-14s : %d/%d, 5년 위험 %.1f%%  | 확정악화 %d건, HR %.2f\n",
                           cl, z$events, z$n, z$risk, y$events, y$HR))
}
cat("\n=== 완료 ===\n")
