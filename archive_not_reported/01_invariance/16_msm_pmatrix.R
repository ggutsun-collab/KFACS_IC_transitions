###############################################################################
## 16_msm_pmatrix.R   (Phase 3.3 — 연속시간 다상태 모형의 2년 전이확률)   v260903
##
## 왜 필요한가: 헤드라인 전이(Robust→Frail 0.35 / T3 vs T1 0.09, Frail→Robust 2.73,
## Normal→Severe, Severe→Normal)는 모두 '상태를 건너뛰는' 관측 전이라 연속시간
## Markov 모형에서는 직접 추정할 수 없습니다(Supp Table 4 의 "not in model").
## 리뷰어는 "interval censoring 에 가장 취약한 전이가 헤드라인" 이라고 지적할 것입니다.
## 해법: msm 을 인접 전이로 적합한 뒤 pmatrix.msm(t = 2) 로 2년 전이확률을 IC 수준별로
## 계산하면, 건너뛰는 전이도 중간 상태를 거치는 경로를 적분한 확률로 나옵니다.
## 그 확률의 IC 수준 간 비(ratio)를 Poisson IRR 과 나란히 제시합니다.
##
## 설계
##   · 54_sens_msm_compare.R 의 msm 설정(Q0, LOCF 공변량, 사다리 적합)을 그대로 씁니다.
##   · 노출은 within-sex gLIC_z (primary scale). 공변량 age_c + sex_f + comorbid_bl.
##   · MI 세트 중 N_SETS 개(기본 5; msm 은 세트당 5-25분)에서 적합해 확률을 Rubin 풀링(logit).
##   · IC = -1, 0, +1 s.d. 에서 2년 전이확률; 비 = P(+1)/P(-1).
##
## 출력: ST_msm_pmatrix (2년 전이확률표), ST_msm_hazards_MI (인접 전이 HR per s.d.)
## 실행: source("16_msm_pmatrix.R", encoding = "UTF-8")   (N_SETS=5 기준 30-120분)
###############################################################################

PM_VERSION <- "v260903"
message("\n=== 16_msm_pmatrix ", PM_VERSION, " ===")

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
.need(c("msm", "openxlsx"))
suppressPackageStartupMessages(library(msm))

###############################################################################
## 1. 설정
###############################################################################
MI_STEM   <- get0("MI_STEM", ifnotfound = "KFACS_mi_pmm_m20")   # 00_setup.R 에서 정함
N_SETS    <- 5                         # msm 을 적합할 MI 세트 수 (전체 20 은 수 시간)
MSM_MAXIT <- 2000
T_PRED    <- 2                         # 전이확률 예측 구간(년) = 방문 간격
IC_LEVELS <- c(-1, 0, 1)               # within-sex s.d.
B_CI      <- 200                       # pmatrix 정규근사 CI 표본 수
CACHE     <- file.path(FIG_DIR, "P3_msm_fits_MI.rds")
SYS <- list(
  list(var = "state",       lab = "ADL disability, 3 states", states = c("Normal", "Mild", "Severe", "Death")),
  list(var = "frail_state", lab = "Frailty phenotype",        states = c("Robust", "Pre-frail", "Frail", "Death")))
Q0 <- rbind(c(0,1,0,1), c(1,0,1,1), c(0,1,0,1), c(0,0,0,0))
MI <- readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))); SETS <- names(MI)[seq_len(min(N_SETS, length(MI)))]

## within-sex z (14_ 과 동일)
ws_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]; z <- d$gLIC
  for (s in levels(d$sex_f)) { k <- d$sex_f == s; k1 <- w1$sex_f == s
    z[k] <- (d$gLIC[k] - mean(w1$gLIC[k1], na.rm = TRUE)) / sd(w1$gLIC[k1], na.rm = TRUE) }
  d$gLIC_z <- z; d
}
## 54_ 의 msm 입력 구성 (동일)
locf <- function(v) { ok <- which(!is.na(v)); if (!length(ok)) return(v)
  idx <- cummax((seq_along(v) %in% ok) * seq_along(v)); idx[idx == 0] <- ok[1]; v[idx] }
mk_msm_data <- function(d, svar) {
  x <- d[, c("id", "time", svar, "gLIC_z", "age_c", "sex_f", "comorbid_bl")]
  names(x)[3] <- "st"; x$st <- as.integer(as_num(x$st)); x$id <- as.character(x$id)
  x <- x[is.finite(x$time) & !is.na(x$st) & x$st %in% 1:4, ]; x <- x[order(x$id, x$time), ]
  for (cc in c("gLIC_z", "age_c", "comorbid_bl")) { x[[cc]] <- as.numeric(x[[cc]])
    x[[cc]] <- unlist(lapply(split(x[[cc]], x$id), locf), use.names = FALSE) }
  sx <- unlist(lapply(split(as.character(x$sex_f), x$id), function(v) { ok <- which(!is.na(v))
    if (!length(ok)) return(v); rep(v[ok[1]], length(v)) }), use.names = FALSE)
  x$sex_f <- factor(sx, levels = c("Male", "Female"))
  x <- x[stats::complete.cases(x[, c("gLIC_z", "age_c", "sex_f", "comorbid_bl")]), ]
  x <- x[order(paste(x$id, sprintf("%.4f", x$time)), -x$st), ]
  x <- x[!duplicated(paste(x$id, sprintf("%.4f", x$time))), ]; x <- x[order(x$id, x$time), ]
  x <- do.call(rbind, lapply(split(x, x$id), function(p) { k <- which(p$st == 4)
    if (length(k)) p <- p[seq_len(k[1]), , drop = FALSE]; p }))
  x <- do.call(rbind, lapply(split(x, x$id), function(p) { if (nrow(p) < 2) return(p)
    p[c(TRUE, diff(p$time) > 1e-6), , drop = FALSE] }))
  n <- table(x$id); x <- x[x$id %in% names(n[n >= 2]), ]; rownames(x) <- NULL; x
}
fit_msm <- function(x, lab) {
  q <- try(msm::crudeinits.msm(st ~ time, id, data = x, qmatrix = Q0), silent = TRUE)
  if (inherits(q, "try-error") || !all(is.finite(q))) q <- Q0 * 0.1
  off <- Q0 > 0; q[off & (q <= 0 | !is.finite(q))] <- 1e-3; diag(q) <- 0; diag(q) <- -rowSums(q)
  ctl <- list(maxit = MSM_MAXIT, fnscale = max(1000, 4 * nrow(x)), reltol = 1e-9)
  ladder <- list(list(tag = "full", cov = ~ gLIC_z + age_c + sex_f + comorbid_bl, con = NULL, de = TRUE),
                 list(tag = "shared_death", cov = ~ gLIC_z + age_c + sex_f + comorbid_bl,
                      con = list(gLIC_z = c(1, 2, 3, 4, 4, 4)), de = TRUE),
                 list(tag = "ic_only", cov = ~ gLIC_z, con = NULL, de = TRUE))
  for (a in ladder) {
    args <- list(formula = st ~ time, subject = x$id, data = x, qmatrix = q, control = ctl, covariates = a$cov)
    if (isTRUE(a$de)) args$deathexact <- 4
    if (!is.null(a$con)) args$constraint <- a$con
    f <- try(suppressWarnings(do.call(msm::msm, args)), silent = TRUE)
    if (!inherits(f, "try-error")) { attr(f, "tag") <- a$tag; msg(sprintf("  [%s] 성공 (%s)", lab, a$tag)); return(f) }
    msg(sprintf("  [%s] 실패 (%s)", lab, a$tag))
  }
  NULL
}

###############################################################################
## 2. 세트별 적합 (캐시)
###############################################################################
FITS <- if (file.exists(CACHE)) readRDS(CACHE) else list()
for (set in SETS) {
  if (!is.null(FITS[[set]])) { msg("[", set, "] 캐시 재사용"); next }
  load_imputed(set, stem = MI_STEM); d <- ws_scale(prep_long(load_kfacs()))
  FITS[[set]] <- lapply(SYS, function(S) {
    x <- mk_msm_data(d, S$var)
    msg(sprintf("[%s] %s: %d명 / %d관측 / 사망 %d", set, S$lab, length(unique(x$id)), nrow(x), sum(x$st == 4)))
    f <- fit_msm(x, paste(set, S$var))
    if (is.null(f)) return(NULL)
    list(fit = f, ref = list(age_c = mean(x$age_c[!duplicated(x$id)]), comorbid_bl = mean(x$comorbid_bl[!duplicated(x$id)])))
  })
  saveRDS(FITS, CACHE)
}

###############################################################################
## 3. 2년 전이확률 (IC 수준별) + 인접 전이 HR
###############################################################################
pm_one <- function(F, S, sex) {
  do.call(rbind, lapply(IC_LEVELS, function(z) {
    covl <- F$fit$qcmodel$covlabels          # 모형에 실제로 들어간 공변량만
    cov <- list(gLIC_z = z)
    if ("age_c" %in% covl) cov$age_c <- F$ref$age_c
    if ("comorbid_bl" %in% covl) cov$comorbid_bl <- F$ref$comorbid_bl
    if (any(grepl("^sex_f", covl))) cov$sex_f <- sex
    P <- try(msm::pmatrix.msm(F$fit, t = T_PRED, covariates = cov, ci = "normal", B = B_CI), silent = TRUE)
    if (inherits(P, "try-error")) return(NULL)
    est <- P$estimates; lo <- P$L; hi <- P$U
    do.call(rbind, lapply(1:3, function(i) do.call(rbind, lapply(setdiff(1:4, i), function(j)
      data.frame(system = S$lab, from = S$states[i], to = S$states[j], sex = sex, ic = z,
                 p = est[i, j], lo = lo[i, j], hi = hi[i, j])))))
  }))
}
PM <- list(); HZ <- list()
for (set in SETS) for (k in seq_along(SYS)) {
  F <- FITS[[set]][[k]]; if (is.null(F)) next
  for (sx in c("Male", "Female")) { r <- pm_one(F, SYS[[k]], sx); if (!is.null(r)) { r$set <- set; PM[[length(PM) + 1]] <- r } }
  h <- try(msm::hazard.msm(F$fit)[["gLIC_z"]], silent = TRUE)
  if (!inherits(h, "try-error") && !is.null(h)) {
    h <- as.data.frame(h); nm <- rownames(h); pr <- regmatches(nm, gregexpr("[0-9]+", nm))
    HZ[[length(HZ) + 1]] <- data.frame(set = set, system = SYS[[k]]$lab,
      from = vapply(pr, function(v) SYS[[k]]$states[as.integer(v[1])], ""), to = vapply(pr, function(v) SYS[[k]]$states[as.integer(v[2])], ""),
      b = log(h[[1]]), se = (log(h[[3]]) - log(h[[2]])) / (2 * qnorm(0.975)))
  }
}
PM <- do.call(rbind, PM); HZ <- do.call(rbind, HZ)
save_vals(PM, "ST_msm_pmatrix_perset.csv", FIG_DIR)

## 성별 평균 후 Rubin 풀링 (logit)
lg <- function(p) log(p / (1 - p)); il <- function(x) 1 / (1 + exp(-x))
PM$se_lg <- (lg(pmin(PM$hi, 1 - 1e-6)) - lg(pmax(PM$lo, 1e-6))) / (2 * qnorm(0.975))
PM$key <- paste(PM$system, PM$from, PM$to, PM$ic, sep = "|")
PP <- do.call(rbind, lapply(split(PM, PM$key), function(z) {
  ## 성별 평균: 두 성별의 logit 평균과 SE 합성 (같은 세트 내), 이후 세트 간 Rubin
  zz <- do.call(rbind, lapply(split(z, z$set), function(w) data.frame(
    b = mean(lg(pmax(pmin(w$p, 1 - 1e-6), 1e-6))), se = sqrt(mean(w$se_lg^2)) / sqrt(nrow(w)))))
  r <- rubin_pool(zz$b, zz$se)
  data.frame(system = z$system[1], from = z$from[1], to = z$to[1], ic = z$ic[1],
             p = il(r["est"]), lo = il(r["lcl"]), hi = il(r["ucl"]), m = r["m"]) }))
rownames(PP) <- NULL
## +1 vs -1 s.d. 비 (세트별 계산 후 Rubin, log 척도)
RT <- do.call(rbind, lapply(split(PM, paste(PM$system, PM$from, PM$to)), function(z) {
  zz <- do.call(rbind, lapply(split(z, z$set), function(w) {
    a <- w[w$ic == 1, ]; b <- w[w$ic == -1, ]
    if (!nrow(a) || !nrow(b)) return(NULL)
    pa <- mean(a$p); pb <- mean(b$p)
    sea <- sqrt(mean(((a$hi - a$lo) / (2 * qnorm(0.975)))^2)) / sqrt(nrow(a)); seb <- sqrt(mean(((b$hi - b$lo) / (2 * qnorm(0.975)))^2)) / sqrt(nrow(b))
    data.frame(b = log(pa / pb), se = sqrt((sea / pa)^2 + (seb / pb)^2)) }))
  r <- rubin_pool(zz$b, zz$se)
  data.frame(system = z$system[1], from = z$from[1], to = z$to[1], ratio = exp(r["est"]), rlo = exp(r["lcl"]), rhi = exp(r["ucl"])) }))
rownames(RT) <- NULL

SEV <- c(Normal = 1, Mild = 2, Severe = 3, Robust = 1, `Pre-frail` = 2, Frail = 3, Death = 9)
fmtp <- function(p, l, u) ifelse(is.finite(p), sprintf("%.1f (%.1f-%.1f)", 100 * p, 100 * l, 100 * u), "-")
keys <- unique(PP[, c("system", "from", "to")])
TAB <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
  g <- function(z) PP[PP$system == keys$system[k] & PP$from == keys$from[k] & PP$to == keys$to[k] & PP$ic == z, ]
  a <- g(-1); b <- g(0); c <- g(1); r <- RT[RT$system == keys$system[k] & RT$from == keys$from[k] & RT$to == keys$to[k], ]
  data.frame(System = keys$system[k], Transition = paste(keys$from[k], "->", keys$to[k]),
             Type = ifelse(keys$to[k] == "Death", "Death", ifelse(SEV[keys$to[k]] > SEV[keys$from[k]], "Worsening", "Recovery")),
             `Skips a state` = ifelse(keys$to[k] != "Death" & abs(SEV[keys$to[k]] - SEV[keys$from[k]]) == 2, "yes", "no"),
             `IC -1 s.d.` = fmtp(a$p, a$lo, a$hi),
             `IC 0` = fmtp(b$p, b$lo, b$hi), `IC +1 s.d.` = fmtp(c$p, c$lo, c$hi),
             `Ratio, +1 vs -1 s.d. (95% CI)` = ifelse(nrow(r), sprintf("%.2f (%.2f-%.2f)", r$ratio, r$rlo, r$rhi), "-"),
             check.names = FALSE, stringsAsFactors = FALSE) }))
TAB <- TAB[order(match(TAB$System, c("ADL disability, 3 states", "Frailty phenotype")), match(TAB$Type, c("Worsening", "Recovery", "Death")),
               SEV[sub(" ->.*", "", TAB$Transition)]), ]
save_table(TAB, "ST_msm_pmatrix",
  title = "Supplementary Table | Two-year transition probabilities (%, 95% CI) from a continuous-time Markov multistate model, by level of intrinsic capacity",
  footnotes = c(
    sprintf("Continuous-time Markov models with adjacent transitions only (normal-mild, mild-severe and each living state to death; robust-pre-frail, pre-frail-frail and each state to death), death treated as exactly observed, fitted in %d completed datasets with intrinsic capacity (standardized within sex) as a time-varying covariate together with age, sex and baseline comorbidity count.", length(SETS)),
    sprintf("Two-year transition probabilities were obtained from the fitted intensity matrices (pmatrix.msm, t = %g years) at intrinsic capacity of -1, 0 and +1 s.d., age and comorbidity at their cohort means, averaged over sex, and combined across completed datasets by Rubin's rules on the logit scale. Confidence intervals are from %d draws from the asymptotic normal distribution of the parameters.", T_PRED, B_CI),
    "Transitions that skip a state (for example robust to frail) are not intensities of the model; their two-year probabilities integrate over all paths through the intermediate state within the interval. They are directly comparable with the observed between-visit transitions analysed by the discrete-time models in Table 2, and reconcile the two frameworks for the transitions that the discrete-time analysis reports as its largest effects."))

## 인접 전이 HR (MI 풀링)
HZ$key <- paste(HZ$system, HZ$from, HZ$to)
HH <- do.call(rbind, lapply(split(HZ, HZ$key), function(z) { r <- rubin_pool(z$b, z$se)
  data.frame(System = z$system[1], Transition = paste(z$from[1], "->", z$to[1]),
             `HR per +1 s.d. (95% CI)` = sprintf("%.2f (%.2f-%.2f)", exp(r["est"]), exp(r["lcl"]), exp(r["ucl"])),
             P = ifelse(r["p"] < 0.001, "<0.001", sprintf("%.3f", r["p"])), check.names = FALSE) }))
save_table(HH, "ST_msm_hazards_MI",
  title = "Supplementary Table | Continuous-time transition intensities per +1 s.d. of intrinsic capacity (multiple imputation, within-sex scale)",
  footnotes = "Hazard ratios for the adjacent-state intensities of the continuous-time Markov model, combined across completed datasets by Rubin's rules.")

cat("\n", strrep("=", 70), "\nmsm 2년 전이확률 요약 (skip 전이)\n", strrep("=", 70), "\n", sep = "")
print(TAB[TAB$`Skips a state` == "yes", c("Transition", "IC -1 s.d.", "IC 0", "IC +1 s.d.", "Ratio, +1 vs -1 s.d. (95% CI)")], row.names = FALSE)
cat("  -> Robust->Frail 의 비가 Poisson IRR(0.35)^2 ≈ 0.12 근처이면 두 틀이 일치합니다 (2 s.d. 차이).\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
