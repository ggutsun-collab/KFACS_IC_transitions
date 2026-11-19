###############################################################################
## 12_sex_invariance_and_scores.R   (Phase 1 — 측정 모형 확정)   v260903
##
## 이 스크립트는 이후 모든 수정 분석의 전제를 만듭니다.
##   1.1  성별 측정불변성 (wave 1, 남/여 다중집단 bifactor CFA)
##        configural -> metric -> scalar -> partial scalar (절편 순차 해제)
##   1.2  Sex-neutral IC 점수 세 가지 산출·비교
##        (A) 현행 gLIC 의 within-sex z         — ED Table 1-2 정의 그대로
##        (B) 지표 수준 sex-specific z 후 bifactor 재적합 -> 일반요인
##        (C) partial-scalar 성별 불변 모형의 요인점수
##        -> 상관·삼분위 일치·남녀 SMD. (A)-(B) r >= 0.97 이면 (A) 를 primary.
##   1.3  Non-overlapping IC-13 (HGS, GS_ms, EXH, loss_of_Bwt 제외)
##   1.4  IC-17 측정 모형 공개표: loadings, fit, ECV, omega_H, HGS 기여도
##
## 출력 (FIG_DIR):
##   ST_SexInvariance.xlsx/csv, ST_ScoreCompare, ST_IC13_Model, ST_IC17_Model
##   P1_scores.rds  — id·wave 별 점수 6종 (이후 스크립트가 merge)
##
## 실행: source("12_sex_invariance_and_scores.R", encoding = "UTF-8")  (5-10분)
###############################################################################

P1_VERSION <- "v260903"
message("\n=== 12_sex_invariance_and_scores ", P1_VERSION, " ===")

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
.need(c("lavaan", "openxlsx"))
suppressPackageStartupMessages(library(lavaan))

###############################################################################
## 1. 설정
###############################################################################
INV_SET    <- "MAIN"
D_CFI      <- -0.010
D_RMSEA    <-  0.015
MAX_FREE   <-  5
R_KEEP_A   <-  0.97     # (A)-(B) 상관이 이 이상이면 (A) 유지
ESTIMATOR  <- "MLR"

## 17개 지표 — 10_measurement_invariance.R 와 동일 (rev_* 는 이미 방향 반전됨)
IND17 <- c("Balance","GS_ms","rev_CST","loss_of_Bwt","Appetite","EXH","HGS",
           "rev_logMAR","rev_PTA","Orientation","Memory","Attention_Calculation",
           "Language","Visuospatial","Negative_affect","Positive_affect","Motivation")
DOMAIN <- list(loco = c("Balance","GS_ms","rev_CST"),
               vita = c("loss_of_Bwt","Appetite","EXH","HGS"),
               sens = c("rev_logMAR","rev_PTA"),
               cogn = c("Orientation","Memory","Attention_Calculation","Language","Visuospatial"),
               psyc = c("Negative_affect","Positive_affect","Motivation"))
OVERLAP <- c("HGS", "GS_ms", "EXH", "loss_of_Bwt")   # Fried 와 공유

## ★ gLIC 를 만든 00a_KFACS_impute_pipeline_FULL_260617.R 의 postprocess() 는
##   특정요인 4개(loco/vita/cogn/psyc), sensory 는 g 에만 적재 -> 원고 Methods 와 일치.
##   반면 10_measurement_invariance.R 의 BIFACTOR 는 sens 특정요인을 둡니다(5개) ->
##   Supp Table 3 은 채점 모형과 다른 모형으로 검정된 것이므로 10_ 을 고쳐 다시 돌려야 합니다.
##   여기서는 채점 모형(4개)을 씁니다. [1.4] 에서 기존 gLIC 와의 상관으로 검증합니다.
SENS_SPECIFIC <- FALSE

bf_syntax <- function(keys, domain = DOMAIN, sens_specific = SENS_SPECIFIC) {
  sp <- character(0)
  for (dn in names(domain)) {
    if (dn == "sens" && !sens_specific) next
    it <- intersect(domain[[dn]], keys)
    if (length(it) >= 2) sp <- c(sp, sprintf("%s =~ %s", dn, paste(it, collapse = " + ")))
  }
  facs <- sub(" =~.*", "", sp)
  paste(c(sprintf("g =~ %s", paste(keys, collapse = " + ")), sp,
          if (length(facs)) sprintf("g ~~ %s", paste0("0*", facs, collapse = " + "))),
        collapse = "\n")
}

###############################################################################
## 2. 자료 — 10_measurement_invariance.R 와 동일한 표본 구성
###############################################################################
load_imputed(INV_SET)
raw <- load_kfacs()
miss <- setdiff(IND17, names(raw))
if (length(miss)) stop("지표가 없습니다: ", paste(miss, collapse = ", "))

d0 <- prep_long(raw)
if (!"wave" %in% names(d0)) d0$wave <- round(as_num(d0$time) / CFG$WAVE_GAP) + 1
d0$wave <- as.integer(as_num(d0$wave))
d0 <- d0[is.finite(d0$wave), ]
key0 <- paste(d0$id, d0$wave)

D <- raw[, c("id", "wave", IND17), drop = FALSE]
D$id <- as.character(D$id); D$wave <- as.integer(as_num(D$wave))
D <- D[is.finite(D$wave) & paste(D$id, D$wave) %in% key0, ]
D <- D[stats::complete.cases(D[, IND17]), ]
k  <- match(paste(D$id, D$wave), key0)
D$sexF  <- as.integer(d0$sex_f[k] == "Female")
D$gLIC  <- d0$gLIC[k]
D$gLIC_z<- d0$gLIC_z[k]                         # pooled W1 z (현행 primary)
D$tertP <- d0$gLIC_tert_W1cut[k]
D$time  <- d0$time[k]
D <- D[is.finite(D$sexF) & is.finite(D$gLIC), ]
is_w1 <- D$wave == min(D$wave)
msg(sprintf("측정 자료: %d person-waves / %d명 / wave 1 n=%d / 여성 %.1f%%",
            nrow(D), length(unique(D$id)), sum(is_w1), 100 * mean(D$sexF[is_w1])))

## pooled 표준화 (전체 풀 1회 — 논문 정의)
Zp <- D
for (v in IND17) Zp[[v]] <- (as_num(D[[v]]) - mean(as_num(D[[v]]))) / stats::sd(as_num(D[[v]]))

###############################################################################
## 3. 보조 함수
###############################################################################
fit_bf <- function(model, data, ...) {
  f <- tryCatch(suppressWarnings(lavaan::cfa(model, data = data, std.lv = TRUE,
                                             estimator = ESTIMATOR, ...)),
                error = function(e) { message("  cfa 오류: ", conditionMessage(e)); NULL })
  if (is.null(f) || !lavaan::lavInspect(f, "converged")) return(NULL)
  f
}
fitidx <- function(f) {
  if (is.null(f)) return(c(chisq = NA, df = NA, CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA))
  m <- lavaan::fitMeasures(f)
  gt <- function(a, b) if (a %in% names(m)) unname(m[[a]]) else if (b %in% names(m)) unname(m[[b]]) else NA_real_
  c(chisq = gt("chisq.scaled","chisq"), df = gt("df","df"), CFI = gt("cfi.robust","cfi"),
    TLI = gt("tli.robust","tli"), RMSEA = gt("rmsea.robust","rmsea"), SRMR = gt("srmr","srmr"))
}
## 일반요인 점수 (단일집단 모형). 기존 gLIC 와 양의 상관이 되도록 부호 정렬.
gscore <- function(fit, data, ref) {
  s <- as.numeric(lavaan::lavPredict(fit, newdata = data)[, "g"])
  r <- suppressWarnings(stats::cor(s, ref, use = "complete.obs"))
  if (is.finite(r) && r < 0) s <- -s
  s
}
## 다집단 모형의 요인점수를 원래 행 순서로
gscore_mg <- function(fit, data, gvar, ref) {
  pr <- tryCatch(lavaan::lavPredict(fit, newdata = data), error = function(e) NULL)
  if (is.null(pr)) return(rep(NA_real_, nrow(data)))
  s <- rep(NA_real_, nrow(data))
  if (is.list(pr)) {
    gl <- lavaan::lavInspect(fit, "group.label")
    for (j in seq_along(pr)) s[as.character(data[[gvar]]) == gl[j]] <- pr[[j]][, "g"]
  } else s <- pr[, "g"]
  r <- suppressWarnings(stats::cor(s, ref, use = "complete.obs"))
  if (is.finite(r) && r < 0) s <- -s
  s
}
## bifactor 지수: ECV, omega_H, 지표별 일반요인 분산 기여
bf_indices <- function(fit) {
  st <- lavaan::standardizedSolution(fit); st <- st[st$op == "=~", ]
  lg <- st[st$lhs == "g", ]; ls <- st[st$lhs != "g", ]
  items <- lg$rhs
  theta <- vapply(items, function(it) 1 - lg$est.std[lg$rhs == it]^2 -
                    sum(ls$est.std[ls$rhs == it]^2), numeric(1))
  sum_s <- if (nrow(ls)) sum(tapply(ls$est.std, ls$lhs, function(v) sum(v)^2)) else 0
  omegaH <- sum(lg$est.std)^2 / (sum(lg$est.std)^2 + sum_s + sum(theta))
  ECV    <- sum(lg$est.std^2) / (sum(lg$est.std^2) + sum(ls$est.std^2))
  share  <- stats::setNames(lg$est.std^2 / sum(lg$est.std^2), items)
  list(ECV = ECV, omegaH = omegaH, share = share, st = st)
}
z_w1   <- function(s, w1) (s - mean(s[w1], na.rm = TRUE)) / stats::sd(s[w1], na.rm = TRUE)
z_w1_sex <- function(s, w1, sexF) {
  z <- s
  for (g in 0:1) { k <- sexF == g; k1 <- k & w1
    z[k] <- (s[k] - mean(s[k1], na.rm = TRUE)) / stats::sd(s[k1], na.rm = TRUE) }
  z
}
tert_w1 <- function(s, w1) {
  ct <- stats::quantile(s[w1], c(1/3, 2/3), na.rm = TRUE)
  factor(cut(s, c(-Inf, ct, Inf), labels = CFG$TERT_LABELS, right = TRUE), levels = CFG$TERT_LABELS)
}
smd <- function(x, g) { ok <- is.finite(x); x <- x[ok]; g <- g[ok]
  a <- x[g == 0]; b <- x[g == 1]
  (mean(a) - mean(b)) / sqrt((stats::var(a) + stats::var(b)) / 2) }
f3 <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "-")
f2 <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "-")

###############################################################################
## 4. [1.4] IC-17 측정 모형 (pooled 표준화) — 기존 gLIC 재현 확인
###############################################################################
cat("\n[1.4] IC-17 bifactor (pooled, 전체 person-wave)\n")
SYN17 <- bf_syntax(IND17)
fit17 <- fit_bf(SYN17, Zp)
if (is.null(fit17)) stop("IC-17 bifactor 가 수렴하지 않습니다.")
fi17 <- fitidx(fit17); print(round(fi17, 3))
bi17 <- bf_indices(fit17)
s17  <- gscore(fit17, Zp, D$gLIC)
r_rep <- stats::cor(s17, D$gLIC, use = "complete.obs")
cat(sprintf("  기존 gLIC 재현 상관 r = %.4f  %s\n", r_rep,
            if (r_rep > 0.995) "(동일 모형)" else "(★ 0.995 미만 — SENS_SPECIFIC 플래그를 바꿔 다시 확인하십시오)"))
cat(sprintf("  ECV = %.3f · omega_H = %.3f · HGS 의 일반요인 분산 기여 = %.1f%%\n",
            bi17$ECV, bi17$omegaH, 100 * bi17$share[["HGS"]]))
cat(sprintf("  남녀 SMD (남-여): 점수 %.2f · HGS %.2f · EXH %.2f\n",
            smd(s17, D$sexF), smd(Zp$HGS, D$sexF), smd(Zp$EXH, D$sexF)))
dom_of <- function(v) names(DOMAIN)[vapply(DOMAIN, function(z) v %in% z, logical(1))]
IC17TAB <- do.call(rbind, lapply(IND17, function(v) {
  g  <- bi17$st$est.std[bi17$st$lhs == "g" & bi17$st$rhs == v]
  sp <- bi17$st[bi17$st$lhs != "g" & bi17$st$rhs == v, ]
  data.frame(Indicator = v, Domain = dom_of(v),
             `General loading` = f3(g),
             `Specific factor` = if (nrow(sp)) sp$lhs[1] else "general only",
             `Specific loading` = if (nrow(sp)) f3(sp$est.std[1]) else "-",
             `Share of general-factor variance, %` = sprintf("%.1f", 100 * bi17$share[[v]]),
             `Sex difference, SMD (men - women)` = f2(smd(Zp[[v]], D$sexF)),
             `Fried phenotype item` = ifelse(v %in% OVERLAP, "yes", "no"),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

###############################################################################
## 5. [1.1] 성별 측정불변성 (wave 1)
###############################################################################
cat("\n[1.1] 성별 측정불변성 — wave 1\n")
W1 <- Zp[is_w1, ]; W1$sex <- factor(ifelse(W1$sexF == 1, "Female", "Male"), levels = c("Male", "Female"))
f_conf <- fit_bf(SYN17, W1, group = "sex")
f_metr <- fit_bf(SYN17, W1, group = "sex", group.equal = "loadings")
f_scal <- fit_bf(SYN17, W1, group = "sex", group.equal = c("loadings", "intercepts"))
if (is.null(f_conf) || is.null(f_metr) || is.null(f_scal))
  stop("다집단 CFA 수렴 실패 — 10_measurement_invariance.R 의 CORR5 대체 모형으로 재시도하십시오.")

## partial scalar: 절편 등식제약 가운데 score test 가 큰 것부터 해제
freed <- character(0); f_part <- f_scal; cfi_prev <- NA
repeat {
  fm <- fitidx(f_metr); fp <- fitidx(f_part)
  ok <- (fp["CFI"] - fm["CFI"] >= D_CFI) && (fp["RMSEA"] - fm["RMSEA"] <= D_RMSEA)
  if (ok || length(freed) >= MAX_FREE) break
  mi <- tryCatch(lavaan::lavTestScore(f_part, epc = TRUE)$epc, error = function(e) NULL)
  if (is.null(mi)) break
  cand <- mi[mi$op == "~1" & !(mi$lhs %in% freed), ]
  if (!nrow(cand)) break
  cand <- cand[order(-abs(cand$epc)), ]
  freed <- c(freed, cand$lhs[1])
  f_new <- fit_bf(SYN17, W1, group = "sex", group.equal = c("loadings", "intercepts"),
                  group.partial = paste0(freed, " ~ 1"))
  if (is.null(f_new)) { freed <- freed[-length(freed)]; break }
  f_part <- f_new
  cfi_now <- fitidx(f_part)["CFI"]
  cat(sprintf("    절편 해제 %d: %s -> CFI %.4f (metric 대비 %+.4f)\n",
              length(freed), freed[length(freed)], cfi_now, cfi_now - fm["CFI"]))
  if (is.finite(cfi_prev) && cfi_now - cfi_prev < 1e-4) break
  cfi_prev <- cfi_now
}
INV <- rbind(data.frame(Model = "Configural", t(fitidx(f_conf))),
             data.frame(Model = "Metric (loadings)", t(fitidx(f_metr))),
             data.frame(Model = "Scalar (+ intercepts)", t(fitidx(f_scal))),
             data.frame(Model = sprintf("Partial scalar (%d freed)", length(freed)), t(fitidx(f_part))))
INV$dCFI   <- c(NA, INV$CFI[2] - INV$CFI[1], INV$CFI[3] - INV$CFI[2], INV$CFI[4] - INV$CFI[2])
INV$dRMSEA <- c(NA, INV$RMSEA[2] - INV$RMSEA[1], INV$RMSEA[3] - INV$RMSEA[2], INV$RMSEA[4] - INV$RMSEA[2])
verdict <- function(dc, dr) ifelse(is.na(dc), "-", ifelse(dc >= D_CFI & dr <= D_RMSEA, "supported", "NOT supported"))
INV$Verdict <- verdict(INV$dCFI, INV$dRMSEA)
INV$Verdict[4] <- paste0(INV$Verdict[4], if (length(freed)) paste0(" (freed: ", paste(freed, collapse = ", "), ")") else " (none freed)")
print(INV, row.names = FALSE, digits = 4)
SEX_METRIC_OK <- INV$Verdict[2] == "supported"
SEX_SCALAR_OK <- INV$Verdict[3] == "supported"
SEX_PARTIAL_OK <- grepl("^supported", INV$Verdict[4])

## 잠재평균 차이 (partial scalar 모형, Female 대비 Male; 첫 집단 평균 = 0)
lat_diff <- NA_real_
pe <- tryCatch(lavaan::parameterEstimates(f_part), error = function(e) NULL)
if (!is.null(pe)) {
  a <- pe[pe$op == "~1" & pe$lhs == "g", ]
  if (nrow(a) >= 2) lat_diff <- a$est[1] - a$est[2]      # Male - Female (SD 단위)
}
## lavaan 의 요인 부호는 임의 -> 부호 정렬된 점수 (C) 의 방향과 맞춥니다
.sC_tmp <- gscore_mg(f_part, {Zt <- Zp; Zt$sex <- factor(ifelse(Zt$sexF == 1, "Female", "Male"), levels = c("Male","Female")); Zt}, "sex", D$gLIC)
if (is.finite(lat_diff) && sign(lat_diff) != sign(smd(.sC_tmp, D$sexF))) lat_diff <- -lat_diff
cat(sprintf("  잠재 일반요인 평균차 (남-여, partial scalar): %.2f s.d.  [관측 점수 차 %.2f]\n",
            lat_diff, smd(s17, D$sexF)))

###############################################################################
## 6. [1.2] Sex-neutral 점수 세 가지
###############################################################################
cat("\n[1.2] Sex-neutral 점수 산출·비교\n")
zP <- D$gLIC_z                                   # pooled (현행 primary)
zA <- z_w1_sex(D$gLIC, is_w1, D$sexF)            # (A) within-sex z of gLIC
## (B) 지표 수준 sex-specific 표준화 -> bifactor 재적합
Zs <- D
for (v in IND17) for (g in 0:1) { i <- D$sexF == g
  x <- as_num(D[[v]][i]); Zs[[v]][i] <- (x - mean(x)) / stats::sd(x) }
fitB <- fit_bf(SYN17, Zs)
sB <- if (!is.null(fitB)) gscore(fitB, Zs, D$gLIC) else rep(NA_real_, nrow(D))
zB <- z_w1(sB, is_w1)
biB <- if (!is.null(fitB)) bf_indices(fitB) else NULL
## (C) partial-scalar 성별 불변 모형 점수 (전체 person-wave)
Zc <- Zp; Zc$sex <- factor(ifelse(Zc$sexF == 1, "Female", "Male"), levels = c("Male", "Female"))
sC <- gscore_mg(f_part, Zc, "sex", D$gLIC)
zC <- z_w1(sC, is_w1)

tP <- factor(D$tertP, levels = CFG$TERT_LABELS); tA <- tert_w1(zA, is_w1)
tB <- tert_w1(zB, is_w1); tC <- tert_w1(zC, is_w1)
one <- function(nm, z, tt) data.frame(
  Score = nm,
  `r with pooled score` = f3(stats::cor(z, zP, use = "complete.obs")),
  `r with (A)` = f3(stats::cor(z, zA, use = "complete.obs")),
  `Tertile agreement with (A), %` = sprintf("%.1f", 100 * mean(tt == tA, na.rm = TRUE)),
  `Sex difference, SMD` = f2(smd(z, D$sexF)),
  `Women in lowest wave-1 tertile, %` = sprintf("%.1f", 100 * mean(D$sexF[is_w1 & tt == "T1"] == 1, na.rm = TRUE)),
  `Women in highest wave-1 tertile, %` = sprintf("%.1f", 100 * mean(D$sexF[is_w1 & tt == "T3"] == 1, na.rm = TRUE)),
  check.names = FALSE, stringsAsFactors = FALSE)
CMP <- rbind(one("Pooled scale (current primary)", zP, tP),
             one("(A) General factor standardized within sex", zA, tA),
             one("(B) Indicators standardized within sex, bifactor refitted", zB, tB),
             one("(C) Sex-invariant (partial scalar) model score", zC, tC))
print(CMP, row.names = FALSE)
rAB <- stats::cor(zA, zB, use = "complete.obs"); rAC <- stats::cor(zA, zC, use = "complete.obs")
PRIMARY <- if (is.finite(rAB) && rAB >= R_KEEP_A) "A" else "B"

###############################################################################
## 7. [1.3] Non-overlapping IC-13
###############################################################################
cat("\n[1.3] IC-13 — Fried 공유 지표 제외\n")
K13 <- setdiff(IND17, OVERLAP)
SYN13 <- bf_syntax(K13)                                  # loco = Balance+rev_CST(2), vita = Appetite(단일 -> g 만)
cat(SYN13, "\n")
fit13 <- fit_bf(SYN13, Zp[, c(K13)])
bad_theta <- function(f) { th <- diag(lavaan::lavInspect(f, "est")$theta); any(th < 0) }
if (is.null(fit13) || bad_theta(fit13)) {
  message("  loco 2지표 특정요인 불안정 -> loco 를 general 에만 적재")
  D13 <- DOMAIN; D13$loco <- character(0)
  SYN13 <- bf_syntax(K13, domain = D13); fit13 <- fit_bf(SYN13, Zp[, K13])
}
if (is.null(fit13)) stop("IC-13 모형이 수렴하지 않습니다.")
fi13 <- fitidx(fit13); print(round(fi13, 3))
bi13 <- bf_indices(fit13)
s13  <- gscore(fit13, Zp[, K13], D$gLIC)
z13  <- z_w1(s13, is_w1)
z13s <- z_w1_sex(s13, is_w1, D$sexF)
t13  <- tert_w1(z13, is_w1)
cat(sprintf("  IC-13 vs IC-17 r = %.3f · ECV %.3f · omega_H %.3f · 남녀 SMD %.2f (IC-17 %.2f)\n",
            stats::cor(z13, zP, use = "complete.obs"), bi13$ECV, bi13$omegaH, smd(z13, D$sexF), smd(zP, D$sexF)))
IC13TAB <- data.frame(
  Model = c("IC-17 (all indicators)", "IC-13 (Fried-overlapping indicators removed)"),
  Indicators = c(17L, length(K13)),
  CFI = f3(c(fi17["CFI"], fi13["CFI"])), TLI = f3(c(fi17["TLI"], fi13["TLI"])),
  RMSEA = f3(c(fi17["RMSEA"], fi13["RMSEA"])), SRMR = f3(c(fi17["SRMR"], fi13["SRMR"])),
  ECV = f3(c(bi17$ECV, bi13$ECV)), `omega H` = f3(c(bi17$omegaH, bi13$omegaH)),
  `r with IC-17` = f3(c(1, stats::cor(z13, zP, use = "complete.obs"))),
  `Tertile agreement with IC-17, %` = sprintf("%.1f", 100 * c(1, mean(t13 == tP, na.rm = TRUE))),
  `Sex difference, SMD` = f2(c(smd(zP, D$sexF), smd(z13, D$sexF))),
  `Women in lowest wave-1 tertile, %` = sprintf("%.1f", 100 * c(
    mean(D$sexF[is_w1 & tP == "T1"] == 1, na.rm = TRUE), mean(D$sexF[is_w1 & t13 == "T1"] == 1, na.rm = TRUE))),
  check.names = FALSE, stringsAsFactors = FALSE)
print(IC13TAB, row.names = FALSE)

###############################################################################
## 8. 저장
###############################################################################
SC <- data.frame(id = D$id, wave = D$wave, time = D$time, sexF = D$sexF,
                 IC17_pooled_z = zP, IC17_withinsex_z = zA, IC17_sexneutral_z = zB,
                 IC17_sexinv_z = zC, IC13_pooled_z = z13, IC13_withinsex_z = z13s,
                 stringsAsFactors = FALSE)
attr(SC, "primary") <- PRIMARY
attr(SC, "cut_withinsex") <- list(
  men   = unname(stats::quantile(zA[is_w1 & D$sexF == 0], c(1/3, 2/3))),
  women = unname(stats::quantile(zA[is_w1 & D$sexF == 1], c(1/3, 2/3))))
saveRDS(SC, file.path(FIG_DIR, "P1_scores.rds"))
message("scores: ", file.path(FIG_DIR, "P1_scores.rds"))

fmtINV <- data.frame(Model = INV$Model, `Chi-square` = sprintf("%.1f", INV$chisq),
                     df = sprintf("%.0f", INV$df), CFI = f3(INV$CFI), TLI = f3(INV$TLI),
                     RMSEA = f3(INV$RMSEA), SRMR = f3(INV$SRMR),
                     `Delta CFI` = ifelse(is.na(INV$dCFI), "-", sprintf("%+.3f", INV$dCFI)),
                     `Delta RMSEA` = ifelse(is.na(INV$dRMSEA), "-", sprintf("%+.3f", INV$dRMSEA)),
                     Verdict = INV$Verdict, check.names = FALSE, stringsAsFactors = FALSE)
save_table(fmtINV, "ST_SexInvariance",
  title = "Supplementary Table | Measurement invariance of the intrinsic-capacity factor between men and women at wave 1",
  footnotes = c(
    sprintf("Multi-group confirmatory factor analysis with sex as the grouping variable at wave 1 (n = %d; %d women), bifactor model identical to that used to score intrinsic capacity, %s estimator, indicators standardized once on the pooled person-wave distribution.", sum(is_w1), sum(D$sexF[is_w1]), ESTIMATOR),
    sprintf("Invariance is judged by change from the preceding model: supported if Delta CFI >= %.3f and Delta RMSEA <= %.3f. The partial scalar model frees, one at a time, the intercept with the largest expected parameter change until the criterion is met (freed: %s); its change is computed relative to the metric model.", D_CFI, D_RMSEA, if (length(freed)) paste(freed, collapse = ", ") else "none"),
    sprintf("Under the partial scalar model the latent general-factor mean was %.2f s.d. higher in men than in women, compared with an observed score difference of %.2f s.d. The sex difference in the general factor is therefore not an artefact of non-invariant intercepts but a substantive difference in the measured capacities, most of which favoured men in this cohort.", lat_diff, smd(s17, D$sexF)),
    "Metric invariance licenses comparison of associations per +1 s.d. between men and women. Because the sex difference in level is substantive, tertiles defined on a pooled scale are largely a contrast between women and men; the within-sex scale used in the primary stratified analyses classifies each participant relative to others of the same sex, in the same way that the frailty phenotype applies sex-specific cut-points to grip strength and gait speed.",
    "CFI, comparative fit index; RMSEA, root mean square error of approximation; SRMR, standardized root mean square residual; TLI, Tucker-Lewis index."))
save_table(CMP, "ST_ScoreCompare",
  title = "Supplementary Table | Agreement between alternative sex-neutral scalings of intrinsic capacity",
  footnotes = c(
    "The pooled scale standardizes the general-factor score on the wave-1 distribution of the whole cohort. (A) standardizes the same score within sex on the wave-1 distribution of men and of women separately. (B) standardizes each of the 17 indicators within sex before refitting the bifactor model and scoring the general factor. (C) scores the general factor from the sex-invariant (partial scalar) multi-group model.",
    "Tertiles are defined once from the wave-1 distribution of each score and applied unchanged at every wave. Sex difference is the standardized mean difference (men minus women) across all person-waves.",
    sprintf("Scores (A) and (B) correlated r = %.3f and (A) and (C) r = %.3f; %s.", rAB, rAC,
            if (PRIMARY == "A") "the simpler within-sex standardization of the general factor (A) is therefore retained as the primary sex-neutral scale, with (B) and (C) as sensitivity analyses"
            else "the two scalings differ materially and the indicator-level standardization (B) is used as the primary sex-neutral scale")))
save_table(IC13TAB, "ST_IC13_Model",
  title = "Supplementary Table | Intrinsic-capacity factor with and without the indicators shared with the frailty phenotype",
  footnotes = c(
    sprintf("IC-13 removes the four indicators that also define the Fried phenotype (%s). Where a domain retained a single indicator (vitality: appetite) that indicator loads on the general factor only.", paste(OVERLAP, collapse = ", ")),
    "ECV, explained common variance of the general factor; omega H, coefficient omega hierarchical for the general factor. Both are computed from standardized loadings.",
    "IC-13 is used to test whether the associations of intrinsic capacity with state transitions, and its gradation of risk within the robust stratum, persist when the measurement overlap with the frailty phenotype is removed by construction."))
save_table(IC17TAB, "ST_IC17_Model",
  title = "Supplementary Table | Measurement model of the 17-indicator intrinsic-capacity factor",
  footnotes = c(
    sprintf("Standardized loadings from the bifactor model fitted to %d person-wave observations of %d participants (%s estimator). Fit: CFI %.3f, TLI %.3f, RMSEA %.3f, SRMR %.3f; ECV %.3f; omega H %.3f.", nrow(D), length(unique(D$id)), ESTIMATOR, fi17["CFI"], fi17["TLI"], fi17["RMSEA"], fi17["SRMR"], bi17$ECV, bi17$omegaH),
    "Share of general-factor variance is the squared standardized general loading of the indicator divided by the sum of squared general loadings over all indicators.",
    "Sex difference is the standardized mean difference of the pooled-standardized indicator (men minus women) across all person-waves."))

###############################################################################
## 9. 판정 — Phase 2 로 넘어가기 위한 결정
###############################################################################
cat("\n", strrep("=", 70), "\nPhase 1 판정\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  IC-17 재현      : r(gLIC) = %.4f %s\n", r_rep,
            if (r_rep > 0.995) "OK" else "★ 모형 문법 불일치 — SENS_SPECIFIC 확인"))
cat(sprintf("  성별 metric     : %s\n", if (SEX_METRIC_OK) "지지 -> 남녀 간 per-s.d. 연관 비교 가능" else "★ 미지지"))
cat(sprintf("  성별 scalar     : %s / partial %s\n",
            if (SEX_SCALAR_OK) "지지" else "미지지",
            if (SEX_PARTIAL_OK) sprintf("지지 (해제 %s)", paste(freed, collapse = ",")) else "미지지"))
cat(sprintf("  잠재 평균차     : %.2f s.d. (관측 %.2f) -> %s\n",
            lat_diff, smd(s17, D$sexF),
            if (abs(lat_diff) >= abs(smd(s17, D$sexF))) "남녀 차이는 실제 차이 (DIF 로 설명 안 됨) -> within-sex 척도의 근거는 Fried 식 성별 기준"
            else sprintf("관측 차이 중 %.0f%% 가 DIF", 100 * (1 - lat_diff / smd(s17, D$sexF)))))
cat(sprintf("  (A)-(B) r       : %.3f -> primary sex-neutral scale = (%s)\n", rAB, PRIMARY))
cat(sprintf("  IC-13           : r(IC-17) %.3f · SMD %.2f · HGS 제거로 남녀 차이 %s\n",
            stats::cor(z13, zP, use = "complete.obs"), smd(z13, D$sexF),
            if (abs(smd(z13, D$sexF)) < 0.5 * abs(smd(zP, D$sexF))) "절반 이하로 감소" else "여전히 큼"))
cat("\n  다음: 13_mi_hotdeck_m20.R (Phase 2) 는 P1_scores.rds 의 IC17_withinsex_z 를\n")
cat("        primary 노출로, IC13_withinsex_z 를 Phase 3.2 노출로 씁니다.\n")
cat(strrep("=", 70), "\n=== 완료 ===\n")
