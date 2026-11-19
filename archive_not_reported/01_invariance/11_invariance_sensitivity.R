###############################################################################
## KF_ST_InvarSensitivity_v260801.R
## Supplementary Data | 측정불변성 민감도 —
##   부분 scalar 불변모형에서 뽑은 gLIC 로 주요 전이분석을 통째로 재실행하고,
##   본 분석 결과와 나란히 놓습니다.
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
##  KF_ST_Invariance 실행 결과(2026-07-31):
##    metric  지지     (ΔCFI -0.005)  -> "per +1 SD" 비교는 정당
##    scalar  미지지   (ΔCFI -0.014)  -> "수준" 비교의 전제가 형식적으로 깨짐
##    그러나 두 점수의 상관 r = 0.997, 삼분위 일치율 95.6%
##
##  즉 통계적 검정은 깨졌지만 실질적 영향은 없어 보입니다. 리뷰어에게
##  "상관이 높으니 괜찮다"고 **주장**하는 것과, 실제로 다시 돌려서
##  "추정치가 바뀌지 않는다"고 **보여주는** 것은 무게가 다릅니다.
##  이 스크립트가 후자를 만듭니다.
##
## ── 하는 일 ──────────────────────────────────────────────────────────────
##  1) ST_Invariance_gLIC_invariant.csv 를 읽어 raw$gLIC 를 교체
##  2) prep_long 을 다시 태워 gLIC_z / W1 삼분위를 새 점수 기준으로 재산출
##  3) Table 2 / Table 4(M0, M1) 의 전이별 IRR 을 재추정
##  4) 원본과 나란히 놓고 최대·중앙 차이를 보고
##
## ※ 12번(KF_ST_Invariance_v260801.R)을 먼저 돌려야 합니다.
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
## ── 불변모형 점수 읽기 ──────────────────────────────────────────────────
INV_FILE <- file.path(FIG_DIR, "ST_Invariance_gLIC_invariant.csv")
if (!file.exists(INV_FILE))
  stop("불변모형 점수 파일이 없습니다: ", INV_FILE,
       "\n  -> KF_ST_Invariance_v260801.R 을 먼저 실행하십시오.")
INV <- utils::read.csv(INV_FILE, stringsAsFactors = FALSE)
if (!all(c("id", "wave", "gLIC_invariant") %in% names(INV)))
  stop("불변모형 점수 파일의 열 구성이 예상과 다릅니다: ",
       paste(names(INV), collapse = ", "))

load_imputed("MAIN")
raw <- load_kfacs()

k  <- match(paste(raw$id, as.integer(as_num(raw$wave))),
            paste(INV$id, as.integer(as_num(INV$wave))))
hit <- mean(!is.na(k))
msg(sprintf("불변모형 점수 결합률 %.1f%% (%d / %d행)", 100 * hit, sum(!is.na(k)), nrow(raw)))
if (hit < 0.5)
  stop("결합률이 너무 낮습니다. id / wave 형식을 확인하십시오.")

raw_inv <- raw
raw_inv$gLIC <- INV$gLIC_invariant[k]      # 미결합 행(사망행 등)은 NA -> 노출로 안 쓰임

## ── 두 버전의 long 자료 ─────────────────────────────────────────────────
d_org <- prep_long(raw)
d_inv <- prep_long(raw_inv)

## 재산출된 삼분위가 얼마나 일치하는가 (본문에 쓸 수치)
kk  <- match(paste(d_inv$id, round(d_inv$time, 3)),
             paste(d_org$id, round(d_org$time, 3)))
agr <- mean(as.character(d_inv$gLIC_tert_W1cut) ==
            as.character(d_org$gLIC_tert_W1cut[kk]), na.rm = TRUE)
rz  <- suppressWarnings(stats::cor(d_inv$gLIC_z, d_org$gLIC_z[kk], use = "complete.obs"))
msg(sprintf("재산출 후: gLIC_z 상관 %.4f | W1 삼분위 일치율 %.1f%%", rz, 100 * agr))

## ── 전이분석 재실행 ─────────────────────────────────────────────────────
blocks <- list(c("state",        "ADL disability, 3 states"),
               c("frailty_3cat", "Frailty phenotype"),
               c("adl_state",    "Any ADL disability (binary)"))

grab_set <- function(dd, tag) {
  do.call(rbind, lapply(blocks, function(bk) {
    iv <- prep_iv_adj(dd, bk[1], "rolling")
    if (is.null(iv)) return(NULL)
    tl <- default_transitions(bk[1], iv)
    r  <- rbind(run_adj_set(iv, "M0", CFG$ADJ,               "gLIC_z", tl),
                run_adj_set(iv, "M1", c(CFG$ADJ, "frail_f"), "gLIC_z", tl))
    r$sysname <- bk[2]; r$src <- tag; r
  }))
}
A <- grab_set(d_org, "original")
B <- grab_set(d_inv, "invariant")
if (is.null(A) || is.null(B)) stop("전이분석 결과가 비어 있습니다.")

key <- function(x) paste(x$sysname, x$transition, x$model)
i <- match(key(A), key(B))
W <- data.frame(System = A$sysname, Transition = A$transition, Model = A$model,
                Events = A$n_events,
                irr_o = A$IRR, lcl_o = A$lcl, ucl_o = A$ucl, p_o = A$p,
                irr_i = B$IRR[i], lcl_i = B$lcl[i], ucl_i = B$ucl[i], p_i = B$p[i],
                stringsAsFactors = FALSE)
W$dlog <- log(W$irr_i) - log(W$irr_o)
W$dpct <- 100 * (exp(W$dlog) - 1)
## 결론이 뒤집혔는가 (유의성 방향 변화)
sg <- function(l, u) ifelse(is.finite(l) & is.finite(u), ifelse(l > 1, "+", ifelse(u < 1, "-", "0")), NA)
W$sig_o <- sg(W$lcl_o, W$ucl_o); W$sig_i <- sg(W$lcl_i, W$ucl_i)
W$flip  <- !is.na(W$sig_o) & !is.na(W$sig_i) & W$sig_o != W$sig_i

W <- W[order(match(W$System, vapply(blocks, `[`, character(1), 2)), W$Model), ]
save_vals(W, "ST_InvarSensitivity_values_long.csv", FIG_DIR)

TAB <- data.frame(
  System                    = W$System,
  Transition                = W$Transition,
  Model                     = W$Model,
  Events                    = W$Events,
  `Main analysis (95% CI)`  = fmt_irr(W$irr_o, W$lcl_o, W$ucl_o),
  `P`                       = fmt_p(W$p_o),
  `Invariance-based score (95% CI)` = fmt_irr(W$irr_i, W$lcl_i, W$ucl_i),
  `P `                      = fmt_p(W$p_i),
  `Difference in IRR`       = ifelse(is.finite(W$dpct), sprintf("%+.1f%%", W$dpct), "-"),
  `Conclusion changed`      = ifelse(is.na(W$flip), "-", ifelse(W$flip, "yes", "no")),
  check.names = FALSE, stringsAsFactors = FALSE)
TAB$System[duplicated(TAB$System)] <- ""

n_cmp  <- sum(is.finite(W$dpct))
mx     <- max(abs(W$dpct), na.rm = TRUE)
md     <- stats::median(abs(W$dpct), na.rm = TRUE)
n_flip <- sum(W$flip, na.rm = TRUE)

save_table(TAB, "ST_InvarSensitivity",
  title = "Supplementary Data | Sensitivity of the transition analyses to the measurement-invariance model",
  footnotes = c(
    "The intrinsic-capacity score was re-estimated from the partial scalar invariance model and substituted for the score used in the main analyses; the wave-1 standardisation and tertile cut-points were recomputed from the new score, and every transition model was refitted.",
    sprintf("The two scores correlated r = %.4f and agreed on the wave-1 tertile for %.1f%% of person-waves.", rz, 100 * agr),
    sprintf("Across %d transition-by-model estimates the median absolute change in the incidence rate ratio was %.1f%% and the largest was %.1f%%; %d estimate%s changed direction or statistical significance.",
            n_cmp, md, mx, n_flip, if (n_flip == 1) "" else "s"),
    "Models are as in the concurrent-frailty supplementary table: M0 adjusted for age, sex, baseline comorbidity count, education, household income and area of residence; M1 additionally for the frailty phenotype at the start of the same interval.",
    "CI, confidence interval; IRR, incidence rate ratio."))

cat("\n=== 측정불변성 민감도 완료 ===\n")
cat(sprintf("  gLIC_z 상관 %.4f | W1 삼분위 일치율 %.1f%%\n", rz, 100 * agr))
cat(sprintf("  비교 추정치 %d개 | IRR 변화 중앙값 %.1f%% | 최대 %.1f%%\n", n_cmp, md, mx))
if (n_flip == 0) {
  cat("  ★ 방향/유의성이 바뀐 추정치가 하나도 없습니다.\n")
  cat("     -> 측정 비불변성은 결론에 영향을 주지 않습니다. 본 분석 유지하고\n")
  cat("        이 표를 Supplementary 에 넣어 리뷰어 지적을 선제 차단하십시오.\n")
} else {
  cat(sprintf("  ※ %d개 추정치의 방향/유의성이 바뀌었습니다:\n", n_flip))
  f <- W[which(W$flip), ]
  for (i in seq_len(nrow(f)))
    cat(sprintf("      %-28s %s  %.2f -> %.2f\n", f$Transition[i], f$Model[i],
                f$irr_o[i], f$irr_i[i]))
  cat("     -> 해당 전이는 본문에서 신중하게 기술하십시오.\n")
}
