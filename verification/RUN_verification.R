# =============================================================================
# verification/RUN_verification.R
#
# 본 파이프라인(99_MASTER_run_all.R)이 끝난 뒤에 실행합니다.
# d 와 iv 가 워크스페이스에 있어야 합니다.
#
# 이 폴더의 목적은 새 추정치를 만드는 것이 아니라, 분석 프레임이 논문에 보고된
# 값과 실제로 일치하는지 확인하는 것입니다. 개정 과정에서 워크스페이스에 남아
# 있던 중간 프레임이 코호트와 쇠약 유병률(robust 1,082 대 1,349)과 사망자 수
# (556 대 500)에서 어긋난 채로 추정에 쓰인 일이 있었고, 그 결과는 재현되지
# 않았습니다. 아래 검사는 그 유형의 오류를 추정 이전에 잡아냅니다.
#
# 기준값: n 3,011 / robust 1,349 / pre-frail 1,416 / frail 246 / 사망 500 /
#         구간 11,571 / person-years 22,730 / ADL 전이 969·141·361·651·74·97·21·30·42
# =============================================================================
## 자기 위치를 스스로 찾습니다 — 저장소 최상위에서 실행해도, verification/ 안에서
## 실행해도 동작합니다.
.vdir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- try(sys.frame(i)$ofile, silent = TRUE)
    if (!inherits(of, "try-error") && !is.null(of) && nzchar(of)) {
      d <- dirname(normalizePath(of, winslash = "/", mustWork = FALSE)); break }
  }
  if (is.null(d)) d <- if (dir.exists("verification/01_build")) "verification" else "."
  d
})

STEPS <- c(
  "01_build/01_rebuild_wave1_frame.R",
  "01_build/02_reconstruct_followup_time.R",
  "01_build/03_lock_death_time_rule.R",
  "01_build/04_lock_adl_state_variable.R",
  "02_verify/01_reconcile_manuscript_benchmarks.R",
  "02_verify/02_adjudicate_frailty_definition.R",
  "02_verify/03_sex_gap_decomposition.R",
  "03_analysis/01_robust_stratum_sexstd.R",
  "03_analysis/02_tables_1_2_3.R",
  "03_analysis/03_adl_transitions_recovery.R",
  "03_analysis/04_standardized_absolute_risks.R",
  "03_analysis/05_tertile_contrasts_and_recovery.R")
for (s in STEPS) {
  message("\n==== ", s, " ====")
  p <- file.path(.vdir, s)
  if (!file.exists(p)) { warning("missing: ", p); next }
  try(source(p, encoding = "UTF-8"), silent = FALSE)
}
message("\n검증 완료. 각 단계의 '<<< 불일치' 표시를 확인하십시오.")
