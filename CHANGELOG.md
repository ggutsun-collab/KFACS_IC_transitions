# Changelog

## v1.1.0 — 2026-11-17 (revised manuscript)

Manuscript retitled *Intrinsic capacity as a marker of resilience across the frailty
spectrum*; the analysis is focused on two questions (deterioration and recovery; outcomes
within frailty states).

- **Baseline covariates imputed.** Education, household income and area of residence
  (missing for 2, 192 and 16 participants) are imputed once within each of the 20
  completed datasets (`01_invariance/13b_`) and applied at load time
  (`R/04_cov_imputation.R`). Previously they entered the models as an `Unknown` level.
  `04_sensitivity/64_` compares the two codings (Supplementary Table 6, missing-value category column).
- **New scripts:** `13b_`, `21_` (from the former `19_coauthor_additions.R`, without the
  horizon tables), `34_` (frail recovery by number of criteria), `35_` (decline after
  hospitalization or a fall), `62_` (additional confounders, transitions and within-state
  models in one script), `64_`; shared code in `R/05_functions_state_models.R`.
- **`R/00_setup.R`:** new input `RAW_FILE` (full pre-imputation master; read only by `35_`).
- **Driver:** run list reduced to the display items of the revised paper.
- **Moved to `archive_not_reported/`:** `10_`, `11_`, `15_`, `16_`, `17_`, `32_`, `33_`,
  `60b_`, `50_`–`61_`.
- `.gitignore` included in the repository.

## v1.0.0 — 2026-09-08

Code as first submitted (*Intrinsic capacity grades risk and recovery within frailty
phenotype categories*). Zenodo DOI 10.5281/zenodo.22654943.

### 1.1.0 (2026-11-17, fix c)
- Education covariate: `edu_high_bl` holds the original 8-level code; `prep_long()` now
  defines `edu_bl_f` (">= middle school") as categories 5-8. Previously only category 1
  (4.6% of participants) was coded as the upper level, so all adjusted models and the
  education subgroup panel used a mislabelled contrast. All analyses rerun.
- 21_: Supplementary Table 1 education row uses the same definition.
- 34_: death counted as a competing event; cumulative-incidence formatting fixed.
- 35_: sensitivity analysis counting deaths before the closing visit as decline.
