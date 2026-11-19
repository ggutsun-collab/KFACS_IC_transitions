# MANIFEST (v1.1.0)

Every script in the reported pipeline, in run order, with the display item it produces.
Paths are relative to the repository root. Numeric prefixes are stable. Every analysis
script begins with the same block: it locates `R/_bootstrap.R` and calls `kf_init()`,
which loads the helpers once. All paths are set in `R/00_setup.R` and nowhere else.

## Helpers (`R/`, loaded once by `kf_init()`)

| File | Lines | Purpose |
|---|---|---|
| `R/_bootstrap.R` | 119 | Locates the helper folder from the calling script's location and loads `00`–`05`. |
| `R/00_setup.R` | 89 | Paths and run options: `IMP_DIR`, `IMP_STEM` (kNN set), `MI_STEM` (mice set), `PRE_FILE`, `RAW_FILE`, `FIG_DIR`. |
| `R/01_functions_common.R` | 853 | `prep_long()`, interval construction, transition-specific Poisson models with participant-clustered SE, Rubin pooling, `CFG`. |
| `R/02_theme_and_output.R` | 536 | Figure standards, colour scheme, save helpers, `load_imputed(set, stem)` (applies `04` to MI sets). |
| `R/03_functions_incremental.R` | 153 | Start-of-interval covariate attachment, adjustment-set sweeps, log-scale attenuation. |
| `R/04_cov_imputation.R` | 55 | Replaces missing education, income and area in each MI set with the values from `13b_`; `options(kf.cov_impute = FALSE)` switches it off. |
| `R/05_functions_state_models.R` | 144 | Wave-1 person frame (5-year death and composite worsening), 18-transition intervals, within-state Cox per s.d., Rubin's rules. Reproduces Table 2 (11,571 intervals) and the Table 3 per-s.d. estimates. |

## Stage 0 — data construction (`00_measurement/`, run once, not in the driver)

| Script | Produces |
|---|---|
| `00b_add_state_to_preimp_260618.R` | Pre-imputation longitudinal master with the four-state disability variable (`PRE_FILE`). Input to 13_ and 18b_. |
| `00a_impute_pipeline_260617.R` | Single-imputation (kNN) four-set file. Input to 12_. |

## Stage A — measurement model and imputation

| Script | Lines | Produces |
|---|---|---|
| `01_invariance/12_sex_invariance_and_scores.R` | 420 | Supplementary Table 2 (bifactor model, 17 indicators) and Supplementary Table 3 (invariance between men and women); within-sex scoring rule. |
| `01_invariance/13_mi_pmm_m20.R` | 236 | mice PMM, m = 20, on the indicators, ADL/IADL and frailty items; bifactor scoring per completed set. |
| `01_invariance/13b_impute_baseline_covariates.R` | 121 | Education, household income and area of residence imputed once within each completed set (`KFACS_cov_imputed_m20.rds`). |
| `01_invariance/18b_longitudinal_cfa_pairs.R` | 226 | Adjacent-wave longitudinal CFA. |
| `01_invariance/18c_longitudinal_cfa_tables.R` | 102 | Supplementary Table 4 (longitudinal invariance). |

## Stage B — main tables

| Script | Lines | Produces |
|---|---|---|
| `02_tables/20_table1_baseline_MI.R` | 234 | Table 1. |
| `01_invariance/14_mi_pool_tables.R` | 272 | Table 2 (18 transitions), Table 3 (frailty state x IC tertile cells; robust per-s.d.). |
| `02_tables/21_st1_bysex_and_st6ab_sex_MI.R` | 308 | Supplementary Table 1 (baseline by sex), Supplementary Table 6, + Sex column (sex as covariate). |

## Stage C — figures and recovery

| Script | Lines | Produces |
|---|---|---|
| `03_figures/30_figure1_transition_rates_MI.R` | 619 | Figure 1. |
| `03_figures/31_figure2_recovery_MI.R` | 263 | Table 4 (pre-frail to robust; mild disability to normal) and Supplementary Fig. 1 (frail median split). The former Figure 2 it also draws is not in the paper. |
| `03_figures/34_frail_recovery_by_criteria_MI.R` | 105 | Table 4 (frail row) and Supplementary Table 5 (by number of frailty criteria). |
| `03_figures/35_edfig1_hospitalization_fall_MI.R` | 137 | Extended Data Fig. 1 and Supplementary Table 7 (decline after hospitalization or a fall). |
| `03_figures/43_edfig1_subgroups_MI.R` | 485 | Supplementary Fig. 2 (subgroups, set MI01). |
| `04_sensitivity/57_edfig2_wave_stability_MI.R` | 358 | Supplementary Fig. 3 (study waves, set MI01). |

## Stage D — sensitivity

| Script | Lines | Produces |
|---|---|---|
| `04_sensitivity/62_sens_additional_confounders.R` | 68 | Supplementary Table 6, + Four covariates column (living alone, polypharmacy, smoking, physical activity added). Its primary column gives the Table 3 per-s.d. rows for the robust and pre-frail states. |
| `04_sensitivity/64_sens_covariate_missing_category.R` | 65 | Supplementary Table 6, Missing-value category column (imputed covariates versus missing-category coding). |

## Stage E — reporting

| Script | Produces |
|---|---|
| `05_reporting/90_collect_captions.R` | Assembles the per-figure caption files. |
| `05_reporting/91_sf_cleanup.R` | Output-folder cleanup (not run by the driver). |
| `99_MASTER_run_all.R` | Driver: runs A to E, skips 12_/13_/13b_ when their outputs exist, continues past failures, logs the console. |

## `archive_not_reported/`

Scripts of v1.0.0 for analyses the revised paper no longer reports, kept unchanged so the
earlier results can be regenerated: `10_`, `11_` (wave multi-group invariance),
`15_` (IC-13), `16_` (msm), `17_` (level versus change), `32_` (joint model),
`33_` (former Figure 3), `60b_` (death-time bracketing) and the kNN sensitivity
programme `50_`–`61_`. They are not run by the driver.

## `verification/` and `legacy_single_imputation/`

Unchanged from v1.0.0; see README §4.
