# MANIFEST

Every script in the reported pipeline, in run order, with the display item it produces.
Paths are relative to the repository root. Only `99_MASTER_run_all.R` and the
documentation files sit at the root. Numeric prefixes are stable, so a script referred to
as "script 14" in correspondence is `01_invariance/14_mi_pool_tables.R`.

Every analysis script begins with the same block: it locates `R/_bootstrap.R` and calls
`kf_init()`, which loads the four helpers once. All paths are set in `R/00_setup.R` and
nowhere else.

## Helpers (`R/`, loaded once by `kf_init()`)

| File | Lines | Purpose |
|---|---|---|
| `R/_bootstrap.R` | 117 | Locates the helper folder from the calling script's own location and loads `00_setup` – `03_functions_incremental`. |
| `R/00_setup.R` | 81 | Paths and run options — the only file containing a path. Resolves `KFACS_ROOT`; defines `IMP_DIR`, `IMP_STEM` (kNN set), `MI_STEM` (mice set), `PRE_FILE` (pre-imputation master), `FIG_DIR`. Checks the inputs exist before anything runs. |
| `R/01_functions_common.R` | 853 | Data preparation (`prep_long`), person-interval construction, transition-specific Poisson models with participant-clustered robust SE, Rubin pooling, save helpers, `CFG` (adjustment set, minimum cell-event rule). |
| `R/02_theme_and_output.R` | 533 | Figure standards for the Nature format, the single semantic colour scheme, table/figure save helpers, imputed-data loader (`load_imputed(set, stem)`). |
| `R/03_functions_incremental.R` | 153 | Incremental-value helpers: start-of-interval covariate attachment, adjustment-set sweeps, log-scale attenuation. |

## Stage 0 — data construction (`00_measurement/`, run once, not in the driver)

| Script | Lines | Produces |
|---|---|---|
| `00b_add_state_to_preimp_260618.R` | 102 | Pre-imputation longitudinal master with the four-state disability variable (`data/KFACS_master_FINAL_preimput_260618_state.xlsx`). Input to 13_, 10_, 15_, 18b_. |
| `00a_impute_pipeline_260617.R` | 234 | Single-imputation (kNN) four-set file with bifactor-scored `gLIC`. Input to 12_ and to the sensitivity programme 50–61. |

## Stage A — measurement model and multiple imputation

| Script | Lines | Produces |
|---|---|---|
| `01_invariance/12_sex_invariance_and_scores.R` | 420 | Sex multi-group bifactor: configural/metric/scalar, partial-scalar search, latent mean difference, within-sex scoring rule. Justifies the within-sex capacity scale used everywhere downstream. Writes `P1_scores.rds`. |
| `01_invariance/13_mi_pmm_m20.R` | 236 | mice PMM, m = 20, wide format, on the 17 indicators, ADL/IADL and frailty items; bifactor scoring per completed set. Writes `IMP_DIR/KFACS_mi_pmm_m20.rds` and a diagnostics workbook. |
| `01_invariance/10_measurement_invariance_v3.R` | 350 | Wave multi-group CFA on observed data with FIML (four specific factors, matching the scoring model). Supplementary table on measurement invariance. |
| `01_invariance/11_invariance_sensitivity.R` | 160 | Sensitivity of the reported associations to the invariance-based factor score. |
| `01_invariance/18b_longitudinal_cfa_pairs.R` | 226 | Adjacent-wave-pair longitudinal CFA (within-person dependence through factor and residual covariances); latent mean decline per interval. Writes `LongCFA_pairs_fits.rds`. |
| `01_invariance/18c_longitudinal_cfa_tables.R` | 102 | Supplementary table on longitudinal invariance from the 18b fits. |

## Stage B — main analysis (m = 20, within-sex scale)

| Script | Lines | Produces |
|---|---|---|
| `01_invariance/14_mi_pool_tables.R` | 272 | Table 2 (capacity and all observed transitions), Table 3 (five-year worsening, death, recovery by frailty stratum × capacity tertile), Extended Data Table 1; per-set values `MI_Table2_perset.csv`, `MI_Table3_perset.csv`. |
| `02_tables/20_table1_baseline_MI.R` | 234 | Table 1, baseline characteristics by within-sex capacity tertile (observed denominators for income, area, alcohol). |
| `01_invariance/15_ic13_headtohead.R` | 328 | IC-13 (Fried-overlap indicators removed) retention of IC-17 effects; within-robust head-to-head of IC-13 against grip strength + gait speed (ΔC, IDI, bootstrap CI). Writes the IC-13 anchor. |
| `01_invariance/16_msm_pmatrix.R` | 217 | Continuous-time Markov model (msm) on `N_SETS` completed sets; two-year transition probabilities at −1/0/+1 s.d. |
| `01_invariance/17_level_vs_change.R` | 210 | A: discrete level-versus-change decomposition (supplementary table). B: JMbayes2 value + slope — reported as non-identifiable in the paper; kept so the claim is reproducible. |

## Stage C — figures

| Script | Lines | Produces |
|---|---|---|
| `03_figures/30_figure1_transition_rates_MI.R` | 619 | Figure 1, adjusted transition rates at −1/0/+1 s.d., common covariate profile; open markers where a tertile had no events. |
| `03_figures/31_figure2_recovery_MI.R` | 263 | Figure 2, recovery (pre-frail → robust; mild disability → normal) by within-sex tertile; Supplementary Figure (frail median split). |
| `03_figures/33_figure3_strata_MI.R` | 208 | Figure 3, mortality and composite worsening across the frailty × capacity cross-classification (uses `MI_Table3_perset.csv` from 14_; recomputes if absent). |
| `03_figures/32_figure4_JM_MI.R` | 314 | Figure 4, dynamic survival predictions from the joint model (set MI01). The slowest step. |
| `03_figures/43_edfig1_subgroups_MI.R` | 485 | Extended Data Fig. 1, consistency of the capacity association across subgroups (set MI01). |
| `04_sensitivity/57_edfig2_wave_stability_MI.R` | 358 | Extended Data Fig. 2, stability across assessment waves (set MI01). |

## Stage D — sensitivity programme (single-imputation kNN set, stated in the paper)

| Script | Lines | Produces |
|---|---|---|
| `04_sensitivity/50_sens_lead_time.R` | 171 | Reverse causation: exclusion of each participant's first one or two intervals. |
| `04_sensitivity/51_sens_evalue.R` | 173 | E-values for the transition analyses. |
| `04_sensitivity/52_sens_exposure_definition.R` | 116 | Rolling landmark versus wave-1-fixed exposure. |
| `04_sensitivity/53_sens_discrimination_ci.R` | 220 | Bootstrap CI for the change in C-index and AUC. |
| `04_sensitivity/54_sens_msm_compare.R` | 261 | Continuous-time msm against the discrete-time analysis. |
| `04_sensitivity/55_sens_complete_case.R` | 209 | Complete-case analysis. |
| `04_sensitivity/56_sens_adjustment_reduced.R` | 201 | Reduced adjustment set. |
| `04_sensitivity/58_sens_mi_pooling.R` | 183 | Pooling across the four single-imputation sets (MAIN / Seq / NoFrailty / MNAR). |
| `04_sensitivity/60_sens_deathtime_ph.R` | 230 | Alternative interval-censored death-time assignment; proportional-hazards diagnostics. |
| `04_sensitivity/61_sens_observed_ic.R` | 143 | Restriction to intervals whose capacity rests entirely on observed indicators. |
| `04_sensitivity/59_sens_summary.R` | 213 | Summary verdict across 50–61. Runs last. |

## Stage E — reporting

| Script | Produces |
|---|---|
| `05_reporting/90_collect_captions.R` | Assembles the per-figure caption files into one document. |
| `05_reporting/91_sf_cleanup.R` | Output-folder cleanup. |
| `99_MASTER_run_all.R` | Driver. Runs A → E in manuscript order, skips 12_/13_ when their outputs already exist, continues past failures, logs the console, flags scripts that produce no output. |

## `verification/`

Frame construction with locked rules and benchmark checks that stop the run on mismatch
(participants 3,011; wave-1 frailty 1,349 / 1,416 / 246; deaths 500; intervals 11,571;
person-years 22,730; ADL transition counts). See README §4. Run separately with
`verification/RUN_verification.R`; console logs of the reported run are in `verification/logs/`.

## `legacy_single_imputation/`

Pre-revision versions (single kNN set, pooled scale) superseded by the `_MI` scripts above.
Not run by the driver; see the README there for the mapping.
