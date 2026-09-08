# legacy_single_imputation — superseded scripts, kept for provenance

These scripts produced the pre-revision analysis (single kNN-imputed data set,
capacity standardized on the pooled wave-1 distribution). The submitted manuscript
reports the multiple-imputation, within-sex analysis in `01_invariance/`,
`02_tables/`, `03_figures/` and `04_sensitivity/*_MI.R`; nothing in the paper is
produced from this folder, and `99_MASTER_run_all.R` does not run it.

They are kept because the sensitivity programme in `04_sensitivity/50–61` still runs
on the same kNN data set, and because the revision history (pooled-scale → within-sex,
single imputation → m = 20) is easier to audit with the earlier versions in place.

| Script | Superseded by |
|---|---|
| `20_table1_baseline.R` | `02_tables/20_table1_baseline_MI.R` |
| `21_table2_transitions.R`, `22_table3_frailty_ic.R` | `01_invariance/14_mi_pool_tables.R` |
| `23_table4_attenuation.R` | no direct MI counterpart — the before/after-frailty-adjustment contrast is no longer a stand-alone table in the submitted version [AUTHOR ACTION: confirm] |
| `30_figure1_transition_rates.R` | `03_figures/30_figure1_transition_rates_MI.R` |
| `31_figure2_dynamic_prediction.R`, `32_figure2_alternative.R` | `03_figures/32_figure4_JM_MI.R` |
| `33_figure3_capacity_stratification.R` | `03_figures/33_figure3_strata_MI.R` |
| `43_suppfig4_subgroups.R` | `03_figures/43_edfig1_subgroups_MI.R` |
| `57_sens_wave_stability.R` | `04_sensitivity/57_edfig2_wave_stability_MI.R` |
| `18_longitudinal_cfa_wide.R` (85-variable five-wave CFA; Heywood cases, > 90 min per model) | `01_invariance/18b_longitudinal_cfa_pairs.R` |

These scripts still use the shared helpers in `R/` and run if sourced individually.
