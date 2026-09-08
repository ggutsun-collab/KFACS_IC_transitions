# verification/logs

Console output of the benchmark-verification runs (`verification/01_build`, `02_verify`,
`03_analysis`) on the single-imputation working frame. They document that the frame
reproduces the design-fixed counts (n = 3,011; robust 1,349 / pre-frail 1,416 / frail 246;
500 deaths; 11,571 intervals; 22,730 person-years; the nine ADL transition counts) and how
the death-time rule (R1) and the disability-state variable (`state_lab`) were locked.

| Log | Script | What it shows |
|---|---|---|
| `rebuild_log.txt` | `01_build/01_rebuild_wave1_frame.R` | All five headline counts reproduced |
| `find_state_log.txt` | `01_build/04_lock_adl_state_variable.R` | `state_lab` reproduces 9/9 ADL transition counts |
| `death_time_log.txt` | `01_build/03_lock_death_time_rule.R` | Rule R1 reproduces the reported deaths |
| `sexgap_log.txt` | `02_verify/03_sex_gap_decomposition.R` | Indicator-level sex differences behind the within-sex scale |
| `TABLES_log.txt`, `FINAL_log.txt` | `03_analysis/` | Pooled-scale versus within-sex contrasts on the working frame |

Estimates printed in `TABLES_log.txt` and `FINAL_log.txt` come from the single-imputation
frame on which the verification was run; the values reported in the paper are from the
20 multiply imputed datasets (`01_invariance/14_`, `03_figures/`) and differ accordingly.
Logs of exploratory or failed diagnostic runs are not kept.
