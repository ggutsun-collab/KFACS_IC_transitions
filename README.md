# Intrinsic capacity as a marker of resilience across the frailty spectrum — KFACS

Analysis code for *Intrinsic capacity as a marker of resilience across the frailty
spectrum* (Korean Frailty and Aging Cohort Study, n = 3,011, five biennial waves,
500 deaths).

**Archive DOI:** https://doi.org/10.5281/zenodo.24358027 (v1.1.0, revised manuscript);
v1.0.0, the code as first submitted, is https://doi.org/10.5281/zenodo.22654943 (see `CHANGELOG.md`)
**Repository:** https://github.com/ggutsun-collab/KFACS_IC_transitions

---

## 1. Run it

```r
# R >= 4.3.  Open the repository root as the working directory, point KFACS_ROOT at the
# project folder that contains data/ and result/, then source one file.
setwd("path/to/KFACS_IC_transitions")
Sys.setenv(KFACS_ROOT = "D:/path/to/Nat_Aging")
source("99_MASTER_run_all.R")
```

`99_MASTER_run_all.R` runs the whole pipeline in manuscript order, continues past a
failing script, writes the full console to `RUN_log_<timestamp>.txt`, and ends with a
summary table. A script that finishes without error but produces no output file is
flagged as suspect, so a silent failure cannot pass unnoticed.

Three steps are run once and skipped automatically when their output already exists: the
sex-invariance scoring (`12_`, writes `P1_scores.rds`), the multiple imputation
(`13_`, writes `KFACS_mi_pmm_m20.rds`, 30–90 min) and the imputation of the baseline
socioeconomic covariates (`13b_`, writes `KFACS_cov_imputed_m20.rds`). Set
`FORCE_REBUILD <- TRUE` to redo them. The remaining run takes about 1–2 hours. To skip
any script add it to `SKIP`.

Any script can also be run on its own. Each begins with the same block, which finds
`R/_bootstrap.R` from the script's own location and calls `kf_init()`; that loads the four
helper files once and does nothing if they are already loaded. Nothing else needs to be
sourced by hand.

Paths are resolved in `R/00_setup.R` and nowhere else, in this order: the environment
variable `KFACS_ROOT`; a `data/` directory found at `.`, `..` or `../..`; and finally
`PROJECT_ROOT_FALLBACK` at the top of that file. `KFACS_IMP_DIR`, `KFACS_IMP_STEM`,
`KFACS_MI_STEM`, `KFACS_PRE_FILE`, `KFACS_RAW_FILE` (full pre-imputation master with the
hospitalization and fall items, read only by `35_`) and `KFACS_OUT_DIR` override the
individual inputs and the output folder if the layout differs.

## 2. Layout

```
99_MASTER_run_all.R         driver — the only script at the root

R/                          helpers, loaded once by kf_init()
  _bootstrap.R              locates the helper folder and loads 00-05
  00_setup.R                paths — the only file with a path in it
  01_functions_common.R     data preparation, interval construction, transition models
  02_theme_and_output.R     figure standards, colour scheme, save helpers, load_imputed()
  03_functions_incremental.R  covariate attachment and attenuation helpers
  04_cov_imputation.R       applies the imputed baseline covariates when an MI set is loaded
  05_functions_state_models.R person frame, intervals, within-state Cox, Rubin (13b_, 34_, 35_, 62_, 64_)

00_measurement/             00b_ pre-imputation master; 00a_ single-imputation (kNN) set   (run first, once)
01_invariance/              12_ sex invariance and within-sex scoring; 13_ mice PMM m = 20;
                            13b_ baseline covariate imputation; 18b_ 18c_ longitudinal invariance;
                            14_ Tables 2-3
02_tables/                  20_ Table 1; 21_ Supplementary Table 1 and Supplementary Table 6 (+ Sex column)
03_figures/                 30_ Figure 1; 31_ Table 4 and Supplementary Fig. 1;
                            34_ Table 4 (frail row) and Supplementary Table 5;
                            35_ Extended Data Fig. 1 and Supplementary Table 7;
                            43_ Supplementary Fig. 2
04_sensitivity/             57_ Supplementary Fig. 3; 62_ Supplementary Table 6 (+ Four covariates column); 64_ Supplementary Table 6 (Missing-value category column)
05_reporting/               90_ caption collection; 91_ output-folder cleanup

verification/               benchmark verification (see §4)
archive_not_reported/       analyses of v1.0.0 that the revised paper no longer reports
                            (IC-13, msm, level vs change, joint model, Figure 3, kNN sensitivity programme)
legacy_single_imputation/   pre-revision scripts, not part of the reported analysis
```

`MANIFEST.md` lists every script with the display item it produces.
Numeric prefixes are stable, so a script referred to as "script 14" is
`01_invariance/14_mi_pool_tables.R`.

## 3. Design rules carried in the code

These are recorded here because each one was learned from an error that reached a
result before it was caught.

**The primary analysis is multiply imputed and the capacity scale is within-sex.**
`13_` imputes the 17 indicators, ADL/IADL and frailty items with mice PMM (m = 20, wide
format). Education, household income and residential area (missing for 2, 192 and 16
participants) are imputed once within each completed dataset by `13b_` and applied by
`load_imputed()` through `R/04_cov_imputation.R`, so every table and figure uses the same
imputed values. The general factor is standardized within sex because the
sex difference in the latent mean is about 1.5 s.d. and multi-group invariance holds
(`12_`); the analogy is the sex-specific cut-points of the Fried phenotype. Every table
and figure in the paper pools 20 completed sets by Rubin's rules unless the caption says
otherwise (the joint model, the subgroup figure and the wave-stability figure use set
MI01; msm uses `N_SETS = 5`).

**Death is an absorbing state, and the death time is not on the visit grid.**
`prep_long()` truncates imputed rows after death, blanks the capacity measurements on the
death row, and places the death at `followup_years` rather than at the nominal visit
time. Fourteen participants have a grid time at or after their death time; without this
step the `duration > 0` filter silently drops those deaths.

**Only the general factor is used.** The bifactor model has one general factor and four
specific factors (loco, vita, cogn, psyc), with the two sensory indicators loading on the
general factor only. `gLIC` is the sole exposure; domain scores are not used, and the
invariance tests (`10_`, `12_`, `18b_`) use the same specification as the scoring model.

**Indicators are standardized on the pooled person-wave distribution, and the tertile
cut-points are fixed at wave 1.** A rolling-landmark analysis therefore compares waves
with the same ruler; standardizing within wave would remove between-wave change by
construction.

**Time-varying comorbidity is not used in the adjustment set.** Its missingness is
correlated with death: a complete-case analysis using the time-varying value loses 146
of 500 deaths (29%). Age is reconstructed deterministically and comorbidity is fixed at
its wave-1 value.

**Missing socioeconomic values are imputed, never folded into the reference.** Folding
them in asserts that non-responders resemble the least advantaged group. v1.0.0 carried
them as an explicit `Unknown` level; v1.1.0 imputes them (`13b_`), and `64_` shows that
the two codings change the reported log ratios by at most 1.8%. `options(kf.cov_impute =
FALSE)` restores the `Unknown` coding for a run.

**Sparse cells are reported as not estimable rather than estimated.** A tertile
contributing fewer than five events makes the contrast diverge, and a single event
produced an IRR of 9.42 (95% CI 1.04–85.69) that would have been printed as a finding.
`CFG$MIN_CELL_EVENTS` governs this; Figure 1 shows such positions as open markers.

**Attenuation is computed on the log scale.** `(IRR_adj − 1)/(IRR_base − 1)` explodes
when the base IRR is near 1; it produced −1041% and +1948% in one run. Attenuation is
also suppressed when the base association is not significant.

**Colour is defined once, in `02_theme_and_output.R`.** Every figure draws from the same
five semantic anchors and the colour footnote is generated from the same object.

**Figures carry no title and no caption inside the image.** Legends live in the
manuscript; `save_na()` writes each caption to a separate `_caption.txt`, and
`90_collect_captions.R` assembles them.

## 4. `verification/` — benchmark checks

Every frame-construction step there ends in an explicit comparison against a count that
is fixed by the study design and reported in the paper, and the analysis refuses to
proceed if a comparison fails.

| Quantity | Value |
|---|---|
| Participants | 3,011 |
| Wave-1 frailty phenotype | robust 1,349 / pre-frail 1,416 / frail 246 |
| Deaths | 500 |
| Between-visit intervals | 11,571 |
| At-risk person-years in the transition models | 22,730 (total follow-up 22,785) |
| ADL transition counts | 969, 141, 361, 651, 74, 97, 21, 30, 42 |

This exists because a derived frame that had been carried forward in a working session
disagreed with the cohort on frailty prevalence (robust 1,082 versus 1,349) and on deaths
(556 versus 500) and produced estimates that did not reproduce. Two rules were recovered
in the process and are locked in `verification/01_build/`:

- the death time is the `time` value on the row where `wave` is missing; `death_event`
  is a person-level constant replicated across every person-wave row and cannot be used
  to locate the event;
- the disability state is `state_lab`, not `adl_3cat` — the latter carries 2,625 missing
  person-wave values and does not reproduce the transition counts.

Run it after the main pipeline, with `d` and `iv` in the workspace:

```r
source("verification/RUN_verification.R")
```

`verification/logs/` holds the console output of the run that produced the reported
numbers, so the checks can be inspected without access to the restricted data.

## 5. Data availability

KFACS data contain identifiable information about participants and are available from
the KFACS Steering Committee (contact: Miji Kim, mijiak@khu.ac.kr) on reasonable request,
subject to approval by the Institutional Review Board of Kyung Hee University Medical
Center. Requests are typically answered within 8 weeks. Source data for Fig. 1 and
Extended Data Fig. 1 are provided with the paper.

`.gitignore` blocks `.rds`, `.RData`, `.csv`, `.xlsx` and other data formats. Confirm
with `git status` before the first push that no participant-level file is staged.

## 6. Known gaps

- **kNN set provenance.** `R/00_setup.R` reads `KFACS_imputed_4sets_LONGITUDINAL_260730`,
  the longitudinal re-run of the single-imputation pipeline; `00a_impute_pipeline_260617.R`
  is the version of that pipeline retained here. In v1.1.0 the kNN set is an input only to
  the sex-invariance step (`12_`); the multiply imputed sets used for every reported
  estimate are produced from the pre-imputation master by `13_`.
- **Rebuilt scripts.** The scripts for Supplementary Table 5 (`34_`) and for Extended Data
  Fig. 1 / Supplementary Table 7 (`35_`) were rewritten for v1.1.0 from the definitions in
  the Methods, because the working versions used for the first submission were not kept.
  Each prints the earlier published values as a reproduction check.
- `legacy_single_imputation/23_table4_attenuation.R` has no stand-alone counterpart in
  the submitted version.
- `SESSION_INFO.txt` (package versions of the reported run) is added with the first
  release; versions are also stated in the Methods.
