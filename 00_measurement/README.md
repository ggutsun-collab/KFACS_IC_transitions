# 00_measurement — 대체 전 master 와 단일대체(kNN) 세트의 구축

한 번만 실행하는 단계입니다. `99_MASTER_run_all.R` 에는 들어 있지 않습니다.
원자료 폴더를 지정한 뒤 `00b` → `00a` 순으로 실행합니다.

```r
Sys.setenv(KFACS_ROOT    = "D:/path/to/Nat_Aging")        # data/ 가 있는 프로젝트 루트 (00b)
Sys.setenv(KFACS_RAW_DIR = "D:/path/to/KF_DATA/build")     # 원자료 (00a)
Sys.setenv(KFACS_IMP_OUT = "D:/path/to/imputation_out")    # 00a 출력
source("00_measurement/00b_add_state_to_preimp_260618.R")
source("00_measurement/00a_impute_pipeline_260617.R")
```

| 파일 | 하는 일 | 산출물 |
|---|---|---|
| `00b_add_state_to_preimp_260618.R` | 원본 추적조사 파일의 FUP 코드를 pre-imputation master 에 붙이고 4-상태 변수를 만듭니다. Death 는 `death_wave` 에서만, Severe 는 ADL ≥ 3 또는 시설·병원 입소(FUP 5·6). | `data/KFACS_master_FINAL_preimput_260618_state.xlsx` — **주 분석 MI(13_)와 측정모형(10_, 15_, 18b_)의 입력** |
| `00a_impute_pipeline_260617.R` | 단일대체 4-세트(MAIN / Seq / NoFrailty / MNAR)를 만들고 bifactor 를 적합해 `gLIC` 를 산출합니다. IADL 을 먼저 대체한 뒤 ADL 을 대체하는 순서가 핵심입니다. | `KFACS_imputed_4sets_FINAL_260617` — 12_(성별 불변성)와 04_sensitivity/50–61 의 입력 |

## 주 분석의 대체는 여기가 아니라 `01_invariance/13_mi_pmm_m20.R` 입니다

투고본(rev9)의 모든 표·그림은 mice PMM m = 20 (wide format) 세트
`IMP_DIR/KFACS_mi_pmm_m20.rds` 를 씁니다. 그 세트는 `00b` 의 산출물을 입력으로
`13_` 이 만듭니다. `00a` 의 kNN 세트는 성별 불변성 검정의 기준 점수와
민감도 분석(원고에 명시)에만 쓰입니다.

## 측정모형

일반요인 하나와 특수요인 네 개(loco, vita, cogn, psyc)이며, 감각 2지표
(`rev_logMAR`, `rev_PTA`)는 일반요인에만 적재합니다. 분석에 쓰는 값은
일반요인 점수 `gLIC` 하나뿐이고, within-sex z 로 척도화합니다(12_ 의 근거).

## [AUTHOR ACTION] — 공개 전 반드시 해결

1. `R/00_setup.R` 은 `KFACS_imputed_4sets_LONGITUDINAL_260730` 을 읽는데, 여기
   있는 `00a` 는 `..._FINAL_260617` 을 만듭니다. 260730 종단 재대체를 만든
   스크립트를 `00c_impute_longitudinal_260730.R` 로 추가하거나, 두 세트가 같은
   파이프라인의 재실행이라면 그 사실을 이 README 에 적으십시오.
2. `IMP_DIR/KFACS_gLIC_anchor_fiml.rds` 를 만드는 스크립트가 저장소에 없습니다.
   현재 `13_` 은 MI 세트마다 bifactor 를 다시 적합해 점수를 냅니다(complete cases, MLR).
   원고 Methods 가 "단일 FIML anchor 로 채점" 이라면 그 버전의 `13_` 을 넣고,
   아니라면 Methods 를 현재 코드에 맞추십시오. 둘이 다르면 심사에서 지적됩니다.
