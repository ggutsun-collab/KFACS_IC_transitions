###############################################################################
## KF_SupplFig1_Subgroups_v260801.R
## Supplementary Figure 4 | 하위군 일관성 —
##   Table 4 의 증분효과(프레일티 보정 후 IC 의 전이별 IRR)가
##   성별·연령 하위군에서도 같은 방향으로 유지되는가.
##
## ── 왜 필요한가 ──────────────────────────────────────────────────────────
##  "IC 효과가 특정 하위군(예: 여성, 초고령)에만 있는 것 아니냐"는 지적은
##  거의 반드시 나옵니다. 층별로 나눠 보여주고 상호작용 검정을 함께 실으면
##  그 지적이 닫힙니다.
##
## ── 설계 ─────────────────────────────────────────────────────────────────
##  · 모형은 Table 4 의 M1 과 동일 (연령·성별·기저 동반질환 + 구간시작 프레일티)
##    단, 층화변수 자체는 층 내 상수가 되므로 자동으로 빠집니다(성별 층화 시 sex_f).
##  · 층: 성별(남/여), 연령(구간 시작 시점 중앙값 기준 2분)
##  · 상호작용 P: 층화하지 않은 전체 자료에 gLIC_z x 층 항을 넣어 검정
##    (층별 CI 가 겹치는지 눈으로 보는 것보다 정확합니다)
##  · 사건이 적은 층은 계수가 발산하므로 MIN_CELL_EVENTS 미만이면 그리지 않습니다
##
## ── 색 ───────────────────────────────────────────────────────────────────
##  ★ 하위군은 '좋음/나쁨' 축이 아니므로 중증도 팔레트(파랑-회색-주황-주홍)를
##    쓰지 않습니다. 남성=주황 이면 "남성이 나쁨"으로 읽힙니다.
##    KF_theme 의 pal_sub() (검정/보라)를 씁니다.
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
if (!exists("THEME_VERSION") || THEME_VERSION < "2026-08-01")
  stop("02_theme_and_output.R (2026-08-01) 이 필요합니다 — pal_sub() 가 없습니다.")
if (!exists("ADJ_VERSION")) stop("03_functions_incremental.R 을 R/ 폴더에 두십시오 (kf_init 이 적재합니다).")
###############################################################################
## ★ 층별 추정치의 신뢰 하한 (2026-07-31, 검증 실행에서 발견)
##
##  fit_transition_pois() 의 기본 하한은 5건인데, 이는 삼분위 같은 범주형
##  노출 기준입니다. 층으로 쪼개면 연속형 노출에서도 계수가 발산합니다.
##  실제로 Severe -> Mild (Age < 79) 가 8건에서 IRR 36.09 (3.21-406) 로
##  튀었습니다. 이런 값이 forest 에 실리면 x축이 망가지고 오독을 만듭니다.
##
##  두 가지 기준으로 걸러 냅니다.
##    (1) 층 내 사건 < SUB_MIN_EV
##    (2) CI 상/하한 비가 SUB_MAX_CIR 초과 (사건수와 무관하게 발산 감지)
###############################################################################
##  ★ 기준 재조정 (2026-07-31, 실행 결과 확인 후)
##    20건 기준은 너무 엄격했습니다. 실행값에서 아래처럼 CI 가 충분히 좁은데도
##    사건수만으로 빠진 추정치가 9개 있었습니다.
##      Severe -> Normal (Age >= 79)  n=19  IRR 2.17 (1.48-3.16)  CI비 2.1
##      Mild -> Severe   (Male)       n=19  IRR 0.63 (0.44-0.90)  CI비 2.1
##    사건수는 정밀도의 대리지표일 뿐이고, 실제로 걸러야 하는 것은 '발산'입니다.
##    CI 비 기준을 조이고(20->10) 사건수 기준을 완화합니다(20->15).
##    -> 표시 58/70 에서 64/70 로 늘고, 발산 사례(n=5, CI비 12.1)는 그대로 제외됩니다.
SUB_MIN_EV  <- 15     # 층 내 최소 사건수
SUB_MAX_CIR <- 10     # ucl/lcl 허용 상한 (발산 감지의 주 기준)

load_imputed("MAIN")
d <- prep_long(load_kfacs())

## 이 그림에 넣을 체계 — 이분형 Any-ADL 은 3범주와 중복이라 제외(값표에는 남김)
## ★ 체계 라벨에 a/b 를 붙이지 않습니다. patchwork 가 패널마다 a/b/c/d 태그를
##   붙이므로, 여기에 또 a/b 가 있으면 "패널 a 안에 a 와 b" 가 되어 읽을 수 없습니다.
blocks <- list(c("state",        "ADL disability"),
               c("frailty_3cat", "Frailty phenotype"))

## ── 층 정의 ────────────────────────────────────────────────────────────
LV_ORDER <- NULL      # add_groups() 가 층의 표시 순서를 여기에 기록합니다
###############################################################################
## ★ 연령 절단 (2026-07-31)
##   중앙값(79세)은 자료에서 나온 값이라 원고에서 임의로 보입니다.
##   80세는 관례적 절단(oldest-old)이면서, 이 자료에서 균형도 더 낫습니다:
##     79세 -> 48.5% / 51.5%   80세 -> 56.1% / 43.9%
##   희소 전이의 사건 균형은 오히려 80세가 유리합니다
##     Mild->Severe   79: 19/103   80: 25/97
##     Normal->Death  79: 71/259   80: 97/233
##
## ★ 주의: 이 층은 **시간가변**입니다. 구간 시작 시점 연령으로 나누므로
##   한 사람이 초반에는 <80, 후반에는 >=80 층에 들어갑니다. 성별(시간고정)과
##   성격이 다르므로 각주에 "person-time before vs from age 80" 으로 씁니다.
##   ("80세 미만인 사람 vs 이상인 사람" 으로 읽히면 안 됩니다.)
###############################################################################
AGE_CUT    <- 80     # 연령 층 절단(세). NA 로 두면 중앙값을 씁니다.
COMORB_CUT <- 2      # 동반질환 층: < CUT vs >= CUT

add_groups <- function(iv) {
  ## 성별 — 시간고정
  iv$sex_g <- factor(as.character(iv$sex_f), levels = c("Male", "Female"))

  ## 연령 — 시간가변 (구간 시작 시점)
  av <- as_num(iv$bage) + as_num(iv$t0)
  cut_age <- if (is.na(AGE_CUT)) stats::median(av, na.rm = TRUE) else AGE_CUT
  bal <- mean(av < cut_age, na.rm = TRUE)
  if (is.finite(bal) && (bal < 0.15 || bal > 0.85))
    warning(sprintf("연령 절단 %.0f세가 %.0f%%/%.0f%% 로 치우쳤습니다 -> AGE_CUT 재검토",
                    cut_age, 100 * bal, 100 * (1 - bal)))
  iv$age_g <- factor(ifelse(av < cut_age, "younger", "older"),
                     levels = c("younger", "older"),
                     labels = c(sprintf("Age < %.0f y", cut_age),
                                sprintf("Age >= %.0f y", cut_age)))

  ## 동반질환 — 기저 고정. "IC 는 결국 질병부담 아니냐"에 대한 답.
  cb <- as_num(iv$comorbid_bl)
  iv$comorb_g <- factor(ifelse(cb < COMORB_CUT, "lo", "hi"), levels = c("lo", "hi"),
                        labels = c(sprintf("%d-%d", 0, COMORB_CUT - 1),
                                   sprintf(">= %d", COMORB_CUT)))

  ## 사회경제 — 기저(wave 1) 고정. "IC-결과 연관이 사회계층에 따라 다른가"에
  ##   대한 답. Unknown 수준은 층화에서 제외합니다(다른 패널에는 그대로 포함).
  ##   (v260813: 쇠약 층화는 Table 3 과 중복되어 제거하고 SES 로 교체)
  .g2 <- function(v, lv, lab) { x <- as.character(v); x[x == "Unknown"] <- NA
    factor(x, levels = lv, labels = lab) }
  iv$edu_g    <- .g2(iv$edu_bl_f, c("PrimaryOrLess", "MidSchoolPlus"),
                     c("Primary or less", "Middle school+"))
  iv$income_g <- .g2(iv$income_bl_f, c("LowerTwo", "TopTertile"),
                     c("Lower two", "Top tertile"))
  .ar <- as.character(iv$area_bl_f); .ar[.ar == "Unknown"] <- NA
  iv$area_g <- factor(ifelse(is.na(.ar), NA,
                             ifelse(.ar == "Rural", "Rural", "Metro or urban")),
                      levels = c("Metro or urban", "Rural"))

  ## ★ 층의 표시 순서를 전역에 기록합니다.
  ##   by_stratum() 이 stratum 을 문자로 저장해 factor 수준이 사라지므로,
  ##   그림에서 순서가 자료 등장순(사실상 무작위)이 되어 패널마다 색 배정이
  ##   달라집니다. 기준(건강한 쪽)이 항상 먼저 오도록 여기서 고정합니다.
  LV_ORDER <<- list(sex_g    = levels(iv$sex_g),
                    age_g    = levels(iv$age_g),
                    comorb_g = levels(iv$comorb_g),
                    edu_g    = levels(iv$edu_g),
                    income_g = levels(iv$income_g),
                    area_g   = levels(iv$area_g))
  attr(iv, "cut_age") <- cut_age
  iv
}

## ── 층별 추정 ──────────────────────────────────────────────────────────
by_stratum <- function(iv, gvar, tl) {
  lv <- levels(iv[[gvar]])
  do.call(rbind, lapply(lv, function(l) {
    sub <- iv[!is.na(iv[[gvar]]) & iv[[gvar]] == l, , drop = FALSE]
    if (!nrow(sub)) return(NULL)
    r <- run_adj_set(sub, "M1", c(CFG$ADJ, "frail_f"), "gLIC_z", tl)
    if (is.null(r)) return(NULL)
    r$gvar <- gvar; r$stratum <- l
    r
  }))
}

## ── 상호작용 검정 ──────────────────────────────────────────────────────
## 층별 CI 가 겹치는지 눈으로 보는 것은 검정이 아닙니다. 전체 자료에
## gLIC_z x 층 항을 넣어 직접 검정합니다.
inter_p <- function(iv, from, to, gvar) {
  x <- iv[iv$from_lab == from & is.finite(iv$dur) & iv$dur > 0 & !is.na(iv[[gvar]]), ]
  if (!nrow(x)) return(NA_real_)
  x$y <- as.integer(x$to_lab == to)
  ev <- tapply(x$y, droplevels(factor(x[[gvar]])), sum)
  if (length(ev) < 2 || any(ev < CFG$MIN_CELL_EVENTS, na.rm = TRUE)) return(NA_real_)
  rhs <- c(CFG$ADJ, "frail_f")
  rhs <- rhs[vapply(rhs, function(v) {
    vv <- all.vars(stats::as.formula(paste("~", v)))[1]
    vv %in% names(x) && length(unique(stats::na.omit(x[[vv]]))) > 1
  }, logical(1))]
  f <- stats::as.formula(paste0("y ~ gLIC_z * ", gvar,
        if (length(rhs)) paste0(" + ", paste(rhs, collapse = " + ")) else "",
        " + offset(log(dur))"))
  m <- try(stats::glm(f, family = stats::poisson(), data = x), silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  ct <- try(.robust_ct(m, x$id), silent = TRUE)
  if (inherits(ct, "try-error") || is.null(ct)) return(NA_real_)
  k <- grep(paste0("^gLIC_z:", gvar), rownames(ct))
  if (!length(k)) return(NA_real_)
  suppressWarnings(min(ct[k, ncol(ct)], na.rm = TRUE))
}

## ── 실행 ────────────────────────────────────────────────────────────────
GV   <- c("sex_g", "age_g", "comorb_g", "edu_g", "income_g", "area_g")
GLAB <- c(sex_g    = "By sex",
          age_g    = "By age (person-time)",
          comorb_g = "By baseline comorbidity",
          edu_g    = "By education",
          income_g = "By household income",
          area_g   = "By area of residence")
GV_FOR <- function(sv) GV

RES <- do.call(rbind, lapply(blocks, function(bk) {
  iv <- prep_iv_adj(d, bk[1], "rolling")
  if (is.null(iv)) return(NULL)
  iv <- add_groups(iv)
  tl <- default_transitions(bk[1], iv)
  gv <- GV_FOR(bk[1])
  out <- do.call(rbind, lapply(gv, function(g) by_stratum(iv, g, tl)))
  if (is.null(out)) return(NULL)
  ip <- do.call(rbind, lapply(gv, function(g) do.call(rbind, lapply(tl, function(tr)
    data.frame(gvar = g, transition = paste(tr[1], "->", tr[2]),
               p_int = inter_p(iv, tr[1], tr[2], g), stringsAsFactors = FALSE)))))
  out <- merge(out, ip, by = c("gvar", "transition"), all.x = TRUE)
  out$sysname <- bk[2]
  out$cut_age <- attr(iv, "cut_age")
  out
}))
if (is.null(RES) || !nrow(RES)) stop("하위군 추정 결과가 비어 있습니다.")

## 신뢰할 수 없는 층별 추정치 제거 (사유는 각주에 기록)
RES$cir  <- RES$ucl / RES$lcl
RES$drop_ev  <- !is.na(RES$n_events) & RES$n_events < SUB_MIN_EV
RES$drop_cir <- is.finite(RES$cir) & RES$cir > SUB_MAX_CIR
RES$shown <- is.finite(RES$IRR) & !RES$drop_ev & !RES$drop_cir
n_drop_ev  <- sum(RES$drop_ev  & is.finite(RES$IRR))
n_drop_cir <- sum(RES$drop_cir & !RES$drop_ev)
msg(sprintf("층별 추정치 %d개 중 표시 %d개 (사건<%d: %d개, CI 발산: %d개)",
            sum(is.finite(RES$IRR)), sum(RES$shown), SUB_MIN_EV, n_drop_ev, n_drop_cir))

RES$kind <- kind_of(sub(" ->.*$", "", RES$transition), sub("^.*-> ", "", RES$transition))
RES$kind <- factor(RES$kind, levels = c("Worsening", "Recovery", "Death"))
save_vals(RES[, c("sysname","gvar","stratum","transition","kind","n_events","pyears",
                  "IRR","lcl","ucl","p","p_int","shown","note")],
          "SF1_subgroup_values.csv", FIG_DIR)

## 상호작용이 유의한 전이 (본문/각주용)
ints <- unique(RES[is.finite(RES$p_int) & RES$p_int < 0.05,
                   c("sysname","gvar","transition","p_int")])
## ★ 2026-07-31 수정: RES$row 열을 없앤 뒤에도 여기서 참조하고 있었습니다
##   ("정의하지 않은 열들이 선택되었습니다"). 전이 식별은 체계+전이로 합니다.
n_int_test <- nrow(unique(RES[is.finite(RES$p_int), c("gvar", "sysname", "transition")]))

int_txt <- if (nrow(ints)) {
  paste0("Interaction with the subgroup variable reached P<0.05 for ",
         paste(sprintf("%s (%s, P=%.3f)", ints$transition,
                       c(sex_g = "sex", age_g = "age", comorb_g = "comorbidity",
                         edu_g = "education", income_g = "income",
                         area_g = "area")[ints$gvar], ints$p_int),
               collapse = "; "),
         sprintf("; %d interaction tests were performed in total, so isolated results should be read as hypothesis-generating.", n_int_test))
} else {
  sprintf("No interaction reached P<0.05 across %d tests.", n_int_test)
}

## ── 그림 ────────────────────────────────────────────────────────────────
## 행 = 전이 (체계별 facet), 색 = 층. 라벨은 전이만 쓰고 체계는 facet 제목으로.
## 다중검정: 검정 건수로 Bonferroni 임계값을 계산해 넘는 것만 dagger 표시
## ★ n_int_test 는 위 블록에서 먼저 계산되어야 합니다(순서 주의).
BONF  <- 0.05 / max(1, n_int_test)
.pfmt <- function(p) ifelse(p < 0.001, "P<0.001", sprintf("P=%.3f", p))
PLOT <- RES[RES$shown, ]
if (!nrow(PLOT)) stop("표시할 층별 추정치가 없습니다. SUB_MIN_EV 를 낮춰 보십시오.")
tr_ord <- unique(PLOT[order(match(PLOT$sysname, vapply(blocks, `[`, character(1), 2)),
                           PLOT$kind, PLOT$ord), "transition"])
PLOT$sysf <- factor(PLOT$sysname, levels = vapply(blocks, `[`, character(1), 2))
## ── 체계 라벨: 우측 세로 스트립 대신 세로축 위의 구획 헤더 (v260813) ──────
##  "ADL disability, 3 states" -> "Disability", "Frailty phenotype" -> "Frailty"
##  각 패널 제목(By sex 등) 아래, 해당 전이 블록 바로 위 축 라벨 자리에 굵게.
SYS_SHORT <- function(s) {
  s <- as.character(s)
  s[grepl("^ADL", s)]     <- "Disability"
  s[grepl("^Frailty", s)] <- "Frailty"
  s[grepl("^Any ADL", s)] <- "Any disability"
  s
}
.disp <- character(0)                       # 표시 순서 (위 -> 아래)
for (bs in vapply(blocks, `[`, character(1), 2)) {
  trs <- tr_ord[PLOT$sysname[match(tr_ord, PLOT$transition)] == bs]
  trs <- unique(PLOT$transition[PLOT$sysname == bs])
  trs <- tr_ord[tr_ord %in% trs]
  .disp <- c(.disp, paste0("HDR_", SYS_SHORT(bs)), arrow_lab(trs))
}
Y_LEVELS <- rev(.disp)                      # ggplot y 는 아래에서 위로
## 라벨 함수: 헤더는 굵게(expression), 나머지는 일반 문자
##  (리스트 대신 expression 벡터를 반환 — ggplot2 버전에 무관하게 동작)
Y_LABFUN <- function(x) {
  out <- vector("expression", length(x))
  for (i in seq_along(x)) {
    l <- as.character(x[i])
    out[[i]] <- if (startsWith(l, "HDR_"))
      as.expression(bquote(bold(.(sub("^HDR_", "", l)))))[[1]]
    else as.expression(l)[[1]]
  }
  out
}
PLOT$trf <- factor(arrow_lab(PLOT$transition), levels = Y_LEVELS)

mk <- function(g) {
  x <- PLOT[PLOT$gvar == g, ]
  if (!nrow(x)) return(NULL)
  lv <- if (!is.null(LV_ORDER) && !is.null(LV_ORDER[[g]])) LV_ORDER[[g]] else
        sort(unique(as.character(x$stratum)))
  lv <- lv[lv %in% unique(as.character(x$stratum))]
  x$stratum <- factor(as.character(x$stratum), levels = lv)
  ## ★ droplevels 하지 않습니다. 두 패널(성별/연령)이 같은 행 순서를 가져야
  ##   나란히 놓았을 때 같은 전이가 같은 높이에 옵니다. 한쪽에서만 제거된
  ##   추정치는 그 자리를 빈 채로 둡니다(그것 자체가 정보입니다).
  ## ★ 상호작용 P 를 오른쪽에 함께 찍습니다.
  ##   "층별 CI 가 겹치니 일관되다"는 눈대중보다 P 가 옆에 있는 편이
  ##   리뷰어에게 훨씬 강합니다. 행마다 값이 하나이므로 중복 제거해서 붙입니다.
  pl <- unique(x[, c("trf", "p_int")])
  pl <- pl[!duplicated(pl$trf), ]
  pl$txt <- ifelse(!is.finite(pl$p_int), "",
             ifelse(pl$p_int < BONF, sprintf("%s\u2020", .pfmt(pl$p_int)), .pfmt(pl$p_int)))

  pd <- position_dodge(width = 0.55)
  ggplot(x, aes(IRR, trf, colour = stratum)) +
    geom_text(data = pl, aes(x = Inf, y = trf, label = txt), inherit.aes = FALSE,
              hjust = 1.03, size = BASE_PT / .pt * 0.85, colour = "grey25") +
    scale_y_discrete(limits = Y_LEVELS, labels = Y_LABFUN) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = LINE_PT,
               colour = SEV_COL[["aux"]]) +
    geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0, linewidth = 0.4,
                   position = pd, na.rm = TRUE) +
    geom_point(size = 1.4, position = pd, na.rm = TRUE) +
    scale_colour_manual(values = pal_sub(lv), name = NULL) +
    scale_x_continuous(trans = "log", breaks = c(0.25, 0.5, 1, 2, 4),
                       labels = c("0.25","0.50","1.00","2.00","4.00"),
                       expand = expansion(mult = c(0.05, 0.58))) +
    coord_cartesian(clip = "off") +
    labs(x = "IRR per +1 s.d. of IC (log scale)", y = NULL, title = GLAB[[g]]) +
    theme_na() +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top",
          axis.text.y = element_text(hjust = 1))
}
###############################################################################
## 그림 렌더 — 2단계 안전장치 (v260813b)
##  1차: 세로축 헤더를 굵은 글꼴(expression)로.  실패하면
##  2차: 헤더를 대문자 일반 문자로 (어떤 ggplot2 버전에서도 동작).
##  둘 다 실패하면 FIG_DIR/SF1_ERROR.txt 에 오류 전문을 기록하고 멈춥니다.
###############################################################################
.Y_LAB_EXPR <- function(x) {
  out <- vector("expression", length(x))
  for (i in seq_along(x)) {
    l <- as.character(x[i])
    out[[i]] <- if (startsWith(l, "HDR_"))
      as.expression(bquote(bold(.(sub("^HDR_", "", l)))))[[1]]
    else as.expression(l)[[1]]
  }
  out
}
.Y_LAB_PLAIN <- function(x) {
  l <- as.character(x)
  ifelse(startsWith(l, "HDR_"), toupper(sub("^HDR_", "", l)), l)
}

.render_sf1 <- function(use_expr) {
  Y_LABFUN <<- if (use_expr) .Y_LAB_EXPR else .Y_LAB_PLAIN
  message("[SF1] 단계 1/4: 패널 생성")
  PN <- lapply(GV, mk); names(PN) <- GV
  PN <- PN[!vapply(PN, is.null, logical(1))]
  if (length(PN) < 2) stop("그림을 만들 층별 추정치가 부족합니다.")
  message("[SF1] 단계 2/4: 패널 ", length(PN), "개 조립")
  fig <- if (length(PN) >= 6) (PN[[1]] | PN[[2]] | PN[[3]]) / (PN[[4]] | PN[[5]] | PN[[6]]) else
         if (length(PN) >= 4) (PN[[1]] | PN[[2]]) / (PN[[3]] | PN[[4]]) else
         if (length(PN) == 3) (PN[[1]] | PN[[2]] | PN[[3]]) else (PN[[1]] | PN[[2]])
  fig <- fig + plot_layout(guides = "keep")
  message("[SF1] 단계 3/4: 제목·각주 부착")
  fig <- fig + annot_na(
    width_mm = NA_WIDTH$double, tag = TRUE,
    title = "Consistency of the effect of intrinsic capacity across subgroups",
    caption = paste0(
      "Incidence rate ratios per +1 s.d. of intrinsic capacity, from the same models as the concurrent-frailty supplementary table ",
      "(M1: adjusted for age, sex, baseline comorbidity count, education, household income and area of residence, plus the frailty phenotype at ",
      "the start of the interval); the stratifying variable is constant within stratum and is ",
      "therefore dropped from the covariate set. Sex, baseline comorbidity, education, household income and area of ",
      "residence are fixed characteristics, whereas the age strata are time-varying: they divide ",
      "person-time before and from age ", sprintf("%.0f", RES$cut_age[1]), " years, so a ",
      "participant contributes to both strata as they age. The socioeconomic panels address ",
      "whether the association between capacity and the transitions differs across social strata; ",
      "participants with an unknown value are omitted from that panel only. ",
      "Estimates are omitted where a stratum contributed fewer than ", SUB_MIN_EV,
      " events or where the confidence interval spanned more than a ", SUB_MAX_CIR,
      "-fold range, because the coefficient then diverges rather than converging (",
      n_drop_ev + n_drop_cir, " of ", sum(is.finite(RES$IRR)), " stratum-specific estimates). ",
      "The number to the right of each row is the P value for the interaction between intrinsic ",
      "capacity and the subgroup variable; a dagger marks values below the Bonferroni threshold ",
      "for the ", n_int_test, " interaction tests performed (P<", sprintf("%.4f", BONF), "). ",
      "The direction of the association was identical in both strata for every pair of ",
      "displayed estimates. ", int_txt,
      " Subgroup colours are deliberately outside the severity palette used ",
      "elsewhere, because the subgroup variables are not ordered from better to worse."))
  message("[SF1] 단계 4/4: 저장 (여기서 실제 렌더가 일어납니다)")
  save_na(fig, "SupplFig1_Subgroups", width_mm = NA_WIDTH$double,
          height_mm = if (length(PN) >= 4) 235 else 170, dir = FIG_DIR)
  invisible(TRUE)
}

.err_log <- function(e, tag) {
  fp <- file.path(FIG_DIR, "SF1_ERROR.txt")
  cl <- tryCatch(paste(deparse(conditionCall(e)), collapse = " "),
                 error = function(x) "(호출 없음)")
  vs <- tryCatch(paste(sprintf("%s %s", c("ggplot2","patchwork","scales","R"),
                 c(as.character(utils::packageVersion("ggplot2")),
                   as.character(utils::packageVersion("patchwork")),
                   as.character(utils::packageVersion("scales")),
                   paste(R.version$major, R.version$minor, sep = "."))),
                 collapse = " | "), error = function(x) "")
  writeLines(c(paste0("[", tag, "] ", format(Sys.time())),
               paste("메시지 :", conditionMessage(e)),
               paste("실패 호출:", cl),
               paste("버전    :", vs)), fp)
  message("[SF1] ", tag, " 렌더 실패: ", conditionMessage(e),
          "\n  실패 호출: ", cl, "\n  오류 전문: ", fp)
}
ok <- tryCatch(.render_sf1(TRUE), error = function(e) { .err_log(e, "expression 라벨"); FALSE })
if (!isTRUE(ok)) {
  message("[SF1] 일반 문자 라벨로 재시도합니다 (헤더는 대문자로 표시).")
  ok <- tryCatch(.render_sf1(FALSE), error = function(e) { .err_log(e, "일반 라벨"); FALSE })
}
if (!isTRUE(ok))
  stop("SF1 그림 렌더에 실패했습니다. FIG_DIR 의 SF1_ERROR.txt 를 보내주십시오.")

cat("\n=== Supplementary Figure 4 완료 ===\n")
cat(sprintf("  층별 추정치 %d개 중 표시 %d개 | 상호작용 검정 %d건 | P<0.05 %d건\n",
            sum(is.finite(RES$IRR)), sum(RES$shown), n_int_test, nrow(ints)))

## ★ 가장 중요한 요약: 방향이 뒤집힌 전이가 있는가
##
##   세 기준으로 나누어 봅니다 (2026-07-31 실행 결과 반영).
##   사건 5~7건짜리 불안정 추정치(CI 가 1 을 크게 포함)가 부호만 뒤집히는 것을
##   "반전"으로 세면 오해를 부릅니다. 실행값에서 전체 61쌍 중 3쌍이 반전으로
##   잡혔는데, 3쌍 모두 한쪽이 표시 기준 미달(n=5-7)이고 전부 비유의였습니다.
##   표시된 쌍만 보면 반전 0, 양쪽 다 유의한 쌍만 보면 역시 반전 0 입니다.
##   원고에는 '표시된 쌍' 기준을 쓰십시오.
.ok  <- is.finite(RES$IRR)
.key <- paste(RES$gvar, RES$sysname, RES$transition)
dir_chk <- do.call(rbind, lapply(split(RES[.ok, ], .key[.ok]), function(z) {
  if (nrow(z) != 2) return(NULL)
  data.frame(gvar = z$gvar[1], sysname = z$sysname[1], transition = z$transition[1],
             same_dir   = all(z$IRR > 1) || all(z$IRR < 1),
             both_shown = all(z$shown),
             both_sig   = all(is.finite(z$lcl) & is.finite(z$ucl) & (z$lcl > 1 | z$ucl < 1)),
             stringsAsFactors = FALSE)
}))
n_rev_shown <- NA_integer_
if (!is.null(dir_chk)) {
  .rep <- function(lab, sel) {
    n <- sum(sel); rv <- sum(sel & !dir_chk$same_dir)
    cat(sprintf("    %-20s %2d쌍 | 방향 동일 %2d | 반전 %d\n", lab, n, n - rv, rv))
    rv
  }
  cat("  방향 일관성\n")
  r_all   <- .rep("전체",           rep(TRUE, nrow(dir_chk)))
  r_shown <- .rep("양쪽 다 표시된 쌍", dir_chk$both_shown)
  r_sig   <- .rep("양쪽 다 유의한 쌍", dir_chk$both_sig)
  n_rev_shown <- r_shown
  if (r_shown == 0)
    cat("  ★ 그림에 표시된 모든 쌍에서 방향이 동일합니다.\n",
        "    -> 상호작용이 유의해도 '크기 차이'일 뿐 '불일치'가 아닙니다.\n", sep = "")
  if (r_all > r_shown)
    cat(sprintf("    (전체 기준 반전 %d쌍은 모두 표시 기준 미달인 불안정 추정치입니다)\n",
                r_all - r_shown))
}

if (nrow(ints)) {
  cat(sprintf("\n  Bonferroni 임계값 P<%.5f (검정 %d건) 을 넘는 상호작용:\n", BONF, n_int_test))
  bi <- ints[ints$p_int < BONF, , drop = FALSE]
  if (nrow(bi)) print(bi, row.names = FALSE) else
    cat("    없음 -- P<0.05 인 ", nrow(ints), "건은 모두 다중검정으로 설명됩니다.\n", sep = "")
  cat("\n  (참고) P<0.05 전체:\n"); print(ints, row.names = FALSE)
} else {
  cat("  ★ 성별·연령 상호작용이 유의한 전이가 없습니다.\n")
}
