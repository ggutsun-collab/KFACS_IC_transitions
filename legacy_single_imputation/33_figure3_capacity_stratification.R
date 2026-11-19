# =============================================================================
# Figure 3 — Capacity scaled within sex stratifies five-year risk within every
# frailty stratum.  KFACS.
#
# 색·서체는 02_theme_and_output.R 의 규칙을 그대로 따릅니다(Fig 1·2와 동일):
#   T1(낮음) #E69F00 · T2(중간) #595959 · T3(높음) #0072B2
#   사망 #D55E00 · 악화 #E69F00 · 보조선 #A6A6A6 · 본문 7pt · 패널라벨 8pt
# 값은 03_analysis/02_tables_1_2_3.R 의 table3_ev.csv / table3_e_comp.csv 에서 옵니다.
# 그림 안에는 제목도 해설문도 두지 않습니다. 전부 원고 legend 가 담당합니다.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
})

OUT3 <- if (exists("FIG_DIR")) FIG_DIR else "."

BASE_PT <- 7; LAB_PT <- 8; LINE_PT <- 0.35
COL_T   <- c("T1" = "#E69F00", "T2" = "#595959", "T3" = "#0072B2")   # 능력 삼분위
COL_DEA <- "#D55E00"                                                  # 사망
COL_WOR <- "#E69F00"                                                  # 악화
COL_AUX <- "#A6A6A6"
FRL <- c("Robust", "Pre-frail", "Frail")
TRL <- c("T1\nlowest", "T2", "T3\nhighest")

theme_na <- function(base = BASE_PT) {
  theme_classic(base_size = base) +
    theme(text = element_text(colour = "black"),
          axis.text = element_text(size = base, colour = "black"),
          axis.title = element_text(size = base),
          axis.line = element_line(linewidth = LINE_PT, colour = "black"),
          axis.ticks = element_line(linewidth = LINE_PT, colour = "black"),
          axis.ticks.length = unit(1.2, "mm"),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = base, face = "bold", hjust = 0),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = base),
          legend.title = element_text(size = base),
          legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(1, "mm"),
          plot.title = element_text(size = base, face = "plain", hjust = 0),
          plot.tag = element_text(size = LAB_PT, face = "bold"),
          plot.tag.position = c(0, 1),
          plot.margin = margin(2, 2, 2, 2, "mm"))
}

dat <- tribble(
  ~fr,        ~tert, ~n,  ~d_ev, ~d_risk, ~d_lo, ~d_hi, ~w_ev, ~w_risk,
  "Robust",    "T1", 104,  12,    11.5,    5.2,  17.5,   25,    24.0,
  "Robust",    "T2", 429,  18,     4.2,    2.3,   6.1,   72,    16.8,
  "Robust",    "T3", 816,  21,     2.6,    1.5,   3.7,  117,    14.3,
  "Pre-frail", "T1", 670,  66,     9.9,    7.6,  12.1,  169,    25.2,
  "Pre-frail", "T2", 560,  24,     4.3,    2.6,   5.9,  108,    19.3,
  "Pre-frail", "T3", 186,   6,     3.2,    0.7,   5.7,   29,    15.6,
  "Frail",     "T1", 231,  36,    15.6,   10.8,  20.1,   83,    35.9,
  "Frail",     "T2",  14,   2,    14.3,    0.0,  30.8,    4,    28.6,
  "Frail",     "T3",   1,   0,      NA,     NA,    NA,    0,      NA
) %>% mutate(fr = factor(fr, FRL),
             tert = factor(tert, c("T1", "T2", "T3")),
             small = n < 20)

## ── a, b : 위험 격자 ────────────────────────────────────────────────────
heat <- function(v, ev, hi_col, lim, ttl) {
  d <- dat %>% mutate(val = .data[[v]], ev = .data[[ev]])
  ggplot(d, aes(tert, fr, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(val), "NE", sprintf("%.1f", val)),
                  colour = !is.na(val) & val > lim[2] * 0.62),
              size = 2.5, fontface = "plain", vjust = -0.15, show.legend = FALSE) +
    geom_text(aes(label = ifelse(is.na(val), "n = 1", sprintf("%d/%d", ev, n)),
                  colour = !is.na(val) & val > lim[2] * 0.62),
              size = 1.8, vjust = 1.6, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradient(low = "#f7f7f7", high = hi_col, limits = lim,
                        na.value = "#eeeeee", name = "%",
                        guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                barheight = unit(22, "mm"),
                                                ticks.colour = "white")) +
    scale_x_discrete(labels = TRL) +
    scale_y_discrete(limits = rev(FRL)) +
    labs(x = "Intrinsic capacity tertile (within sex)", y = NULL, title = ttl) +
    theme_na() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = "right", legend.title = element_text(size = BASE_PT - 1))
}
pa <- heat("d_risk", "d_ev", COL_DEA, c(0, 20), "Five-year mortality")
pb <- heat("w_risk", "w_ev", COL_WOR, c(0, 40), "Five-year composite worsening")

## ── c : 9개 칸을 한 축에. 쇠약 층은 facet 으로 묶고(Fig 1 과 같은 방식),
##        색은 Fig 1·2 와 동일하게 능력 삼분위를 뜻합니다.
dc <- dat %>%
  mutate(lab = ifelse(is.na(d_risk), "", sprintf("%.1f%%", d_risk)),
         lab_x = pmin(d_hi + 1.0, 40),
         ne = is.na(d_risk),
         tert_lab = factor(tert, levels = c("T3", "T2", "T1"),
                           labels = c("T3", "T2", "T1")))

pc <- ggplot(dc, aes(d_risk, tert_lab, colour = tert)) +
  geom_vline(xintercept = 11.5, linetype = "22", colour = COL_AUX, linewidth = LINE_PT) +
  geom_errorbarh(aes(xmin = d_lo, xmax = d_hi), height = 0, linewidth = 0.45) +
  geom_point(aes(shape = small), size = 1.9, fill = "white", stroke = 0.5) +
  geom_text(aes(x = lab_x, label = lab), hjust = 0, size = 2.1, colour = "black",
            na.rm = TRUE) +
  geom_text(data = dplyr::filter(dc, ne),
            aes(x = 0.6, y = tert_lab, label = "not estimable (n = 1)"),
            hjust = 0, size = 2.0, colour = COL_AUX, fontface = "italic",
            inherit.aes = FALSE, show.legend = FALSE) +
  facet_grid(fr ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  scale_colour_manual(values = COL_T, name = "Intrinsic capacity tertile (within sex)",
                      labels = c("T1 (lowest)", "T2", "T3 (highest)"),
                      breaks = c("T1", "T2", "T3")) +
  scale_x_continuous(limits = c(0, 46), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Five-year mortality (%)", y = NULL,
       title = "Five-year mortality in all nine cells, grouped by frailty stratum") +
  theme_na() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1,
                                         margin = margin(r = 2, unit = "mm")),
        panel.spacing.y = unit(1.5, "mm"),
        legend.position = "top", legend.justification = "left")

fig <- (pa | pb) / pc + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = LAB_PT, face = "bold"))

ggsave(file.path(OUT3, "Fig3_capacity_stratification.pdf"), fig,
       width = 180, height = 140, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT3, "Fig3_capacity_stratification.png"), fig,
       width = 180, height = 140, units = "mm", dpi = 450)
message("saved: ", file.path(OUT3, "Fig3_capacity_stratification.pdf/.png"))
