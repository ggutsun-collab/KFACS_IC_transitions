###############################################################################
## 35_edfig1_hospitalization_fall_MI.R                                  v261117
##
## Extended Data Fig. 1 and Supplementary Table 7. Resilience to a defined
## stressor: functional decline over a between-visit interval in which a
## hospitalization or a fall was reported at the closing visit.
##   - events: self-reported for the 12 months before each visit (hospi == 1;
##     fall_ex == 1; fall with fracture = fall plus any fall_frac1-7 == 1),
##     read from RAW_FILE (R/00_setup.R); not imputed (missing or code 9 ->
##     interval excluded from that event's analysis)
##   - an event reported at visit k+1 is assigned to the interval k -> k+1;
##     only participants alive at visit k+1 contribute
##   - any decline: worse frailty or worse ADL state at k+1;
##     major decline: transition to frailty or to a worse ADL state
##   - exposure: within-sex IC at visit k, per +1 s.d.
##   - risk ratios: Poisson regression (log link) with participant-clustered
##     robust variance, adjusted for age at visit k, baseline comorbidity,
##     education, household income, area of residence and wave
##   - groups: robust at k with the event; robust at k without the event;
##     all states at k with the event. Interaction: ratio of the IC risk
##     ratios with versus without the event among robust participants.
## Reproduction check (manuscript): robust, any decline after hospitalization
##   0.60 (0.51-0.72); after a fall 0.61 (0.53-0.70); P interaction >= 0.09.
###############################################################################

message("\n=== 35_edfig1_hospitalization_fall_MI v261117b ===")
## -- helpers (R/_bootstrap.R -> kf_init(); all paths come from R/00_setup.R) --
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
suppressPackageStartupMessages(library(survival))
.need(c("ggplot2"))
suppressPackageStartupMessages(library(ggplot2))

if (!exists("RAW_FILE")) RAW_FILE <- Sys.getenv("KFACS_RAW_FILE",
  unset = file.path(PROJECT_ROOT, "data", "KFACS_master_FINAL_preimput_260616.rds"))
if (!file.exists(RAW_FILE)) stop("RAW_FILE not found: ", RAW_FILE,
  "\n  -> set Sys.setenv(KFACS_RAW_FILE = '<path to KFACS_master_FINAL_preimput_260616.rds>')")
raw <- readRDS(RAW_FILE)
need_cols <- c("id", "wave", "hospi", "fall_ex", paste0("fall_frac", 1:7))
miss <- setdiff(need_cols, names(raw)); if (length(miss)) stop("RAW_FILE lacks: ", paste(miss, collapse = ", "))
EVT <- raw[, need_cols]; rm(raw); invisible(gc())
EVT$id <- as.character(EVT$id); EVT$wave <- as.integer(as_num(EVT$wave))
yn <- function(x) { x <- as_num(x); ifelse(x %in% 1, 1L, ifelse(x %in% 2, 0L, NA_integer_)) }
EVT$hosp <- yn(EVT$hospi); EVT$fall <- yn(EVT$fall_ex)
fr <- rowSums(sapply(paste0("fall_frac", 1:7), function(v) as_num(EVT[[v]]) %in% 1)) > 0
EVT$fracture <- ifelse(is.na(EVT$fall), NA_integer_, as.integer(EVT$fall == 1 & fr))
EVENTS <- c(hosp = "Hospitalization", fall = "Fall", fracture = "Fall with fracture")

SETS <- names(readRDS(file.path(IMP_DIR, paste0(MI_STEM, ".rds"))))
COV  <- c(SM_COV_PRIMARY, "wave_f")
fit_rr <- function(dd, term = "z", form_extra = NULL) {
  dd <- droplevels(dd); if (sum(dd$y) < 5) return(c(b = NA, v = NA, ev = sum(dd$y), n = nrow(dd)))
  m <- try(glm(reformulate(c(form_extra %||% "z", .sm_keep(dd, COV)), response = "y"),
               family = poisson(link = "log"), data = dd, control = list(maxit = 100)), silent = TRUE)
  if (inherits(m, "try-error") || !term %in% names(coef(m)) || is.na(coef(m)[term]))
    return(c(b = NA, v = NA, ev = sum(dd$y), n = nrow(dd)))
  c(b = unname(coef(m)[term]), v = .sm_vcov_cl(m, dd$id)[term, term], ev = sum(dd$y), n = nrow(dd))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## deaths within this many years after the opening visit close the interval
DEATH_WINDOW <- 2.5
SENS_EVENTS  <- c("hosp", "fall")
R <- list(); NDEC <- list()
for (set in SETS) {
  load_imputed(set, stem = MI_STEM)
  d <- ws_scale(prep_long(load_kfacs()))
  K <- sm_person_frame(d)
  d$st <- as.integer(as_num(d$state)); d$fst <- as.integer(as_num(d$frail_state))
  v <- d[!is.na(d$wave) & d$st %in% 1:3 & d$fst %in% 1:3, c("id", "wave", "st", "fst", "gLIC_z", "age_c")]
  v$wave <- as.integer(v$wave)
  nx <- v[, c("id", "wave", "st", "fst")]; names(nx)[3:4] <- c("st1", "fst1"); nx$wave <- nx$wave - 1L
  iv <- merge(v[is.finite(v$gLIC_z), ], nx, by = c("id", "wave"))
  ## events are reported at the closing visit (wave + 1)
  e1 <- EVT[, c("id", "wave", names(EVENTS))]; e1$wave <- e1$wave - 1L
  iv <- merge(iv, e1, by = c("id", "wave"), all.x = TRUE)
  iv <- merge(iv, K[, c("id", "comorbid_bl", "edu_bl_f", "income_bl_f", "area_bl_f")], by = "id")
  iv$z <- iv$gLIC_z; iv$wave_f <- factor(iv$wave)
  iv$any   <- as.integer(iv$fst1 > iv$fst | iv$st1 > iv$st)
  iv$major <- as.integer((iv$fst1 == 3 & iv$fst < 3) | iv$st1 > iv$st)
  ## -- sensitivity: intervals ending in death (no closing interview, so the
  ##    event is unknown). Death counts as any and major decline; the event is
  ##    set to 1 (scenario A) or 0 (scenario B) for these intervals.
  al <- d[d$st %in% 1:3 & d$fst %in% 1:3 & !is.na(d$wave) & is.finite(d$time), c("id", "wave", "time")]
  al <- al[order(al$id, -al$time), ]; last <- al[!duplicated(al$id), ]
  isd <- (d$st %in% 4 | d$fst %in% 4) & is.finite(d$time)
  td <- tapply(d$time[isd], d$id[isd], min)
  last$tdeath <- as.numeric(td[as.character(last$id)])
  dk <- last[is.finite(last$tdeath) & last$tdeath > last$time & last$tdeath - last$time <= DEATH_WINDOW, c("id", "wave")]
  dk$wave <- as.integer(dk$wave)
  di <- merge(v[is.finite(v$gLIC_z), ], dk, by = c("id", "wave"))
  di <- merge(di, K[, c("id", "comorbid_bl", "edu_bl_f", "income_bl_f", "area_bl_f")], by = "id")
  di$z <- di$gLIC_z; di$wave_f <- factor(di$wave, levels = levels(iv$wave_f)); di$any <- 1L; di$major <- 1L
  for (e in SENS_EVENTS) for (o in c("any", "major")) {
    x <- iv[!is.na(iv[[e]]), ]; x$y <- x[[o]]; x$ev <- x[[e]]
    keep <- c("id", "wave_f", "z", "y", "ev", "fst", COV[COV != "wave_f"])
    for (sc in c("A", "B")) {
      dd <- di; dd$y <- 1L; dd$ev <- if (sc == "A") 1L else 0L
      xs <- rbind(x[, keep], dd[, keep]); xs <- xs[!is.na(xs$wave_f), ]
      gname <- if (sc == "A") "sensA_robust_event" else "sensB_robust_no_event"
      r <- fit_rr(xs[xs$fst == 1 & xs$ev == (sc == "A"), ])
      R[[length(R) + 1]] <- data.frame(set = set, event = e, outcome = o, group = gname, b = r["b"], v = r["v"], ev = r["ev"], n = r["n"])
      r <- fit_rr(xs[xs$fst == 1, ], term = "z:ev", form_extra = c("z", "ev", "z:ev"))
      R[[length(R) + 1]] <- data.frame(set = set, event = e, outcome = o, group = paste0("sens", sc, "_interaction"),
                                       b = r["b"], v = r["v"], ev = NA, n = r["n"])
    }
  }
  NDEC[[set]] <- c(all = nrow(di), robust = sum(di$fst == 1))
  for (e in names(EVENTS)) for (o in c("any", "major")) {
    x <- iv[!is.na(iv[[e]]), ]; x$y <- x[[o]]; x$ev <- x[[e]]
    grp <- list(robust_event = x[x$fst == 1 & x$ev == 1, ], robust_no_event = x[x$fst == 1 & x$ev == 0, ],
                all_event = x[x$ev == 1, ])
    for (g in names(grp)) { r <- fit_rr(grp[[g]])
      R[[length(R) + 1]] <- data.frame(set = set, event = e, outcome = o, group = g, b = r["b"], v = r["v"], ev = r["ev"], n = r["n"]) }
    xr <- x[x$fst == 1, ]; r <- fit_rr(xr, term = "z:ev", form_extra = c("z", "ev", "z:ev"))
    R[[length(R) + 1]] <- data.frame(set = set, event = e, outcome = o, group = "interaction", b = r["b"], v = r["v"], ev = NA, n = r["n"])
  }
  msg("[", set, "] ", nrow(iv), " intervals")
}
R <- do.call(rbind, R)
GL <- c(robust_event = "Robust before the event", robust_no_event = "Robust, no event", all_event = "All states, event",
        interaction = "Interaction (event vs no event, robust)",
        sensA_robust_event = "Deaths as decline, event assumed: robust with event",
        sensA_interaction  = "Deaths as decline, event assumed: interaction",
        sensB_robust_no_event = "Deaths as decline, no event assumed: robust without event",
        sensB_interaction  = "Deaths as decline, no event assumed: interaction")
OUT <- do.call(rbind, lapply(split(R, paste(R$event, R$outcome, R$group)), function(z) {
  p <- sm_rubin(z$b, z$v); m <- sum(is.finite(z$b)); bb <- stats::var(z$b, na.rm = TRUE); ub <- mean(z$v, na.rm = TRUE)
  data.frame(event = z$event[1], outcome = z$outcome[1], group = z$group[1],
             Event = EVENTS[z$event[1]], Outcome = ifelse(z$outcome[1] == "any", "Any decline", "Major decline"),
             Group = GL[z$group[1]], Intervals = median(z$n), Declines = median(z$ev),
             rr = exp(p["est"]), lo = exp(p["lo"]), hi = exp(p["hi"]), RR = sm_fmt(p), P = fmt_p(p["p"]),
             FMI = round((1 + 1/m) * bb / (ub + (1 + 1/m) * bb), 2), row.names = NULL, stringsAsFactors = FALSE)
}))
OUT <- OUT[order(match(OUT$event, names(EVENTS)), OUT$outcome, match(OUT$group, names(GL))), ]
TAB <- OUT[, c("Event", "Outcome", "Group", "Intervals", "Declines", "RR", "P", "FMI")]
names(TAB)[6] <- "RR per +1 s.d. (95% CI)"
cat("\n=== Supplementary Table 7 ===\n"); print(TAB, row.names = FALSE)
NDEC <- do.call(rbind, NDEC)
cat("\nIntervals ending in death within", DEATH_WINDOW, "years (median over sets): all",
    median(NDEC[, "all"]), "; robust at opening visit", median(NDEC[, "robust"]), "\n")
save_table(TAB, "ST7_Hospitalization_Fall_MI",
  title = "Supplementary Table 7 | Intrinsic capacity and functional decline after hospitalization or a fall",
  footnotes = c("Hospitalization and falls in the 12 months before each visit were self-reported and were not imputed; an event was assigned to the between-visit interval closed by the visit at which it was reported, and only participants alive at that visit contribute. Fall with fracture is a fall with any reported fracture.",
                "Any decline is a worse frailty state or a worse ADL disability state at the closing visit than at the opening visit; major decline is a transition to frailty or to a worse ADL disability state.",
                "Risk ratios per +1 s.d. of intrinsic capacity standardized within sex and measured at the opening visit, from Poisson models (log link) adjusted for age at the opening visit, baseline comorbidity count, education, household income, area of residence and wave, with standard errors clustered by participant, combined across 20 completed datasets by Rubin's rules. The interaction row gives the ratio of risk ratios with versus without the event among participants robust at the opening visit.",
                sprintf("Rows labelled 'Deaths as decline' are a sensitivity analysis that adds the intervals ending in death before the closing visit (median %s intervals, %s of them in participants robust at the opening visit). Death is counted as both any and major decline; because no closing interview took place, hospitalization or a fall in these intervals was unknown and was assumed present for all (event assumed) or absent for all (no event assumed).",
                        format(median(NDEC[, "all"])), format(median(NDEC[, "robust"]))),
                "ADL, activities of daily living; FMI, fraction of missing information; RR, risk ratio."))
save_vals(OUT, "EDFig1_Hospitalization_Fall_values.csv", FIG_DIR)

## Extended Data Fig. 1
P <- OUT[OUT$event %in% c("hosp", "fall") & OUT$group %in% names(GL)[1:3], ]
P$Group <- factor(P$Group, levels = rev(GL[1:3]))
P$Event <- factor(P$Event, levels = EVENTS[1:2])
P$lab <- paste0("[", formatC(floor(P$Declines + 0.5), format = "d", big.mark = ","), "]")   # half up, as in Supplementary Table 7
## Style follows Figure 1 (theme_na, Okabe-Ito severity palette, bold panel letters):
##   blue = robust, no event (best case); orange = robust, with event (stressor);
##   grey = all frailty states, with event.
COLS <- c("Robust before the event" = SEV_COL[["bad"]], "Robust, no event" = SEV_COL[["good"]],
          "All states, event" = SEV_COL[["mid"]])
YLAB <- c("Robust before the event" = "Robust, with event", "Robust, no event" = "Robust, no event",
          "All states, event" = "All states, with event")
one_panel <- function(ev, ttl) {
  ggplot(P[P$Event == ev, ], aes(x = rr, y = Group, colour = Group)) +
    geom_vline(xintercept = 1, linetype = "22", colour = SEV_COL[["aux"]], linewidth = LINE_PT) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.45) +
    geom_point(size = 1.6) +
    geom_text(aes(x = hi, label = lab), hjust = -0.25, size = BASE_PT / ggplot2::.pt, colour = "black") +
    scale_x_log10(breaks = c(0.4, 0.6, 0.8, 1, 1.25), labels = c("0.4", "0.6", "0.8", "1.0", "1.25"),
                  limits = c(0.33, 1.6)) +
    scale_y_discrete(labels = YLAB) +
    scale_colour_manual(values = COLS, guide = "none") +
    facet_wrap(~ Outcome, nrow = 1) +
    labs(title = ttl, x = "Risk ratio per +1 s.d. of IC (log scale)", y = NULL) +
    theme_na() + theme(plot.title = element_text(size = BASE_PT, face = "bold"),
                       panel.spacing.x = unit(6, "mm"))
}
if (requireNamespace("patchwork", quietly = TRUE)) {
  g <- patchwork::wrap_plots(one_panel(EVENTS[1], "Hospitalization"), one_panel(EVENTS[2], "Fall"), ncol = 1) +
       patchwork::plot_annotation(tag_levels = "a") &
       theme(plot.tag = element_text(size = LAB_PT, face = "bold"))
} else {
  g <- ggplot(P, aes(x = rr, y = Group, colour = Group)) + geom_point() + facet_grid(Event ~ Outcome) + theme_na()
}
save_na(g, "EDFig1_Hospitalization_Fall", width_mm = NA_WIDTH$double, height_mm = 110, dir = FIG_DIR)
