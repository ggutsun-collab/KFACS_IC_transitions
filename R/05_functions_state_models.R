###############################################################################
## 05_functions_state_models.R - shared code for 13b_, 62_ and 64_
##
## Builds, from the long data returned by prep_long(), (i) one row per
## participant at wave 1 with 5-year death and composite-worsening times, and
## (ii) the between-visit intervals for the 18 transitions. Both reproduce the
## published Table 2 (11,571 intervals) and Table 3 per-s.d. estimates.
###############################################################################

SM_VERSION <- "2026-11-17"

## within-sex z of the general factor, standardized on wave 1
ws_scale <- function(d) {
  w1 <- d[!duplicated(d$id), ]; z <- d$gLIC
  for (s in levels(d$sex_f)) {
    k <- d$sex_f == s; k1 <- w1$sex_f == s
    z[k] <- (d$gLIC[k] - mean(w1$gLIC[k1], na.rm = TRUE)) / sd(w1$gLIC[k1], na.rm = TRUE)
  }
  d$gLIC_z <- z; d
}

## additional wave-1 covariates (living alone, polypharmacy, smoking, physical activity)
SM_EXTRA <- c(alone = "living_alone", poly = "polypharmacy",
              smoke = "smoking_current", pa = "pa_met_min_wk")
.sm_cat <- function(x, tertiles = FALSE) {
  if (tertiles) {
    x <- suppressWarnings(as.numeric(x))
    q <- quantile(x, c(1/3, 2/3), na.rm = TRUE, names = FALSE)
    x <- paste0("T", cut(x, c(-Inf, unique(q), Inf), labels = FALSE, include.lowest = TRUE))
    x[x == "TNA"] <- NA
  } else x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  f <- factor(x); tb <- table(f[f != "Unknown"])
  if (length(tb)) f <- relevel(f, ref = names(tb)[which.max(tb)])
  f
}

## one row per participant at wave 1
sm_person_frame <- function(d) {
  d$st  <- as.integer(as_num(d$state))         # ADL 1-3, 4 = death
  d$fst <- as.integer(as_num(d$frail_state))   # frailty 1-3, 4 = death
  d <- d[is.finite(d$time), ]; d <- d[order(d$id, d$time), ]
  vis <- d[!is.na(d$st) & d$st %in% 1:3 & !is.na(d$fst) & d$fst %in% 1:3, ]
  w1  <- vis[!duplicated(vis$id), ]
  dth <- d[(!is.na(d$st) & d$st == 4) | (!is.na(d$fst) & d$fst == 4), ]
  td  <- tapply(dth$time, dth$id, min); tl <- tapply(vis$time, vis$id, max)
  K <- data.frame(id = w1$id, z = w1$gLIC_z, adl0 = w1$st,
                  fr = factor(w1$fst, 1:3, c("Robust", "Pre-frail", "Frail")),
                  age_c = as_num(w1$age_c), sex_f = w1$sex_f,
                  comorbid_bl = as_num(w1$comorbid_bl),
                  edu_bl_f = w1$edu_bl_f, income_bl_f = w1$income_bl_f, area_bl_f = w1$area_bl_f,
                  stringsAsFactors = FALSE)
  K$comorbid_bl[is.na(K$comorbid_bl)] <- median(K$comorbid_bl, na.rm = TRUE)
  K$ev   <- as.integer(K$id %in% names(td))
  K$t_fu <- ifelse(K$ev == 1, td[K$id], tl[K$id])
  later <- vis[vis$time > ave(vis$time, vis$id, FUN = min), ]
  later <- merge(later[, c("id", "time", "st")], K[, c("id", "adl0")], by = "id")
  wr <- later[later$st > later$adl0, ]; tw <- tapply(wr$time, wr$id, min)
  K$t_comp <- pmin(ifelse(K$id %in% names(tw), tw[K$id], Inf), K$t_fu)
  K$e_comp <- as.integer(K$id %in% names(tw) | K$ev == 1)
  for (k in names(SM_EXTRA)) {
    v <- SM_EXTRA[[k]]
    K[[k]] <- .sm_cat(if (v %in% names(w1)) w1[[v]] else NA, tertiles = (k == "pa"))
  }
  K[is.finite(K$z) & is.finite(K$t_fu) & K$t_fu > 0, ]
}

## between-visit intervals for the 18 transitions
sm_intervals <- function(d, K) {
  d$st  <- as.integer(as_num(d$state)); d$fst <- as.integer(as_num(d$frail_state))
  d <- d[is.finite(d$time), ]; d <- d[order(d$id, d$time), ]
  nx <- function(x) c(x[-1], NA)
  same <- d$id == nx(d$id)
  iv <- data.frame(id = d$id, st = d$st, fst = d$fst, z = d$gLIC_z, age_c = as_num(d$age_c),
                   pt = nx(d$time) - d$time, st1 = nx(d$st), fst1 = nx(d$fst),
                   stringsAsFactors = FALSE)[which(same), ]
  iv <- iv[iv$st %in% 1:3 & iv$fst %in% 1:3 & is.finite(iv$z) & is.finite(iv$pt) & iv$pt > 0, ]
  iv$st1[iv$fst1 %in% 4] <- 4; iv$fst1[iv$st1 %in% 4] <- 4
  iv <- iv[!is.na(iv$st1) & !is.na(iv$fst1), ]
  merge(iv, K[, setdiff(names(K), c("z", "age_c", "adl0", "fr", "ev", "t_fu", "t_comp", "e_comp"))],
        by = "id")
}

SM_TR <- local({
  LAB <- list(st = c("Normal", "Mild", "Severe"), fst = c("Robust", "Pre-frail", "Frail"))
  SYS <- c(st = "Disability", fst = "Frailty")
  tr <- do.call(rbind, lapply(names(LAB), function(s) do.call(rbind, lapply(1:3, function(o)
    data.frame(sys = s, from = o, to = c(setdiff(1:3, o), 4), stringsAsFactors = FALSE)))))
  tr$label <- sprintf("%s: %s -> %s", SYS[tr$sys], mapply(function(s, o) LAB[[s]][o], tr$sys, tr$from),
                      mapply(function(s, t) if (t == 4) "Death" else LAB[[s]][t], tr$sys, tr$to))
  tr
})

## cluster-robust (HC0, G/(G-1)) variance for a Poisson glm
.sm_vcov_cl <- function(m, cl) {
  X <- model.matrix(m); X <- X[, !is.na(coef(m)), drop = FALSE]
  r <- m$y - fitted(m)
  B <- solve(crossprod(X * sqrt(fitted(m))))
  U <- rowsum(X * r, cl); G <- nrow(U)
  (G / (G - 1)) * B %*% crossprod(U) %*% B
}
.sm_keep <- function(dd, cov)
  cov[vapply(cov, function(k) !(is.factor(dd[[k]]) && nlevels(droplevels(dd[[k]])) < 2), TRUE)]

## log rate ratio per s.d. of IC and its variance, for one transition
sm_fit_transition <- function(iv, j, cov) {
  s <- SM_TR$sys[j]; dd <- droplevels(iv[iv[[s]] == SM_TR$from[j], ])
  dd$y <- as.integer(dd[[paste0(s, "1")]] == SM_TR$to[j])
  if (sum(dd$y) < 3) return(c(b = NA, v = NA, ev = sum(dd$y)))
  m <- try(glm(reformulate(c("z", .sm_keep(dd, cov), "offset(log(pt))"), response = "y"),
               family = poisson, data = dd), silent = TRUE)
  if (inherits(m, "try-error")) return(c(b = NA, v = NA, ev = sum(dd$y)))
  c(b = unname(coef(m)["z"]), v = .sm_vcov_cl(m, dd$id)["z", "z"], ev = sum(dd$y))
}

## log hazard ratio per s.d. of IC within a wave-1 frailty state, 5-year horizon
sm_fit_within <- function(K, state, outcome, cov, horizon = 5) {
  x <- K[K$fr == state, ]
  tv <- if (outcome == "death") "t_fu" else "t_comp"; ev <- if (outcome == "death") "ev" else "e_comp"
  x$.t <- pmin(x[[tv]], horizon); x$.e <- ifelse(x[[tv]] > horizon, 0L, x[[ev]])
  x <- droplevels(x[x$.t > 0, ])
  m <- try(survival::coxph(as.formula(paste("survival::Surv(.t, .e) ~ z +",
                                            paste(.sm_keep(x, cov), collapse = " + "))), data = x),
           silent = TRUE)
  if (inherits(m, "try-error")) return(c(b = NA, v = NA, ev = sum(x$.e), n = nrow(x)))
  c(b = unname(coef(m)["z"]), v = vcov(m)["z", "z"], ev = sum(x$.e), n = nrow(x))
}

## Rubin's rules on log scale -> est, lo, hi, p
sm_rubin <- function(b, v) {
  ok <- is.finite(b) & is.finite(v); b <- b[ok]; v <- v[ok]; m <- length(b)
  if (!m) return(c(est = NA, lo = NA, hi = NA, p = NA))
  qb <- mean(b); ub <- mean(v); bb <- if (m > 1) var(b) else 0
  tv <- ub + (1 + 1/m) * bb
  df <- if (bb > 0) (m - 1) * (1 + ub / ((1 + 1/m) * bb))^2 else Inf
  se <- sqrt(tv); tq <- qt(0.975, df)
  c(est = qb, lo = qb - tq * se, hi = qb + tq * se, p = 2 * pt(-abs(qb / se), df))
}
sm_fmt <- function(r) sprintf("%.2f (%.2f-%.2f)", exp(r["est"]), exp(r["lo"]), exp(r["hi"]))

SM_COV_PRIMARY <- c("age_c", "comorbid_bl", "edu_bl_f", "income_bl_f", "area_bl_f")
SM_COV_EXTRA   <- c(SM_COV_PRIMARY, "alone", "poly", "smoke", "pa")
SM_WITHIN <- list(c("Robust", "death"), c("Robust", "worsen"),
                  c("Pre-frail", "death"), c("Pre-frail", "worsen"))
