rm(list=ls())

# ==============================================================================
# NSABP B-14 APPLICATION CODE (OSQC: IW, DR; PSQC: PS) with Cox vs RSF censoring models
# - Censoring models:
#     (1) Cox PH model for censoring
#     (2) randomForestSRC for censoring
# - Outcome regression + pi-model: stratified Cox (by trt)
# - Propensity: logistic regression (even though randomized), with truncation
# - Estimators:
#     IW, DR  -> targets OSQC
#     PS  -> targets PSQC
# - Inference: nonparametric bootstrap, Wald SE + 95% Wald CI
# ==============================================================================

# Load libraries
library(survival)
library(doRNG)
library(doParallel)
library(foreach)
library(parallel)
library(randomForestSRC)
library(dplyr)
library(ggplot2)

set.seed(123)

# Load and clean data (This dataset is not publicaly available and will not be provided)
b14 = read.csv("nsabp.b14.data.csv")[,-1]  

# KM analysis by treatment
b14$trt_f <- factor(b14$trt, levels = c(0, 1), labels = c("Placebo", "Tamoxifen"))
fit_km <- survfit(Surv(time, event) ~ trt_f, data = b14)

# Ensure numeric types
b14 <- within(b14, {
  trt   <- as.integer(trt)
  event <- as.integer(event)
  age   <- as.numeric(age)
  tsize <- as.numeric(tsize)
  time  <- as.numeric(time)
})

# Ensure numeric types
b14 <- within(b14, {
  trt   <- as.integer(trt)
  event <- as.integer(event)
  age   <- as.numeric(age)
  tsize <- as.numeric(tsize)
  time  <- as.numeric(time)
})

# Helper: Extract baseline hazard from Cox model
get_bh_lookup <- function(cox_model) {
  bh <- basehaz(cox_model, centered = FALSE)
  if ("strata" %in% colnames(bh)) {
    strata_vals <- unique(bh$strata)
    bh_list <- lapply(strata_vals, function(s) bh[bh$strata == s, , drop=FALSE])
    names(bh_list) <- strata_vals
    for (nm in names(bh_list)) bh_list[[nm]] <- bh_list[[nm]][order(bh_list[[nm]]$time), , drop=FALSE]
    return(list(strata = TRUE, bh_list = bh_list))
  } else {
    bh <- bh[order(bh$time), ]
    return(list(strata = FALSE, t = bh$time, h = bh$hazard))
  }
}

# Helper: Predict survival from Cox model
cox_surv_at <- function(bh_lookup, lp, times, strata_key = NULL) {
  n <- length(times)
  out <- rep(NA_real_, n)
  if (!bh_lookup$strata) {
    idx <- findInterval(times, bh_lookup$t)
    H0 <- rep(0, n)
    ok <- idx > 0
    H0[ok] <- bh_lookup$h[idx[ok]]
    out <- exp(-H0 * exp(lp))
    return(out)
  }
  bh_list <- bh_lookup$bh_list
  for (s in unique(strata_key)) {
    rows <- which(strata_key == s)
    bh_s <- bh_list[[s]]
    idx <- findInterval(times[rows], bh_s$time)
    H0 <- rep(0, length(rows))
    ok <- idx > 0
    H0[ok] <- bh_s$hazard[idx[ok]]
    out[rows] <- exp(-H0 * exp(lp[rows]))
  }
  out
}

# Helper: Predict survival from Random Survival Forest
rsf_surv_at <- function(rsf_fit, newdata, times) {
  pred <- predict(rsf_fit, newdata = newdata)
  S_mat <- pred$survival
  tgrid <- pred$time.interest
  n <- nrow(S_mat)
  m <- ncol(S_mat)
  tt <- pmin(as.numeric(times), max(tgrid))
  col_idx <- findInterval(tt, tgrid, rightmost.closed = TRUE)
  col_idx[col_idx < 1] <- 1
  col_idx[col_idx > m] <- m
  out <- rep(NA_real_, n)
  ok <- which(is.finite(col_idx) & !is.na(col_idx))
  out[ok] <- S_mat[cbind(ok, col_idx[ok])]
  return(out)
}

# Helper: Robust root solver for theta
solve_root_theta <- function(f, upper_init = 5, upper_max = 20, n_grid = 40, tol = 1e-4) {
  upper <- upper_init
  for (iter in 1:10) {
    grid <- seq(0, upper, length.out = n_grid)
    vals <- suppressWarnings(sapply(grid, f))
    ok <- is.finite(vals)
    if (sum(ok) >= 2) {
      grid2 <- grid[ok]; vals2 <- vals[ok]
      sgn <- sign(vals2)
      idx <- which(diff(sgn) != 0)
      if (length(idx) > 0) {
        return(tryCatch(uniroot(f, c(grid2[idx[1]], grid2[idx[1] + 1]), tol = tol)$root,
                        error = function(e) lo))
      }
    }
    if (upper >= upper_max) {
      f0 <- suppressWarnings(f(0))
      fu <- suppressWarnings(f(upper))
      if (!is.finite(f0)) f0 <- Inf
      if (!is.finite(fu)) fu <- Inf
      return(ifelse(abs(f0) <= abs(fu), 0, upper))
    }
    upper <- min(upper * 1.8, upper_max)
  }
  0
}

# Core function: Estimate IW, DR, and PS contrasts
estimate_one <- function(dat, t0, tau, censor_method = c("cox", "rsf"), trunc_G = 0.01) {
  censor_method <- match.arg(censor_method)
  dat <- dat[complete.cases(dat[, c("trt","age","tsize","time","event")]), ]
  dat$trt <- as.integer(dat$trt)
  dat$event <- as.integer(dat$event)
  n <- nrow(dat)
  
  # Propensity score model
  f_ps <- as.formula(trt ~ age + tsize)
  mod_ps <- glm(f_ps, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(mod_ps, type = "response"))
  ps_hat <- pmin(pmax(ps_hat, 0.01), 0.99)
  
  # Censoring models
  dat$Dc <- 1L - dat$event
  if (censor_method == "cox") {
    f_c <- as.formula(Surv(time, Dc) ~ trt + age + tsize)
    mod_c <- coxph(f_c, data = dat, x = TRUE)
    bh_c <- get_bh_lookup(mod_c)
    lp_c <- predict(mod_c, newdata = dat, type = "lp")
    Ghat_time <- cox_surv_at(bh_c, lp_c, dat$time)
    Ghat_t0   <- cox_surv_at(bh_c, lp_c, rep(t0, n))
  } else {
    rsf_c <- rfsrc(Surv(time, Dc) ~ trt + age + tsize, data = dat, ntree = 200, nodesize = 15, forest = TRUE)
    Ghat_time <- rsf_surv_at(rsf_c, newdata = dat, times = dat$time)
    Ghat_t0   <- rsf_surv_at(rsf_c, newdata = dat, times = rep(t0, n))
  }
  w_cen_time <- 1 / pmax(Ghat_time, trunc_G)
  w_cen_t0   <- 1 / pmax(Ghat_t0, trunc_G)
  
  # Outcome model (Stratified Cox)
  f_out <- as.formula(Surv(time, event) ~ strata(trt) + age + tsize)
  mod_out <- coxph(f_out, data = dat, x = TRUE)
  bh_out <- get_bh_lookup(mod_out)
  X_out <- model.matrix(~ age + tsize, data = dat)
  if ("(Intercept)" %in% colnames(X_out)) X_out <- X_out[, -1, drop = FALSE]
  lp_out <- as.vector(X_out %*% coef(mod_out))
  
  strata_key_0 <- NULL
  strata_key_1 <- NULL
  if (bh_out$strata) {
    strata_names <- names(bh_out$bh_list)
    s_num <- suppressWarnings(as.numeric(gsub(".*=|\\D", "", strata_names)))
    ord <- order(s_num, na.last = TRUE)
    strata_key_0 <- rep(strata_names[ord][1], n)
    strata_key_1 <- rep(strata_names[ord][2], n)
  }
  
  # Pi-model for PS weights
  dat_pi <- dat
  dat_pi$time_pi  <- pmin(dat_pi$time, t0)
  dat_pi$event_pi <- as.integer(dat_pi$time <= t0 & dat_pi$event == 1)
  f_pi <- as.formula(Surv(time_pi, event_pi) ~ strata(trt) + age + tsize)
  mod_pi <- coxph(f_pi, data = dat_pi, x = TRUE)
  bh_pi <- get_bh_lookup(mod_pi)
  lp_pi <- as.vector(X_out %*% coef(mod_pi))
  strata_names_pi <- names(bh_pi$bh_list)
  s_num_pi <- suppressWarnings(as.numeric(gsub(".*=|\\D", "", strata_names_pi)))
  ord_pi <- order(s_num_pi, na.last = TRUE)
  S_pi_t0_A0 <- cox_surv_at(bh_pi, lp_pi, rep(t0, n), strata_key = rep(strata_names_pi[ord_pi][1], n))
  S_pi_t0_A1 <- cox_surv_at(bh_pi, lp_pi, rep(t0, n), strata_key = rep(strata_names_pi[ord_pi][2], n))
  
  # Principal score weights
  ratio <- exp(log(pmax(S_pi_t0_A0, 1e-12)) - log(pmax(S_pi_t0_A1, 1e-12)))
  w_sel <- rep(1, n)
  w_sel[(dat$trt == 1) & (dat$time > t0)] <- ratio[(dat$trt == 1) & (dat$time > t0)]
  
  w_trt <- function(a) if (a == 1) 1 / ps_hat else 1 / (1 - ps_hat)
  
  # Estimating functions: IW, DR, PS
  ee_iw <- function(theta, a) {
    wa <- w_trt(a)
    Ia <- as.integer(dat$trt == a)
    num <- Ia * wa * as.integer(dat$time > t0 & dat$time <= (t0 + theta) & dat$event == 1) * w_cen_time
    den <- Ia * wa * as.integer(dat$time > t0) * w_cen_t0
    mean(num) - tau * mean(den)
  }
  
  ee_dr <- function(theta, a) {
    if (a == 1) {
      S_theta <- cox_surv_at(bh_out, lp_out, times = rep(t0 + theta, n), strata_key = strata_key_1)
      S_base <- S_pi_t0_A1
      w  <- 1 / ps_hat
      aw <- (dat$trt - ps_hat) / ps_hat
    } else {
      S_theta <- cox_surv_at(bh_out, lp_out, times = rep(t0 + theta, n), strata_key = strata_key_0)
      S_base <- S_pi_t0_A0
      w  <- 1 / (1 - ps_hat)
      aw <- ((1 - dat$trt) - (1 - ps_hat)) / (1 - ps_hat)
    }
    m_ax <- (S_base - S_theta) - tau * S_base
    idx_n <- (dat$trt == a & dat$time > t0 & dat$time <= (t0 + theta) & dat$event == 1)
    idx_d <- (dat$trt == a & dat$time > t0)
    mean(w * (as.integer(idx_n) * w_cen_time - tau * as.integer(idx_d) * w_cen_t0) - aw * m_ax)
  }
  
  ee_ps <- function(theta, a) {
    Ia <- as.integer(dat$trt == a)
    wa <- w_trt(a)
    w_extra <- rep(1, n)
    if (a == 1) w_extra <- w_sel
    num <- Ia * wa * w_extra * as.integer(dat$time > t0 & dat$time <= (t0 + theta) & dat$event == 1) * w_cen_time
    den <- Ia * wa * w_extra * as.integer(dat$time > t0) * w_cen_t0
    mean(num) - tau * mean(den)
  }
  
  resid_all <- dat$time[dat$time > t0] - t0
  upper_init <- ifelse(length(resid_all) >= 10, max(quantile(resid_all, 0.95, na.rm=TRUE), 2), 5)
  
  q1_iw <- solve_root_theta(function(th) ee_iw(th, 1), upper_init = upper_init)
  q0_iw <- solve_root_theta(function(th) ee_iw(th, 0), upper_init = upper_init)
  q1_dr <- solve_root_theta(function(th) ee_dr(th, 1), upper_init = upper_init)
  q0_dr <- solve_root_theta(function(th) ee_dr(th, 0), upper_init = upper_init)
  q1_ps <- solve_root_theta(function(th) ee_ps(th, 1), upper_init = upper_init)
  q0_ps <- solve_root_theta(function(th) ee_ps(th, 0), upper_init = upper_init)
  
  c(IW = q1_iw - q0_iw, DR = q1_dr - q0_dr, PS = q1_ps - q0_ps)
}

# Main function: Wrapper for landmark analysis with bootstrap
analyze_B14 <- function(dat, landmarks = c(3, 5), taus = c(0.3), B = 1000, 
                        censor_methods = c("cox", "rsf"), trunc_G = 0.01) {
  set.seed(123)
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  out <- list()
  for (cm in censor_methods) {
    for (t0 in landmarks) {
      for (tau in taus) {
        theta_hat <- estimate_one(dat, t0 = t0, tau = tau, censor_method = cm, trunc_G = trunc_G)
        n <- nrow(dat)
        boot_mat <- foreach(b = 1:B, .combine = rbind, .packages = c("survival", "randomForestSRC"),
                            .export = c("estimate_one", "get_bh_lookup", "cox_surv_at", "rsf_surv_at", "solve_root_theta")) %dorng% {
                              idx <- sample.int(n, n, replace = TRUE)
                              dat_b <- dat[idx, , drop = FALSE]
                              if (length(unique(dat_b$trt)) < 2) return(rep(NA_real_, length(theta_hat)))
                              tryCatch(estimate_one(dat_b, t0 = t0, tau = tau, censor_method = cm, trunc_G = trunc_G),
                                       error = function(e) rep(NA_real_, length(theta_hat)))
                            }
        colnames(boot_mat) <- names(theta_hat)
        se <- apply(boot_mat, 2, sd, na.rm = TRUE)
        out[[paste(cm, t0, tau, sep = "_")]] <- data.frame(
          censor_model = cm, t0 = t0, tau = tau, estimator = names(theta_hat),
          estimate = as.numeric(theta_hat), boot_se = as.numeric(se),
          ci_lo = as.numeric(theta_hat - 1.96 * se), ci_hi = as.numeric(theta_hat + 1.96 * se),
          B_eff = colSums(is.finite(boot_mat))
        )
        cat(sprintf("[DONE] censor=%s | t0=%g | tau=%g\n", cm, t0, tau))
      }
    }
  }
  do.call(rbind, out)
}

# Run analysis
results <- analyze_B14(b14, landmarks = c(1,2,3,4,5), taus = c(0.3), B = 1000, censor_methods = c("cox", "rsf"))
