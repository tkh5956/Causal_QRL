rm(list=ls())

# ==============================================================================
# SUPPORT APPLICATION CODE (OSQC: IW, DR) with Cox censoring model
# - Censoring models: Cox PH model for censoring
# - Outcome regression + pi-model: stratified Cox (by trt)
# - Propensity: logistic regression with truncation
# - Estimators (PS will not be utilized due to violation of monotonicity:
#     IW, DR  -> targets OSQC
# - Inference: nonparametric bootstrap, Wald SE + 95% Wald CI
# ==============================================================================
# ---------------------------------------------------------------------------
# Estimation of Residual Lifetime Quantiles: SUPPORT RHC Study
# ---------------------------------------------------------------------------

rhc = read.csv("https://hbiostat.org/data/repo/rhc.csv", header = TRUE)[,-1] # Remove index

### Missing value processing
# Expand disease categories (cat1/cat2) into dummy variables
rhc <- rhc %>% 
  mutate(val = 1) %>%
  tidyr::pivot_wider(names_from = cat1, values_from = val, values_fill = 0, names_prefix = "cat_") 

rhc <- rhc %>% mutate(
  cat_Cirrhosis = if_else((!is.na(cat2))&(cat2=="Cirrhosis"),1,cat_Cirrhosis),
  `cat_Colon Cancer` = if_else((!is.na(cat2))&(cat2=="Colon Cancer"),1,`cat_Colon Cancer`),
  cat_Coma = if_else((!is.na(cat2))&(cat2=="Coma"),1,cat_Coma),
  `cat_Lung Cancer` = if_else((!is.na(cat2))&(cat2=="Lung Cancer"),1,`cat_Lung Cancer`),
  `cat_MOSF w/Malignancy` = if_else((!is.na(cat2))&(cat2=="MOSF w/Malignancy"),1,`cat_MOSF w/Malignancy`),
  `cat_MOSF w/Sepsis` = if_else((!is.na(cat2))&(cat2=="MOSF w/Sepsis"),1,`cat_MOSF w/Sepsis`)
) %>% dplyr::select(-cat2)

# Handle numerical missingness: ADL and Urine Output (Mean-fill + Indicator)
rhc <- rhc %>%
  mutate(
    ADL_miss = ifelse(is.na(adld3p), 1, 0),
    ADL_filled = ifelse(is.na(adld3p), mean(adld3p, na.rm = TRUE), adld3p),
    Urine_miss = ifelse(is.na(urin1), 1, 0),
    Urine_filled = ifelse(is.na(urin1), mean(urin1, na.rm = TRUE), urin1)
  ) %>% dplyr::select(-adld3p, -urin1)

# Set treatment and mortality indicators
rhc$treatment = as.numeric(rhc$swang1=="RHC"); rhc = rhc %>% dplyr::select(-swang1)
rhc$deathn = as.numeric(rhc$death=="Yes")

# Clean column names
names(rhc) <- gsub(" ", "_", names(rhc))
names(rhc) <- gsub("w/", "", names(rhc))

#-----------------------------
# 0) Utility helpers
#-----------------------------

# Extract baseline hazard for Cox models (supports stratified models)
get_bh_lookup <- function(cox_model) {
  bh <- basehaz(cox_model, centered = FALSE)
  
  if ("strata" %in% colnames(bh)) {
    slev <- unique(bh$strata)
    s0 <- slev[grep("(=|:)0\\b", slev)][1]
    s1 <- slev[grep("(=|:)1\\b", slev)][1]
    
    if (is.na(s0) || is.na(s1)) { s0 <- slev[1]; s1 <- slev[2] }
    
    bh0 <- bh[bh$strata == s0, , drop = FALSE]
    bh1 <- bh[bh$strata == s1, , drop = FALSE]
    bh0 <- bh0[order(bh0$time), , drop = FALSE]
    bh1 <- bh1[order(bh1$time), , drop = FALSE]
    
    return(list(strata = TRUE, t0 = bh0$time, h0 = bh0$hazard, t1 = bh1$time, h1 = bh1$hazard))
  } else {
    bh <- bh[order(bh$time), , drop = FALSE]
    return(list(strata = FALSE, t = bh$time, h = bh$hazard))
  }
}

# Manually compute survival probabilities from Cox nuisance models
calc_surv_manual <- function(bh_lookup, lp, times, treatment_arm = NULL) {
  times <- pmax(times, 0)
  cumhaz <- numeric(length(times))
  
  if (bh_lookup$strata) {
    if (is.null(treatment_arm)) stop("treatment_arm required for stratified hazards.")
    idx0 <- which(treatment_arm == 0L); idx1 <- which(treatment_arm == 1L)
    
    if (length(idx0) > 0) {
      k0 <- findInterval(times[idx0], bh_lookup$t0)
      v0 <- numeric(length(idx0)); m0 <- k0 > 0; v0[m0] <- bh_lookup$h0[k0[m0]]
      cumhaz[idx0] <- v0
    }
    if (length(idx1) > 0) {
      k1 <- findInterval(times[idx1], bh_lookup$t1)
      v1 <- numeric(length(idx1)); m1 <- k1 > 0; v1[m1] <- bh_lookup$h1[k1[m1]]
      cumhaz[idx1] <- v1
    }
  } else {
    k <- findInterval(times, bh_lookup$t); m <- k > 0; cumhaz[m] <- bh_lookup$h[k[m]]
  }
  exp(-cumhaz * exp(lp))
}

# Robust grid-based root finder for non-smooth quantile equations
solve_fast <- function(f, upper_init, upper_max, n_grid = 50) {
  upper <- upper_init
  repeat {
    grid <- seq(0, upper, length.out = n_grid)
    vals <- suppressWarnings(sapply(grid, f)); ok <- is.finite(vals)
    
    if (sum(ok) >= 2) {
      grid_ok <- grid[ok]; vals_ok <- vals[ok]
      idx <- which(diff(sign(vals_ok)) != 0)
      if (length(idx) > 0) {
        return(tryCatch(uniroot(f, c(grid_ok[idx[1]], grid_ok[idx[1] + 1]), tol = 1e-4)$root,
                        error = function(e) grid_ok[idx[1]]))
      }
    }
    if (upper >= upper_max) return(ifelse(all(vals[is.finite(vals)] > 0), 0, upper_max))
    upper <- min(upper * 1.8, upper_max)
  }
}

#-----------------------------
# 1) Data prep for RHC
#-----------------------------

# Process dates and follow-up times for landmarking
prep_rhc <- function(rhc, horizon_days = 180) {
  rhc <- as.data.frame(rhc)
  rhc$sadmdte <- as.Date(rhc$sadmdte); rhc$dthdte <- as.Date(rhc$dthdte); rhc$lstctdte <- as.Date(rhc$lstctdte)
  rhc = rhc %>% mutate(lstctdte = if_else((!is.na(dthdte))&(lstctdte<dthdte), dthdte, lstctdte))
  
  t_death <- as.numeric(rhc$dthdte - rhc$sadmdte); t_cens <- as.numeric(rhc$lstctdte - rhc$sadmdte)
  t_cens[is.na(t_cens)] <- horizon_days
  
  rhc$Y <- pmin(t_death, t_cens, horizon_days, na.rm = TRUE)
  rhc$Delta <- ifelse(!is.na(t_death) & (t_death <= t_cens) & (t_death <= horizon_days), 1L, 0L)
  rhc$treatment <- as.integer(rhc$treatment)
  rhc <- rhc[is.finite(rhc$Y) & rhc$Y >= 0, ]
  
  char_cols <- sapply(rhc, is.character); rhc[char_cols] <- lapply(rhc[char_cols], factor)
  rhc
}

#-----------------------------
# 2) Main estimator function
#-----------------------------

# Calculates IW and DR residual-life quantiles
rhc_landmark_quantile <- function(rhc, t0_days = 3, tau = 0.5, horizon_days = 180, trunc_g = 0.001, caliper_ps = 1e-6, verbose = TRUE, plot = FALSE) {
  theta_max <- horizon_days - t0_days; rhc <- prep_rhc(rhc, horizon_days = horizon_days)
  
  covars <- c("age", "sex", "race", "edu", "income", "ninsclas", colnames(rhc)[grep("cat_",colnames(rhc))],
              "resp","card","neuro","gastr","renal","meta","hema","seps","trauma","ortho",
              "ADL_miss","ADL_filled","das2d3pc","dnr1","ca","surv2md1","aps1","scoma1",
              "wtkilo1","temp1","meanbp1","resp1","hrt1","pafi1","paco21","ph1","wblc1","hema1",
              "sod1","pot1","crea1","bili1","alb1","Urine_miss","Urine_filled",
              "cardiohx","chfhx","dementhx","psychhx","chrpulhx","renalhx","liverhx","gibledhx",
              "malighx","immunhx","transhx","amihx")
  
  covars <- covars[covars %in% names(rhc)]
  dat <- rhc[complete.cases(rhc[, unique(c("treatment", "Y", "Delta", covars))]), ]
  
  # (1) Propensity model
  f_ps <- as.formula(paste("treatment ~", paste(covars, collapse = " + ")))
  mod_ps <- glm(f_ps, data = dat, family = binomial())
  ps_hat <- pmin(pmax(predict(mod_ps, type = "response"), caliper_ps), 1 - caliper_ps)
  
  # (2) Censoring model
  f_c <- as.formula(paste("Surv(Y, 1-Delta) ~ treatment +", paste(covars, collapse = " + ")))
  mod_c <- coxph(f_c, data = dat, x = TRUE)
  bh_c <- get_bh_lookup(mod_c); lp_c <- predict(mod_c, type = "lp")
  
  ghat_Y  <- calc_surv_manual(bh_c, lp_c, dat$Y, treatment_arm = NULL)
  ghat_t0 <- calc_surv_manual(bh_c, lp_c, rep(t0_days, nrow(dat)), treatment_arm = NULL)
  w_cen_Y <- 1 / pmax(ghat_Y, trunc_g); w_cen_t0 <- 1 / pmax(ghat_t0, trunc_g)
  
  # (3) Outcome model
  f_out <- as.formula(paste("Surv(Y, Delta) ~ strata(treatment) +", paste(covars, collapse = " + ")))
  mod_out <- coxph(f_out, data = dat, x = TRUE)
  bh_out <- get_bh_lookup(mod_out)
  
  X_out <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), data = dat)
  if ("(Intercept)" %in% colnames(X_out)) X_out <- X_out[, -1, drop = FALSE]
  lp_out <- as.vector(X_out %*% coef(mod_out))
  
  S_t0_A1 <- calc_surv_manual(bh_out, lp_out, rep(t0_days, nrow(dat)), treatment_arm = rep(1L, nrow(dat)))
  S_t0_A0 <- calc_surv_manual(bh_out, lp_out, rep(t0_days, nrow(dat)), treatment_arm = rep(0L, nrow(dat)))
  
  # Estimating equations
  ee_iw <- function(theta, arm) {
    idx <- which(dat$treatment == arm); w_trt <- if (arm == 1L) 1/ps_hat[idx] else 1/(1-ps_hat[idx])
    num <- w_trt * (dat$Y[idx] > t0_days & dat$Y[idx] <= t0_days + theta & dat$Delta[idx] == 1L) * w_cen_Y[idx]
    den <- w_trt * (dat$Y[idx] > t0_days) * w_cen_t0[idx]
    mean(num) - tau * mean(den)
  }
  
  ee_dr <- function(theta, arm) {
    S_theta <- calc_surv_manual(bh_out, lp_out, rep(t0_days + theta, nrow(dat)), treatment_arm = rep(arm, nrow(dat)))
    S_base <- if (arm == 1L) S_t0_A1 else S_t0_A0
    or_term <- (S_base - S_theta) - tau * S_base
    w <- if (arm == 1L) 1/ps_hat else 1/(1-ps_hat)
    aw <- if (arm == 1L) (dat$treatment - ps_hat)/ps_hat else ((1-dat$treatment) - (1-ps_hat))/(1-ps_hat)
    idx_n <- (dat$treatment == arm & dat$Y > t0_days & dat$Y <= t0_days + theta & dat$Delta == 1L)
    idx_d <- (dat$treatment == arm & dat$Y > t0_days)
    mean(w * (idx_n * w_cen_Y - tau * idx_d * w_cen_t0) - aw * or_term)
  }
  
  # Solve
  q1_iw <- solve_fast(function(th) ee_iw(th, 1L), theta_max, theta_max)
  q0_iw <- solve_fast(function(th) ee_iw(th, 0L), theta_max, theta_max)
  q1_dr  <- solve_fast(function(th) ee_dr(th, 1L), theta_max, theta_max)
  q0_dr  <- solve_fast(function(th) ee_dr(th, 0L), theta_max, theta_max)
  
  list(t0_days = t0_days, tau = tau, q1 = c(IW = q1_iw, DR = q1_dr), q0 = c(IW = q0_iw, DR = q0_dr),
       contrast = c(IW = q1_iw - q0_iw, DR = q1_dr - q0_dr))
}

# 4) Bootstrap Inference
rhc_landmark_quantile_boot <- function(rhc, t0_days = 3, tau = 0.5, horizon_days = 180, B = 500, ci_level = 0.95, n_cores = max(1, parallel::detectCores() - 1)) {
  
  fit0 <- rhc_landmark_quantile(rhc, t0_days, tau, horizon_days, verbose = FALSE)
  theta0 <- c(DR = unname(fit0$contrast["DR"]), IW = unname(fit0$contrast["IW"]))
  n <- nrow(rhc)
  
  cl <- parallel::makeCluster(n_cores); doParallel::registerDoParallel(cl)
  boot_mat <- foreach::foreach(b = 1:B, .combine = rbind, .packages = c("survival", "dplyr", "tidyr"),
                               .export = c("rhc_landmark_quantile", "prep_rhc", "get_bh_lookup", "calc_surv_manual", "solve_fast")) %dopar% {
                                 set.seed(123 + b)
                                 rhc_b <- rhc[sample.int(n, size = n, replace = TRUE), ]
                                 fb <- tryCatch(rhc_landmark_quantile(rhc_b, t0_days, tau, horizon_days, verbose = FALSE), error = function(e) NULL)
                                 if (is.null(fb)) return(c(DR = NA, IW = NA))
                                 c(DR = unname(fb$contrast["DR"]), IW = unname(fb$contrast["IW"]))
                               }
  parallel::stopCluster(cl)
  
  cat("\n==== Results: t0 =", t0_days, ", tau =", tau, "====\n")
  for (m in c("DR", "IW")) {
    vals <- boot_mat[, m]; vals <- vals[is.finite(vals)]; se <- sd(vals)
    cat(m, "Estimator:\n")
    cat("  Point Est:", round(theta0[m], 4), "\n")
    cat("  Std. Error:", round(se, 4), "\n")
    cat("  Wald 95% CI: [", round(theta0[m] - 1.96*se, 4), ",", round(theta0[m] + 1.96*se, 4), "]\n\n")
  }
}

# Execution for t0 = 3, 7, 14
rhc_landmark_quantile_boot(rhc, t0_days = 3, tau = 0.3, B = 1000)
rhc_landmark_quantile_boot(rhc, t0_days = 7, tau = 0.3, B = 1000)
rhc_landmark_quantile_boot(rhc, t0_days = 14, tau = 0.3, B = 1000)
