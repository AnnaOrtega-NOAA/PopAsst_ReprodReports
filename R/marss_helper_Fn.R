### MARSS Helper Functions for Population Trend Analysis
### PIFSC 2026 Update

run_turtle_marss <- function(data_raw, years_to_include) {
  # Prep: Rows = Sites, Cols = Years
  dat_subset <- data_raw[data_raw$Season %in% years_to_include, ]
  marss_dat  <- t(log(dat_subset[, -1] + 0.1)) # 0.1 offset for zeros
  colnames(marss_dat) <- dat_subset$Season
  
  # # Structure: Shared U (growth), Shared Q (process), Shared R (obs error)
  # model.list <- list(
  #   Z = matrix(1, nrow(marss_dat), 1), 
  #   U = "equal", Q = "diagonal and equal", R = "diagonal and equal",
  #   B = "identity", A = "scaling", x0 = "unequal"
  # )
  
  # --- CHANGE TO INDEPENDENT STRUCTURE ---
  model.list <- list(
    Z = "identity",           # Each beach is its own hidden state
    U = "unequal",            # Each beach gets its own growth rate (U)
    Q = "diagonal and equal", # Keep variance shared for stability, or "diagonal and unequal"
    R = "diagonal and equal", # Observation error
    B = "identity", 
    A = "zero",               # No scaling needed when Z is identity
    x0 = "unequal"
  )
  
  fit <- MARSS(marss_dat, model = model.list, method = "kem", silent = TRUE)
  return(MARSSparamCIs(fit))
}

# get_marss_stats <- function(fit, spp_name, remig) {
#   # Extracts
#   u_est <- as.numeric(coef(fit, type="matrix")$U)
#   u_low <- as.numeric(fit$par.low$U)
#   u_up  <- as.numeric(fit$par.up$U)
#   
#   # Abundance
#   last_yr <- ncol(fit$model$data)
#   n_nesters <- sum(exp(fit$states[, last_yr]))
#   total_pop <- n_nesters * remig
#   
#   data.frame(
#     Metric = c("Annual Growth Rate (%)", "Final Year Nesters", "Total Reproductive Pop", "Process Var (Q)"),
#     Median = c(round((exp(u_est)-1)*100, 2), round(n_nesters, 0), round(total_pop, 0), round(as.numeric(coef(fit, type="matrix")$Q), 4)),
#     LCI_95 = c(round((exp(u_low)-1)*100, 2), NA, NA, NA),
#     UCI_95 = c(round((exp(u_up)-1)*100, 2), NA, NA, NA)
#   )
# }

# --- CHANGE TO INDEPENDENT STRUCTURE ---
get_marss_stats <- function(fit, spp_name, remig, period_label) {
  u_estimates  <- coef(fit, type="matrix")$U
  site_names   <- rownames(fit$model$data)
  growth_rates <- round((exp(u_estimates) - 1) * 100, 2)
  
  # Get Confidence Intervals for U
  u_low <- fit$par.low$U
  u_up  <- fit$par.up$U
  lci_rates <- round((exp(u_low) - 1) * 100, 2)
  uci_rates <- round((exp(u_up) - 1) * 100, 2)
  
  last_yr_idx    <- ncol(fit$model$data)
  site_abundance <- exp(fit$states[, last_yr_idx])
  
  data.frame(
    Species = spp_name,
    Site = site_names,
    Period = period_label,
    Growth_Rate_Pct = as.vector(growth_rates),
    LCI_95 = as.vector(lci_rates),
    UCI_95 = as.vector(uci_rates),
    End_Year_Abundance = round(as.vector(site_abundance), 0),
    Total_Repro_Pop = round(as.vector(site_abundance * remig), 0)
  )
}

plot_marss_trend <- function(fit, spp_name, fig_path) {
  years <- as.numeric(colnames(fit$model$data))
  obs   <- colSums(exp(fit$model$data), na.rm=TRUE)
  pred  <- colSums(exp(fit$states))
  
  png(paste0(fig_path, spp_name, "_MARSS_Trend.png"), width=1800, height=1200, res=300)
  plot(years, obs, pch=16, col="black", xlab="Season", ylab="Annual Nesters",
       main=paste(spp_name, "Recent Trend"))
  lines(years, pred, col="blue", lwd=2)
  legend("topright", legend=c("Observed", "MARSS Fit"), pch=c(16, NA), lty=c(NA, 1), col=c("black", "blue"))
  dev.off()
}

run_loo_years <- function(data_raw) {
  years <- data_raw$Season
  results <- data.frame()
  for(y in years) {
    dat_loo <- data_raw
    dat_loo[dat_loo$Season == y, -1] <- NA
    fit <- run_turtle_marss(dat_loo, years)
    results <- rbind(results, data.frame(Removed_Year = y, 
                                         Growth_Rate = round((exp(as.numeric(coef(fit, type="matrix")$U))-1)*100, 2)))
  }
  return(results)
}