library(MARSS)

run_turtle_marss <- function(data_raw, years_to_include) {
  dat_subset <- data_raw[data_raw$Season %in% years_to_include, ]
  marss_dat  <- t(log(dat_subset[, -1] + 0.01)) 
  colnames(marss_dat) <- dat_subset$Season
  
  model.list <- list(Z = "identity", U = "unequal", Q = "diagonal and equal", 
                     R = "diagonal and equal", B = "identity", A = "zero", x0 = "unequal")
  
  fit <- MARSS(marss_dat, model = model.list, method = "kem", 
               control = list(maxit = 5000, allow.degen = TRUE), silent = TRUE)
  return(MARSSparamCIs(fit))
}

get_marss_stats <- function(fit, spp_name, remig, period_label) {
  u_estimates <- coef(fit, type="matrix")$U
  site_names <- rownames(fit$model$data)
  growth_rates <- round((exp(u_estimates) - 1) * 100, 2)
  u_low <- fit$par.low$U; u_up <- fit$par.up$U
  lci_rates <- round((exp(u_low) - 1) * 100, 2)
  uci_rates <- round((exp(u_up) - 1) * 100, 2)
  last_yr_idx <- ncol(fit$model$data)
  site_abundance <- exp(fit$states[, last_yr_idx])
  
  data.frame(Species = spp_name, Site = site_names, Period = period_label,
             Growth_Pct = as.vector(growth_rates), LCI_95 = as.vector(lci_rates), 
             UCI_95 = as.vector(uci_rates), Ann_Nesters = round(as.vector(site_abundance), 0),
             Total_Pop = round(as.vector(site_abundance * remig), 0))
}

plot_marss_trend <- function(fit, spp_name, fig_path) {
  # --- STEP 1: ISOLATE DATA ---
  d_mat  <- fit$model$data
  s_mat  <- fit$states
  se_mat <- MARSSkfss(fit)$xtT.se
  
  # Ensure we extract the time axis correctly from the original data
  # If colnames are NULL, we fallback to a sequence based on the data width
  time_axis <- colnames(d_mat)
  if(is.null(time_axis)) {
    time_axis <- 1:ncol(d_mat)
  } else {
    time_axis <- as.numeric(time_axis)
  }
  
  n_rows    <- nrow(d_mat)
  site_labs <- rownames(d_mat)
  if(is.null(site_labs)) site_labs <- paste("Site", 1:n_rows)
  
  tag       <- ifelse(length(time_axis) > 15, "Full", "Recent")
  
  # --- STEP 2: PLOT ---
  png(paste0(fig_path, spp_name, "_Grouped_Trends_", tag, ".png"), 
      width=2400, height=800 * n_rows, res=300)
  
  par(mfrow=c(n_rows, 1), mar=c(3,4,2,1), oma=c(2,0,2,0))
  
  for(i in 1:n_rows){
    # Use drop=FALSE or explicit indexing to ensure we get a vector
    y_obs  <- as.vector(exp(d_mat[i, ]))
    y_fit  <- as.vector(exp(s_mat[i, ]))
    y_up   <- as.vector(exp(s_mat[i, ] + 1.96 * se_mat[i, ]))
    y_low  <- as.vector(exp(s_mat[i, ] - 1.96 * se_mat[i, ]))
    
    # DEBUGGING CRITICAL CHECK
    if(length(time_axis) != length(y_obs)){
      stop(paste0("Dimension Mismatch in ", spp_name, ": Year axis has ", 
                  length(time_axis), " points, but data has ", length(y_obs)))
    }
    
    plot(x = time_axis, y = y_obs, pch=16, 
         main=paste(spp_name, "-", site_labs[i]), 
         xlab="", ylab="Nesters", las=1, 
         ylim=c(0, max(c(y_obs, y_up), na.rm=T) * 1.2))
    
    # Polygon requires non-NA values to draw correctly
    # Filter out NAs for the polygon to avoid 'hole' errors
    valid_idx <- which(!is.na(y_up) & !is.na(y_low))
    polygon(x = c(time_axis[valid_idx], rev(time_axis[valid_idx])), 
            y = c(y_up[valid_idx], rev(y_low[valid_idx])), 
            col=rgb(0, 0, 1, 0.15), border=NA)
    
    lines(x = time_axis, y = y_fit, col="blue", lwd=2.5)
  }
  dev.off()
}


run_loo_years <- function(data_raw) {
  y_loo <- data_raw$Season; results_df <- data.frame()
  for(y in y_loo) {
    dat_loo <- data_raw; dat_loo[dat_loo$Season == y, -1] <- NA
    fit_loo <- try(run_turtle_marss(dat_loo, y_loo), silent=TRUE)
    if(!inherits(fit_loo, "try-error")){
      # LOO impact usually reported on first site trend
      u_val <- as.numeric(coef(fit_loo, type="matrix")$U[1]) 
      results_df <- rbind(results_df, data.frame(Removed_Year = y, Growth_Impact = round((exp(u_val)-1)*100, 2)))
    }
  }
  return(results_df)
}