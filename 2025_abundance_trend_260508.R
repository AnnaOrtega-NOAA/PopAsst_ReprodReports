#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   PIFSC Sea Turtle DPS Assessment: Master Trend Script (Integrated)
#   Compiler: Anna Ortega (AO EDIT May 2026)
#   Original Framework: Boyd et al. (2016); Revised: Zach Siders (2022)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(jagsUI) # Load jagsUI for Bayesian modeling
library(dplyr) # Load dplyr for data manipulation
library(tidyr) # Load tidyr for reshaping data (long to wide)
library(MARSS) # Load MARSS for Frequentist state-space trends
library(ggplot2) # Load ggplot2 for visualization
library(magrittr) # Load magrittr for the pipe operator (%>%)

# --- 1. DIRECTORY SETUP ---
main.folder <- "C:/Users/anna.ortega/Desktop/PopAsst_RR/" # Set root folder path
code.path    <- paste0(main.folder, "R/") # Path for R scripts and .txt models
data.path    <- paste0(main.folder, "data/") # Path for .csv data files, hidden on Github for now
out.base     <- paste0(main.folder, "output/") # Base output directory, hidden on GitHub for now
fig.path     <- paste0(out.base, "figures/") # Folder for plot images
table.path   <- paste0(out.base, "tables/") # Folder for result tables

if(!dir.exists(fig.path)) dir.create(fig.path, recursive = TRUE) # AO EDIT May 2026: Recursive check ensures folders exist before plotting
if(!dir.exists(table.path)) dir.create(table.path, recursive = TRUE) # AO EDIT May 2026: Ensure table directory exists

source(paste0(code.path, "take_helper_Fn.R")) # Load data extraction helper (links to Oct 2025 data)
source(paste0(code.path, "marss_helper_Fn_midend.R")) # Load Frequentist MARSS execution helper

col2rgbA <- function(color, alpha = 0.3) { # AO EDIT May 2026: Transparency helper for Ghost plots
  rgb_vals <- col2rgb(color) # Extract RGB components from color string
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha * 255, maxColorValue = 255) # Return hex code with alpha
}

# --- 2. SPECIES LOOP ---
species_to_run <- c("Leatherback", "Loggerhead") # AO EDIT May 2026: Multi-species loop replaces manual scenario switching

for(spp in species_to_run){ # Begin loop for each species
  cat("\n--- PROCESSING:", spp, "---\n") # Print species progress to console
  
  ### 2.2 DATA PREPARATION
  if(spp == "Leatherback"){ # Leatherback-specific data preparation
    data.jags.JM <- data.extract("JM", 2001, 2025, file.path = data.path) # Extract Jamursba-Medi data (2001-2025)
    data.jags.W  <- data.extract("W", 2006, 2025, file.path = data.path) # Extract Wermon data (2006-2025)
    
    thedata <- data.frame(Season = 2001:2025, # AO EDIT May 2026: Automate asynchronous alignment of site data
                          JeenYessa = rowSums(exp(data.jags.JM$jags.data2$y), na.rm=T), # Sum JM monthly counts to annual, renamed to Jeen Yessa
                          JeenSyuab = c(rep(NA, 5), rowSums(exp(data.jags.W$jags.data2$y), na.rm=T))) # Renamed W to JS; Pad Jeenn Syuab with NAs for 2001-2005
    remig <- 3.06; clutch <- 5.5 # Biological constants for Leatherbacks
    raw_title <- "Leatherback Imputed Nest Counts" # Define title for raw plot
    
    df_monthly <- data.frame(Year = rep(2001:2025, each=12), Month = rep(1:12, 25), # AO EDIT May 2026: Monthly df for cohort splitting
                             JM = as.vector(t(exp(data.jags.JM$jags.data2$y))), # Flatten JM nesting matrix to vector
                             W = c(rep(NA, 60), as.vector(t(exp(data.jags.W$jags.data2$y))))) # Flatten W nesting matrix with 5-year pad
    
    thedata_midend <- df_monthly %>% # AO EDIT May 2026: Categorize into MidYear and EndYear nesting cohorts
      mutate(Cohort = ifelse(Month %in% 4:9, "MidYear", "EndYear")) %>% # Cohort split: MidYear (Apr-Sept) vs EndYear (Oct-Mar)
      group_by(Year, Cohort) %>% # Group monthly data by Year and Cohort
      summarise(Total = sum(JM + W, na.rm=T), .groups="drop") %>% # Sum all nests per cohort
      pivot_wider(names_from = Cohort, values_from = Total) %>% # Pivot to wide format (columns = cohorts)
      rename(Season = Year) # Rename Year to Season for model compatibility
    
  } else { # Loggerhead-specific data preparation
    cc_raw <- read.csv(paste0(data.path, "Yakushima_data_for_BiOp_2025.csv")) %>% filter(!is.na(Year)) # Read Loggerhead Yakushima data
    thedata <- data.frame(Season = cc_raw$Year, Maehama = cc_raw$Maehama.Beach, # Build Loggerhead beach dataframe
                          Inakahama = cc_raw$Inakahama.Beach, Yotsuse = cc_raw$Yakushima.Yotsusehama) # Map 3 unique beach streams
    remig <- 3.3; clutch <- 4.6 # Biological constants for Loggerheads
    raw_title <- "Loggerhead Raw Nest Counts" # Define title for raw plot
  }
  
  # --- 3. DIAGNOSTIC PLOTS ---
  png(paste0(fig.path, spp, "_1_Raw_Imputed_Nests.png"), width=2000, height=1200, res=300) # Open PNG device for raw plot
  matplot(thedata$Season, thedata[,-1], type="b", pch=16, main=raw_title, xlab="Year", ylab="Nests") # Plot raw nest counts
  legend("topright", legend=colnames(thedata)[-1], col=1:(ncol(thedata)-1), lty=1, pch=16, bty="n") # Add beach legend
  dev.off() # Close PNG device
  
  log_data_mat <- as.matrix(thedata[,-1] / clutch) # AO EDIT May 2026: Convert nest counts to female counts
  log_data_mat <- log(log_data_mat) # Calculate natural log of female counts
  log_data_mat[is.infinite(log_data_mat) | is.nan(log_data_mat)] <- NA # AO EDIT May 2026: Handle log-zeros to prevent crashes
  
  png(paste0(fig.path, spp, "_2_Log_Annual_Females.png"), width=2000, height=1200, res=300) # Open PNG device for log plot
  matplot(thedata$Season, log_data_mat, type="l", lwd=2, main=paste(spp, "Log Annual Females"), ylab="Log(Females)") # Plot log females
  dev.off() # Close PNG device
  
  # --- 4. MARSS TRENDS (FREQUENTIST) ---
  fit_full <- run_turtle_marss(thedata, unique(thedata$Season)) # Fit MARSS model to site data
  plot_marss_trend(fit_full, spp, fig.path) # Generate site-level MARSS trend plots
  
  if(spp == "Leatherback"){ # AO EDIT May 2026: Execute cohort-specific MARSS for Leatherbacks
    fit_me <- run_turtle_marss(thedata_midend, unique(thedata_midend$Season)) # Fit MARSS to Mid/End Year cohort data
    plot_marss_trend(fit_me, "Leatherback_MidEnd", fig.path) # Save cohort MARSS trend plots
  }
  
  # --- 5. BAYESIAN ANALYSIS (JAGS) ---
  jags.data <- list(Y = t(log_data_mat), n.yrs = nrow(thedata), n.timeseries = ncol(log_data_mat), # Prepare JAGS data list
                    Z = matrix(1, nrow = ncol(log_data_mat), ncol = 1), a_mean = 0, a_sd = 4, # Define Z matrix and priors
                    u_mean = 0, u_sd = 0.5, q_alpha = 0.01, q_beta = 0.01, r_alpha = 0.01, r_beta = 0.01, # Set hyper-priors
                    x0_mean = mean(log_data_mat, na.rm=T), x0_sd = 10) # Set initial state priors
  
  jags.model <- jags(jags.data, inits = NULL, parameters.to.save = c("A", "U", "Q", "X0", "X"), # Initialize JAGS model
                     model.file = paste0(code.path, "singleUQ.txt"), # Reference standard singleUQ model
                     n.chains = 3, n.burnin = 10000, n.thin = 10, n.iter = 50000, parallel = TRUE) # Run parallel MCMC
  
  for (i in 1:jags.data$n.timeseries){ # Generate Bayesian site fits
    png(paste0(fig.path, spp, "_3_Bayes_Fit_", colnames(thedata)[i+1], ".png"), width = 1500, height = 1100, res = 300) # Open beach PNG
    beach_posts <- jags.model$sims.list$X + jags.model$sims.list$A[,i] # Combine hidden state and site scaling
    beach_q <- apply(beach_posts, 2, quantile, probs = c(0.025, 0.5, 0.975)) # Calculate 95% Credible Intervals
    plot(thedata$Season, log_data_mat[, i], type="n", ylim=range(beach_q, na.rm=T), main=colnames(thedata)[i+1], ylab="Ln Females") # AO EDIT 2026: Axes first
    polygon(c(thedata$Season, rev(thedata$Season)), c(beach_q[1, ], rev(beach_q[3, ])), col = "grey90", border = FALSE) # Shaded CI behind points
    lines(thedata$Season, beach_q[2, ], col = "blue", lwd = 2) # Blue median line
    points(thedata$Season, log_data_mat[, i], pch = 19) # AO EDIT 2026: Overlay points on top of CI
    dev.off() # Close PNG
  }
  
  # --- 6. INTEGRATED DPS PLOTS (LOG & NATURAL) ---
  nsim <- length(jags.model$sims.list$U); X.thin <- 1:nsim # Set simulation constants
  calc_total <- function(yr_idx) { # Function to sum all site posteriors into population total
    rowSums(sapply(1:jags.data$n.timeseries, function(j) { # Map across all time series j
      exp(jags.model$sims.list$X[X.thin, yr_idx] + jags.model$sims.list$A[X.thin, j]) # Add offset A to state X and exponentiate
    })) # Return summed female count
  } # End calculation function
  X.total_all_yrs <- sapply(1:jags.data$n.yrs, calc_total) # Apply sum to all observed years
  X0_total <- rowSums(sapply(1:jags.data$n.timeseries, function(j) { # Calculate Year T-1 abundance total
    exp(jags.model$sims.list$X0[X.thin] + jags.model$sims.list$A[X.thin, j]) # Transform starting state X0
  })) # Return Year T-1 abundance
  
  ext_yrs <- c(min(thedata$Season)-1, thedata$Season) # Build x-axis from Year T-1
  
  # A. INTEGRATED LOG SCALE
  q_log <- apply(log(cbind(X0_total, X.total_all_yrs)), 2, quantile, probs = c(0.025, 0.5, 0.975)) # Calculate log quantiles
  obs_sum_log <- log(rowSums(thedata[,-1], na.rm = TRUE) / clutch) # Sum observed data to log scale
  png(paste0(fig.path, spp, "_4_DPS_Status_LOG.png"), width = 2100, height = 1386, res = 300) # Open Log PNG
  layout(matrix(1:2, 1, 2), width = c(1, 0.25)); par(mar = c(4, 4, 2, 1)) # Split plot for legend
  plot(thedata$Season, obs_sum_log, type = "n", main=paste(spp, "DPS Status (Log)"), ylab = "log(Annual Nesters)", xlab = "Season", ylim = range(q_log, na.rm=T)) # Axes
  polygon(c(ext_yrs, rev(ext_yrs)), c(q_log[1, ], rev(q_log[3, ])), col = 'grey85', border = FALSE) # Shaded band
  lines(ext_yrs, q_log[2, ], lwd = 3, col = 'gray50') # Median trend
  den.X0 <- density(log(X0_total), adj = 2); den.Nf <- density(log(X.total_all_yrs[, jags.data$n.yrs]), adj = 2) # Calculate densities
  polygon(x = c(rep(ext_yrs[1], 512), rev((den.X0$y/max(den.X0$y))+ext_yrs[1])), y = c(den.X0$x, rev(den.X0$x)), border = FALSE, col = col2rgbA("dodgerblue3", 0.4)) # Blue Ghost
  polygon(x = c(rep(max(thedata$Season), 512), rev(-(den.Nf$y/max(den.Nf$y))+max(thedata$Season))), y = c(den.Nf$x, rev(den.Nf$x)), border = FALSE, col = col2rgbA("darkorchid3", 0.4)) # Purple Ghost
  points(ext_yrs, q_log[2,], pch = 16, col = c('dodgerblue3', rep('red', length(ext_yrs)-2), "darkorchid3")); points(thedata$Season, obs_sum_log, pch = 16, col = "black") # Points
  dev.off() # Close device
  
  # B. INTEGRATED NATURAL SCALE
  q_nat <- apply(cbind(X0_total, X.total_all_yrs), 2, quantile, probs = c(0.025, 0.5, 0.975)) # AO EDIT 2026: Calculate natural scale quantiles
  obs_sum_nat <- rowSums(thedata[,-1], na.rm = TRUE) / clutch # Observed natural counts
  png(paste0(fig.path, spp, "_5_DPS_Status_NATURAL.png"), width = 2100, height = 1386, res = 300) # Open Natural PNG
  layout(matrix(1:2, 1, 2), width = c(1, 0.25)); par(mar = c(4, 4, 2, 1)) # Plot layout
  plot(thedata$Season, obs_sum_nat, type = "n", main=paste(spp, "DPS Status (Natural)"), ylab = "Annual Females", xlab = "Season", ylim = c(0, max(q_nat, na.rm=T)*1.1)) # Natural axes
  polygon(c(ext_yrs, rev(ext_yrs)), c(q_nat[1, ], rev(q_nat[3, ])), col = 'grey85', border = FALSE) # Natural shaded CI
  lines(ext_yrs, q_nat[2, ], lwd = 3, col = 'gray50') # Natural median trend
  den.X0n <- density(X0_total, adj = 2); den.Nfn <- density(X.total_all_yrs[, jags.data$n.yrs], adj = 2) # AO EDIT 2026: Natural scale densities
  polygon(x = c(rep(ext_yrs[1], 512), rev((den.X0n$y/max(den.X0n$y))+ext_yrs[1])), y = c(den.X0n$x, rev(den.X0n$x)), border = FALSE, col = col2rgbA("dodgerblue3", 0.4)) # Natural Blue Ghost
  polygon(x = c(rep(max(thedata$Season), 512), rev(-(den.Nfn$y/max(den.Nfn$y))+max(thedata$Season))), y = c(den.Nfn$x, rev(den.Nfn$x)), border = FALSE, col = col2rgbA("darkorchid3", 0.4)) # Natural Purple Ghost
  points(ext_yrs, q_nat[2,], pch = 16, col = c('dodgerblue3', rep('red', length(ext_yrs)-2), "darkorchid3")); points(thedata$Season, obs_sum_nat, pch = 16, col = "black") # Natural points
  dev.off() # Close device
  
  # --- 7. THRESHOLD TABLES (DPS RISK) ---
  yrf <- 100; N_proj <- matrix(NA, nrow = nsim, ncol = yrf + 1); N_proj[, 1] <- X.total_all_yrs[, jags.data$n.yrs] # AO EDIT 2026: Setup risk matrix
  for (t in 2:(yrf + 1)) { N_proj[, t] <- N_proj[, t-1] * exp(jags.model$sims.list$U + rnorm(nsim, 0, sqrt(jags.model$sims.list$Q))) } # Stochastic steps
  N_curr <- median(N_proj[, 1]); thresh_lvls <- c(0.5, 0.25, 0.125); thresh_names <- c("50% Decline", "25% Decline", "12.5% Decline") # Set risk levels
  results <- data.frame(); for(k in 1:3){ # AO EDIT 2026: Calculate risk crossings
    target <- N_curr * thresh_lvls[k] # Numerical abundance target
    hit_yrs <- apply(N_proj, 1, function(x) { y <- which(x <= target); if(length(y) > 0) return(max(thedata$Season) + min(y) - 1) else return(NA) }) # First year crossing
    results <- rbind(results, data.frame(Threshold = thresh_names[k], Median_Year = median(hit_yrs, na.rm=T), Prob_by_2125 = mean(!is.na(hit_yrs)) * 100)) # Crossing stats
  } # End calculation
  write.csv(results, paste0(table.path, spp, "_DPS_Thresholds.csv"), row.names=F) # Export final risk CSV
} # End Species Loop