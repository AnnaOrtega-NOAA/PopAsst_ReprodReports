#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   PIFSC Hawaiian DSLL Sea Turtle Trend Update 2026
#   Focus: Imputation (JAGS) and Trend Analysis (MARSS)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(mvtnorm); library(truncnorm); library(doParallel); library(foreach)
library(abind); library(jagsUI); library(coda); library(dplyr)
library(reshape2); library(tidyr); library(loo); library(magrittr); library(MARSS)

# --- DIRECTORY SETUP ---
main.folder <- "C:/Users/anna.ortega/Desktop/PopAsst_RR/"
code.path   <- paste0(main.folder, "R/")
data.path   <- paste0(main.folder, "data/")
out.base    <- paste0(main.folder, "output/")
imput.path  <- paste0(out.base, "imputation/")
trend.path  <- paste0(out.base, "trend/")
fig.path    <- paste0(out.base, "figures/")
table.path  <- paste0(out.base, "tables/")

# Create missing folders
dirs <- c(imput.path, trend.path, fig.path, table.path)
lapply(dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# --- LOAD HELPERS ---
source(paste0(code.path, "take_helper_Fn.R"))
source(paste0(code.path, "marss_helper_Fn.R")) # NEW HELPER
redo <- TRUE

###-----------------------------------------------------
#   1. IMPUTATION (Leatherbacks - DC)
###-----------------------------------------------------
check.prerun <- file.exists(paste0(imput.path, "N_imput.Rdata"))
if(!check.prerun | redo){
  
  all.years <- 2001:2025
  n.years <- length(all.years)
  
  data.jags.JM <- data.extract(location = "JM", year.begin = 2001, year.end = 2025, file.path = data.path)
  data.jags.W  <- data.extract(location = "W", year.begin = 2006, year.end = 2025, file.path = data.path)
  
  # FIX: Matrix Padding to avoid sub-multiple warnings
  y_jm_vec <- as.vector(t(data.jags.JM$jags.data2$y))
  y_w_vec  <- as.vector(t(data.jags.W$jags.data2$y))
  
  # Pad Winter data with NAs for years 2001-2005
  padding <- rep(NA, (n.years * 12) - length(y_w_vec))
  y_w_aligned <- c(padding, y_w_vec)
  
  y <- cbind(y_jm_vec, y_w_aligned)
  
  jags.data <- list(y = y, m = rep(1:12, times = n.years), n.steps = nrow(y), 
                    n.months = 12, pi = pi, period = c(12, 6), 
                    n.timeseries = 2, n.years = n.years)
  
  jags.params <- c("c", "beta.cos", "beta.sin", 'sigma.X', "sigma.y", "N", "y", "X", "deviance")
  
  jm <- jags(jags.data, inits = NULL, parameters.to.save= jags.params, 
             model.file = paste0(code.path,'model_norm_norm_Four_imputation.txt'), 
             n.chains = 5, n.burnin = 50000, n.thin = 5, n.iter = 100000, parallel=T)
  
  Ns.stats.JM <- data.frame(time = all.years, low = as.vector(t(jm$q2.5$N[,1])), 
                            median = as.vector(t(jm$q50$N[,1])), high = as.vector(t(jm$q97.5$N[,1])), location = "Jeen Yessa")
  Ns.stats.W  <- data.frame(time = all.years, low = as.vector(t(jm$q2.5$N[,2])), 
                            median = as.vector(t(jm$q50$N[,2])), high = as.vector(t(jm$q97.5$N[,2])), location = "Jeen Syuab")
  save(Ns.stats.JM, Ns.stats.W, file=paste0(imput.path, "N_imput.Rdata"))
} else {
  load(paste0(imput.path, "N_imput.Rdata"), verbose=T)
}

###-----------------------------------------------------
#   2. TREND ANALYSIS (MARSS VERSION)
###-----------------------------------------------------
species_to_run <- c("Leatherback", "Loggerhead")

for(spp in species_to_run){
  
  if(spp == "Leatherback"){
    # Prep from Imputation output
    thedata <- merge(Ns.stats.JM[,c("time","median")], Ns.stats.W[,c("time","median")], by="time", all=T) %>%
      mutate(across(2:3, exp)) # Back-transform from log-space for MARSS helper
    colnames(thedata) <- c("Season", "Jeen_Yessa", "Jeen_Syuab")
    remig <- 3.06
  } else {
    # Prep Loggerhead (Yakushima)
    cc_raw <- read.csv(paste0(data.path, "Yakushima_data_for_BiOp_2025.csv"))
    thedata <- data.frame(Season = cc_raw$Year, Maehama = cc_raw$Maehama.Beach, 
                          Inakahama = cc_raw$Inakahama.Beach, Yotsuse = cc_raw$Yakushima.Yotsusehama)
    remig <- 3.3
  }
  
  # Define Window (Last 3 remigration intervals)
  window_size <- round(remig * 3)
  recent_years <- (2025 - window_size):2025
  
  cat("\n--- STARTING COMPARATIVE MARSS ANALYSIS FOR:", spp, "---\n")
  
  # Run the two timeframes
  fit_full   <- run_turtle_marss(thedata, unique(thedata$Season))
  fit_recent <- run_turtle_marss(thedata, recent_years)
  
  # Extract statistics for both
  stats_full   <- get_marss_stats(fit_full, spp, remig, "Historical (2001-2025)")
  stats_recent <- get_marss_stats(fit_recent, spp, remig, paste0("Recent (Last ", window_size, "yrs)"))
  
  # Combine into one Master Table
  master_report <- rbind(stats_full, stats_recent)
  
  # Print to console for immediate review
  print(master_report)
  
  # Save the comparative report
  write.csv(master_report, file = paste0(table.path, "Master_Trend_Report_", spp, ".csv"), row.names = F)
  
  # 2. Plotting (We will plot the Full trend to see the long-term context)
  plot_marss_trend(fit_full, spp, fig.path)
  # 
  # # Run Full and Recent Models
  # fit_full   <- run_turtle_marss(thedata, unique(thedata$Season))
  # fit_recent <- run_turtle_marss(thedata, recent_years)
  # 
  # 
  # # 1. Statistics Summary Table
  # stats_tab <- get_marss_stats(fit_recent, spp, remig)
  # write.csv(stats_tab, file = paste0(table.path, "Summary_", spp, "_Recent.csv"), row.names = F)
  # print(stats_tab)
  
  # 2. Plotting Trend Lines
  plot_marss_trend(fit_recent, spp, fig.path)
  
  # 3. Influence Analysis (LOO)
  cat("Calculating Influence (LOO) for", spp, "...\n")
  infl <- run_loo_years(thedata)
  write.csv(infl, file = paste0(table.path, "Influence_LOO_", spp, ".csv"), row.names = F)
}