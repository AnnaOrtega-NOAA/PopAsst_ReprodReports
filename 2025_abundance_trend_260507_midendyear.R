#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   PIFSC Hawaiian DSLL Sea Turtle Trend Update 2026
#   Revised: Seasonal Leatherback Grouping & MARSS
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(mvtnorm); library(truncnorm); library(doParallel); library(foreach)
library(abind); library(jagsUI); library(coda); library(dplyr)
library(reshape2); library(tidyr); library(loo); library(magrittr); library(MARSS)

# --- DIRECTORY SETUP ---
main.folder <- "C:/Users/anna.ortega/Desktop/PopAsst_RR/"
code.path    <- paste0(main.folder, "R/")
data.path    <- paste0(main.folder, "data/")
out.base     <- paste0(main.folder, "output/")
imput.path   <- paste0(out.base, "imputation/")
trend.path   <- paste0(out.base, "trend/")
fig.path     <- paste0(out.base, "figures/")
table.path   <- paste0(out.base, "tables/")

dirs <- c(imput.path, trend.path, fig.path, table.path)
lapply(dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

source(paste0(code.path, "take_helper_Fn.R"))
source(paste0(code.path, "marss_helper_Fn_midend.R")) 
redo <- TRUE

###-----------------------------------------------------
#   1. IMPUTATION & GROUPED SPLITTING (Leatherbacks)
###-----------------------------------------------------
if(!file.exists(paste0(imput.path, "N_imput_grouped.Rdata")) | redo){
  all.years <- 2001:2025
  n.years <- length(all.years)
  
  data.jags.JM <- data.extract(location = "JM", year.begin = 2001, year.end = 2025, file.path = data.path)
  data.jags.W  <- data.extract(location = "W", year.begin = 2006, year.end = 2025, file.path = data.path)
  
  y_jm_vec <- as.vector(t(data.jags.JM$jags.data2$y))
  y_w_vec  <- as.vector(t(data.jags.W$jags.data2$y))
  padding  <- rep(NA, (n.years * 12) - length(y_w_vec))
  y_w_aligned <- c(padding, y_w_vec)
  y <- cbind(y_jm_vec, y_w_aligned)
  
  jags.data <- list(y = y, m = rep(1:12, times = n.years), n.steps = nrow(y), 
                    n.months = 12, pi = pi, period = c(12, 6), n.timeseries = 2, n.years = n.years)
  jm <- jags(jags.data, inits = NULL, parameters.to.save= c("y"), 
             model.file = paste0(code.path,'model_norm_norm_Four_imputation.txt'), 
             n.chains = 5, n.burnin = 50000, n.thin = 5, n.iter = 100000, parallel=T)
  
  y_monthly <- jm$mean$y 
  df_monthly <- data.frame(Year = rep(all.years, each=12), Month = rep(1:12, times=n.years),
                           JM = exp(y_monthly[,1]), W = exp(y_monthly[,2]))
  
  site_totals <- df_monthly %>% group_by(Year) %>% summarise(JM = sum(JM), W = sum(W), .groups = "drop")
  cohort_totals <- df_monthly %>%
    mutate(Cohort = ifelse(Month %in% 4:9, "MidYear", "EndYear")) %>%
    group_by(Year, Cohort) %>%
    summarise(Total_Nests = sum(JM + W), .groups = "drop") %>%
    pivot_wider(names_from = Cohort, values_from = Total_Nests)
  
  thedata.leathers.grouped <- merge(site_totals, cohort_totals, by="Year") %>% rename(Season = Year)
  save(thedata.leathers.grouped, file=paste0(imput.path, "N_imput_grouped.Rdata"))
} else { load(paste0(imput.path, "N_imput_grouped.Rdata")) }

###-----------------------------------------------------
#   2. TREND ANALYSIS (MARSS)
###-----------------------------------------------------
species_to_run <- c("Leatherback", "Loggerhead")

for(spp in species_to_run){
  if(spp == "Leatherback"){
    thedata <- thedata.leathers.grouped 
    remig <- 3.06
  } else {
    # CLEAN DATA: Removes empty rows at the bottom of the CSV
    cc_raw <- read.csv(paste0(data.path, "Yakushima_data_for_BiOp_2025.csv")) %>% 
      filter(!is.na(Year)) 
    
    thedata <- data.frame(Season = cc_raw$Year, 
                          Maehama = cc_raw$Maehama.Beach, 
                          Inakahama = cc_raw$Inakahama.Beach, 
                          Yotsuse = cc_raw$Yakushima.Yotsusehama)
    remig <- 3.3
  }
  
  window_size <- round(remig * 3)
  # Calculate recent window based on the actual maximum year in the CURRENT data
  recent_years <- (max(thedata$Season) - window_size):max(thedata$Season)
  
  cat("\n--- STARTING MARSS ANALYSIS FOR:", spp, "---\n")
  fit_full   <- run_turtle_marss(thedata, unique(thedata$Season))
  fit_recent <- run_turtle_marss(thedata, recent_years)
  
  master_report <- rbind(get_marss_stats(fit_full, spp, remig, "Historical"),
                         get_marss_stats(fit_recent, spp, remig, "Recent"))
  
  write.csv(master_report, file = paste0(table.path, "Master_Report_", spp, ".csv"), row.names = F)
  print(master_report)
  
  # PLOTTING
  plot_marss_trend(fit_full, spp, fig.path)
  plot_marss_trend(fit_recent, spp, fig.path) # Added so you get both plots!
  
  cat("Calculating Influence (LOO) for", spp, "...\n")
  infl <- run_loo_years(thedata)
  write.csv(infl, file = paste0(table.path, "Influence_LOO_", spp, ".csv"), row.names = F)
}
