#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   PIFSC Hawaiian DSLL Sea Turtle Trend Update 2026
#   Focus: Imputation and Trend Analysis 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(mvtnorm); library(truncnorm); library(doParallel); library(foreach)
library(abind); library(jagsUI); library(coda); library(dplyr)
library(reshape2); library(tidyr); library(loo); library(cmdstanr); library(magrittr)

# --- REVISED DIRECTORY SETUP (AO EDIT 2026) ---
main.folder <- "C:/Users/anna.ortega/Desktop/PopAsst_RR/"

# Inputs 
code.path   <- paste0(main.folder, "R/")
data.path   <- paste0(main.folder, "data/")

# Outputs (Organized by stage) 
out.base    <- paste0(main.folder, "output/")
imput.path  <- paste0(out.base, "imputation/")
trend.path  <- paste0(out.base, "trend/")
fig.path    <- paste0(out.base, "figures/")
table.path  <- paste0(out.base, "tables/")
take.path   <- paste0(out.base, "take/") 

# Create missing folders 
dirs <- c(imput.path, trend.path, fig.path, table.path, take.path)
lapply(dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# --- LOAD HELPERS ---
source(paste0(code.path, "take_helper_Fn.R"))
redo <- TRUE

###-----------------------------------------------------
#   1. IMPUTATION (Leatherbacks - DC)
###-----------------------------------------------------
check.prerun <- file.exists(paste0(imput.path, "N_imput.Rdata"))
if(!check.prerun | redo){
  
  period.JM <- 12
  period.W <- 6
  year.begin.JM <- 2001
  year.end <- 2025                       
  all.years <- year.begin.JM:year.end
  n.years <- length(all.years)           
  
  data.jags.JM <- data.extract(location = "JM", year.begin = year.begin.JM, year.end = year.end, file.path = data.path)
  data.jags.W  <- data.extract(location = "W", year.begin = 2006, year.end = year.end, file.path = data.path)
  
  # Matrix Padding
  y.W <- rbind(array(data = NA, dim = c(nrow(data.jags.JM$jags.data2$y) - nrow(data.jags.W$jags.data2$y), ncol(data.jags.JM$jags.data2$y))), data.jags.W$jags.data2$y)
  y <- cbind(as.vector(t(data.jags.JM$jags.data2$y)), as.vector(t(y.W)))
  
  y.raw <- as.data.frame(apply(y,2,exp))
  colnames(y.raw) <- c("JM","W")
  y.raw$year.frac <- rep(all.years, each=12) + seq(1,12)/12 
  save(y.raw, file=paste0(data.path,"DC_raw_nests_month.Rdata"))
  
  jags.data <- list(y = y, m = rep(1:12, times = n.years), n.steps = nrow(y), n.months = 12, pi = pi, period = c(period.JM, period.W), n.timeseries = 2, n.years = n.years)
  jags.params <- c("c", "beta.cos", "beta.sin", 'sigma.X', "sigma.y", "N", "y", "X", "deviance")
  
  jm <- jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file = paste0(code.path,'model_norm_norm_Four_imputation.txt'), 
             n.chains = 5, n.burnin = 50000, n.thin = 5, n.iter = 100000, parallel=T)
  
  Ns.stats.JM <- data.frame(time = all.years, low = as.vector(t(jm$q2.5$N[,1])), median = as.vector(t(jm$q50$N[,1])), high = as.vector(t(jm$q97.5$N[,1])), location = "Jeen Yessa")
  Ns.stats.W  <- data.frame(time = all.years, low = as.vector(t(jm$q2.5$N[,2])), median = as.vector(t(jm$q50$N[,2])), high = as.vector(t(jm$q97.5$N[,2])), location = "Jeen Syuab")
  save(Ns.stats.JM, Ns.stats.W, file=paste0(imput.path, "N_imput.Rdata"))
} else {
  load(paste0(imput.path, "N_imput.Rdata"), verbose=T)
}

###-----------------------------------------------------
#   2. HISTORICAL ANE 
###-----------------------------------------------------
# Setting up a dummy table in case ANE calculation fails or data are missing
ANE.tab.pm.long <- data.frame(y = 2001:2025, sr = 0)

try({
  if(file.exists(paste0(data.path,"DSLL_interactions.csv"))){
    td.dat.raw <- read.csv(paste0(data.path,"DSLL_interactions.csv"), stringsAsFactors=FALSE)
    colnames(td.dat.raw)[1:6] <- c("Year", "Spp","Status",'Ryder', "M.low", "M.high")
    
    DC.VBGF <- list(Linf = 142.7, K = 0.2262, tknot=-0.17, Amat = 16.1, Pj = 0.81, Pa = 0.893, PF = 0.73)
    DC.VBGF$max_age <- vbgf_bc(prop=0.99, VBGF=DC.VBGF)
    
    ANE.full <- do.call('rbind', mapply(function(Year, Age, M.mu, max_year) with(DC.VBGF, {
      y <- seq(Year, max_year)
      safe_age_end <- max(Age, DC.VBGF$max_age)
      if(Age >= DC.VBGF$max_age){
        v <- rep(DC.VBGF$Linf * 0.99, length(y))
      } else {
        age_seq <- seq(Age, safe_age_end, by = 1)
        y2mat <- Year + length(age_seq) - 1
        if(y2mat < max_year){
          v <- c(DC.VBGF$Linf * (1 - exp(-DC.VBGF$K * (age_seq - DC.VBGF$tknot))), rep(DC.VBGF$Linf * 0.99, max_year - y2mat))
        } else {
          v <- DC.VBGF$Linf * (1 - exp(-DC.VBGF$K * (seq(Age, length.out = length(y), by = 1) - DC.VBGF$tknot)))
        }
      }
      p <- maturity(v, casewise=TRUE, p=0.99, VBGF=DC.VBGF)
      s <- cumprod((1-p) * 0.81 + p * 0.893)
      sr <- s * p * 0.73 * M.mu * (1/3)
      return(as.data.frame(cbind(y,v,p,s,sr)))
    }), Year=td.dat.raw$Year, Age=15, M.mu=0.5, max_year=2025, SIMPLIFY=FALSE))
    
    if(!is.null(ANE.full) && nrow(ANE.full) > 0) {
      ANE.tab.pm.long <- aggregate(sr~y, data=ANE.full, FUN=sum)
    }
  }
}, silent = FALSE)

###-----------------------------------------------------
#   3. TREND ANALYSIS (Consolidated Loop)
###-----------------------------------------------------
species_to_run <- c("Leatherback", "Loggerhead")

for(spp in species_to_run){
  
  if(spp == "Leatherback"){
    # JM and W data preparation
    JM.dat <- Ns.stats.JM %>% mutate(across(2:4, exp), Females_median = median/5.5, Season = time)
    W.dat  <- Ns.stats.W %>% filter(time %in% 2006:2025) %>% mutate(across(2:4, exp), Females_median = median/5.5, Season = time)
    
    # Combined dataframe for projection script
    thedata.leathers <- merge(JM.dat[,c("Season","Females_median")], W.dat[,c("Season","Females_median")], by="Season", all=T)
    
    scenario <- "Dc_JM&W_MEDIAN_sUQ"
    remig <- 3.06; remigLB <- 3.06
    clutch.freq <- 5.5; clutch.freqLB <- 5.5
    
  } else if(spp == "Loggerhead"){
    # Yakushima data preparation
    cc_raw <- read.csv(paste0(data.path, "Yakushima_data_for_BiOp_2025.csv"))
    
    # Column mapping fixed to match CC beach structure
    thedata.loggers <- data.frame(
      Season    = cc_raw$Year, 
      Maehama   = cc_raw$Maehama.Beach / 4.6, 
      Inakahama = cc_raw$Inakahama.Beach / 4.6,
      Yotsuse   = cc_raw$Yakushima.Yotsusehama / 4.6
    )
    
    scenario <- "Cc_Yakushima_sUQ" 
    remig <- 3.3; remigLH <- 3.3
    clutch.freq <- 4.6; clutch.freqLH <- 4.6
  }
  
  save.dir <- paste0(trend.path, scenario)
  if(!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
  
  cat("\n--- STARTING TREND ANALYSIS FOR:", spp, "(", scenario, ") ---\n")
  source(paste0(code.path, "singleUQ_indeptUQs_PROJECTIONS_20260506.R"))
  
  # Summary calculation block
  n_2025_posts <- posts.out[, "N_fym0"]
  n_curr_posts <- remig * (posts.out[, "N_fym0"] + posts.out[, "N_fym1"] + 
                             posts.out[, "N_fym2"] + posts.out[, "N_fym3"]) / 4
  
  get_summary <- function(x) quantile(x, probs = c(0.5, 0.025, 0.975))
  u_stats      <- get_summary(jags.model$sims.list$U)
  q_stats      <- get_summary(jags.model$sims.list$Q)
  n_2025_stats <- get_summary(n_2025_posts)
  n_curr_stats <- get_summary(n_curr_posts)
  
  master_summary <- data.frame(
    Metric = c("Annual Growth Rate (%)", "Process Variance (Q)", "2025 Nesters", "Total Abundance"),
    Median = c(round((exp(u_stats[1]) - 1) * 100, 2), round(q_stats[1], 4), round(n_2025_stats[1], 0), round(n_curr_stats[1], 0)),
    LCI_95 = c(round((exp(u_stats[2]) - 1) * 100, 2), round(q_stats[2], 4), round(n_2025_stats[2], 0), round(n_curr_stats[2], 0)),
    UCI_95 = c(round((exp(u_stats[3]) - 1) * 100, 2), round(q_stats[3], 4), round(n_2025_stats[3], 0), round(n_curr_stats[3], 0)),
    Unit = c("% growth/yr", "Variance", "Annual Females", "Total Reproductive Females")
  )
  
  print(master_summary, row.names = FALSE)
  write.csv(master_summary, file = paste0(table.path, "Final_Results_Table_", scenario, ".csv"), row.names = FALSE)
}