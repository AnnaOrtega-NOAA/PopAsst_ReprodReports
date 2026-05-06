### Bayesian state-space model for estimating long-term trend in nest count time series data
### Revised May 6, 2026 by Anna Ortega (AO EDIT 2026)

require(ggplot2)

# Helper function for plot transparency
col2rgbA <- function(color, alpha = 0.3) {
  rgb_vals <- col2rgb(color)
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], 
      alpha = alpha * 255, maxColorValue = 255)
}

##==============================================================================
# 1. SCENARIO SETUP (Loggerhead vs Leatherback)
##==============================================================================
file.tag <- paste0(scenario, "_", Sys.Date())

if (scenario == "Cc_Yakushima_sUQ") {
  # Loggerhead logic: Select the 3 beaches defined in main script  
  thedata <- thedata.loggers[, c("Season", "Maehama", "Inakahama", "Yotsuse")] 
  pop.name <- c("Maehama", "Inakahama", "Yotsuse")
  pop.name.combo <- "N. Pacific Loggerheads - Yakushima"
} else {
  # Leatherback logic: Map JM and W from the merge  
  thedata <- thedata.leathers
  colnames(thedata) <- c("Season", "Jeen Yessa", "Jeen Syuab")
  pop.name <- c("Jeen Yessa", "Jeen Syaub")
  pop.name.combo <- "Western Pacific Leatherback Turtles"
}

data.cols <- 2:ncol(thedata)
data.rows <- 1:nrow(thedata)

# Save the exact data used for this run  
write.csv(thedata[data.rows, c(1, data.cols)], 
          file = paste0(table.path, file.tag, "_0_data_used.csv"), 
          quote = FALSE, row.names = FALSE)

##==============================================================================
# 2. DIAGNOSTIC DATA PLOTS
##==============================================================================
for (i in 1:length(pop.name)) {
  # Raw Nest Counts 
  png(filename = paste0(fig.path, file.tag, "_0.", i, "_data_rawNests.png"), width = 650, height = 575)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(thedata[,1], clutch.freq * thedata[, i+1], pch = 19, type = "b",
       xlab = "Season", ylab = "Nest Counts", main = pop.name[i], 
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
  
  # Log-transformed Annual Females 
  png(filename = paste0(fig.path, file.tag, "_1.", i, "_data_lnFemales.png"), width = 650, height = 575)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(thedata[,1], log(thedata[, i+1]), pch = 19, type = "b",
       xlab = "Season", ylab = "Ln(Annual Females)", main = pop.name[i],
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
}

##==============================================================================
# 3. TRANSFORM DATA & JAGS EXECUTION
##==============================================================================
log.data <- log(thedata[data.rows, data.cols])
data.mat <- t(log.data) 
n.yrs <- ncol(data.mat)
n.timeseries <- nrow(data.mat)

Y <- data.mat
Z <- matrix(1, nrow = n.timeseries, ncol = 1) 

jags.data <- list(
  Y = Y, n.yrs = n.yrs, n.timeseries = n.timeseries, Z = Z,
  a_mean = 0, a_sd = 4, u_mean = 0, u_sd = 0.5,
  q_alpha = 0.01, q_beta = 0.01, r_alpha = 0.01, r_beta = 0.01,
  x0_mean = mean(data.mat, na.rm = TRUE), x0_sd = 10
)

# Run JAGS model using singleUQ.txt [cite: 5, 8]
jags.model <- jags(
  jags.data, inits = NULL, 
  parameters.to.save = c("A", "U", "Q", "R", "X0", "X"),
  model.file = paste0(code.path, "singleUQ.txt"), 
  n.chains = 2, n.burnin = 50000, n.thin = 50, n.iter = 100000, parallel = TRUE
)

##==============================================================================
# 4. POSTERIOR EXTRACTION & COMBINED ABUNDANCE
##==============================================================================
nsim <- 10000
X.thin <- round(seq(1, length(jags.model$sims.list$U), length.out = nsim))

calc_total <- function(yr_idx) {
  rowSums(sapply(1:n.timeseries, function(j) {
    exp(jags.model$sims.list$X[X.thin, yr_idx] + jags.model$sims.list$A[X.thin, j])
  }))
}

X.fym0 <- calc_total(n.yrs)
X.fym1 <- calc_total(n.yrs - 1)
X.fym2 <- calc_total(n.yrs - 2)
X.fym3 <- calc_total(n.yrs - 3)

posts.out <- cbind(
  "U" = jags.model$sims.list$U[X.thin], 
  "Q" = jags.model$sims.list$Q[X.thin], 
  "N_fym0" = X.fym0, "N_fym1" = X.fym1, "N_fym2" = X.fym2, "N_fym3" = X.fym3
)

write.csv(posts.out, paste0(table.path, "Posteriors_", file.tag, ".csv"), row.names = FALSE)

##==============================================================================
# 5. INTEGRATED TREND PLOT (Density Polygons)
##==============================================================================
raw_obs_sum <- log(rowSums(thedata[data.rows, 2:ncol(thedata)], na.rm = TRUE))
raw_obs_sum[!is.finite(raw_obs_sum)] <- NA 

X.total_all_yrs <- sapply(1:n.yrs, calc_total)

# Calculate X0 separately to extend plot to T-1
X0.total <- rowSums(sapply(1:n.timeseries, function(j) {
  exp(jags.model$sims.list$X0[X.thin] + jags.model$sims.list$A[X.thin, j])
}))

# Combine X0 and X for the full credible interval band
X.full_matrix <- cbind(X0.total, X.total_all_yrs)
X.q <- apply(log(X.full_matrix), 2, quantile, probs = c(0.025, 0.5, 0.975))

# Define yrange for plot clipping
yrange <- range(c(raw_obs_sum, X.q), na.rm = TRUE)

# Densities for Ghosts
den.X0 <- density(log(X0.total), adj = 2)
den.X0$y2 <- (den.X0$y / max(den.X0$y))
q.X0 <- quantile(log(X0.total), probs = c(0.025, 0.975))
xid <- sapply(q.X0, function(x) which.min(abs(x - den.X0$x)))

den.N0 <- density(log(X.total_all_yrs[, n.yrs]), adj = 2)
den.N0$y2 <- (den.N0$y / max(den.N0$y))
q.N0 <- quantile(log(X.total_all_yrs[, n.yrs]), probs = c(0.025, 0.975))
xid.n0 <- sapply(q.N0, function(x) which.min(abs(x - den.N0$x)))

png(filename = paste0(fig.path, file.tag, "_model_fit_U.png"), width = 7, height = 4.62, units = "in", res = 300)
layout(matrix(1:2, 1, 2), width = c(1, 0.23))
par(mar = c(4, 4, 1, 1))

# Initial Plot
plot(thedata$Season, raw_obs_sum, pch = 16, type = "n", las = 1, 
     ylab = "log(Annual Nesters)", xlab = "Season", ylim = range(pretty(yrange)), 
     xlim = c(min(thedata$Season) - 1.1, max(thedata$Season) + 0.25), xaxs = 'i')

# Extended CI Band (includes Year 0)
extended_years <- c(min(thedata$Season) - 1, thedata$Season)
polygon(c(extended_years, rev(extended_years)), c(X.q[1, ], rev(X.q[3, ])), col = 'grey85', border = FALSE)

# Median Trend Line
lines(extended_years, X.q[2, ], lwd = 3, col = 'gray50')

# Density Polygons
polygon(x = c(rep((thedata[1,1]-1), length(xid[1]:xid[2])), rev(den.X0$y2[xid[1]:xid[2]]+(thedata[1,1]-1))), 
        y = c(den.X0$x[xid[1]:xid[2]], rev(den.X0$x[xid[1]:xid[2]])), 
        lwd = 2, border = FALSE, col = col2rgbA("dodgerblue3", 0.3))

polygon(x = c(rep(thedata[nrow(thedata),1], length(xid.n0[1]:xid.n0[2])), rev(-den.N0$y2[xid.n0[1]:xid.n0[2]]+thedata[nrow(thedata),1])), 
        y = c(den.N0$x[xid.n0[1]:xid.n0[2]], rev(den.N0$x[xid.n0[1]:xid.n0[2]])), 
        lwd = 2, border = FALSE, col = col2rgbA("darkorchid3", 0.3))

# Highlighted Median Points
points(extended_years, X.q[2,], pch = 16, 
       col = c('dodgerblue3', rep('red', length(extended_years) - 2), "darkorchid3"))

# Observation Points
points(thedata$Season, raw_obs_sum, pch = 16, col = "black")

# Restoration of Legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", 
       legend = c(expression(sum(N[list(obs,j)], j, "")), 
                  expression(sum(T[j] + a[j], j, "")), 
                  "Median r", "95% r", 
                  expression(paste(T[0], " (95%CI)")), 
                  expression(paste(N[final], " (95% CI)"))), 
       pch = c(16, 16, NA, 15, NA, NA), 
       lwd = c(NA, NA, 3, NA, 2, 2), 
       pt.cex = c(1, 1, NA, 3, NA, NA), 
       col = c("black", "red", "gray50", "gray85", "dodgerblue3", "darkorchid3"), 
       bty = "n", y.intersp = 0.9, xpd = NA)
dev.off()

##==============================================================================
# 6. INDIVIDUAL BEACH MODEL FIT PLOTS
##==============================================================================
if(n.timeseries > 1){
  mygray <- rgb(red=190, green=190, blue=190, alpha=200, maxColorValue=255)
  for (i in 1:n.timeseries){
    png(filename = paste0(fig.path, file.tag, "_6.", i, "_model_fit_", pop.name[i], ".png"), width = 650, height = 575)
    beach_A <- jags.model$sims.list$A[X.thin, i]
    beach_q <- apply(jags.model$sims.list$X[X.thin, ] + beach_A, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    plot(thedata$Season, log.data[, i], pch = 19, ylim = range(c(log.data[, i], beach_q), na.rm = TRUE),
         xlab = "Year", ylab = "Ln(Annual Females)", main = pop.name[i])
    polygon(c(thedata$Season, rev(thedata$Season)), c(beach_q[1, ], rev(beach_q[3, ])), col = mygray, border = FALSE)
    lines(thedata$Season, beach_q[2, ], col = "blue", lwd = 2)
    points(thedata$Season, log.data[, i], pch = 16)
    dev.off()
  }
}

##==============================================================================
# 7. 100-YEAR PROJECTIONS
##==============================================================================
yrf <- 100 
N_proj <- matrix(NA, nrow = nsim, ncol = yrf + 1)
N_proj[, 1] <- posts.out[, "N_fym0"]

for (t in 2:(yrf + 1)) {
  N_proj[, t] <- N_proj[, t-1] * exp(posts.out[, "U"] + rnorm(nsim, 0, sqrt(posts.out[, "Q"])))
}

N_med_proj <- apply(N_proj, 2, median)
proj_years <- seq(max(thedata$Season), max(thedata$Season) + yrf)

png(filename = paste0(fig.path, file.tag, "_7_Projections_100yr.png"), width = 7, height = 5, res = 300, units = "in")
plot(proj_years, N_med_proj, type = "l", lwd = 2, col = "blue", 
     ylim = c(0, quantile(N_proj[, yrf+1], 0.95)), 
     ylab = "Projected Total Females", xlab = "Year", main = "100-Year Abundance Projection")
polygon(c(proj_years, rev(proj_years)), 
        c(apply(N_proj, 2, quantile, 0.025), rev(apply(N_proj, 2, quantile, 0.975))), 
        col = rgb(0, 0, 1, 0.2), border = NA)
abline(h = N_med_proj[1], lty = 2, col = "red")
dev.off()

saveRDS(jags.model, paste0(trend.path, file.tag, ".rds"))