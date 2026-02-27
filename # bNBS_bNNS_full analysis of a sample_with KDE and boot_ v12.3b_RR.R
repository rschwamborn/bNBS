# The bNBS and bNNS size spectra -------
# Version 12.3b, with a simple function at the end of the script  -------------
# histogram, bootstrap and KDE (kernel density estimation)
# Example: analysis of a zooplankton sample
# Uses Robust Regression (the "rlm( )" function in the MASS package)

# NNSS: normalized numbers size spectrum (Abundance vs body mass, normalized, but not rescaled)
# bNNS (new!): backtransformed and normalized numbers  spectrum (Abundance vs body mass, normalized, AND rescaled)
# NBSS: normalized biomass size spectrum (Biomass vs body mass, normalized, but not rescaled)
# bNBS (new!): backtransformed  and normalized numbers  spectrum (Biomass vs body mass, normalized, AND rescaled)
# KDE:  kernel density estimation

# BWIA: bin-width-inflated Abundance (distorted by non-linear binning, not normalized)
# BWIB: bin-width-inflated Biomass (distorted by non-linear binning, not normalized)

# FUNCTION
# A  convenient function is  provided at the end of the script:
# function: "fun.bNNS_and_bNBS_plots_summary_table()" 
# gives a  summary results table and two plots (bNNS_and_bNBS)

# The function 
# runs a two-step bootstrap (for each bin and for the linear regression)
# Produces  two  plots (bNNS + bNBS)
# Uses  max to last non- zero bin selection
# Shows  95% CI envelope on regression envelope 
# Show KDE high-resolution density distribution line (for bNNs only)
# Computes bootstrap slope + intercept CIs
# Computes totals (total Abundance and total Biomass, for all data, and for the linear range only)
# Returns a  summary results table

# Examples with bootstrap and KDE (kernel density estimation)

# CONTENT: ######### ---------

# 0. Preparations ---------

# 1. Input data and vectors --------------
# 1 .0 Example of a tropical zooplankton sample ----------------
#  1.1 Low variability simulated data (perfect power law)
#  1.2 High variability simulated data 
#  1.3 Very high Variability simulated data 
# 1.4 Standard log2-binning vector (breaks for binning) ---------------

# 2. bNNS  --------------------

# 3. bNBS -----------
# 3.1 bNBS (simple, with common regression line, no bootstrapping) -------------
# 3.2 bNBS (bootstrapping the histogram, with 95% CI error bars) -------------
# 3.3 bNBS (bootstrapping the histogram, with with regression line and bootstrapped regression envelope) -------------

# 4. The FUNCTION ---------------

# 5. Test the FUNCTION ---------------


##########################

# 0. Preparations

# clean memory (caution...)-----------
# rm(list = ls())
# gc()

library(scales)
library(MASS)

##########################

# 1. Input data and vectors --------------

# 1.0 Example of a tropical zooplankton sample ----------------

#  import example dataset N = 7088 individuals in the sample  ----------

# Zooplankton caught with a 300 micron mesh net in the Tropical Atlantic off Fernando de Noronha Island, Brazil (ABRACOS Project)  

# Dataset: "no_Artef zz_Merged_ zz_Infl_ abr1_noronha_st01_300_d1_d1_1_dat1"
# N = 7088 individuals in the sample
# sampled volume: 100 cubic meters (hypothetical)
# dataset (csv file) can be obtained at:

url <- "https://raw.githubusercontent.com/rschwamborn/bNBS/refs/heads/main/carbon_no_Artef%20zz_Merged_%20zz_Infl_%20abr1_noronha_st01_300.csv"

# df <- read.csv(url)

carbon_g_ind_1_zoopldataset.table <- read.csv(url, sep = ";")
# View(carbon_g_ind_1_zoopldataset.table)

data.carbon_g_ind_1 <- carbon_g_ind_1_zoopldataset.table$g_carbon
data.carbon_g_ind_1 <- as.numeric(data.carbon_g_ind_1)

# Sampling unit (sampled volume, filtered volume)

# Filtered volume --------
S = 100 # 100 cubic metres = filtered volume (example)

# Total number of individuals in the sample --------------
N_total = length(data.carbon_g_ind_1)# N = 7088, indiv. per sample
N_total

# Abundance (ind.  m-3) -------------
A_total = length(data.carbon_g_ind_1) / S # A =  N/S , indiv. per cubic meter
A_total # 70.88 total Abundance (ind.  m-3)

# Biomass (g.  m-3) ------------
B_total = sum(data.carbon_g_ind_1) / S # B =  sum(M)/S, g. per cubic meter
B_total # 0.00342 total Biomass (g.  m-3)


# Simulated body mass data (power-law distributed) --------------

#  1.1 low Variability simulated data (perfect power law)

set.seed(123)
n <- 100000
alpha <- 2
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)
body_mass.pw <- x 
vec1_2 <- x

length(body_mass.pw) # 100,000 ind.
sum(body_mass.pw) # 90.679 g

# Total number of individuals in the sample --------------
N_total = length(body_mass.pw)# N = 100,000 indiv. per sample
N_total

# Abundance (ind.  m-3) -------------
A_total = length(body_mass.pw) / S # A =  N/S , indiv. per cubic meter
A_total # 1000 total Abundance (ind.  m-3)

# Biomass (g.  m-3) ------------
B_total = sum(body_mass.pw) / S # B =  sum(M)/S, g. per cubic meter
B_total # 0.90679 total Biomass (g.  m-3)


#  1.2 High variability Data ----------------

#  High variability Data -------------------
fun_high.variab.makes_powerlaw_vector<-  function(alpha, n=100000, sd = 0.0002, seed =  123) {
  
  set.seed(seed)
  n <- n
    # 
  alpha <- alpha
  xmin <- 1e-4   # small cutoff to avoid divergence
  
  u <- runif(n)
  
  x <- (abs(rnorm(sd = sd, n = n)))  + (( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))  )
  
  x[100:6000] <- 0.1
  x[4500:15000] <- 0.2
  x[10800:20000] <- 0.01
  
  ; x
  
}

vec1_2.high.variab <- fun_high.variab.makes_powerlaw_vector(2)

length(vec1_2.high.variab) # 100,000 ind.
sum(vec1_2.high.variab) # 1877.934

# Total number of individuals in the sample --------------
N_total = length(vec1_2.high.variab)# N = 100,000 indiv. per sample
N_total

# Abundance (ind.  m-3) -------------
A_total = length(vec1_2.high.variab) / S # A =  N/S , indiv. per cubic meter
A_total # 1000 total Abundance (ind.  m-3)

# Biomass (g.  m-3) ------------
B_total = sum(vec1_2.high.variab) / S # B =  sum(M)/S, g. per cubic meter
B_total #18.7793 total Biomass (g.  m-3)

#  1.3 Very high variability Data ----------------

fun_very.high.variab.makes_powerlaw_vector<-  function(alpha, n=100000, sd = 0.004, seed =  123) {
  
  set.seed(seed)
  n <- n
  
  alpha <- alpha
  xmin <- 1e-4   # small cutoff to avoid divergence
  
  u <- runif(n)
  
  x <- (abs(rnorm(sd = sd, n = n)))  + (( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))  )
  
  x[300:3000] <- 0.005
  x[4500:30000] <- 0.03
  x[10500:50000] <- 0.01
  x[50000:60000] <- 0.5
  
  
  ; x
  
}


vec1_2.very.high.variab <- fun_very.high.variab.makes_powerlaw_vector(2)

length(vec1_2.very.high.variab) # 100,000 ind.
sum(vec1_2.very.high.variab) # 5760.534 g

# Total number of individuals in the sample --------------
N_total = length(vec1_2.very.high.variab)# N = 100,000 indiv. per sample
N_total

# Abundance (ind.  m-3) -------------
A_total = length(vec1_2.very.high.variab) / S # A =  N/S , indiv. per cubic meter
A_total # 1000 total Abundance (ind.  m-3)

# Biomass (g.  m-3) ------------
B_total = sum(vec1_2.very.high.variab) / S # B =  sum(M)/S, g. per cubic meter
B_total # 57.6053 total Biomass (g.  m-3)



# 1.4 Standard log2-binning vector (breaks for binning) ---------------
# length of the vector: 80 elements 
# n <- 80
half <- floor(80 / 2)
exponents <- seq(-half, half - 1)   # symmetric around 0
log2vec <-   2^exponents  # the standard log2-binning vector, these are the breaks!
breaks <- log2vec


##########################

# Select dataset for analysis ---------------

 body_mass.pw <- vec1_2 # perfect log-log-linear power law data, b_bNBS = -1, b_bNNS = -2, 

# body_mass.pw <- vec1_2 # perfect log-log-linear power law data
# body_mass.pw <- vec1_2.high.variab # high variab. power law data
# body_mass.pw <- vec1_2.very.high.variab # very high variab. power law data
# body_mass.pw <- data.carbon_g_ind_1 #  zooplankton data

##########################



# 2. bNNS  --------------------
# ("back-transformed normalized numbers spectrum") ------
# non-linear binning (log2) ---------------

# BWIA: bin-width-inflated Abundance ----------------
# BWIA suffers the non-linear binning distortion effect -------------
# The slope of the Abundance-body mass spectrum shifts from b = -2 (OK) to  b= -1 (distorted by BWIA)

# ============================================================
# 2. bNNS  ( selection: maximum  to first non-zero  bin)
# ============================================================


# ---- Histogram ----
hpw2 <- hist(
  body_mass.pw,
  breaks = log2vec,
  plot = FALSE
)

w2    <- diff(hpw2$breaks)
mids2 <- hpw2$mids
R2    <- diff(range(hpw2$breaks))
k2    <- w2 / R2

BWIA  <- hpw2$counts  # BWIA ("Bin-width inflated abundance": distorted abundance spectrum. non-normalized) ------
NNSS  <- BWIA / w2    # NNSS ("normalized numbers size spectrum: "w-normalized" abundance spectrum) ------
knNSS <- BWIA / k2    # kNSS ("k-normalized numbers size  spectrum": "k-normalized" abundance spectrum) ------

# Rescaling to original dimensions with F-prime:   --------

F_p <- sum(BWIA / S, na.rm = TRUE) /
  sum(knNSS, na.rm = TRUE)

bNNS <- knNSS * F_p   # bNSS ("back-transformed normalized numbers  spectrum" ---------
#                            (i.e., the "k-normalized" and "F-scaled" abundance spectrum)


# ============================================================
# Data SELECTION from max to last non-zero bin (AS IN SECTION 3)
# ============================================================

select_linear <- function(y_original) {
  
  valid <- which(is.finite(y_original) & y_original > 0)
  
  if (length(valid) < 3) return(integer(0))
  
  peak <- valid[which.max(y_original[valid])]
  
  sel <- peak
  i <- peak + 1
  
  while (i <= length(y_original)) {
    
    if (!is.finite(y_original[i]) || y_original[i] <= 0)
      break
    
    sel <- c(sel, i)
    i <- i + 1
  }
  
  return(sel)
}

selected_idx <- select_linear(bNNS)


# ============================================================
# MAIN REGRESSION (only selected bins)
# ============================================================

x_sel <- log10(mids2[selected_idx])
y_sel <- log10(bNNS[selected_idx])

fit_main <- rlm(y_sel ~ x_sel)


# ============================================================
# BOOTSTRAP 
# ============================================================

B <- 500
nbins <- length(mids2)

boot_bNNS   <- matrix(NA, B, nbins)
boot_slope  <- numeric(B)
boot_interc <- numeric(B)

set.seed(123)

for (i in 1:B) {
  
  sample_i <- sample(body_mass.pw, replace = TRUE)
  
  h_i <- hist(sample_i,
              breaks = log2vec,
              plot = FALSE)
  
  NNSS_i  <- h_i$counts / w2
  knNSS_i <- h_i$counts / k2
  
  Fp_i <- sum(h_i$counts / S, na.rm = TRUE) /
    sum(knNSS_i, na.rm = TRUE)
  
  bNNS_i <- knNSS_i * Fp_i
  
  boot_bNNS[i, ] <- bNNS_i
  
  # ---- Data SELECTION ----
  sel_i <- select_linear(bNNS_i)
  
  if (length(sel_i) > 3) {
    
    fit_i <- rlm(log10(bNNS_i[sel_i]) ~
                  log10(mids2[sel_i]))
    
    boot_slope[i]  <- coef(fit_i)[2]
    boot_interc[i] <- coef(fit_i)[1]
  }
}

# ---- Histogram 95% CI ----
lower_ci <- apply(boot_bNNS, 2, quantile, 0.025, na.rm = TRUE)
upper_ci <- apply(boot_bNNS, 2, quantile, 0.975, na.rm = TRUE)

# ---- Regression envelope ----
x_seq <- seq(min(x_sel), max(x_sel), length.out = 200)

pred_mat <- sapply(1:B, function(i) {
  boot_interc[i] + boot_slope[i] * x_seq
})

env_lower <- apply(pred_mat, 1, quantile, 0.025, na.rm = TRUE)
env_upper <- apply(pred_mat, 1, quantile, 0.975, na.rm = TRUE)
env_mid   <- coef(fit_main)[1] +
  coef(fit_main)[2] * x_seq


# ============================================================
# KDE (kernel density estimation), test different bandwidths (bw)
# ============================================================

log_mass <- log10(body_mass.pw) # kernel density estimation in log space (not log-log)

# Test different bandwidths (bw) --------------

# kde <- density(log_mass,  bw = "SJ")
# kde <- density(log_mass,  bw = "0.2")
# kde <- density(log_mass,  bw = "0.25")
# kde <- density(log_mass,  bw = "nrd0")

 kde <- density(log_mass) # The default, bw = "nrd0". Silverman's ‘rule of thumb’, Silverman (1986).  "bw.nrd0" implements a rule-of-thumb for choosing the bandwidth of a Gaussian kernel density estimator. It defaults to 0.9 times the minimum of the standard deviation and the interquartile range divided by 1.34 times the sample size to the negative one-fifth power (= Silverman's ‘rule of thumb’, Silverman (1986, page 48, eqn (3.31))) unless the quartiles coincide when a positive result will be guaranteed

mass_kde <- 10^kde$x
kde_mass_density <- kde$y / (mass_kde * log10(10))
kde_abund <- kde_mass_density * S

valid_kde <- kde_abund > 0
x_kde <- log10(mass_kde[valid_kde])
y_kde <- log10(kde_abund[valid_kde])


# ============================================================
# PLOT
# ============================================================

valid_plot <- bNNS > 0 & lower_ci > 0 & upper_ci > 0

x_plot <- log10(mids2[valid_plot])
y_plot <- log10(bNNS[valid_plot])
y_lower <- log10(lower_ci[valid_plot])
y_upper <- log10(upper_ci[valid_plot])

ylim_all <- range(c(y_lower,
                    y_upper,
                    env_lower,
                    env_upper,
                    y_kde),
                  na.rm = TRUE)

plot(x_plot, y_plot,
     pch = 16,
     col = "darkgreen",
     ylim = ylim_all,
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Abundance (ind. m-3)",
     main = "bNNS with bootstrap, regression envelope and KDE")

# Histogram CI
arrows(x_plot,
       y_lower,
       x_plot,
       y_upper,
       angle = 90,
       code = 3,
       length = 0.04,
       col = "grey40")

# Regression envelope
polygon(c(x_seq, rev(x_seq)),
        c(env_lower, rev(env_upper)),
        col = "darkorange",
        border = NA)

# Regression line
lines(x_seq, env_mid,
      col = "black",
      lwd = 1.4)

# KDE
lines(x_kde, y_kde,
      col = "red",
      lwd = 2)

points(x_plot, y_plot,
       pch = 16,
       col = "darkgreen")


##########################


# 3. bNBS -----------
# bNBS (new!): backtransformed  and normalized numbers  spectrum (Biomass vs body mass, normalized, AND rescaled)

# 2.1 bNBS (simple, with common regression line, no bootstrapping) -------------

#  Biomass is directly binned (not through abundance) ---------------

# non-linear binning, log2 ----------------

binnedB2 <- cut(body_mass.pw, breaks = log2vec, include.lowest = TRUE)

h4 <- hist(body_mass.pw, breaks = log2vec, plot = F) # to obtain mids, etc
mids2 <- h4$mids

biomass_per_bin2 <- tapply(body_mass.pw, binnedB2, sum)
biomass_per_bin2

# Biomass per bin
biomass_per_bin2b <- biomass_per_bin2
BWIB <- biomass_per_bin2b/S

sum(BWIB, na.rm = T)
BWIB_total <- sum(BWIB, na.rm = T)

B_total =sum(biomass_per_bin2, na.rm = T) # g m-3 Total Biomass,   OK!
B_total
# 0.00342 total Biomass (g.  m-3), for zooplankton  sample example OK, ind
# 0.9067911 for Biomass (g.  m-3)simulated power-law data , 90.67 g per sample

sum(biomass_per_bin2b,na.rm = T)
# 0.00342 total Biomass (g.  m-3), OK

biomass_per_bin2 <- biomass_per_bin2b


# Biomass (g.  m-3)
B_total = sum(body_mass.pw) / 100 # B =  sum(M)/S, g. per cubic meter
B_total # 0.00342 total Biomass (g.  m-3)


plot(biomass_per_bin2b ~ mids2)

plot(log10(biomass_per_bin2b) ~ log10(mids2), main  = "BWIB, bin-width inflated biomass")  # flat slpe b = 0, as expected for BWIB!
#abline ( h = -8)
# biomass_per_bin2 ins BWIB! flat dostributin because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log10(biomass_per_bin2[selected_idx]) ~ log10(mids2[selected_idx]), col = "navyblue", pch = 16)

lm.BWIB <- rlm(log10(biomass_per_bin2[selected_idx]) ~ log10(mids2[selected_idx]) )
summary(lm.BWIB) 
abline(lm.BWIB)
# nearly flat BWIB- body mass spectrum , b = 0 !

# "For a common power-law distribution, the use of 
# geometrically increasing bin sizes will lead to a flat (slope = 0) , 
# unrealistic distribution of the BWIB,
# while the NBSS has a realistic power-law shape, with a linearly 
# downtrending  slope on a  log-log-scaled  plot, that represents well the original 
# biomass-body-mass relationship of the natural ecosystem." Schwamborn, 2026 

# biomass_per_bin2 is BWIB! flat distributin because of distorting effect of non-linear binning!
#  BWIB = Bin-width inflated Biomass

# the common NBSS ------------
NBSS <- as.vector(  BWIB/w2)
sum(NBSS, na.rm = T) #  OK

plot(log10(NBSS) ~ log10(mids2), main = "NBSS", 
     xlab = "body mass, log10 (g ind.-1)",
     ylab= "NBSS ( g m-3 / g ind.-i)") 

points(log10(NBSS[selected_idx]) ~ log10(mids2[selected_idx]), col = "navyblue", pch = 16)

lm.NBSS <- rlm(log10(NBSS[selected_idx]) ~ log10(mids2[selected_idx]))
summary(lm.NBSS) # slope = -1 , OK
abline(lm.NBSS)

# All NBSS data  are then transformed back (re-dimensionalized back and rescaled back)  into the original scale, dimension, and units of biomass, i.e., into  “back-transformed normalized biomass” (bNB), where


# bNBi = NBSSi * D

# Since D has units of B * M (e.g.,  g² m-3 ind.-1), this calculation backtransforms the units of the NBSS into the original biomass units, within the  bNBS (e.g., g m-3). 


# D =  BWIB_total  / NBSS_total

BWIB <- biomass_per_bin2   # BWIB: "bin-width-inflated biomass" (not yet normalized)

B_total  #0.9067

BWIB_total = sum(BWIB,na.rm = T)/S # 0.003427 = total Biomass, OK!
BWIB_total # 0.9067911 , should  be equal to B_total,  in g m-3

# Check  if  B_total is equal to  BWIB_total
if (isTRUE(all.equal(B_total, BWIB_total))) {
  message("OK: B_total and BWIB_total are equal.")
} else {
  stop("Error: B_total and BWIB_total are not equal.")
}

NBSS_total = sum(NBSS,na.rm = T)
NBSS_total

# rescale with D ----------
# After normalizing (dividing BWIB by w2, where NBSS = BWIB / w2) to obtain the common NBSS, rescaling is needed. 
#This can be done by using the w-scaled dimensional correction factor D, where
# D = BWIB_total / NBSS_total 

D = BWIB_total / NBSS_total
D

bNB = NBSS * D   #  # bNb ("back-transformed normalized biomass" ---------
#                            (i.e., the "w-normalized" and "D-scaled" biomass spectrum)


bNB

sum(bNB, na.rm = T) # 0.00342 g m-3 = total Biomass, OK!

bNB_total <- sum(bNB, na.rm = T) # should be equal to B_total

# Check  if  B_total is equal to  BWIB_total
if (isTRUE(all.equal(B_total, bNB_total))) {
  message("OK: B_total and bNB_total are equal.")
} else {
  stop("Error: B_total and bNB_total are not equal.")
}


# Make a plot for bNBS  vs Body mmass (the bNBS plot) -----------

plot(log10(bNB) ~ log10(mids2), 
     main = "bNBS plot",
     ylab ="Biomass (bNB), log10(g m-3)" ,
     xlab = "Body mass, log10(g ind.-1)") 


points(log10(bNB[selected_idx]) ~ log10(mids2[selected_idx]), col = "navyblue", pch = 16)

lm.bNB <- rlm(log10(bNB[selected_idx]) ~ log10(mids2[selected_idx]))
summary(lm.bNB)
abline(lm.bNB)

# slope = approx b = -1, OK!


# 2.2 bNBS (bootstrapping the histogram, with 95% CI error bars) -------------

# 2.2 bNBS (bootstrapped histogram, 95% CI) -------------

B <- 500
nbins <- length(mids2)

boot_bNB <- matrix(NA, nrow = B, ncol = nbins)

set.seed(123)

for (i in 1:B) {
  
  # --- resample individuals ---
  sample_i <- sample(body_mass.pw, replace = TRUE)
  
  # --- biomass per bin ---
  binned_i <- cut(sample_i, breaks = log2vec, include.lowest = TRUE)
  biomass_i <- tapply(sample_i, binned_i, sum)
  
  biomass_i[is.na(biomass_i)] <- 0
  
  # --- BWIB ---
  BWIB_i <- biomass_i
  
  BWIB_total_i <- sum(BWIB_i, na.rm = TRUE) / S
  
  # --- NBSS ---
  NBSS_i <- BWIB_i / w2
  
  NBSS_total_i <- sum(NBSS_i, na.rm = TRUE)
  
  # --- rescale factor ---
  D_i <- BWIB_total_i / NBSS_total_i
  
  # --- backtransformed normalized biomass ---
  bNB_i <- NBSS_i * D_i
  
  boot_bNB[i, ] <- bNB_i
}

# ---- 95% CI in linear space ----
lower_ci <- apply(boot_bNB, 2, quantile, 0.025, na.rm = TRUE)
upper_ci <- apply(boot_bNB, 2, quantile, 0.975, na.rm = TRUE)

# ---- Plot in log10 space ----
valid <- bNB > 0 & lower_ci > 0 & upper_ci > 0

x_plot <- log10(mids2[valid])
y_plot <- log10(bNB[valid])
y_lower <- log10(lower_ci[valid])
y_upper <- log10(upper_ci[valid])

ylim_all <- range(c(y_lower, y_upper), na.rm = TRUE)

plot(x_plot, y_plot,
     pch = 16,
     col = "darkgreen",
     ylim = ylim_all,
     main = "bNBS with bootstrap 95% CI",
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Biomass (g m-3)")

# Error bars
arrows(x_plot,
       y_lower,
       x_plot,
       y_upper,
       angle = 90,
       code = 3,
       length = 0.05,
       col = "grey40")

points(x_plot, y_plot,
       pch = 16,
       col = "darkgreen")


# 2. bNBS (bootstrapping the histogram, with with regression line and bootstrapped regression envelope) -------------

# 2.3a no range selection (uses ALL Data) bNBS (bootstrapping the histogram, with with regression line and bootstrapped regression envelope) -------------

B <- 500
nbins <- length(mids2)

boot_bNB   <- matrix(NA, B, nbins)
boot_slope <- numeric(B)
boot_int   <- numeric(B)

set.seed(123)

for (i in 1:B) {
  
  # --- resample individuals ---
  sample_i <- sample(body_mass.pw, replace = TRUE)
  
  # --- biomass per bin ---
  binned_i <- cut(sample_i, breaks = log2vec, include.lowest = TRUE)
  biomass_i <- tapply(sample_i, binned_i, sum)
  biomass_i[is.na(biomass_i)] <- 0
  
  # --- BWIB ---
  BWIB_i <- biomass_i
  
  BWIB_total_i <- sum(BWIB_i, na.rm = TRUE) / S
  
  # --- NBSS ---
  NBSS_i <- BWIB_i / w2
  NBSS_total_i <- sum(NBSS_i, na.rm = TRUE)
  
  # --- rescale ---
  D_i <- BWIB_total_i / NBSS_total_i
  bNB_i <- NBSS_i * D_i
  
  boot_bNB[i, ] <- bNB_i
  
  # --- regression (log-log) ---
  valid_i <- bNB_i > 0
  
  if (sum(valid_i) > 3) {
    fit_i <- rlm(log10(bNB_i[valid_i]) ~ log10(mids2[valid_i]))
    boot_slope[i] <- coef(fit_i)[2]
    boot_int[i]   <- coef(fit_i)[1]
  }
}

# ---- Histogram 95% CI ----
lower_ci <- apply(boot_bNB, 2, quantile, 0.025, na.rm = TRUE)
upper_ci <- apply(boot_bNB, 2, quantile, 0.975, na.rm = TRUE)

# ---- Main regression (original data) ----

# ---- Main regression (original data) ----
valid <- is.finite(bNB) & bNB > 0 &
  is.finite(mids2) & mids2 > 0

x <- log10(mids2[valid])
y <- log10(bNB[valid])

# remove non-finite log results explicitly
finite_idx <- is.finite(x) & is.finite(y)

x <- x[finite_idx]
y <- y[finite_idx]

fit_main <- rlm(y ~ x)


# ---- Regression envelope ----
x_seq <- seq(min(x), max(x), length.out = 200)

pred_mat <- sapply(1:B, function(i) {
  boot_int[i] + boot_slope[i] * x_seq
})

env_lower <- apply(pred_mat, 1, quantile, 0.025, na.rm = TRUE)
env_upper <- apply(pred_mat, 1, quantile, 0.975, na.rm = TRUE)
env_mid   <- coef(fit_main)[1] + coef(fit_main)[2] * x_seq

# ---- Y limits ----
y_lower_plot <- log10(lower_ci[valid])
y_upper_plot <- log10(upper_ci[valid])

ylim_all <- range(c(y_lower_plot,
                    y_upper_plot,
                    env_lower,
                    env_upper),
                  na.rm = TRUE)

# ======================
#         PLOT
# ======================

plot(x, y,
     pch = 16,
     col = "darkgreen",
     ylim = ylim_all,
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Biomass (g m-3)",
     main = "bNBS with bootstrap and regression envelope")

# Histogram CI
arrows(x,
       y_lower_plot,
       x,
       y_upper_plot,
       angle = 90,
       code = 3,
       length = 0.04,
       col = "grey40")

# Regression envelope
polygon(c(x_seq, rev(x_seq)),
        c(env_lower, rev(env_upper)),
        col = rgb(0,0,0,0.15),
        border = NA)

# Regression line
lines(x_seq, env_mid,
      col = "black",
      lwd = 1.4)

points(x, y,
       pch = 16,
       col = "darkgreen")


# 2.3c bNBS ,  plots all data, bnBS with range selection, but plots ALL selcted data ------------------------
#  bNBS (bootstrap + regression + envelope, dynamic range selection) ---------
# FINAL bNBS plot for Paper ----------------
# (bootstrap + regression + envelope, dynamic range selection) ---------

B <- 500
nbins <- length(mids2)

boot_bNB   <- matrix(NA, B, nbins)
boot_slope <- numeric(B)
boot_int   <- numeric(B)

set.seed(123)

for (i in 1:B) {
  
  # --- resample individuals ---
  sample_i <- sample(body_mass.pw, replace = TRUE)
  
  # --- biomass per bin ---
  binned_i <- cut(sample_i, breaks = log2vec, include.lowest = TRUE)
  biomass_i <- tapply(sample_i, binned_i, sum)
  biomass_i[is.na(biomass_i)] <- 0
  
  BWIB_i <- biomass_i
  BWIB_total_i <- sum(BWIB_i) / S
  
  NBSS_i <- BWIB_i / w2
  NBSS_total_i <- sum(NBSS_i)
  
  D_i <- BWIB_total_i / NBSS_total_i
  bNB_i <- NBSS_i * D_i
  
  boot_bNB[i, ] <- bNB_i
  
  # ============================
  #   DYNAMIC RANGE SELECTION
  # ============================
  
  y_i <- log10(bNB_i)
  x_i <- log10(mids2)
  
  # must be finite
  finite_idx <- is.finite(y_i) & is.finite(x_i)
  
  if (sum(finite_idx) > 5) {
    
    y_i <- y_i[finite_idx]
    x_i <- x_i[finite_idx]
    
    # find maximum
    max_idx_i <- which.max(y_i)
    
    # last finite bin
    last_non_na_i <- max(which(is.finite(y_i)))
    
    if (last_non_na_i > max_idx_i) {
      
      selected_idx_i <- max_idx_i:last_non_na_i
      
      if (length(selected_idx_i) > 3) {
        
        fit_i <- rlm(y_i[selected_idx_i] ~ x_i[selected_idx_i])
        
        boot_slope[i] <- coef(fit_i)[2]
        boot_int[i]   <- coef(fit_i)[1]
      }
    }
  }
}

# ============================
#   MAIN DATA REGRESSION
# ============================

y <- log10(bNB)
x <- log10(mids2)

finite_idx <- is.finite(y) & is.finite(x)

y <- y[finite_idx]
x <- x[finite_idx]

max_idx <- which.max(y)
last_non_na <- max(which(is.finite(y)))

selected_idx <- max_idx:last_non_na

fit_main <- rlm(y[selected_idx] ~ x[selected_idx])

# ============================
#   REGRESSION ENVELOPE
# ============================

x_seq <- seq(min(x[selected_idx]),
             max(x[selected_idx]),
             length.out = 200)

pred_mat <- sapply(1:B, function(i) {
  boot_int[i] + boot_slope[i] * x_seq
})

env_lower <- apply(pred_mat, 1, quantile, 0.025, na.rm = TRUE)
env_upper <- apply(pred_mat, 1, quantile, 0.975, na.rm = TRUE)
env_mid   <- coef(fit_main)[1] + coef(fit_main)[2] * x_seq

# ============================
#   PLOT
# ============================

ylim_all <- range(c(y[selected_idx],
                    env_lower,
                    env_upper),
                  na.rm = TRUE)

plot(x[selected_idx],
     y[selected_idx],
     pch = 16,
     col = "darkgreen",
     ylim = ylim_all,
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Biomass (g m-3)",
     main = "bNBS with bootstrap regression envelope")

# envelope
polygon(c(x_seq, rev(x_seq)),
        c(env_lower, rev(env_upper)),
        col = rgb(0,0,0,0.15),
        border = NA)

# regression line
lines(x_seq, env_mid, lwd = 2)

points(x[selected_idx],
       y[selected_idx],
       pch = 16,
       col = "darkgreen")

# ============================
#   PREPARE DATA
# ============================

y_all <- log10(bNB)
x_all <- log10(mids2)

finite_idx <- is.finite(y_all) & is.finite(x_all)

x_all <- x_all[finite_idx]
y_all <- y_all[finite_idx]

# dynamic selection (already validated as correct)
max_idx <- which.max(y_all)
last_non_na <- max(which(is.finite(y_all)))
selected_idx <- max_idx:last_non_na

# ============================
#   PLOT LIMITS
# ============================

ylim_all <- range(c(y_all,
                    env_lower,
                    env_upper),
                  na.rm = TRUE)

# ============================
#   PLOT
# ============================

plot(x_all, y_all,
     pch = 16,
     col = "grey70",      # all bins in grey
     ylim = ylim_all,
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Biomass (g m-3)",
     main = "bNBS with bootstrap regression envelope")

# Highlight selected linear section
points(x_all[selected_idx],
       y_all[selected_idx],
       pch = 16,
       col = "darkgreen")

# Regression envelope
polygon(c(x_seq, rev(x_seq)),
        c(env_lower, rev(env_upper)),
        col = rgb(0,0,0,0.15),
        border = NA)

# Regression line
lines(x_seq, env_mid,
      lwd = 2)

# Optional legend
legend("topright",
       legend = c("All bins",
                  "Selected linear range",
                  "Bootstrap envelope",
                  "Regression line"),
       col = c("grey70", "darkgreen", rgb(0,0,0,0.3), "black"),
       pch = c(16,16,15,NA),
       lwd = c(NA,NA,NA,2),
       bty = "n")

# 2.3d bNBS with 95% CI  bars,  plots all data, bnBS with range selection, but plots ALL selcted data ------------------------
#  bNBS (bootstrap + regression + envelope, dynamic range selection) ---------
# FINAL bNBS plot for Paper ----------------
# (bootstrap + regression + envelope, dynamic range selection) ---------
# Selected regression range (dark green)
# 95% bootstrap CI error bars (histogram uncertainty)
# Regression line
# Bootstrap 95% Confidence  envelope around regression line

####################################
# FINAL bnBS PLOT FOR PAPER ----------
####################################

# ============================
#   PREPARE DATA
# ============================

y_all <- log10(bNB)
x_all <- log10(mids2)

finite_idx <- is.finite(y_all) & is.finite(x_all)

x_all <- x_all[finite_idx]
y_all <- y_all[finite_idx]

lower_log <- log10(lower_ci[finite_idx])
upper_log <- log10(upper_ci[finite_idx])

# dynamic regression range
max_idx <- which.max(y_all)
last_non_na <- max(which(is.finite(y_all)))
selected_idx <- max_idx:last_non_na

# ============================
#   PLOT LIMITS
# ============================

ylim_all <- range(c(y_all,
                    lower_log,
                    upper_log,
                    env_lower,
                    env_upper),
                  na.rm = TRUE)

# ============================
#   PLOT
# ============================

plot(x_all, y_all,
     pch = 16,
     col = "grey70",
     ylim = ylim_all,
     xlab = "log10 Body mass (g ind.-1)",
     ylab = "log10 Biomass (g m-3)",
     main = "bNBS with bootstrap CI and regression envelope")

# ---- 95% CI error bars (all bins) ----
arrows(x_all,
       lower_log,
       x_all,
       upper_log,
       angle = 90,
       code = 3,
       length = 0.04,
       col = "grey50")

# ---- Highlight selected regression range ----
points(x_all[selected_idx],
       y_all[selected_idx],
       pch = 16,
       col = "darkgreen")

# ---- Regression envelope ----
polygon(c(x_seq, rev(x_seq)),
        c(env_lower, rev(env_upper)),
        col = rgb(0,0,0,0.15),
        border = NA)

# ---- Regression line ----
lines(x_seq, env_mid,
      lwd = 2)

# ---- Legend ----
legend("topright",
       legend = c("All bins",
                  "Selected linear range",
                  "95% CI (histogram)",
                  "Bootstrap envelope",
                  "Regression line"),
       col = c("grey70", "darkgreen", "grey50",
               rgb(0,0,0,0.3), "black"),
       pch = c(16,16,NA,15,NA),
       lwd = c(NA,NA,1,NA,2),
       bty = "n")


#  4. FUNCTION ---------------

# the complete analysis outputs as a function
# FUNCTION
# A  convenient function is  provided at the end of the script
# "fun.bNNS_and_bNBS_plots_summary_table" 
#       gives a  summary results table and two plots 

# The function 
# runs a two-step bootstrap (for raw data and for the linear regression)
# Produces  two  plots (bNNS + bNBS)
# Uses  "maximum to last non- zero bin" data selection
# Shows 95% confidence regression envelope for each graph 
# The 95% confidence envelopes  are not  parametric. They are generated from the 2.5% and 97.5% quantiles of the bootstrap posterior distributions of the regression parameters.
# Shows a KDE high-resolution density distribution line, for bNNS only (not possible for biomass, for abundance only))
# Computes bootstrap slope + intercept 95% CIs
# Computes totals (total Abundance and total Biomass, for all data, and for the linear range only)
# Returns a  summary results table ("res$summary_table")
# plots contains the "a" and "b" estimates, with bootstrap 95% CI,and R-squared


fun.bNNS_and_bNBS_plots_summary_table  <- function(body_mass,
                                           S = 100,
                                           B = 500,
                                           breaks = NULL){
    
    body_mass <- as.numeric(body_mass)
    body_mass <- body_mass[is.finite(body_mass) & body_mass > 0]
    
    if(length(body_mass) < 10)
      stop("Not enough valid observations.")
    
    if(is.null(breaks)){
      half <- floor(80/2)
      breaks <- 2^seq(-half, half-1)
    }
    
    h <- hist(body_mass, breaks = breaks, plot = FALSE)
    
    mids  <- h$mids
    w     <- diff(h$breaks)
    R     <- diff(range(h$breaks))
    k     <- w/R
    x_log <- log10(mids)
    nbins <- length(mids)
    
    select_linear <- function(y){
      valid <- which(is.finite(y) & y > 0)
      if(length(valid) < 3) return(integer(0))
      peak <- valid[which.max(y[valid])]
      sel <- peak
      i <- peak + 1
      while(i <= length(y)){
        if(!is.finite(y[i]) || y[i] <= 0) break
        sel <- c(sel, i)
        i <- i + 1
      }
      sel
    }
    
    set.seed(123)
    
    # ============================================================
    # ========================== bNNS ============================
    # ============================================================
    
    counts <- h$counts
    kn     <- counts/k
    Fp     <- sum(counts/S) / sum(kn)
    bNNS   <- kn * Fp
    
    selA <- select_linear(bNNS)
    yA   <- log10(bNNS)
    
    fitA <- rlm(yA[selA] ~ x_log[selA])
    aA  <- coef(fitA)[1]
    bA  <- coef(fitA)[2]
  
    
    custom_r2 <- function(model) {
      y_true <- model$model[[1]]      # response variable
      y_pred <- fitted(model)         # predicted values
      
      ss_res <- sum((y_true - y_pred)^2)
      ss_tot <- sum((y_true - mean(y_true))^2)
      
      r2 <- 1 - (ss_res / ss_tot)
      return(r2)
    }
    
    R2A <- custom_r2(fitA)
    
    
     # R2A <- summary(fitA)$r.squared
    
    boot_sA <- boot_iA <- rep(NA, B)
    boot_mat_A <- matrix(NA, B, nbins)
    
    for(i in 1:B){
      samp <- sample(body_mass, replace=TRUE)
      h_i  <- hist(samp, breaks=breaks, plot=FALSE)
      kn_i <- h_i$counts/k
      Fp_i <- sum(h_i$counts/S)/sum(kn_i)
      b_i  <- kn_i*Fp_i
      boot_mat_A[i,] <- b_i
      
      sel_i <- select_linear(b_i)
      if(length(sel_i) > 3){
        fit_i <- rlm(log10(b_i[sel_i]) ~ x_log[sel_i])
        boot_sA[i] <- coef(fit_i)[2]
        boot_iA[i] <- coef(fit_i)[1]
      }
    }
    
    validA <- which(is.finite(boot_sA))
    slopeCI_A <- quantile(boot_sA[validA], c(.025,.975))
    intCI_A   <- quantile(boot_iA[validA], c(.025,.975))
    
    bin_lo_A <- apply(boot_mat_A,2,quantile,.025,na.rm=TRUE)
    bin_hi_A <- apply(boot_mat_A,2,quantile,.975,na.rm=TRUE)
    
    x_seqA <- seq(min(x_log[selA]), max(x_log[selA]), length=200)
    predA  <- sapply(validA,function(i)
      boot_iA[i] + boot_sA[i]*x_seqA)
    envA_lo <- apply(predA,1,quantile,.025)
    envA_hi <- apply(predA,1,quantile,.975)
    
    # KDE
    log_mass <- log10(body_mass)
    kde <- density(log_mass)
    mass_kde <- 10^kde$x
    dens_mass <- kde$y/(mass_kde*log(10))
    kde_abund <- dens_mass*(length(body_mass)/S)
    valid_kde <- kde_abund > 0
    
    # limits
    xlim_A <- range(x_log[bNNS > 0], finite=TRUE)
    ylim_A <- range(c(
      yA[bNNS>0],
      log10(bin_lo_A[bNNS>0]),
      log10(bin_hi_A[bNNS>0]),
      envA_lo, envA_hi,
      log10(kde_abund[valid_kde])
    ), finite=TRUE)
    
    # ============================================================
    # ========================== bNBS ============================
    # ============================================================
    
    bins  <- cut(body_mass,breaks=breaks,include.lowest=TRUE)
    bio   <- tapply(body_mass,bins,sum)
    bio[is.na(bio)] <- 0
    
    NBSS  <- bio/w
    D     <- (sum(bio)/S)/sum(NBSS)
    bNBS  <- NBSS*D
    
    selB <- select_linear(bNBS)
    yB   <- log10(bNBS)
    
    fitB <- rlm(yB[selB] ~ x_log[selB])
    aB  <- coef(fitB)[1]
    bB  <- coef(fitB)[2]
   
    
    
    custom_r2 <- function(model) {
      y_true <- model$model[[1]]      # response variable
      y_pred <- fitted(model)         # predicted values
      
      ss_res <- sum((y_true - y_pred)^2)
      ss_tot <- sum((y_true - mean(y_true))^2)
      
      r2 <- 1 - (ss_res / ss_tot)
      return(r2)
    }
    
    R2B <- custom_r2(fitB)
    
     #R2B <- summary(fitB)$r.squared
    
    boot_sB <- boot_iB <- rep(NA,B)
    boot_mat_B <- matrix(NA,B,nbins)
    
    for(i in 1:B){
      samp <- sample(body_mass, replace=TRUE)
      bins_i <- cut(samp,breaks=breaks,include.lowest=TRUE)
      bio_i  <- tapply(samp,bins_i,sum)
      bio_i[is.na(bio_i)] <- 0
      NBSS_i <- bio_i/w
      D_i    <- (sum(bio_i)/S)/sum(NBSS_i)
      b_i    <- NBSS_i*D_i
      boot_mat_B[i,] <- b_i
      
      sel_i <- select_linear(b_i)
      if(length(sel_i)>3){
        fit_i <- rlm(log10(b_i[sel_i]) ~ x_log[sel_i])
        boot_sB[i] <- coef(fit_i)[2]
        boot_iB[i] <- coef(fit_i)[1]
      }
    }
    
    validB <- which(is.finite(boot_sB))
    slopeCI_B <- quantile(boot_sB[validB],c(.025,.975))
    intCI_B   <- quantile(boot_iB[validB],c(.025,.975))
    
    bin_lo_B <- apply(boot_mat_B,2,quantile,.025,na.rm=TRUE)
    bin_hi_B <- apply(boot_mat_B,2,quantile,.975,na.rm=TRUE)
    
    x_seqB <- seq(min(x_log[selB]), max(x_log[selB]), length=200)
    predB  <- sapply(validB,function(i)
      boot_iB[i] + boot_sB[i]*x_seqB)
    envB_lo <- apply(predB,1,quantile,.025)
    envB_hi <- apply(predB,1,quantile,.975)
    
    xlim_B <- range(x_log[bNBS > 0], finite=TRUE)
    ylim_B <- range(c(
      yB[bNBS>0],
      log10(bin_lo_B[bNBS>0]),
      log10(bin_hi_B[bNBS>0]),
      envB_lo, envB_hi
    ), finite=TRUE)
    
    # ============================================================
    # ======================== PLOTS =============================
    # ============================================================
    
    add_stats <- function(a,aCI,b,bCI,R2){
      legend("topright",
             legend=c(
               sprintf("a (95%% CI): %.3f [%.3f, %.3f]",a,aCI[1],aCI[2]),
               sprintf("b (95%% CI): %.3f [%.3f, %.3f]",b,bCI[1],bCI[2]),
               sprintf("R²: %.3f",R2)
             ),
             bty="n",cex=.9)
    }
    
    # 1) bNNS
    plot(x_log,yA,pch=16,col="grey70",
         xlim=xlim_A,ylim=ylim_A,
         main="bNNS",
         xlab="log10 Body mass",
         ylab="log10 Abundance")
    
    arrows(x_log,log10(bin_lo_A),
           x_log,log10(bin_hi_A),
           angle=90,code=3,length=.03,col="grey50")
    
    polygon(c(x_seqA,rev(x_seqA)),
            c(envA_lo,rev(envA_hi)),
            col=adjustcolor("steelblue",0.4),
            border="steelblue")
    
    lines(log10(mass_kde[valid_kde]),
          log10(kde_abund[valid_kde]),
          col="red",lwd=2)
    
    points(x_log[selA],yA[selA],col="darkgreen",pch=16)
    abline(fitA,lwd=2)
    add_stats(aA,intCI_A,bA,slopeCI_A,R2A)
    
    # 2) bNBS
    plot(x_log,yB,pch=16,col="grey70",
         xlim=xlim_B,ylim=ylim_B,
         main="bNBS",
         xlab="log10 Body mass",
         ylab="log10 Biomass")
    
    arrows(x_log,log10(bin_lo_B),
           x_log,log10(bin_hi_B),
           angle=90,code=3,length=.03,col="grey50")
    
    polygon(c(x_seqB,rev(x_seqB)),
            c(envB_lo,rev(envB_hi)),
            col=adjustcolor("steelblue",0.4),
            border="steelblue")
    
    points(x_log[selB],yB[selB],col="darkgreen",pch=16)
    abline(fitB,lwd=2)
    add_stats(aB,intCI_B,bB,slopeCI_B,R2B)
    
    # 3) two-in-one
    dev.new()
    par(mfrow=c(1,2))
    plot(x_log,yA,pch=16,col="grey70",
         xlim=xlim_A,ylim=ylim_A,
         main="bNNS")
    polygon(c(x_seqA,rev(x_seqA)),
            c(envA_lo,rev(envA_hi)),
            col=adjustcolor("steelblue",0.4),
            border="steelblue")
    abline(fitA,lwd=2)
    
    plot(x_log,yB,pch=16,col="grey70",
         xlim=xlim_B,ylim=ylim_B,
         main="bNBS")
    polygon(c(x_seqB,rev(x_seqB)),
            c(envB_lo,rev(envB_hi)),
            col=adjustcolor("steelblue",0.4),
            border="steelblue")
    abline(fitB,lwd=2)
    par(mfrow=c(1,1))
    
    # ============================================================
    # ====================== SUMMARY TABLE =======================
    # ============================================================
    
    summary_table <- data.frame(
      Metric = c("Slope (b)",
                 "Intercept (a)",
                 "R2",
                 "Total abundance",
                 "Abundance linear",
                 "Total biomass",
                 "Biomass linear"),
      bNNS = c(
        sprintf("%.3f [%.3f, %.3f]", bA, slopeCI_A[1], slopeCI_A[2]),
        sprintf("%.3f [%.3f, %.3f]", aA, intCI_A[1], intCI_A[2]),
        sprintf("%.3f", R2A),
        sprintf("%.2f", sum(counts)),
        sprintf("%.2f", sum(bNNS[selA])),
        NA,
        NA
      ),
      bNBS = c(
        sprintf("%.3f [%.3f, %.3f]", bB, slopeCI_B[1], slopeCI_B[2]),
        sprintf("%.3f [%.3f, %.3f]", aB, intCI_B[1], intCI_B[2]),
        sprintf("%.3f", R2B),
        NA,
        NA,
        sprintf("%.2f", sum(bio)),
        sprintf("%.2f", sum(bNBS[selB]))
      )
    )
    
    result <- list(
      x = x_log,
      y_bNNS = bNNS,
      y_bNBS = bNBS,
      summary_table = summary_table
    )
    
    return(result)
  }
  
# 5. Test the function --------

# 1. perfect log-log-linear power law data ( n= 100,000) ----

res <- fun.bNNS_and_bNBS_plots_summary_table(vec1_2, S = 100, B   = 400, 
                                          breaks = log2vec) 

res$summary_table

# Metric                    bNNS                      bNBS
# 1        Slope (b)       -2.002 [-2.052, -1.973]   -1.006 [-1.044, -0.977]
# 2    Intercept (a)       -4.828 [-4.960, -4.748]   -4.210 [-4.344, -4.103]
# 3               R2       1.000                      1.000
# 4  Total abundance       1000.00                    <NA>
# 5 Abundance linear       600.64                     <NA>
# 6    Total biomass       <NA>                       0.91
# 7   Biomass linear       <NA>                       0.70


# 2. High variability synthetic power law data ( n= 100,000) ----

res <- fun.bNNS_and_bNBS_plots_summary_table(vec1_2.high.variab, S = 100, B   = 400, 
                                             breaks = log2vec) 

res$summary_table



# 3. Very high variability synthetic power law data ( n= 100,000) ----

res <- fun.bNNS_and_bNBS_plots_summary_table(vec1_2.very.high.variab, S = 100, 
                                             B   = 400, 
                                             breaks = log2vec) 

res$summary_table


# 4. zooplankton  data  ----

res <- fun.bNNS_and_bNBS_plots_summary_table(data.carbon_g_ind_1, S = 100, 
                                             B   = 400, 
                                             breaks = log2vec) 

res$summary_table






######################
  # End of Script #
######################

