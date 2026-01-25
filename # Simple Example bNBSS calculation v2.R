#  Testing the bNBSS size spectrum method
# Version   2.1

# bNBSS: backtransformed and normalized numbers size spectrum (Biomass vs body mass, normalized, AND rescaled)

# R script for the paper
# "Is normalized biomass really abundance? 
# Ralf Schwamborn 
# January 2026
# Source : "https://raw.githubusercontent.com/rschwamborn/bNBS


# I.   Import Data
# II.  NNSS: Abundance spectrum
# III. bNNS: the "backtransformed” NNSS (new method)
# IV.  NBSS: Biomass (NBSS: normalized biomass size spectrum)
# V.   bNBS: "back-transformed normalized biomass" (new method)



# I. Import Data, example dataset ----------

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


# S: sampling volume (S = 100 cubic meters, in this simple example)

S = 100

# N: Total number of individual in the sample
N_total = length(data.carbon_g_ind_1)# N = 7088, indiv. per sample
N_total

# A; Abundance (ind.  m-3)
A_total = length(data.carbon_g_ind_1) / 100 # A =  N/S , indiv. per cubic meter
A_total # 70.88 total Abundance (ind.  m-3)

# B: Biomass (g.  m-3)
B_total = sum(data.carbon_g_ind_1) / 100 # B =  sum(M)/S, g. per cubic meter
B_total # 0.00342 total Biomass (g.  m-3)


# II. NNSS: Abundance spectrum  ------------
#        Abundance-body mass spectrum, slope is expected to be close to (- alpha) =  -2

# BINNING --------------

# use the standard log2 binning vector

# Create and use the standard log2-binning vector (bin mids = ...,0.25,0.5,1,2,4,8,…) 
# Create log2 bins, with 80 elements, centered at "1"

#length of the vector: n elements (80 in this case)
n <- 80
half <- floor(n / 2)
exponents <- seq(-half, half - 1)   # symmetric around 0
x <- 2^exponents
log2vec <- x  # the standard log2-binning vector (used as "breaks",  not mids)

hpw2c <- hist(
  data.carbon_g_ind_1,
  breaks = log2vec,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

hpw2c$counts
sum(hpw2c$counts)# total N = 7088, OK

# BWIA: bin-width-inflated abundance
BWIA = hpw2c$counts / S
BWIA
sum(BWIA)  #  total N = 70,88, OK

plot(BWIA ~ hpw2c$mids)
plot(log10(BWIA) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "BWIA using the complete log2 standard binning vector")

# Define useful (linear) size range for linear model
# (from maximum to first empty bin
# Corresponds to the size range, where sampliinng is really reprentative of the communty size spectrum
# = size range of quantitative sampling (!)

points(log10(BWIA[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

summary(lm(log10(BWIA[24:30]) ~ log10(hpw2c$mids[24:30])))

# Numbers to Abundance (divide by S)
# sampled volume "S": 100 cubic meters (hypothetical)
# A = N / S
# binned (not normalized) Abundance = BWIA
# Bin-width-inflated Abundance (BWIA)

#BWIA <- hpw2c$counts /100


# NORMALIZE ------------
# from BWIA to NNSS

# Bin widths ("w")

w2c <- diff(hpw2c$breaks)

# NNSS (divide Ni by wi) , normalized Abundance ------
# Normalized numbers size spectrum


NNSS <- BWIA / w2c

plot(NNSS ~ hpw2c$mids)
plot(log10(NNSS) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "NNSS using the log2 standard binning vector")

points(log10(NNSS[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1A <- lm(log10(NNSS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1A)
# slope : b = -2, perfectly OK for Abundance vs body mass 

# III. bNNS, the "backtransformed” NNSS ----------------

# BACKTRANSFORM with D -----------
# fom NNSS to bNNS
# After normalizing (dividing BWIA by w, where NNSS = BWIA / w) to obain the common NBSS, rescaling is needed. This can be done by using the 
# w-scaled dimensional correction factor D, where
# 
# D =  BWIA_total  / NNSS_total  
# 
# All NNSSi  are then transformed back into the scale of A,
# i.e., into a "backtransformed” NNSS (bNNS), where

# calculate D (dimensional correction factor D) -------

# D =  BWIA_total  / NNSS_total
BWIA_total = sum(BWIA) # 70.88 = total Abundance, OK!
NNSS_total = sum(NNSS)

D_Abund = BWIA_total / NNSS_total

   bNN = NNSS * D_Abund

bNN 

# bNN: backtransformed normalized numbers

sum(bNN) # 70.88 = total Abundance, OK!



plot(bNN ~ hpw2c$mids)
plot(log10(bNN) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "bNNS using the log2 standard binning vector",
     ylab ="Abundance (bNN), log10(ind. m-3)" ,
     xlab = "Body mass, log10(g ind.-1)") 


points(log10(bNN[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1 <- lm(log10(bNN[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2.01, perfectly OK for Abundance vs body mass 
# expected: -2,  for Abundance vs body mass 


# IV. NBSS, Biomass (NBSS: normalized biomass size spectrum) -------------------------

#  Biomass is directly binned (not through abundance) ---------------

# non-linear binning, log2 ----------------

binnedB2 <- cut(data.carbon_g_ind_1, breaks = log2vec, include.lowest = TRUE)

h4 <- hist(data.carbon_g_ind_1, breaks = log2vec, plot = F) # to obtain mids, etc
mids2 <- h4$mids

biomass_per_bin2 <- tapply(data.carbon_g_ind_1, binnedB2, sum)
biomass_per_bin2

# Biomass per bin, divied by sampled volume (100 m3 in this example)
biomass_per_bin2b <- biomass_per_bin2/100


B_total =sum(biomass_per_bin2, na.rm = T) / 100 # g m-3 Total Biomass,   OK!
B_total
# 0.00342 total Biomass (g.  m-3), OK

sum(biomass_per_bin2b,na.rm = T)
# 0.00342 total Biomass (g.  m-3), OK

biomass_per_bin2 <- biomass_per_bin2b


# Biomass (g.  m-3)
B_total = sum(data.carbon_g_ind_1) / 100 # B =  sum(M)/S, g. per cubic meter
B_total # 0.00342 total Biomass (g.  m-3)


plot(biomass_per_bin2b ~ mids2)
plot(log(biomass_per_bin2b) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
#abline ( h = -8)
# biomass_per_bin2 ins BWIB! flat dostributin because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[24:30]) ~ log(mids2[24:30]), col = "navyblue", pch = 16)

lm.BWIB <- lm(log(biomass_per_bin2[24:30]) ~ log(mids2[24:30]))
summary(lm.BWIB) 
abline(lm.BWIB)
# nearly flat BWIB- body mass spectrum , b = 0.04 !

# "For a common power-law distribution, the use of 
# geometrically increasing bin sizes will lead to a flat (slope = 0) , 
# unrealistic distribution of the BWIB,
# while the NBSS has a realistic power-law shape, with a linearly 
# downtrending  slope on a  log-log-scaled  plot, that represents well the original 
# biomass-body-mass relationship of the natural ecosystem." Schwamborn, 2026 

# biomass_per_bin2 is BWIB! flat distributin because of distorting effect of non-linear binning!
#  BWIB = Bin-width inflated Biomass


# calculate the common NBSS (normalized biomass size specptrum)

w2b <- diff(log2vec)

NBSS <- as.vector(  biomass_per_bin2/w2b)
sum(NBSS, na.rm = T) # 106.32, OK

plot(NBSS ~ mids2)
plot(log(NBSS) ~ log(mids2)) 

points(log(NBSS[24:30]) ~ log(mids2[24:30]), col = "navyblue", pch = 16)

lm1b <- lm(log(NBSS[24:30]) ~ log(mids2[24:30]))
summary(lm1b)
abline(lm1b)
        
# slope = -1.041, approx b = -1, OK!

# V. bNBS (new method) "back-transformed normalized biomass" ----------------------

# Example calculation of bNBS
# bNB = "back-transformed normalized biomass", name of variable and numerical biomass data
# bNBS = "back-transformed normalized biomass spectrum", name of plot, model, and method

# BWIB: “bin-with-inflated biomass” 
# The binning process distorts and transforms the original biomass into the binned biomass, or  “bin-with-inflated biomass” (BWIB).
# After normalizing (dividing BWIBi by wi, where NBSSi = BWIBi / wi) to obain the common NBSS, rescaling is needed. This can be done by using the w-scaled dimensional correction factor D, where
# 
# D =  BWIBtotal  / NBSStotal  
# 
# All NBSS data  are then transformed back (re-dimensionalized back and rescaled back)  into the original scale, dimension, and units of biomass, i.e., into  “back-transformed normalized biomass” (bNB), where


# bNBi = NBSSi * D

# Since D has units of B * M (e.g.,  g² m-3 ind.-1), this calculation backtransforms the units of the NBSS into the original biomass units, within the  bNBS (e.g., g m-3). 


# D =  BWIB_total  / NBSS_total

BWIB <- biomass_per_bin2

BWIB_total = sum(BWIB,na.rm = T) # 0.003427 = total Biomass, OK!
BWIB_total

NBSS_total = sum(NBSS,na.rm = T)
NBSS_total

D = BWIB_total / NBSS_total
D

bNB = NBSS * D

bNB

sum(bNB, na.rm = T) # 0.00342 g m-3 = total Biomass, OK!

# Makke a plot for bnBS (the bNBS plot)


plot(log10(bNB) ~ log10(mids2), 
     main = "bNBS plot",
     ylab ="Biomass (bNB), log10(g m-3)" ,
     xlab = "Body mass, log10(g ind.-1)") 


points(log10(bNB[24:30]) ~ log10(mids2[24:30]), col = "navyblue", pch = 16)

lm1g <- lm(log10(bNB[24:30]) ~ log10(mids2[24:30]))
summary(lm1g)
abline(lm1g)

# Intercept =      -7.95
# slope = -1.041, approx b = -1, OK!

sum(bNB, na.rm = T) # 0.00342 g m-3 = total Biomass, OK!



#### End of Script ####


