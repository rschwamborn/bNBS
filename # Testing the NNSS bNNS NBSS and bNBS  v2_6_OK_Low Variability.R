#  Testing the NNSS, bNNSS, NBSS, and bNBSS size spectra
# Version   2.5
# Low Variability (perfectly linear spectra)


# NNSS: normalized numbers size spectrum (Abundance vs body mass, normalized, but not rescaled)
# bNNSS (new!): backtransformed and normalized numbers size spectrum (Abundance vs body mass, normalized, AND rescaled)
# NBSS: normalized biomass size spectrum (Biomass vs body mass, normalized, but not rescaled)
# bNBSS (new!): backtransformed  and normalized numbers size spectrum (Biomass vs body mass, normalized, AND rescaled)


# Testing and comparing 

# original distribution (synthetic power law distribution)
# vs 
# binned biomass distribution (binning from individual body mass data) 
# vs
# binned abundance 
# vs
#  biomass distribution from binned abundance (B = A * M)


# R script for the paper
# "Is normalized biomass really abundance? 
# Ralf Schwamborn 
# January 2026
# Source : "https://raw.githubusercontent.com/rschwamborn/bNBS


# I.1. Uniform distribution
# I.2. Linearly decreasing distribution
# I.3. Power-law Abundance-body-mass distribution
# I.4. Power-law Biomass-body-mass calculated from Abundance distribution
# I.5. Calculate  Biomass per bin directly, without first calculating Abundance
# I.6. Compare input vs output distribution results  
# I.7 As function and loop ... bias tests, input vales vs output, with bias
# 1.8 Test the bNBSS Backtransformed (backtransform = rescaling and redimensionalizing
# I.8.A. k-normalization and rescaling - the kNBS
# I.8.B. w-normalization and F-rescaling of the NBSS, the bNBSS 
# II. Example, using a standardized log2 binning vector - example calculation
# III. Import example in situ data
# IV. Comparison of variability
# V. bNBS (new method)


# clean memory (caution...)-----------
# rm(list = ls())
# gc()


#I. 1. Uniform body_mass distribution ---------

set.seed(1)

# Individual body_masss
n <- 10000
body_mass.u <- runif(n)


h1.u <- hist(body_mass.u)


h1.u <- hist(
  body_mass.u,
  freq = F,
  plot = TRUE,
  col = "lightgray",
  main = "y: Abund., Uniform distrib., linear binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

h1.u$counts

plot(h1.u$counts ~ h1.u$mids)

sum(h1.u$counts) ## sum = 10,000, OK


# Create log2 breaks

x<- c(0.1,2000)

breaks <- 2^(floor(log2(min(x))) : ceiling(log2(max(x))))


non_lin_bins <- b2 <- breaks /2000


h1.u2 <- hist(
  body_mass.u,
  freq = F,
  breaks = non_lin_bins,
  plot = TRUE,
  col = "lightgray",
  main = "y: Abund., Uniform distrib., log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

h1.u2$counts

sum(h1.u2$counts) ## sum = 10,000, OK

plot(h1.u2$counts ~ h1.u2$mids)



#I. 2. Linearly decreasing body mass distribution -----------

#We’ll make many small individuals and few large ones.

set.seed(1)

# Individual body_masss
n <- 1000
body_mass <- runif(n)

# Linearly decreasing probability (more small than large)
prob <- 1 - body_mass
keep <- runif(n) < prob
body_mass <- body_mass[keep]
length(body_mass) # N = 513

hist(body_mass)

#2. Histogram where area of each bar = abundance

#This uses density = FALSE and equal-width bins, so bar area = count.

bins <- seq(0, 1, by = 0.1)

h <- hist(
  body_mass,
  breaks = bins,
  plot = TRUE,
  col = "lightgray",
  main = "Histogram: y: Abund., area of bars = abundance* bin width",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

h$counts
sum(h$counts) # # N = 513, OK

# now use a non-linear binning vector

# Create log2 breaks

x<- c(0.1,2000)

breaks <- 2^(floor(log2(min(x))) : ceiling(log2(max(x))))


non_lin_bins <- b2 <- breaks /2000


h2 <- hist(
  body_mass,
  breaks = non_lin_bins,
  plot = TRUE,
  col = "lightgray",
  main = "Histogram: y: Abund., area of bars = abundance* bin width",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

sum(h2$counts) # # N = 513, OK

# non-linear binning does not change the total! Same scale! Same parameter "N"!


#I. 3. Power-law Abundance - body_mass distributions -----------

set.seed(123)
n <- 100000
alpha <- 2
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)

body_mass.pw <- x 

# x is in (0,1)

summary(x)

# Linear binning -----------

hpw <- hist(
  body_mass.pw,
  plot = TRUE,
  col = "lightgray",
  main = "Histogram: y: Abund., area of bars = abundance* bin width",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw$breaks # the linear binning vector
w <- diff(hpw$breaks) # constant bin with (w) of the linear binning vector
hpw$mids # midpoints (mp)  of the linear binning vector
range(hpw$breaks)# min, max of the binning vector
R <- diff(range(hpw$breaks))# Range (R) of the binning vector

k <- w/R # unitless and dimensionless binning factor "k"
k
sum(k) # sum = 1, OK

plot( hpw$counts ~  hpw$mids)

plot( log(hpw$counts) ~  log(hpw$mids))

hpw$counts

# [1] 99830    84    30    18     9     6     4     4     1     3     2     3     3     1
# [15]     1     0     0     0     1
# Many zeros (empty bins)!

plot( log(hpw$counts) ~  log(hpw$mids)  ,col = "pink",
      main =       "Power-law, linear binning, Abundance")

# Nicely lnear log-log plot, except for the zeros (!) and extremes!

points(log(hpw$counts[2:7]) ~  log(hpw$mids[2:7]), col = "blue" )

lm1 <- lm( log(hpw$counts[2:7]) ~  log(hpw$mids[2:7]) )

abline(lm1)
summary(lm1)

# Intercept a =          -0.87
# slope     b =         -2.06  (approx. b = -2) 

text(x = -1.5, y = 8, " slope b =  -2.06, Intercept a =  -0.87" )

# non-linear binning (log2) ---------------

# BWIA: bin-width-inflated Abundance ----------------
# BWIA suffers the non-linear binning distortion effect -------------
# The slope of the Abundance-body mass spectrum shifts from b = -2 (OK) to  b= -1 (distorted by BWIA)

# Create log2 breaks ------------

x<- c(0.1,2000)

breaks2 <- 2^(floor(log2(min(x))) : ceiling(log2(max(x))))

non_lin_breaks <- breaks2 /2000


hpw2 <- hist(
  body_mass.pw,
  breaks = non_lin_breaks,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw2$breaks # the non-linear binning vector (log2)
w2 <- diff(hpw2$breaks) #  bin withs (w) of the linear binning vector
hpw2$mids # midpoints (mp)  of the linear binning vector
range(hpw2$breaks)# min, max of the binning vector
R2 <- diff(range(hpw2$breaks))# Range (R) of the binning vector

k2 <- w2/R2 # unitless and dimensionless binning factor "k"
k2
sum(k2) # sum = 1, OK

plot( hpw2$counts ~  hpw2$mids)

plot( log(hpw2$counts) ~  log(hpw2$mids))

hpw2$counts
# [1]     0 20033 40105 19939  9993  4930  2502  1231   656   314   160    69    40    17
# [15]    11
# Many zeros (empty bins)!

plot( log(hpw2$counts) ~  log(hpw2$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, Abundance, BWIA")


# Nicely linear log-log plot, except for the zeros (!) and extremes!

points(log(hpw2$counts[3:14]) ~  log(hpw2$mids[3:14]), col = "blue" )

lm1 <- lm( log(hpw2$counts[3:14]) ~  log(hpw2$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         3.8
# slope     b =  -1  ( b was  -2, but now shows the non-linear binning distortion effect) 

text(x = -7, y = 4.5, "b = -1" )
text(x = -6, y = 3, "b was -2, now shows the non-linear binning distortion effect" )
text(x = -6, y = 1.5, "BWIA suffers the non-linear binning distortion effect" )


# normalized abundance, NNSS (divide by w) --------------------
# NNSS: normalized numbers size spectrum  -----------------------

NNSS <- hpw2$counts/w2

plot( log(NNSS) ~  log(hpw2$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(NNSS[3:14]) ~  log(hpw2$mids[3:14]), col = "blue" )

lm1 <- lm( log(NNSS[3:14]) ~  log(hpw2$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         3.8
# slope     b =  -2, NNSS, OK  

text(x = -7, y = 4.5, "b = -2, OK" )


# k-normalized abundance, knNSS (divide by k) --------------------
# knNSS: k-normalized numbers size spectrum ------------------

knNSS <- hpw2$counts/k2

plot( log(knNSS) ~  log(hpw2$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(knNSS[3:14]) ~  log(hpw2$mids[3:14]), col = "blue" )

lm1 <- lm( log(knNSS[3:14]) ~  log(hpw2$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         2.4
# slope     b =  -2, knNSS, OK  

text(x = -7, y = 4.5, "b = -2, OK" )



#I.4. Power-law Biomass - body-mass distributions -----------

# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

sum(body_mass.pw) # 90.67 g m-3  = Total Biomass, Btotal

# calculate Biomass for each bin, from  Abundance A and body mass M

# linear binning ----------

mids1 <- hpw$mids # midpoints (mp)  of the linear binning vector
# total abundance
sum(hpw$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal <- sum(hpw$counts * hpw$mids) # Biomass = 2,524.25 g m-3
B <- hpw$counts * hpw$mids
sum(B)

NBSS <- B/w
sum(NBSS)

plot(B ~ mids1)
plot(log(B) ~ log(mids1))
points(log(B[2:7]) ~ log(mids1[2:7]), col = "navyblue", pch = 16)

lm1 <- lm(log(B[2:7]) ~ log(mids1[2:7]))
summary(lm1) # slope = - 1, original biomass-body mass slope, OK
# intercept = - 0.87, OK
abline (lm1)


# non-linear binning ----------

mids2 <- hpw2$mids # midpoints (mp)  of the non-linear binning vector (log2)
# total abundance
sum(hpw2$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal2 <- sum(hpw2$counts * hpw2$mids) # Biomass =  99.02 g m-3 !!
B2 <- hpw2$counts * hpw2$mids
NBSS2<- B2/w2
sum(NBSS2) # 150000

# I.5. Calculate  biomass per bin directly, without first calculating abundance ---------

# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

# linear binning --------------

binnedB <- cut(body_mass.pw, breaks = hpw$breaks, include.lowest = TRUE)

biomass_per_bin <- tapply(body_mass.pw, binnedB, sum)
biomass_per_bin
sum(biomass_per_bin, na.rm = T) # 90.7 g m-3 Total Biomass


plot(biomass_per_bin ~ mids1)
plot(log(biomass_per_bin) ~ log(mids1))
points(log(biomass_per_bin[2:7]) ~ log(mids1[2:7]), col = "navyblue", pch = 16)

lm1 <- lm(log(biomass_per_bin[2:7]) ~ log(mids1[2:7]))
summary(lm1) 
# slope = -1, original biomass-body mass slope, OK
# intercept = -0.8, OK
abline (lm1)

# non-linear binning ----------------

binnedB2 <- cut(body_mass.pw, breaks = non_lin_breaks, include.lowest = TRUE)

biomass_per_bin2 <- tapply(body_mass.pw, binnedB2, sum)
biomass_per_bin2
sum(biomass_per_bin, na.rm = T) # 90.7 g m-3 Total Biomass,   OK!

plot(biomass_per_bin2 ~ mids2)
plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!

# biomass_per_bin2 ins BWIB! flat dostributin because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)

lm.BWIB <- lm(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]))
summary(lm.BWIB) 
# flat BWIB- body mass spectrum , b = 0 !

# "For a common power-law distribution, the use of 
# geometrically increasing bin sizes will lead to a flat (slope = 0) , 
# unrealistic distribution of the BWIB,
# while the NBSS has a realistic power-law shape, with a linearly 
# downtrending  slope on a  log-log-scaled  plot, that represents well the original 
# biomass-body-mass relationship of the natural ecosystem." Schwamborn, 2026 

# calculate the common NBSS (normalized biomass size specptrum)

NBSS <- biomass_per_bin2/w2
sum(NBSS, na.rm = T)

plot(biomass_per_bin2 ~ mids2)
plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!

# biomass_per_bin2 is BWIB! flat distributin because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)

lm.BWIB <- lm(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]))
summary(lm.BWIB) 
# flat BWIB- body mass spectrum , b = 0 (OK)




# I.6. Compare input vs output distribution results  ------------------------

# Input: power law data  
# output: plots and linear models (a,b, R²) for A (NNSS), B from A, B ,
#          coefficient of variation of residuals , CVR, % 


# How to measure of Variability in the spectrum --------

# Additionally to comparing linear model slopes (“b” values), 
# I also compared the representation of the variability of the data,
# around the linear model. To verify whether the variability
# in the original data is consistently represented in the output data
# models, a non-dimensional coefficient of variation of residuals 
# (CVR, %) was calculated for BWIB, NNSS, and NBSS, for all simulations.
# CVR (%) was calculated as the standard error of the residuals (SER) of the log-log-linear model 
# (fitted with linear regression), divided by the mean “y” value, 
# and multiplied by 100:
#   CVR(%) = (SER / mean) * 100.

# general formula:
# coefficient of variation of residuals , CVR, % :
# CVR <- abs(  sd( resid(lm1) ) / mean(lm1$fitted.values) * 100 )



# Input: power law data  -------------

set.seed(123)
n <- 100000

# Standard Power-Law distribution (NNNS = -2, NBSS = -1): -------------
alpha <- 2
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)


# Standard Power Law distribution: -------------

x12  <- body_mass.pw <- x # biomass b = -1, abundance = b = -2


# Flat Power Law distribution: -------------
# biomass spectrum b = -05, abundance spectrum b = -1.5

alpha <- 1.5
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)

# Flat Power Law distribution:
x05_15  <- body_mass.pw05_15 <- x # biomass b = -0.5, abundance = b = -1.5

# Steep Power Law distribution: -------------

alpha <- 2.5
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)

# Steep Power Law distribution:
x15_25  <- body_mass.pw15_25 <- x 
# biomass b = -1.5, abundance = b = -2.5



# output: plots and linear models (a,b, R²) for A (NNSS), BfromA, B  ------------

# 1.6 A. FLAT spectrum ---------------------------

# Linear binning -----------

hpw05_15 <- hist(
  body_mass.pw05_15,
  plot = TRUE,
  col = "lightgray",
  main = "flat spectrum, b_abund = -1.5",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw05_15$breaks # the linear binning vector
w <- diff(hpw05_15$breaks) # constant bin with (w) of the linear binning vector
hpw05_15$mids # midpoints (mp)  of the linear binning vector
range(hpw05_15$breaks)# min, max of the binning vector
R <- diff(range(hpw05_15$breaks))# Range (R) of the binning vector

k <- w/R # unitless and dimensionless binning factor "k"
k
sum(k) # sum = 1, OK

plot( hpw05_15$counts ~  hpw05_15$mids)

plot( log(hpw05_15$counts) ~  log(hpw05_15$mids))

hpw05_15$counts

# [1] 99830    84    30    18     9     6     4     4     1     3     2     3     3     1
# [15]     1     0     0     0     1
# Many zeros (empty bins)!

plot( log(hpw05_15$counts) ~  log(hpw05_15$mids)  ,col = "pink",
      main =       "Power-law, linear binning, b = -1.5, Abundance")

# Nicely linear log-log plot, except for the zeros (!) and extremes!

points(log(hpw05_15$counts[2:7]) ~  log(hpw05_15$mids[2:7]), col = "blue" )

lm1 <- lm( log(hpw05_15$counts[2:7]) ~  log(hpw05_15$mids[2:7]) )

abline(lm1)
summary(lm1)


# Intercept a =         3.33
# slope     b =         -1.46  (approx. b = -1.5 ,  OK) 

text(x = -1.5, y = 8, " slope b =  -1.45, Intercept a = 3.33, OK" )

# non-linear binning (log2) ---------------

# BWIA: bin-width-inflated Abundance ----------------
# BWIA suffers the non-linear binning distortion effect -------------
# The slope of the Abundance-body mass spectrum shifts from b = -2 (OK) to  b= -1 (distorted by BWIA)

# Create log2 breaks ------------

x<- c(0.1,2000)

breaks2 <- 2^(floor(log2(min(x))) : ceiling(log2(max(x))))

non_lin_breaks <- breaks2 /2000


hpw05_152 <- hist(
  body_mass.pw05_15,
  breaks = non_lin_breaks,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw05_152$breaks # the non-linear binning vector (log2)
w2 <- diff(hpw05_152$breaks) #  bin withs (w) of the linear binning vector
hpw05_152$mids # midpoints (mp)  of the linear binning vector
range(hpw05_152$breaks)# min, max of the binning vector
R2 <- diff(range(hpw05_152$breaks))# Range (R) of the binning vector

k2 <- w2/R2 # unitless and dimensionless binning factor "k"
k2
sum(k2) # sum = 1, OK

plot( hpw05_152$counts ~  hpw05_152$mids)

plot( log(hpw05_152$counts) ~  log(hpw05_152$mids))

hpw05_152$counts
# [1]     0 20033 40105 19939  9993  4930  2502  1231   656   314   160    69    40    17
# [15]    11
# Many zeros (empty bins)!

plot( log(hpw05_152$counts) ~  log(hpw05_152$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, Abundance, BWIA,    OK")


# Nicely linear log-log plot, except for the zeros (!) and extremes!

points(log(hpw05_152$counts[3:14]) ~  log(hpw05_152$mids[3:14]), col = "blue" )

lm1 <- lm( log(hpw05_152$counts[3:14]) ~  log(hpw05_152$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         5.9
# slope     b =  -0.5  ( b was  -2, but now shows the non-linear binning distortion effect) 

 text(x = -7, y = 8.5, "b = -0.5,  OK" )
 text(x = -6, y = 7, "b was -1.5, now shows the non-linear binning distortion effect" )
 text(x = -6, y = 6.5, "BWIA suffers the non-linear binning distortion effect" )

# normalized abundance, NNSS (divide by w) --------------------
# NNSS: normalized numbers size spectrum  -----------------------

NNSS <- hpw05_152$counts/w2

plot( log(NNSS) ~  log(hpw05_152$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(NNSS[3:14]) ~  log(hpw05_152$mids[3:14]), col = "blue" )

lm1 <- lm( log(NNSS[3:14]) ~  log(hpw05_152$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         6.3
# slope     b =       -1.5, NNSS, OK, perfect slope , ZERO BIAS!!!  

text(x = -7, y = 4.5, "b = -1.5, perfect slope, zero bias, OK" )

# calculate percent bias:
x1<- -1.4966
x2<- -1.5
percent_bias <- abs(x1 - x2) / x1  * 100
percent_bias # 0.22 percent
# percent bias in slope "b" is less than 0.3%


# k-normalized abundance, knNSS (divide by k) --------------------
# knNSS: k-normalized numbers size spectrum ------------------

knNSS <- hpw05_152$counts/k2

plot( log(knNSS) ~  log(hpw05_152$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(knNSS[3:14]) ~  log(hpw05_152$mids[3:14]), col = "blue" )

lm1 <- lm( log(knNSS[3:14]) ~  log(hpw05_152$mids[3:14]) )

abline(lm1)
summary(lm1)

# Intercept a =         6.3
# slope     b =  -1.5, knNSS, OK  

text(x = -7, y = 4.5, "b = 1.5, OK" )


# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

sum(body_mass.pw05_15) # 997.3 g m-3  = Total Biomass, Btotal

# calculate Biomass for each bin, from  Abundance A and body mass M

# linear binning ----------

mids1 <- hpw05_15$mids # midpoints (mp)  of the linear binning vector
# total abundance
sum(hpw05_15$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal <- sum(hpw05_15$counts * hpw05_15$mids) # Biomass = 2,524.25 g m-3
B <- hpw05_15$counts * hpw05_15$mids
sum(B)

NBSS <- B/w
sum(NBSS)

plot(B ~ mids1)
plot(log(B) ~ log(mids1))
points(log(B[2:7]) ~ log(mids1[2:7]), col = "navyblue", pch = 16)

lm1 <- lm(log(B[2:7]) ~ log(mids1[2:7]))
summary(lm1) 
# slope b = - 0.46, OK, similar to -0.5 original biomass-body mass slope, OK
# intercept = 3.33
abline (lm1)


# non-linear binning ----------

mids2 <- hpw05_152$mids # midpoints (mp)  of the non-linear binning vector (log2)
# total abundance
sum(hpw05_152$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal2 <- sum(hpw05_152$counts * hpw05_152$mids) # Biomass =  99.02 g m-3 !!
B2 <- hpw05_152$counts * hpw05_152$mids
NBSS2<- B2/w2
sum(NBSS2) # 150000

plot(log(NBSS2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
points(log(NBSS2[3:13]) ~ log(mids2[3:13]), col = "navyblue", pch = 16)

lm1 <- lm(log(NBSS2[3:13]) ~ log(mids2[3:13]))
summary (lm1) # b = -0.5, OK, excellent fit! ZERO bias!!!
text( x = -7,y = 7,"b = -0.5, perfect fit, zero bias!" )
abline(lm1)


# calculate percent bias:
x1<- -0.5
x2<- -0.50365
percent_bias <- abs(x1 - x2) / x1  * 100
percent_bias # 0.73 percent
# percent bias in slope "b" is less than 0.8%



# I.5. Calculate  biomass per bin directly, without first calculating abundance ---------

# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

# linear binning --------------

binnedB <- cut(body_mass.pw05_15, breaks = hpw05_15$breaks, include.lowest = TRUE)

biomass_per_bin <- tapply(body_mass.pw05_15, binnedB, sum)
biomass_per_bin
sum(biomass_per_bin, na.rm = T) # 90.7 g m-3 Total Biomass

plot(biomass_per_bin ~ mids1)
plot(log(biomass_per_bin) ~ log(mids1))
points(log(biomass_per_bin[2:7]) ~ log(mids1[2:7]), col = "navyblue", pch = 16)

lm1 <- lm(log(biomass_per_bin[2:7]) ~ log(mids1[2:7]))
summary(lm1) 
# slope = -0.42, original biomass-body mass slope = -0.5, +_OK, slight bias...
# intercept = 3.38, OK
abline (lm1)

# non-linear binning ----------------

binnedB2 <- cut(body_mass.pw05_15, breaks = non_lin_breaks, include.lowest = TRUE)

biomass_per_bin2 <- tapply(body_mass.pw05_15, binnedB2, sum)
biomass_per_bin2
sum(biomass_per_bin, na.rm = T) # 997.2 g m-3 Total Biomass,   OK!

plot(biomass_per_bin2 ~ mids2)
plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!

# biomass_per_bin2 ins BWIB! 
# flat B distributin gives POSITIVE BWIB distribution (!) because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)

NBSS2b<- biomass_per_bin2/w2
sum(NBSS2b) # 150000

plot(log(NBSS2b) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
points(log(NBSS2b[3:13]) ~ log(mids2[3:13]), col = "navyblue", pch = 16)

lm2b <- lm(log(NBSS2b[3:13]) ~ log(mids2[3:13]))
summary (lm2b) # b = -0.5, OK, excellent fit! ZERO bias!!!
text( x = -7,y = 7,"b = -0.5, perfect fit, zero bias!" )
abline(lm2b)


# 1.6 B. Steep spectrum ---------------------------

# Linear binning -----------

hpw15_25 <- hist(
  body_mass.pw15_25,
  plot = TRUE,
  col = "lightgray",
  main = "steep spectrum, b_abund = -2.5",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw15_25$breaks # the linear binning vector
w <- diff(hpw15_25$breaks) # constant bin with (w) of the linear binning vector
hpw15_25$mids # midpoints (mp)  of the linear binning vector
range(hpw15_25$breaks)# min, max of the binning vector
R <- diff(range(hpw15_25$breaks))# Range (R) of the binning vector

k <- w/R # unitless and dimensionless binning factor "k"
k
sum(k) # sum = 1, OK

plot( hpw15_25$counts ~  hpw15_25$mids)

plot( log(hpw15_25$counts) ~  log(hpw15_25$mids))

hpw15_25$counts

# [1] 99893    68    21    10     2     0     0     2     0     0     1     0     0     2
# [15]     0     0     0     1
#   Many zeros (empty bins)!

plot( log(hpw15_25$counts) ~  log(hpw15_25$mids)  ,col = "red",
      main =       "Power-law, linear binning, b = -2.5, Abundance")

# Nicely linear log-log plot, except for the zeros (!) and extremes!

points(log(hpw15_25$counts[2:4]) ~  log(hpw15_25$mids[2:4]), col = "blue" )

lm1 <- lm( log(hpw15_25$counts[2:4]) ~  log(hpw15_25$mids[2:4]) )

abline(lm1)
summary(lm1)


# Intercept a =         3.33
# slope     b =         -2.2  (approx. b = -2.5 ,  +_ OK) 



# non-linear binning (log2) ---------------

# BWIA: bin-width-inflated Abundance ----------------
# BWIA suffers the non-linear binning distortion effect -------------
# The slope of the Abundance-body mass spectrum shifts from b = -2 (OK) to  b= -1 (distorted by BWIA)

# Create log2 breaks ------------

x<- c(0.1,2000)

breaks2 <- 2^(floor(log2(min(x))) : ceiling(log2(max(x))))

non_lin_breaks <- breaks2 /2000


hpw15_252 <- hist(
  body_mass.pw15_25,
  breaks = non_lin_breaks,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

# Calculating "k"  --------------------

hpw15_252$breaks # the non-linear binning vector (log2)
w2 <- diff(hpw15_252$breaks) #  bin withs (w) of the linear binning vector
hpw15_252$mids # midpoints (mp)  of the linear binning vector
range(hpw15_252$breaks)# min, max of the binning vector
R2 <- diff(range(hpw15_252$breaks))# Range (R) of the binning vector

k2 <- w2/R2 # unitless and dimensionless binning factor "k"
k2
sum(k2) # sum = 1, OK

plot( hpw15_252$counts ~  hpw15_252$mids)

plot( log(hpw15_252$counts) ~  log(hpw15_252$mids))

hpw15_252$counts
# [1]     0 20033 40105 19939  9993  4930  2502  1231   656   314   160    69    40    17
# [15]    11
# Many zeros (empty bins)!

plot( log(hpw15_252$counts) ~  log(hpw15_252$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, Abundance, BWIA,    OK")


# Nicely linear log-log plot, except for the zeros (!) and extremes!

points(log(hpw15_252$counts[3:11]) ~  log(hpw15_252$mids[3:11]), col = "blue" )

lm1 <- lm( log(hpw15_252$counts[3:11]) ~  log(hpw15_252$mids[3:11]) )

abline(lm1)
summary(lm1)

# Intercept a =         2.3
# slope     b =  -1.5 
# ( b was  -2, but now shows the non-linear binning distortion effect) 

text(x = -7, y = 8.5, "b = -1.5,  OK" )
text(x = -6, y = 7, "b was -2.5, now shows the non-linear binning distortion effect" )
text(x = -6, y = 6.5, "BWIA suffers the non-linear binning distortion effect" )

# normalized abundance, NNSS (divide by w) --------------------
# NNSS: normalized numbers size spectrum  -----------------------

NNSS <- hpw15_252$counts/w2

plot( log(NNSS) ~  log(hpw15_252$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(NNSS[3:11]) ~  log(hpw15_252$mids[3:11]), col = "blue" )

lm1 <- lm( log(NNSS[3:11]) ~  log(hpw15_252$mids[3:11]) )

abline(lm1)
summary(lm1)

# Intercept a =         -1.9
# slope     b =       -2.5, NNSS, OK, perfect slope , ZERO BIAS!!!  

text(x = -7, y = 4.5, "b = -2.5, perfect slope, zero bias, OK" )

+ (-2.5 - -2.5304)/ (-2.5304)

x1 <- -2.5
x2 <-  -2.5304

percent_bias <- abs(x1 - x2) / x1  * 100
percent_bias #  1.2 percent
# percent bias in slope "b" is less than 1.3%


# k-normalized abundance, knNSS (divide by k) --------------------
# knNSS: k-normalized numbers size spectrum ------------------

knNSS <- hpw15_252$counts/k2

plot( log(knNSS) ~  log(hpw15_252$mids)  ,col = "pink",
      main =       "Power-law, log2-binning, norm. Abundance, NNSS")

# Nicely linear log-log plot, except for the  lower extreme

points(log(knNSS[3:11]) ~  log(hpw15_252$mids[3:11]), col = "blue" )

lm1 <- lm( log(knNSS[3:11]) ~  log(hpw15_252$mids[3:11]) )

abline(lm1)
summary(lm1)

# Intercept a =         -1.9
# slope     b =  -2.53, knNSS, OK  

text(x = -7, y = 4.5, "b = -2.5, OK" )


# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

sum(body_mass.pw15_25) #  29.6  g m-3  = Total Biomass, Btotal

# calculate Biomass for each bin, from  Abundance A and body mass M

# linear binning ----------

mids1 <- hpw15_25$mids # midpoints (mp)  of the linear binning vector
# total abundance
sum(hpw15_25$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal <- sum(hpw15_25$counts * hpw15_25$mids) # Biomass = 2,524.25 g m-3
B <- hpw15_25$counts * hpw15_25$mids
sum(B)

NBSS <- B/w
sum(NBSS)

plot(B ~ mids1)
plot(log(B) ~ log(mids1))
points(log(B[2:4]) ~ log(mids1[2:4]), col = "navyblue", pch = 16)

lm1 <- lm(log(B[2:4]) ~ log(mids1[2:4]))
summary(lm1) 
# slope b = -1.26, OK, similar to -1.5 original biomass-body mass slope, OK

abline (lm1)


# non-linear binning ----------

mids2 <- hpw15_252$mids # midpoints (mp)  of the non-linear binning vector (log2)
# total abundance
sum(hpw15_252$counts) # N = 100,000 OK
# total biomass B = A * M
Btotal2 <- sum(hpw15_252$counts * hpw15_252$mids) # Biomass =  99.02 g m-3 !!
B2 <- hpw15_252$counts * hpw15_252$mids
NBSS2<- B2/w2
sum(NBSS2) # 150000

plot(log(NBSS2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
points(log(NBSS2[3:13]) ~ log(mids2[3:13]), col = "navyblue", pch = 16)

lm1 <- lm(log(NBSS2[3:13]) ~ log(mids2[3:13]))
summary (lm1) # b = -0.5, OK, excellent fit! ZERO bias!!!
text( x = -7,y = 7,"b = -0.5, perfect fit, zero bias!" )
abline(lm1)




# I.5. Calculate  biomass per bin directly, without first calculating abundance ---------

# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

# linear binning --------------

binnedB <- cut(body_mass.pw15_25, breaks = hpw15_25$breaks, include.lowest = TRUE)

biomass_per_bin <- tapply(body_mass.pw15_25, binnedB, sum)
biomass_per_bin
sum(biomass_per_bin, na.rm = T) # 29.6 g m-3 Total Biomass

plot(biomass_per_bin ~ mids1)
plot(log(biomass_per_bin) ~ log(mids1))
points(log(biomass_per_bin[2:4]) ~ log(mids1[2:4]), col = "navyblue", pch = 16)
# ISSUE: only 3 useful bins when using linear binning ... needs non-linear binning.. 


lm1 <- lm(log(biomass_per_bin[2:4]) ~ log(mids1[2:4]))
summary(lm1) 
# slope = -1.16, original biomass-body mass slope = -1.5, +_OK, slight bias... 
# ISSUE: only 3 useful bins in linear binning ... needs non-linear binning.. 
abline (lm1)

# non-linear binning ----------------

binnedB2 <- cut(body_mass.pw15_25, breaks = non_lin_breaks, include.lowest = TRUE)

biomass_per_bin2 <- tapply(body_mass.pw15_25, binnedB2, sum)
biomass_per_bin2
sum(biomass_per_bin, na.rm = T) # 29.6 g m-3 Total Biomass,   OK!

plot(biomass_per_bin2 ~ mids2)
plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!

# biomass_per_bin2 ins BWIB! 
# steep B distributin gives POSITIVE BWIB distribution (!) because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)

NBSS2b<- biomass_per_bin2/w2
sum(NBSS2b) # 150000

plot(log(NBSS2b) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
points(log(NBSS2b[3:11]) ~ log(mids2[3:11]), col = "navyblue", pch = 16)

lm2b <- lm(log(NBSS2b[3:11]) ~ log(mids2[3:11]))
summary (lm2b) # b = -1.5, OK, excellent fit! ZERO bias!!!
text( x = -7,y = 7,"b = -1.5, perfect fit, zero bias!" )
abline(lm2b)

x1 <-  -1.53993
x2 <-  -1.5

percent_bias <- abs(x1 - x2) / x1  * 100
percent_bias #  2.6 percent

# percent bias in slope "b" was always less than 3%, when using non-lner binning, and at least 5 bins 


# I.6 Functions ... BIAS tests, input vales vs output, with bias ----------

# input : Input Power law vector "b_abundance"  value
# output : Power law vector "b_abundance", and "b_biomass"  values

# input function

  fun_makes_powerlaw_vector<-  function(alpha, n=100000, seed =  123) {

set.seed(seed)
n <- n

# Standard Power-Law distribution (NNNS = -2, NBSS = -1): -------------
alpha <- alpha
xmin <- 1e-4   # small cutoff to avoid divergence

u <- runif(n)

x <- ( (1 - u) * xmin^(1 - alpha) + u * 1^(1 - alpha) )^(1 / (1 - alpha))

summary(x)
range(x)


# Standard Power Law distribution: -------------

x12  <- body_mass.pw <- x # biomass b = -1, abundance = b = -2

; x12

}

vec1_2 <- fun_makes_powerlaw_vector(2)



# Variability (CVR and VSR)  -----------------
# Variability represented in the plot and model...R-squared:  0.96
# CVR <- abs(  sd( resid(lm1) ) / mean(lm1$fitted.values) * 100 )

cvr.fun <- function(lm) { CVR <- round ( abs(  sd( resid(lm) ) / mean(lm$fitted.values) * 100 ), 2)
; CVR }

 # cvr.fun(lm1b)

# CVR: 18.2 percent (direct binning of biomass, N = 7 useful bins)
# R-squared:  0.96 (direct binning of biomass, N = 7 useful bins)

# VSR: Variability-to-slope-ratio
# The CVR-to-slope-ratio (i.e., error / information/ or noise / signal ratio, or 
#         variability-to-slope-ratio
#          Abbrev.: VSR) was then computed and compared across multiple methods (e.g., NNSS and NBSS). 

cvr.VSR.fun <- function(lm) { 
  CVR <- round ( abs(  sd( resid(lm) ) / mean(lm$fitted.values) * 100 ), 2)
  VSR <- abs( CVR / lm$coef[2] )
  ; c(CVR,VSR) }

# cvr.VSR.fun(lm1b)








# output functions -------------

# input_body_mass_data_vector <- vec1_2
#   inputalpha <- 2

#fun_direct_A_b_CVR.VSR.biomass_input_output_bias - Abundance spectrum (NNSS) fom input alpha 
# with slope , slope bias, and variability estimates

fun_direct_A_b_CVR.VSR.biomass_input_output_bias <- function(inputalpha, input_body_mass_data_vector) {
  
  body_mass.pw15_25<- input_body_mass_data_vector
  
  #length of the vector: n elements (80 in this case)
  n <- 80
  half <- floor(n / 2)
  exponents <- seq(-half, half - 1)   # symmetric around 0
  x <- 2^exponents
  log2vec <- x  # the standard log2-binning vector, these are the breaks!
  
  w2c <- diff(log2vec)
  
  log2mids  <- (head(log2vec, -1) + tail(log2vec, -1)) / 2

  hpw2 <- hist(
    body_mass.pw15_25,
    breaks = log2vec,
    plot = F,
    col = "lightgray",
    main = "Abund., power-law, b = -2, log2-binning",
    xlab = "Individual body_mass",
    ylab = "N of individuals (N) = Abundance"
  )
  
  # Calculating "k"  --------------------
  
  hpw2$breaks # the non-linear binning vector (log2)
  w2 <- diff(hpw2$breaks) #  bin withs (w) of the linear binning vector
  hpw2$mids # midpoints (mp)  of the linear binning vector
  

    Abund_per_bin2 <- hpw2$counts 
    
  #biomass_per_bin2
  #sum(biomass_per_bin, na.rm = T) # 29.6 g m-3 Total Biomass,   OK!
  
  # plot(biomass_per_bin2 ~ mids2)
  # plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
  # 
  # biomass_per_bin2 ins BWIB! 
  # steep B distributin gives POSITIVE BWIB distribution (!) because of distorting effect of non-linear binning!
  # BWIB = Bin-width inflated Biomass
  
  # points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)
  
  NNSS2b<- Abund_per_bin2/w2c
  sum(NNSS2b) # 150000
  
  max.pos <-  which.max(NNSS2b)
  end.pos <- max.pos+6
  
  # plot(log(NBSS2b) ~ log(log2mids))  # flat slpe b = 0, as expected for BWIB!
  #   points(log(NBSS2b[max.pos:end.pos ]) ~ log(log2mids[max.pos:end.pos]), col = "navyblue", pch = 16)
  # 
  lm2b <- lm(log(NNSS2b[max.pos:end.pos]) ~ log(log2mids[max.pos:end.pos]))
  # summary (lm2b) # b = -1.5, OK, excellent fit! ZERO bias!!!
  # text( x = -7,y = 7,"b = -1.5, perfect fit, zero bias!" )
  # abline(lm2b)
  # 
  inputalpha <- inputalpha
  
  input_b_abund <- (inputalpha  * -1) 
  
  x1 <-  input_b_abund
  x2 <- b_estimate <-   lm2b$coefficients[2]
  percent_bias <-  abs( abs(x1 - x2) / x1  * 100)
  
  
  res.cvr.vsr <- cvr.VSR.fun(lm2b)
  
  
  res <- list(inputalpha = inputalpha,
              input_b_abund=  input_b_abund,
              b_estimate=  as.numeric( b_estimate), 
              percent_bias= as.numeric(percent_bias),
              CVR <- res.cvr.vsr[1],
              VSR <- res.cvr.vsr[2]  ) 
  
  ;  res }

fun_direct_A_b_CVR.VSR.biomass_input_output_bias(2,vec1_2 )


fun_indirect_B_b_CVR.VSR.biomass_input_output_bias <- function(inputalpha, input_body_mass_data_vector) {
  
  body_mass.pw15_25<- input_body_mass_data_vector
  
  #length of the vector: n elements (80 in this case)
  n <- 80
  half <- floor(n / 2)
  exponents <- seq(-half, half - 1)   # symmetric around 0
  x <- 2^exponents
  log2vec <- x  # the standard log2-binning vector, these are the breaks!
  
  w2c <- diff(log2vec)
  
  log2mids  <- (head(log2vec, -1) + tail(log2vec, -1)) / 2
  
  hpw2 <- hist(
    body_mass.pw15_25,
    breaks = log2vec,
    plot = F,
    col = "lightgray",
    main = "Abund., power-law, b = -2, log2-binning",
    xlab = "Individual body_mass",
    ylab = "N of individuals (N) = Abundance"
  )
  

  hpw2$breaks # the non-linear binning vector (log2)
  w2 <- diff(hpw2$breaks) #  bin withs (w) of the linear binning vector
  hpw2$mids # midpoints (mp)  of the linear binning vector
  
  
  Abund_per_bin2 <- hpw2$counts 

  Biom_per_bin2 <- hpw2$counts * hpw2$mids
  
    
  #biomass_per_bin2
  #sum(biomass_per_bin, na.rm = T) # 29.6 g m-3 Total Biomass,   OK!
  
  # plot(biomass_per_bin2 ~ mids2)
  # plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
  # 
  # biomass_per_bin2 ins BWIB! 
  # steep B distributin gives POSITIVE BWIB distribution (!) because of distorting effect of non-linear binning!
  # BWIB = Bin-width inflated Biomass
  
  # points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)
  
  NBSS2b<- Biom_per_bin2/w2c
  
  max.pos <-  which.max(NBSS2b)
  end.pos <- max.pos+6
  
  # plot(log(NBSS2b) ~ log(log2mids))  # flat slpe b = 0, as expected for BWIB!
  #   points(log(NBSS2b[max.pos:end.pos ]) ~ log(log2mids[max.pos:end.pos]), col = "navyblue", pch = 16)
  # 
  lm2b <- lm(log(NBSS2b[max.pos:end.pos]) ~ log(log2mids[max.pos:end.pos]))
  # summary (lm2b) # b = -1.5, OK, excellent fit! ZERO bias!!!
  # text( x = -7,y = 7,"b = -1.5, perfect fit, zero bias!" )
  # abline(lm2b)
  # 
  inputalpha <- inputalpha
  
  input_b_biomass <- (inputalpha  * -1) +1
  
  x1 <-  input_b_biomass
  x2 <- b_estimate <-   lm2b$coefficients[2]
  percent_bias <-  abs( abs(x1 - x2) / x1  * 100)
  
  
  res.cvr.vsr <- cvr.VSR.fun(lm2b)
  
  
  res <- list(inputalpha = inputalpha,
              input_b_biomass=  input_b_biomass,
              b_estimate=  as.numeric( b_estimate), 
              percent_bias= as.numeric(percent_bias),
              CVR <- res.cvr.vsr[1],
              VSR <- res.cvr.vsr[2]  ) 
  
  ;  res }

fun_indirect_B_b_CVR.VSR.biomass_input_output_bias(2,vec1_2 )


fun_direct_B_b_CVR.VSR.biomass_input_output_bias <- function(inputalpha, input_body_mass_data_vector) {
  
  body_mass.pw15_25<- input_body_mass_data_vector
  
  #length of the vector: n elements (80 in this case)
  n <- 80
  half <- floor(n / 2)
  exponents <- seq(-half, half - 1)   # symmetric around 0
  x <- 2^exponents
  log2vec <- x  # the standard log2-binning vector, these are the breaks!
  
  w2c <- diff(log2vec)
  
  log2mids  <- (head(log2vec, -1) + tail(log2vec, -1)) / 2
  
  binnedB2 <- cut(body_mass.pw15_25, breaks = log2vec, include.lowest = TRUE)
  
  biomass_per_bin2 <- tapply(body_mass.pw15_25, binnedB2, sum)
  #biomass_per_bin2
  #sum(biomass_per_bin, na.rm = T) # 29.6 g m-3 Total Biomass,   OK!
  
  # plot(biomass_per_bin2 ~ mids2)
  # plot(log(biomass_per_bin2) ~ log(mids2))  # flat slpe b = 0, as expected for BWIB!
  # 
  # biomass_per_bin2 ins BWIB! 
  # steep B distributin gives POSITIVE BWIB distribution (!) because of distorting effect of non-linear binning!
  # BWIB = Bin-width inflated Biomass
  
  # points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)
  
  NBSS2b<- biomass_per_bin2/w2c
  sum(NBSS2b) # 150000
  
 max.pos <-  which.max(NBSS2b)
 end.pos <- max.pos+6
 
 # plot(log(NBSS2b) ~ log(log2mids))  # flat slpe b = 0, as expected for BWIB!
 #   points(log(NBSS2b[max.pos:end.pos ]) ~ log(log2mids[max.pos:end.pos]), col = "navyblue", pch = 16)
  # 
  lm2b <- lm(log(NBSS2b[max.pos:end.pos]) ~ log(log2mids[max.pos:end.pos]))
 # summary (lm2b) # b = -1.5, OK, excellent fit! ZERO bias!!!
  # text( x = -7,y = 7,"b = -1.5, perfect fit, zero bias!" )
  # abline(lm2b)
  # 
  inputalpha <- inputalpha
  
  input_b_biomass <- (inputalpha  * -1) +1
  
  x1 <-  input_b_biomass
  x2 <- b_estimate <-   lm2b$coefficients[2]
  percent_bias <-  abs( abs(x1 - x2) / x1  * 100)
  
  
  res.cvr.vsr <- cvr.VSR.fun(lm2b)
  
  
  res <- list(inputalpha = inputalpha,
             input_b_biomass=  input_b_biomass,
             b_estimate=  as.numeric( b_estimate), 
             percent_bias= as.numeric(percent_bias),
             CVR <- res.cvr.vsr[1],
             VSR <- res.cvr.vsr[2]  ) 
  
  ;  res }

fun_direct_B_b_CVR.VSR.biomass_input_output_bias(2,vec1_2 )



# 1.7 LOOPS --------------------------- 
# run loops  with different input values, Abundance and Biomass spectra ------------


# 1.7a prepare  LOOPs --------------

# N simulations per loop ---------------

 n.sim <- 20 #    (very fast setting, for first tests only)
# n.sim <- 200 #  (fast setting, for first tests only)
# n.sim <- 500 #  (this may take some time...)
# n.sim <- 2000 # (caution, this may take some time...)


set.seed(1)


input.alpha.values.vector <- runif(n.sim, min = 1.4 , max  = 2.6)

output.table <- data.frame(inputalpha = rep(NA,n.sim), 
                           input_b  = rep(NA,n.sim),
                           b_estimate  = rep(NA,n.sim),
                          percent_bias  = rep(NA,n.sim),
                          CVR_estimate  = rep(NA,n.sim),
                          SVR_estimate  = rep(NA,n.sim))


output.table1 <- output.table2 <- output.table3 <- output.table

# 1.7b run LOOPs, run as loop with different input values ------------

t1 <- Sys.time() # starting time


# run LOOP 1 ----- 
# input vs Abundance (NNSS) ----------------
fun_direct_A_b_CVR.VSR.biomass_input_output_bias


for (w in 1 : n.sim) { 
  
  alpha.in <-   input.alpha.values.vector[w] 
  
  vec.in <- fun_makes_powerlaw_vector(alpha.in)
  
  #res.out <-  fun_b_biomass_input_output_bias(alpha.in,vec.in )
  
  res.out <-  fun_direct_A_b_CVR.VSR.biomass_input_output_bias(alpha.in,vec.in )
  
  
  
  output.table1[w,1] <-  res.out[1]
  output.table1[w,2] <- res.out[2]
  output.table1[w,3] <- res.out[3]
  output.table1[w,4] <-  res.out[4]
  
  output.table1[w,5] <- res.out[5]
  output.table1[w,6] <-  res.out[6]
  
}

output.table1

summary(output.table1)

hist(output.table1$percent_bias)


# run LOOP 2 ----- 
# input vs indirect B (B = A * M)----------------
# indirect bimass calculation, from Abundance

fun_indirect_B_b_CVR.VSR.biomass_input_output_bias


for (w in 1 : n.sim) { 
  
  alpha.in <-   input.alpha.values.vector[w] 
  
  vec.in <- fun_makes_powerlaw_vector(alpha.in)
  
  #res.out <-  fun_b_biomass_input_output_bias(alpha.in,vec.in )
  
  res.out <-  fun_indirect_B_b_CVR.VSR.biomass_input_output_bias(alpha.in,vec.in )
  
  
  
  output.table2[w,1] <-  res.out[1]
  output.table2[w,2] <- res.out[2]
  output.table2[w,3] <- res.out[3]
  output.table2[w,4] <-  res.out[4]
  
  output.table2[w,5] <- res.out[5]
  output.table2[w,6] <-  res.out[6]
  
}

output.table2

summary(output.table2)

hist(output.table2$percent_bias)


# run LOOP 3 ----- 
# input vs direct B

for (w in 1 : n.sim) { 
  
  alpha.in <-   input.alpha.values.vector[w] 
  
  vec.in <- fun_makes_powerlaw_vector(alpha.in)
  
 #res.out <-  fun_b_biomass_input_output_bias(alpha.in,vec.in )
 
 res.out <-  fun_direct_B_b_CVR.VSR.biomass_input_output_bias(alpha.in,vec.in )
 
 output.table3[w,1] <-  res.out[1]
 output.table3[w,2] <- res.out[2]
 output.table3[w,3] <- res.out[3]
 output.table3[w,4] <-  res.out[4]

  output.table3[w,5] <- res.out[5]
 output.table3[w,6] <-  res.out[6]
 
  }

output.table3

summary(output.table3)

hist(output.table3$percent_bias)

# Bias in b (output vs input b_biomass)
#   is always below  2 percent, mostly below 0.7 percent)

t2 <- Sys.time()
duration <- t2-t1
duration # the time it takes for  3 loops (2.5 minutes, for 2000 simulations)

# 1.7c write outputs to disk ----------------- 
# Save outputs to disk---------------
setwd("~/Papers/0000 - paper NBSS Metod- bias and critique - new method")
write.table(output.table1, file =  "output.table1.csv", sep = ";")
write.table(output.table2, file =  "output.table2.csv", sep = ";")
write.table(output.table3, file =  "output.table3.csv", sep = ";")

output.table1 <-   read.table(file  =  "output.table1.csv", sep = ";")
output.table2 <-   read.table( file =  "output.table2.csv", sep = ";")
output.table3 <-   read.table( file =  "output.table3.csv", sep = ";")

# output.table1 <-  read.table(file  =  "output.table1b.csv", sep = ";") # n = 2000
# output.table2 <-  read.table( file =  "output.table2b.csv", sep = ";") # n = 2000
# output.table3 <-  read.table( file =  "output.table3b.csv", sep = ";") # n = 2000


# 1.7d Compare outputs ----------------
# tests (Permutation and Mann-Whitney tests)

# variability in Abund vs Biomass spectra ---------------

# Mann-Whitney tests -----------
wilcox.test(output.table1$SVR_estimate , output.table2$SVR_estimate)
#p-value < 2.2e-16

median(output.table1$SVR_estimate) # 0.04189656
median (output.table2$SVR_estimate) #  0.152067 , Biomass has 3 times more variabiability than A, OK due to added M effect

wilcox.test(output.table1$SVR_estimate , output.table3$SVR_estimate)
#p-value < 2.2e-16

wilcox.test(output.table2$SVR_estimate , output.table3$SVR_estimate)
# p-value = 0.007702 ( n = 200),  p-value = 1.563e-06 (n = 2000)

wilcox.test(output.table2$CVR_estimate , output.table3$CVR_estimate)
# p-value = 0.08, ( n = 200),p-value = 0.002478 (n = 2000)

wilcox.test(output.table2$percent_bias , output.table3$percent_bias)
# p-value = 0.45, ( n = 200), p-value = 0.2159 ( n = 2000)

# permutation test (coin, median_test and ) ----------

# stack both vectors (table 1 and 2, Abund vs Biomass, SVR variability index) into one data frame for testing
x1 <- output.table1$SVR_estimate
x2 <- output.table2$SVR_estimate
SVR_estimate  <- c(x1, x2)
group <- factor(c(
  rep("table1", length(x1)),  rep("table2", length(x2))  ))
df12 <- data.frame(SVR_estimate = SVR_estimate, group = group)
summary(df12)


# stack both vectors (table 2 and 3, indirect vs direct biomass estimations, SVR variability index) into one data frame for testing
x2 <- output.table2$SVR_estimate
x3 <- output.table3$SVR_estimate
SVR_estimate  <- c(x2, x3)
group <- factor(c(
  rep("table2", length(x2)),  rep("table3", length(x3))  ))
df23 <- data.frame(SVR_estimate = SVR_estimate, group = group)
summary(df23)

median( output.table2$SVR_estimate)
# indirect B (B = A * M) : SVR =  0.152, slightly less variability (linearizing SAA effect?) 

median( output.table3$SVR_estimate) 
# directly binned  B: SVR = 0.166, slightly more variability 

library(coin)
coin::oneway_test(statistic="median", df23$SVR_estimate ~ df23$group)
# p-value = 0.03 (n = 2000)# medians are different between direct and indirect Biomass! 

library(coin)

set.seed(123)
median_test(df23$SVR_estimate ~ df23$group, distribution = "approximate")
# p-value < 1e-04 (n = 2000);  p = 0.0201 (n = 200)
# medians are different (p < 0.0001) between direct and indirect Biomass!

# ?median_test
# ? oneway_test

set.seed(123)
oneway_test(df23$SVR_estimate ~ df23$group,
            statistic = "median",
            distribution = approximate(nresample = 10000))

set.seed(123)
y  <- df23$SVR_estimate
grp <- df23$group
obs <- abs(median(y[grp=="table2"]) - median(y[grp=="table3"]))

perm <- replicate(10000, {
  g <- sample(grp)
  abs(median(y[g=="table2"]) - median(y[g=="table3"]))
})

p_value <- mean(perm >= obs)
p_value
# p = 0.0053

# percent difference in variability (bias in variability, direct  indirect B estimate)

x1 = 0.152 # SVR of table 2 (SVR with indirect B  estimation, B = A * M) 
x2 =  0.166  # SVR of table 3, direct binning of body mass values into biomass
# calculate percent bias:
percent_bias <- abs(x1 - x2) / x1  * 100
percent_bias # 9.2 percent





boxplot(output.table2$SVR_estimate , output.table3$SVR_estimate, main = "Direct vs indirect biomass spectra" ,
        ylab = "SVR , Slope to variabllity ratio",
        xlab = "Indirect (B = A * B)     ;   Direct biomass")

summary(output.table1$SVR_estimate)

summary(output.table2$SVR_estimate)
summary     (output.table3$SVR_estimate)


# table 1 vs 2 ratio
table1vs2ratio <- output.table1$SVR_estimate /output.table2$SVR_estimate
summary(table1vs2ratio)
hist(table1vs2ratio)

# table 1 vs 2 diff
table1vs2diff <- output.table1$SVR_estimate - output.table2$SVR_estimate
summary(table1vs2diff)
hist(table1vs2diff)



# Output table 3 (direct estimation) has significantly higher variabilty! SAA effect!

# I.8 Test the bNBSS Backtransformed (backtransform = rescaling and redimensionalizing ) ----------

# I.8.A. k-normalization and rescaling (!) - the kNBS --------------
#      kNBS:  "k-normalized NBSS"
# check total Biomass (BWIB vs kNBS)


# k-normalized abundance, knNSS (divide by k) --------------------
# knNSS: k-normalized numbers size spectrum ------------------

A_total <-   sum(hpw15_252$counts) # total Abundance = 100.000 ind m-3, OK!
A_total

# k-normalization
kNNSS <- hpw15_252$counts/k2 
# normalize using the untitless and non-dinasional rescling factor k
# k = w / R 
# ( k = bin width / Range)


sum(k2) # 1, OK!

sum(kNNSS) # 925301081 , not equal to total abundance, due to non-linear effects

#now rescaling gives the
# "bNNSS" -  the backtransformed (rescaled, k-normalized)  numbers size spectrum


# rescaling with correction factor "F_prime" F_p

 F_p <-  sum(hpw15_252$counts) / sum(kNNSS)
F_p


#rescale kNNSS with F_p, to obtain bNNSS (backtransformed NNSS):

bNNSS <- kNNSS * F_p

sum(bNNSS) # OK! 1e+05, equals 100,000. OK!

plot( log(kNNSS) ~  log(hpw15_252$mids)  ,col = "darkgrey",
      main =       "Power-law, log2-binning, norm. Abundance, bNNSS")


plot( log(bNNSS) ~  log(hpw15_252$mids)  ,col = "darkgreen",
      main =       "Power-law, log2-binning, norm. Abundance, bNNSS, OK")

sum(bNNSS)

# Nicely linear log-log plot, except for the  lower extreme

points(log(bNNSS[3:11]) ~  log(hpw15_252$mids[3:11]), 
       pch = 16, col = "blue" )

lm1 <- lm( log(bNNSS[3:11]) ~  log(hpw15_252$mids[3:11]) )

abline(lm1)
summary(lm1)

# Intercept a =         -11.03
# slope     b =  -2.53, bNNSS, OK  (original alpha: 2.5, b_biomass: -1.5)


text(x = -7, y = 4.5, "b = -2.5, OK" )


# Calculate total Biomass 
# assume that the sampling unit represents  1 cubic meter

sum(body_mass.pw15_25) #  29.6  g m-3  = Total Biomass, Btotal



# I.8.B. w-normalization and D-rescaling of NBSS, to otain  bNBS ----------

#      bNBS:  "back-transformed NBSS"
# check total Biomass (BWIB vs bNBS)


# II. Example, using a standard log2 binning vector - example calculation -------

# Crate and use the standard log2-binning vector (bin mids = ...,0.25,0.5,1,2,4,8,…) 
# Create log2 bins, with 80 elements, centered at "1"

#length of the vector: n elements (80 in this case)
n <- 80
half <- floor(n / 2)
exponents <- seq(-half, half - 1)   # symmetric around 0
x <- 2^exponents
log2vec <- x  # the standard log2-binning vector


log2(x)
log10(x)
plot(x)


hpw2 <- hist(
  body_mass.pw,
  breaks = non_lin_breaks,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

plot(hpw2$counts ~ hpw2$mids)
plot(log10(hpw2$counts) ~ log10(hpw2$mids))
summary(lm(log10(hpw2$counts[3:12]) ~ log10(hpw2$mids[3:12])))

# use the standard log2 binning vector

hpw2c <- hist(
  body_mass.pw,
  breaks = log2vec,
  plot = TRUE,
  col = "lightgray",
  main = "Abund., power-law, b = -2, log2-binning",
  xlab = "Individual body_mass",
  ylab = "N of individuals (N) = Abundance"
)

plot(hpw2c$counts ~ hpw2c$mids)
plot(log10(hpw2c$counts) ~ log10(hpw2c$mids),
     main = "BWIA using the complete log2 standard binning vector")
points(log10(hpw2c$counts[28:40]) ~ log10(hpw2c$mids[28:40]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm(log10(hpw2c$counts[28:40]) ~ log10(hpw2c$mids[28:40])))

plot(hpw2c$counts ~ hpw2c$mids)
plot(log10(hpw2c$counts) ~ log10(hpw2c$mids),
     xlim = c(-5, 1),
     main = "BWIA using the complete log2 standard binning vector")

points(log10(hpw2c$counts[28:40]) ~ log10(hpw2c$mids[28:40]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm(log10(hpw2c$counts[28:40]) ~ log10(hpw2c$mids[28:40])))


# III. import example dataset ----------

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


# Total number of individual in the sample
N_total = length(data.carbon_g_ind_1)# N = 7088, indiv. per sample
N_total

# Abundance (ind.  m-3)
A_total = length(data.carbon_g_ind_1) / 100 # A =  N/S , indiv. per cubic meter
A_total # 70.88 total Abundance (ind.  m-3)

# Biomass (g.  m-3)
B_total = sum(data.carbon_g_ind_1) / 100 # B =  sum(M)/S, g. per cubic meter
B_total # 0.00342 total Biomass (g.  m-3)



# III.2 Abundance (BWIA, NNSS, bNNS) ------------

# BINNING --------------

# use the standard log2 binning vector

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

plot(hpw2c$counts ~ hpw2c$mids)
plot(log10(hpw2c$counts) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "BWIA using the complete log2 standard binning vector")

# Define useful (linear) size range for linear model
# (from maximum to first empty bin
# Corresponds to the size range, where sampliinng is really reprentative of the communty size spectrum
# = size range of quantitative sampling (!)

points(log10(hpw2c$counts[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm(log10(hpw2c$counts[24:30]) ~ log10(hpw2c$mids[24:30])))


# Numbers to Abundwnce (divide by S)
# sampled volume "S": 100 cubic meters (hypothetical)
# A = N / S
# binned (not normalized) Abundance = BWIA
# Bin-width-nflated Abundance (BWIA)

BWIA <- hpw2c$counts /100


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

   bNNS = NNSS * D_Abund

bNNS 

sum(bNNS) # 70.88 = total Abundance, OK!



plot(bNNS ~ hpw2c$mids)
plot(log10(bNNS) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "bNNSS using the log2 standard binning vector")

points(log10(bNNS[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1 <- lm(log10(bNNS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2, perfectly OK for Abundance vs body mass 

plot(bNNS ~ hpw2c$mids)

plot(log10(bNNS) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "bNNSS using the log2 standard binning vector",
     xlab = "log10 Body mass (log10 (gC ind-1) ",
     ylab = "log10 Abundance (bNNSS, log10(ind. m-3) )"
     )

points(log10(bNNS[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1 <- lm(log10(bNNS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2, perfectly OK for Abundance vs body mass 

sum(bNNS) # total bNNS represents total Abundance (ind. m-3), OK!


# Variability represented in the plot and model...R-squared:  0.96
# CVR <- abs(  sd( resid(lm1) ) / mean(lm1$fitted.values) * 100 )

cvr.fun <- function(lm) { CVR <- round ( abs(  sd( resid(lm) ) / mean(lm$fitted.values) * 100 ), 2)
; CVR }
cvr.fun(lm1)

# CVR: 43.6 percent (direct binning of Abundance, N = 7 useful bins)
# R-squared:  0.99 (direct binning of Abundance, N = 7 useful bins)




# III.2Biomass (BWIB, NBSS, bNBS) ------------
# Calculate Total Biomass per Bin (non-linear binning with log2)

# III.2.a ----------------
#  Biomass from Abundance ---------------
# Calculation: B = A * M 
# or:  B = (N * M) / S ,    S = sampling unit (e.g., 100 m-3)

# non-linear binning, log2 ----------------

hpw3 <- hist(data.carbon_g_ind_1, breaks = log2vec, plot = F) # to obtain mids, etc
mids3 <- hpw3$mids

# B = N * M /S
Btotal3 <- sum(hpw3$counts * hpw3$mids) / 100# Biomass = 0.003608723 g m-3
B3 <- (hpw3$counts * hpw3$mids )/100
sum(B3)

NBSS <- B3/w2c
sum(NBSS) # 106.32, OK

mids2 <- hpw3$mids
plot(NBSS ~ mids2)
plot(log(NBSS) ~ log(mids2)) 

points(log(NBSS[24:30]) ~ log(mids2[24:30]), col = "navyblue", pch = 16)

lm1 <- lm(log(NBSS[24:30]) ~ log(mids2[24:30]))
summary(lm1)
abline(lm1)

# slope = -1.0176, approx b = -1, OK!


# Variability represented in the plot and model...R-squared:  0.96
# CVR <- abs(  sd( resid(lm1) ) / mean(lm1$fitted.values) * 100 )

cvr.fun <- function(lm) { CVR <- round ( abs(  sd( resid(lm) ) / mean(lm$fitted.values) * 100 ), 2)
; CVR }
cvr.fun(lm1)

# CVR: 17.46 percent (biomass from A * M, N = 7 useful bins)
# R-squared:  0.96 (biomass  from A * M, N = 7 useful bins)



# III.2.b ----------------
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



# IV. Comparison of variability (CVR): -----------------------

# Abundance (NNSS)

# (Intercept)                 -3.3431 +- 0.3756  
#   log10(hpw2c$mids[24:30])  -2.0176 +- 0.0920  
# CVR: 43.6 percent ( binning of Abundance, N = 7 useful bins)
# R-squared:  0.99 ( binning of Abundance, N = 7 useful bins)
# Abundance-body-mass slope alpha = -2.02, OK!

# Biomass from Abundance, (NBSS and bNBS)
# (Intercept)        -7.69     +_   0.8649 SE   
#   log(mids2[24:30])  -1.02    +-  0.0920SE
# CVR: 17.46 percent (biomass from A * M, N = 7 useful bins)
# R-squared:  0.96 (biomass from A * M, N = 7 useful bins)
# Biomass-body-mass slope b = -1.02, OK!

# Biomass (NBSS and bNBS) - direct binning
# (Intercept)       -7.98    +- 0.86 SE
#   log(mids2[24:30]) -1.04   +-  0.092 SE   
# CVR: 18.2 percent (direct binning of biomass, N = 7 useful bins)
# R-squared:  0.96 (direct binning of biomass, N = 7 useful bins)
# Biomass-body-mass  slope b  = -1.04, OK!

# V. bNBS (new method) ----------------------

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


