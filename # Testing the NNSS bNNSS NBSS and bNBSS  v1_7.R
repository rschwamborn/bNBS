#  Testing the NNSS, bNNSS, NBSS, and bNBSS size spectra
# Version   1.7

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

# biomass_per_bin2 ins BWIB! flat dostributin because of distorting effect of non-linear binning!
# BWIB = Bin-width inflated Biomass

points(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]), col = "navyblue", pch = 16)

lm.BWIB <- lm(log(biomass_per_bin2[3:7]) ~ log(mids2[3:7]))
summary(lm.BWIB) 
# flat BWIB- body mass spectrum , b = 0 (OK)



# I.6. Compare input vs output distribution results  ------------------------

# Input: power law data  
# output: plots and linear models (a,b, R²) for A (NNSS), BfromA, B  

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


# I.7 As function ... BIAS tests, input vales vs output, with bias ----------

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

# output function


fun_b_biomass_input_output_bias <- function(inputalpha, input_body_mass_data_vector) {

  body_mass.pw15_25<- input_body_mass_data_vector
  
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

inputalpha <- inputalpha

input_b_biomass <- (inputalpha  * -1) +1

x1 <-  input_b_biomass
x2 <- b_estimate <-   lm2b$coefficients[2]
percent_bias <-  abs( abs(x1 - x2) / x1  * 100)

res<- list(inputalpha = inputalpha,
           input_b_biomass=  input_b_biomass,
           b_estimate=  as.numeric( b_estimate), 
           percent_bias= as.numeric(percent_bias) ) 

;  res }

fun_b_biomass_input_output_bias(2,vec1_2 )

# LOOP run as loop with different input values ------------

# N simulations

n.sim <- 200

input.alpha.values.vector <- runif(n.sim, min = 1.4 , max  = 2.6)

output.table <- data.frame(inputalpha = rep(NA,n.sim), 
                           input_b_biomass  = rep(NA,n.sim),
                           b_estimate  = rep(NA,n.sim),
                           percent_bias  = rep(NA,n.sim))
# run LOOP -----

for (w in 1 : n.sim) { 
  
  alpha.in <-   input.alpha.values.vector[w] 
  
  vec.in <- fun_makes_powerlaw_vector(alpha.in)
  
 res.out <-  fun_b_biomass_input_output_bias(alpha.in,vec.in )
  
 output.table[w,1] <-  res.out[1]
 output.table[w,2] <- res.out[2]
 output.table[w,3] <- res.out[3]
 output.table[w,4] <-  res.out[4]
 
  }

output.table

summary(output.table)

# Bias in b (output vs input b_biomass)
#   is always below  2 percent, mostly below 0.7 percent)



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

#now rescaling goiives the
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



# I.8.B. w-normalization and F-rescaling of the NBSS, the bNBSS ----------
#      bNBSS:  "back-transformed NBSS"
# check total Biomass (BWIB vs bNBSS)



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


carbon_g_ind_1_zoopldataset.table <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/carbon_no_Artef zz_Merged_ zz_Infl_ abr1_noronha_st01_300.csv", sep=";")

data.carbon_g_ind_1 <- carbon_g_ind_1_zoopldataset.table$g_carbon


no_Artef.zz_Merged_.zz_Infl_.abr1_plataforma_st27_300_d1_d1_1_da <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/no_Artef zz_Merged_ zz_Infl_ abr1_plataforma_st27_300_d1_d1_1_da.csv", sep=";")



# Total number of individual in the sample
N_total = length(data.carbon_g_ind_1)# N = 7088, indiv. per sample
N_total

# Abundance (ind.  m-3)
A_total = length(data.carbon_g_ind_1) / 100 # A =  N/S  , indiv. per cubic meter
A_total # 70.88 total Abundance (ind.  m-3)



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
summary(lm1 <- lm(log10(NNSS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2, perfectly OK for Abundance vs body mass 

# BACKTRANSFORM -----------
# fom NNSS to bNNSS
# After normalizing (dividing BWIA by w, where NNSS = BWIA / w) to obain the common NBSS, rescaling is needed. This can be done by using the 
# w-scaled dimensional correction factor C, where
# 
# C =  BWIA_total  / NNSS_total  
# 
# All NBSSi  are then transformed back into the scale of B,
# i.e., into a "backtransformed” NBSS (bNBSS), where

# calculate C (dimensional correction factor) -------

# C =  BWIA_total  / NNSS_total
BWIA_total = sum(BWIA) # 70.88 = total Abundance, OK!
NNSS_total = sum(NNSS)

C_Abund = BWIA_total / NNSS_total

   bNNSS = NNSS * C_Abund

bNNSS 

sum(bNNSS) # 70.88 = total Abundance, OK!



plot(bNNSS ~ hpw2c$mids)
plot(log10(bNNSS) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "bNNSS using the log2 standard binning vector")

points(log10(bNNSS[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1 <- lm(log10(bNNSS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2, perfectly OK for Abundance vs body mass 

plot(bNNSS ~ hpw2c$mids)

plot(log10(bNNSS) ~ log10(hpw2c$mids),
     xlim = c(-7, 1),
     main = "bNNSS using the log2 standard binning vector",
     xlab = "log10 Body mass (log10 (gC ind-1) ",
     ylab = "log10 Abundance (bNNSS, log10(ind. m-3) )"
     )

points(log10(bNNSS[24:30]) ~ log10(hpw2c$mids[24:30]), 
       pch = 16, col = "navy")

hpw2c$counts
summary(lm1 <- lm(log10(bNNSS[24:30]) ~ log10(hpw2c$mids[24:30])))
abline(lm1)
# slope : b = -2, perfectly OK for Abundance vs body mass 

sum(bNNSS) # total bNNNS represents total Abundance (ind. m-3), OK!

# Biomass (BWIB, NBSS, bNBSS) ------------
# Calculate Total Biomass per Bin (non-linear binning with log2)

...

# calclate the 



