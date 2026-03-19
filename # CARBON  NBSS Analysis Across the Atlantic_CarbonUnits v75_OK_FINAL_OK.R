# NBSS analysis across the Atlantic and TLs 
# v.74.1 Carbon - medians of ALL samples (unfiltered), no picoplankton in linear models
# Corrected data-rich phytopl. ranges
# Carbon units
# Ralf Schwamborn
# Sept 8, 2025

## clear memory 

  # gc() # caution, clears memory
  # rm(list=ls()) # caution, clears memory
  
  # dev.off() # caution, deletes all plots and restores margins... 


# load libraries
  
  # install.packages("scales")
  # install.packages("ggplot2")
  # install.packages("maps")
  # install.packages("rgl")
  # install.packages("purr")
  # install.packages("broom")
  # install.packages("tidyverse")
  # install.packages("robustbase")
  #  install.packages("pacman")
 
  # load packages with  pacman::p_load (better than "library()")
  pacman::p_load(scales, ggplot2, maps, rgl)
  pacman::p_load(car, colourvalues,robustbase, tidyverse,BSDA,
                 plot3D, relaimpo, wesanderson, complmrob)
     
 
 opar <- par() # save plot parameters

 par()
 
 
# set working directory

setwd("C:/Users/RALF/Documents/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa")



# 1. Import Data  --------------------------------------------------------------

# X_vector_gCind = Bin mids for NBSS, in g C ind.^-1

X_vector_gCind <- read.table("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/Carbon_means_X_vector_v6.csv", quote="\"", comment.char="")

X_vector_gCind <- X_vector_gCind$V1
length(X_vector_gCind) # 55 carbon mids (g C ind.-1)
  
X_vector_gCind
 cat("c(", paste(X_vector_gCind, collapse = ", "), ")", sep = "")

logCvector<- log10(X_vector_gCind)

diff(logCvector)
length(logCvector)

# read ALL data (v6), NBSS in CARBON Units!!! ------------------

Data.All.Atl <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/Stations overview Atlantic size spectraV6_Carbon_units.csv")



# 2. A first quick look at the data -----------------------------------------

summary(Data.All.Atl)

attach(Data.All.Atl)

#summary(slope)   
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -3.2900 -1.2575 -0.9400 -0.8909 -0.5900  2.1000    1307 


# Numbers of samples
# Phytoplankton -------------
#   Carbon units: 545 Phytoplankton samples with data, unfiltered, includes non-linear and low-N samples
# Carbon units: 311 Phytoplankton highly linear samples  (depending on selected size range, this may be 372 samples)

(311 /545 ) # 57 % linear-shaped Phytoplankton samples





# Net-caught zoopl. --------
# Carbon: 169 Net zooplankton samples, unfiltered, no NA's (empty samples) 
#       Carbon: 134 Net zoopl. high-quality linear filtered samples
#    #Biovolume units: 111 Net zoopl. high-quality linear filtered samples

(134 /169 ) # 79 % linear-shaped net zooplankton samples


# UVP Zoopl. -----------
#    Carbon: 739   UVP samples, unfiltered  (linearly and non-linearly shaped ), no NA (empty samples) 
#      Carbon:   NEW:       207  high-quality linear filtered samples, -4.9 to -2 log10gC,  Carbon units
#     Carbon:  271 UVP  high-quality linear filtered samples (all sizes)
#   # Biovolume: 484 UVP high-quality linear filtered samples
 

 
# Micronekt. --------

#dim(meso_fi.NBSSmatrix_no_na)
# 85 useful micronekton samples, with NBSS data


# Mesopelagic fish -----------

#dim(meso_fi.NBSSmatrix_no_na)
# 86 useful samples, with NBSS data

545 + 169 + 739 + 85 + 86

727+ 16+ 739+ 86+ 85 # 1653 samples analyzed

# Total = 1624 samples ( 545 , 169 , 739 , 85 , and 86 samples for    )
# A total  1624 samples wre anlzed in tis stdt,  composed of 545, 169, 739, 86, and 85  (NBSS carbon units) samples
#   Phytoplankton,  Net-caught zooplankton,  UVP zooplankton, Mesopelagic fish, 
# and Micronekton samples, respectively. 
         # apply median by columns
       




Data.All.Atl$target_organisms <- as.factor(Data.All.Atl$target_organisms)


Data.All.Atl$SST_C_insitu <- as.numeric(Data.All.Atl$SST_C_insitu)


summary(Data.All.Atl$target_organisms)
#New:
# invmicronekton mesopelagic fish      micronekton    phytoplankton totalmicronekton 
#     4              390                  95              727               21 
# totmicronekton   totMicronekton      zooplankton 
#     13               11             1579 

#OLD
# detritus + zooplankton       mesopelagic fish            micronekton 
# 741                    356                    119 
# phytoplankton            zooplankton 
# 575                   1580 

# join and cleanup target organisms ----------
Data.All.Atl$Old_target_organisms <- Data.All.Atl$target_organisms
Data.All.Atl$OKtarget_organisms <- as.character(Data.All.Atl$target_organisms)

Data.All.Atl$OKtarget_organisms[ Data.All.Atl$OKtarget_organisms == "totmicronekton"] <- "micronekton"
Data.All.Atl$OKtarget_organisms[ Data.All.Atl$OKtarget_organisms == "totMicronekton"] <- "micronekton"
Data.All.Atl$OKtarget_organisms[ Data.All.Atl$OKtarget_organisms == "totalmicronekton"] <- "micronekton"

Data.All.Atl$OKtarget_organisms[ Data.All.Atl$OKtarget_organisms == "Phytoplankton"] <- "phytoplankton"


Data.All.Atl$target_organisms <- as.factor(Data.All.Atl$OKtarget_organisms)
summary(Data.All.Atl$target_organisms)

# NEW, OK
# invmicronekton mesopelagic fish      micronekton    phytoplankton      zooplankton 
# 4              390              140              727             1579 



# first look at the data

dim(Data.All.Atl)# 2840 rows=lines,  93 columns, ok
length(Data.All.Atl$target_organisms)# 2840 lines, ok
names (Data.All.Atl)
summary (Data.All.Atl)

names(Data.All.Atl[,37:91]) # 55 carbon mids (g C ind.-1)
length( Data.All.Atl[,37:91]  )

Data.All.Atl.NBSSmatrix <- as.matrix(Data.All.Atl[,37:91])
# 55 rows of carbon NBSS values (g C m-3 g-1C)

#first NBSS plots--------------

plot (log10(Data.All.Atl.NBSSmatrix[120,]) ~ log10(X_vector_gCind),
      main = paste( Data.All.Atl$target_organisms[120]  , " , ",  Data.All.Atl$CruiseID[120] )     )

lnumb <- 2221
plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
      main = paste( Data.All.Atl$target_organisms[lnumb]  , ", ",  
                    Data.All.Atl$CruiseID[lnumb], ", ", Data.All.Atl$Gear[lnumb]  ))

Data.All.Atl.NBSSmatrix[lnumb,]
Data.All.Atl[lnumb,]
#View(Data.All.Atl[lnumb,32:95])

lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
      main = paste( Data.All.Atl$target_organisms[lnumb]  , ", ",  
                    Data.All.Atl$CruiseID[lnumb], ", ", Data.All.Atl$Gear[lnumb], 
                    ", line = " , lnumb  ), ylim = c(-20, 13))



# plot ALL NBSS Data (loop) --------------

library(scales)

   Ylab.NBBS <- "log10(gC ind.-1)"
   XlabNBSS  <-  "log10(gC m-3 / gC ind.-1)"
  

lnumb.all <- length(Data.All.Atl$Longitude)
lnumb.all# 2840 data (lines) in C format

lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)

plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
       ylim = c(-20, 13),
      xlab = "log10(gC ind.-1)",
      ylab = "log10(gC m-3 / gC ind.-1)",
            main= "NBSS, All Atlantic, All data, n = 2840, C units")
      

for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){

  points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
  pch = 16, col = alpha("darkgreen", 0.1))
  }

# First Naive linear model (all data): ------------------
# matrix to vector (transpose), TSWA -----------------------  
library (MASS)

x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
y1 <- c ( t(Data.All.Atl.NBSSmatrix))

df1 <- data.frame(x1, y1, ylog = log10(y1))
names(df1)

# -INF to NA!!!

# df1$ylog <- replace(df1$ylog, "-Inf", NA) 

df1[df1=="-Inf"]<-NA

df2<- na.omit(df1)
dim(df1)
dim(df2)

summary(df1)
summary(df2)

mean(df2$y1)

summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.954 (naive lm model)
# Intercept=  -3.377955 (naive lm model)
summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.955  (naive rlm model)
# Intercept =  -3.2915  (naive rlm model)

abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
       lwd = 2.5, lty = 2, col = "darkorange" )


# Naive global model, All Atlantic, All Taxa ALL data, OK

# Slope: 
#     -0.9548    +- 0.0026
# Intercept:
#   -3.2915   +-  0.0191 


abline(a = -3.2915,  b= -0.9548  , xlim = c(-11, -1), ylim = c(0, 14), 
       lwd = 2.5, lty = 2, col = "darkorange" )


# 3. cleanup and data formatting -----------------------------------------------

# 3.1 convert zeros to empty bins ("NA") in NBSS  ------------------------------

# NBSS data (X5.68E.14 to X1.02E.03, columns no. 37 to 91, column 92 is NBSS slope)

summary(Data.All.Atl)

head (Data.All.Atl[,37:92])

# View (Data.All.Atl)


NBSS.raw.only <- Data.All.Atl[,37:92]

NBSS.raw.only[NBSS.raw.only == 0] <- NA

Data.All.Atl[,37:92] <- NBSS.raw.only

summary(Data.All.Atl)

# View (Data.All.Atl) 

# 3.2  Data formats: text to factors or numbers, colours  -------------------------------

# factors
Data.All.Atl$target_organisms <- as.factor(Data.All.Atl$target_organisms)
Data.All.Atl$target_target_org_code_code <- as.factor(Data.All.Atl$target_org_code)
Data.All.Atl$Gear            <- as.factor(Data.All.Atl$Gear)
Data.All.Atl$operation_mode            <- as.factor(Data.All.Atl$operation_mode)
Data.All.Atl$Biomass_determination            <- as.factor(Data.All.Atl$Biomass_determination)
Data.All.Atl$Contact_Person            <- as.factor(Data.All.Atl$Contact_Person)
Data.All.Atl$day_night            <- as.factor(Data.All.Atl$day_night)


# numbers
Data.All.Atl$SST            <- as.numeric(Data.All.Atl$SST)

Data.All.Atl$primprod            <- as.numeric(Data.All.Atl$primprod)
Data.All.Atl$salinity            <- as.numeric(Data.All.Atl$salinity)
Data.All.Atl$pH            <- as.numeric(Data.All.Atl$pH)
Data.All.Atl$chla            <- as.numeric(Data.All.Atl$chla)
Data.All.Atl$nitrate            <- as.numeric(Data.All.Atl$nitrate)
Data.All.Atl$silicate            <- as.numeric(Data.All.Atl$silicate)

Data.All.Atl$log10_1X_chla_insitu <- log10(1+Data.All.Atl$Chlorophyll__mg_m_3_insitu)



# colour codes -----------------------------------------------------------------


summary(Data.All.Atl$target_organisms)
Data.All.Atl$target.org.colour <- as.character( Data.All.Atl$target_organisms)

Data.All.Atl$target_organisms <- as.character( Data.All.Atl$target_organisms)

Data.All.Atl$target_organisms[target_organisms == "phytoplankton"]  <- "phytoplankton"


dim(Data.All.Atl)

# Separate UVP zooplankton from Net zooplankton -------------------------

# create new columns (descriptor) target_organisms_UVP_NET, with Net_zoopl and UVP_zoopl 
Data.All.Atl$target_organisms_UVP_NET <- as.character( Data.All.Atl$target_organisms)

Data.All.Atl$target_organisms_UVP_NET[Data.All.Atl$Gear =="UVP"] <- "UVP_zoopl"

summary(Data.All.Atl$Gear)

Data.All.Atl$target_organisms_UVP_NET<- as.factor(Data.All.Atl$target_organisms_UVP_NET)

summary(Data.All.Atl$target_organisms_UVP_NET)

summary(Data.All.Atl$target_organisms)



attach(Data.All.Atl)
Data.All.Atl$target.org.colour[target_organisms == "phytoplankton"]  <- "forestgreen"
Data.All.Atl$target.org.colour[target_organisms_UVP_NET == "zooplankton"]  <- "dodgerblue"
Data.All.Atl$target.org.colour[target_organisms_UVP_NET == "UVP_zoopl"]  <- "darkmagenta"

Data.All.Atl$target.org.colour[target_organisms == "micronekton"]  <- "darkred"
Data.All.Atl$target.org.colour[target_organisms == "mesopelagic fish"]  <- "darkorange"
Data.All.Atl$target.org.colour[target_organisms == "invmicronekton"]  <- "white"

Data.All.Atl$target.org.colour <- as.factor( Data.All.Atl$target.org.colour)
Data.All.Atl$target_organisms <- as.factor( Data.All.Atl$target_organisms)


summary(Data.All.Atl$target.org.colour)
# darkorange  dodgerblue forestgreen      purple 
# 376        1579         727         119 


summary(Data.All.Atl$target_organisms)
# detritus + zooplankton       mesopelagic fish            micronekton 
#                    741                    356                    119 
#          phytoplankton            zooplankton 
#                    575                   1580 

#New, OK:
# invmicronekton mesopelagic fish      micronekton    phytoplankton      zooplankton 
#       4              390              140 (+4)            727             1579 


1579-741 # 838 net-caught (includes NAs)
741 # 741 UVP zoopl. (incldes NAs)


summary(Data.All.Atl$Gear)
# 
# Aalnetz             Bongo             BT140              IKMT              Meso 
# 93                35                 8                32                31 
# Micro            MT1000            MT1216            MT1600          Multinet 
# 51                12                10               242                85 
# Multinet Midi    Niskin Bottles               RMT               UVP              WP-2 
# 581               727                55               741                33 
# Zooscan. Multinet 
# 104 
# 

# delete empty lines (no NBSS slopes) -----------------------------------------



length(Data.All.Atl$ConsecutiveNumber) # 2840 lines (datasets, ALL Data)

# df[!(is.na(df$start_pc) | df$start_pc==""), ]

#Data.All.Atl <- Data.All.Atl[!(is.na(Data.All.Atl$NBSS_slope) | Data.All.Atl$NBSS_slope ==""), ]


#length(Data.All.Atl$Longitude) # 1488 useful lines (datasets with NBSS slopes, ALL Data)


summary(Data.All.Atl$target_organisms)




# invmicronekton mesopelagic fish      micronekton    phytoplankton      zooplankton 
# 0               30               47              565              846 
# All Data with NBSS slopes (n = 1488)


# 4. Quality checks, plots and further cleanups  -------------------------------


# eliminate "exotic" NBSS (outliers) -------------------------------------
# keep only NBSS slopes between -2 and -0.1 
Data.All.Atl$slope.clean <-    Data.All.Atl$NBSS_slope

# Data.All.Atl$slope.clean[Data.All.Atl$slope.clean < -2]  <- NA
# Data.All.Atl$slope.clean[Data.All.Atl$slope.clean > -0.1]  <- NA

length(Data.All.Atl$slope.clean) # 1488 data points

# 135 exotic data deleted

summary(Data.All.Atl$NBSS_slope)  # # 1494 NBSS data  
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.9720 -1.2651 -0.9755 -0.9604 -0.6874 -0.1002     135 


hist(Data.All.Atl$NBSS_slope)
hist(Data.All.Atl$slope.clean)

# "slope.clean" new data vector with "clean" slopes (no exotic ones), tails were cut off from the distribtion
summary(Data.All.Atl$slope.clean) # "clean" only: 135 exotic data deleted (approx 9% of the data)
hist(Data.All.Atl$slope.clean)





# 4.2 New productivity indices -----------------------

# ...
# ...
# ...
# ...




# 5 Subsets by communities  ------------------------------------------------------

summary(Data.All.Atl$target_organisms)

phytopl.All.Atl <- subset(Data.All.Atl, target_organisms == "phytoplankton" )# 565 data rows
dim(phytopl.All.Atl) # 727  datasets total, but many NAs, 565 datasetss with NBSS slopes  
names(phytopl.All.Atl)


zoopl.All.Atl <- subset(Data.All.Atl, target_organisms == "zooplankton" )# 846 data rows
dim(zoopl.All.Atl)
summary(zoopl.All.Atl$Gear)


mesop_fish.All.Atl <-  subset(Data.All.Atl, target_organisms == "mesopelagic fish" )
# 47 data rows

micronekton.All.Atl <-  subset(Data.All.Atl, target_organisms == "micronekton" )
# 47 data rows


# NBSS plots by communities --------
# plot ALL NBSS Data (loop) --------------

library(scales)

lnumb.all <- length(Data.All.Atl$Longitude)
lnumb.all# 2840 data (lines) in C format


lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)

plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
      ylim = c(-7, 13), col = "white",
      xlim = c(-14, 3),
      xlab = "log10(gC ind.-1)",
      ylab = "log10(gC m-3 / gC ind.-1)",
      main= "NBSS, All Atlantic, All data, n = 2840, C units")


for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
  
  points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
          pch = 16, col = alpha("darkgreen", 0.1))
}

# Naive linear model (all data): ------------------
# matrix to to vector (transpose), TSWA -----------------------  
library (MASS)

x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
y1 <- c ( t(Data.All.Atl.NBSSmatrix))

df1 <- data.frame(x1, y1, ylog = log10(y1))
names(df1)

# -INF to NA!!!

# df1$ylog <- replace(df1$ylog, "-Inf", NA) 

df1[df1=="-Inf"]<-NA

df2<- na.omit(df1)
dim(df1)
dim(df2)

summary(df1)
summary(df2)

mean(df2$y1)

summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.954
summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.955


abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
       lwd = 2.5, lty = 2, col = "darkorange" )


# summary(log10(X_vector_gCind))
# summary(as.vector((log10(Data.All.Atl.NBSSmatrix))))

# Useful functions ----------
# (medians   and bootstrap 95% Conf. intervals)

# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

mean_   <- function(...) mean(..., na.rm=T)
median_ <- function(...) median(..., na.rm=T)

mean_ (c(4,NA, 8, 34)) #15.33, OK
median_ (c(4,NA, 8,34)) #8, OK

# Functions for RR models, for Bootstrap confidence intervals
# Bootstrap confidence intervals for slope and intercept (log10 of y values )

### 1. NOT LOG10 (simple EXAMPLE) --------------

### 2. WITH  LOG10(Y)  --------------

# EXAMPLE DATA 
x = 1:100
y = 2000 + (2* x) + rnorm(100, sd = 50)
plot(y ~ x)

df4 <- data.frame (x = x , y = y) 


### 1. NOT LOG10 (simple EXAMPLE) --------------

# confidence interval for slope:

# Bootstrap function ------------
# creating vectors to store bootstrap regression coefficients
# As fuuncion 
fun.bootstr.rr.slope.SIMPLE <- function (x,y, nruns) {
  
  # B = 2500 # n of runs
  B <- nruns # n of runs
  
  df <- data.frame (x = x , y = y) 
  df <- na.omit(df) # delete rows with NA
  
  boot.out <- rep(NA, B)
  vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
  ### LOOP ###
  for(i in 1:B){ #starting loop
    
    ##creating samples of observation IDs with replacement, of same size    as original sample
    boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
    
    #matching response and explanatory variable values to bootstrap sample   observation IDs
    ysam <- df$y[boot.id]
    xsam <- df$x[boot.id]
    
    #generating bootstrap SLR model for each bootstrap sample
    boot.mod <-rlm( ( ysam) ~  xsam )       
    
    #storing regression coefficient values for each bootstrap SLR
    boot.out[i] <- coef(boot.mod)[2]
  }
  #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
  (medianboot <-median(boot.out))
  # -0.883927
  ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
  # 2.5%      97.5% 
  #   -1.0279427 -0.7463216 
  
  
  ;  boot.CI 
}

fun.bootstr.rr.slope.SIMPLE(df4$x, df4$y, 200)

summary(lm( (df4$y) ~df4$x))

# fun.bootstr.rr.slope.SIMPLE(df4$x, df4$y, 20000)


# Boostrap 95% CI (lm):
# 2.5%        50%      97.5% 
# -0.98 -0.84 -0.71 

# OLSR CI 95% (summary, mean +_ (se*1.96)) ,
mean.slope = -0.84749 ; se.slope =  0.06547
lower.lim = mean.slope - (1.96 * se.slope);   upper.lim = mean.slope + (1.96 * se.slope);   
(lower.lim); (upper.lim)
# -0.98 to -0.72, mean = -0.85

# RLM:
# Boostrap 95% CI (rlm, nruns = 20,000):
# 2.5%          50%         97.5% 
# -1.03       -0.88       -0.74 


# "p" value
permuco::aovperm(  rlm( (df4$y) ~ df4$x )  )





## for INTERCEPT
###  not LOG, SIMPLE   --------------

# confidence interval for INTERCEPT:

# Bootstrap function ------------
# creating vectors to store bootstrap regression coefficients
# As fuuncion 
fun.bootstr.rr.intercept.SIMPLE <- function (x,y, nruns) {
  
  # B = 2500 # n of runs
  B <- nruns # n of runs
  
  df <- data.frame (x = x , y = y) 
  df <- na.omit(df) # delete rows with NA
  
  boot.out <- rep(NA, B)
  vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
  ### LOOP ###
  for(i in 1:B){ #starting loop
    
    ##creating samples of observation IDs with replacement, of same size    as original sample
    boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
    
    #matching response and explanatory variable values to bootstrap sample   observation IDs
    ysam <- df$y[boot.id]
    xsam <- df$x[boot.id]
    
    #generating bootstrap SLR model for each bootstrap sample
    boot.mod <-rlm( ( ysam) ~  xsam )       
    
    #storing regression coefficient values for each bootstrap SLR
    boot.out[i] <- coef(boot.mod)[1]
  }
  #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
  (medianboot <-median(boot.out))
  # -0.883927
  ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
  # 2.5%      97.5% 
  #   -1.0279427 -0.7463216 
  
  
  ;  boot.CI 
}

fun.bootstr.rr.intercept.SIMPLE(df4$x, df4$y, 200)

summary(lm( (df4$y) ~df4$x))

summary(rlm( (df4$y) ~df4$x))

# fun.bootstr.rr.intercept.SIMPLE(df4$x, df4$y, 20000)

# "p" value
permuco::aovperm(  rlm( (df4$y) ~ df4$x )  )




### 2. WITH  LOG10(Y)  --------------
# confidence interval for slope:

# Bootstrap function ------------
# creating vectors to store bootstrap regression coefficients
# As fuuncion 
fun.bootstr.rr.slopeLOGY <- function (x,y, nruns) {
  
  # B = 2500 # n of runs
  B <- nruns # n of runs
  
  df <- data.frame (x = x , y = y) 
  df <- na.omit(df) # delete rows with NA
  
  boot.out <- rep(NA, B)
  vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
  ### LOOP ###
  for(i in 1:B){ #starting loop
    
    ##creating samples of observation IDs with replacement, of same size    as original sample
    boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
    
    #matching response and explanatory variable values to bootstrap sample   observation IDs
    ysam <- df$y[boot.id]
    xsam <- df$x[boot.id]
    
    #generating bootstrap SLR model for each bootstrap sample
    boot.mod <-rlm( log10( ysam) ~  xsam )       
    
    #storing regression coefficient values for each bootstrap SLR
    boot.out[i] <- coef(boot.mod)[2]
  }
  #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
  (medianboot <-median(boot.out))
  # -0.883927
  ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
  # 2.5%      97.5% 
  #   -1.0279427 -0.7463216 
  
  
  ;  boot.CI 
}

fun.bootstr.rr.slopeLOGY(df4$x, df4$y, 200)

summary(lm( log10(df4$y) ~df4$x))

# fun.bootstr.rr.slope(df4$x, df4$y, 20000)


# Boostrap 95% CI (lm):
# 2.5%        50%      97.5% 
# -0.98 -0.84 -0.71 

# OLSR CI 95% (summary, mean +_ (se*1.96)) ,
mean.slope = -0.84749 ; se.slope =  0.06547
lower.lim = mean.slope - (1.96 * se.slope);   upper.lim = mean.slope + (1.96 * se.slope);   
(lower.lim); (upper.lim)
# -0.98 to -0.72, mean = -0.85

# RLM:
# Boostrap 95% CI (rlm, nruns = 20,000):
# 2.5%          50%         97.5% 
# -1.03       -0.88       -0.74 


## for INTERCEPT
###  not LOG, SIMPLE   --------------

# confidence interval for INTERCEPT:

# Bootstrap function ------------
# creating vectors to store bootstrap regression coefficients
# As funcTion 
fun.bootstr.rr.interceptLOGY <- function (x,y, nruns) {
  
  # B = 2500 # n of runs
  B <- nruns # n of runs
  
  df <- data.frame (x = x , y = y) 
  df <- na.omit(df) # delete rows with NA
  
  boot.out <- rep(NA, B)
  vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
  ### LOOP ###
  for(i in 1:B){ #starting loop
    
    ##creating samples of observation IDs with replacement, of same size    as original sample
    boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
    
    #matching response and explanatory variable values to bootstrap sample   observation IDs
    ysam <- df$y[boot.id]
    xsam <- df$x[boot.id]
    
    #generating bootstrap SLR model for each bootstrap sample
    boot.mod <-rlm( log10( ysam) ~  xsam )       
    
    #storing regression coefficient values for each bootstrap SLR
    boot.out[i] <- coef(boot.mod)[1]
  }
  #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
  (medianboot <-median(boot.out))
  # -0.883927
  ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
  # 2.5%      97.5% 
  #   -1.0279427 -0.7463216 
  
  
  ;  boot.CI 
}

fun.bootstr.rr.interceptLOGY(df4$x, df4$y, 200)

summary(lm( log10(df4$y) ~df4$x))

summary(rlm( log10(df4$y) ~df4$x))

# fun.bootstr.rr.intercept(df4$x, df4$y, 200)

# "p" value
permuco::aovperm(  rlm( log10(df4$y) ~ df4$x )  )







# plot ALL NBSS Data, colors by communities, without cutoff lengths (loop) --------------

# NBSS plots by communities --------
# plot ALL NBSS Data (loop) --------------

library(scales)

lnumb <- length(Data.All.Atl$Longitude)
lnumb# 2840 data (lines) in C format

lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)

plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
      ylim = c(-7, 13),col = "white",
      xlim = c(-14, 3),
      xlab = "log10(gC ind.-1)",
      ylab = "log10(gC m-3 / gC ind.-1)",
      main= "NBSS, All Atlantic, All data, n = 2840, C units")


for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
  
  points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
          pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
}

# Naive linear model (all data): ------------------
# matrix to to vector (transpose), TSWA -----------------------  
library (MASS)

x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
y1 <- c ( t(Data.All.Atl.NBSSmatrix))

df1 <- data.frame(x1, y1, ylog = log10(y1))
names(df1)

# -INF to NA!!!

# df1$ylog <- replace(df1$ylog, "-Inf", NA) 

df1[df1=="-Inf"]<-NA

df2<- na.omit(df1)
dim(df1)
dim(df2)

summary(df1)
summary(df2)

mean(df2$y1)

summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.954
summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.955


abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
       lwd = 2.5, lty = 2, col = "darkorange" )


legend( "topright", 
        legend=c("Phytopl.",
                 "Net Zoopl.", "UVP zoopl.", "mesop. fish",  "micronekton"),
        fill =c("forestgreen",
              "dodgerblue", "darkmagenta", "darkred","darkorange" ),
         cex = 1.0,  box.col = "white")


#?legend



# plot ALL NBSS Data with cutoff lengths (loop) --------------

# cutoff for zooplankton ------------
cutoff_zoo <-  c(-4.5 ,-1.8)
# cutoff for UVP ------------
cutoff_zoo <-  c(-1.8)
# cutoff for micronekton and mesop. fish ------------
cutoff_zoo <-  c(-1.9, 1)

# plot ALL NBSS Data, colors by communities, without cutoff lengths (loop) --------------

# NBSS plots by communities --------
# plot ALL NBSS Data (loop) --------------

library(scales)

lnumb <- length(Data.All.Atl$Longitude)
lnumb# 2840 data (lines) in C format

lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)

plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
      ylim = c(-7, 13),col = "white",
      xlim = c(-14, 3),
      xlab = "log10(gC ind.-1)",
      ylab = "log10(gC m-3 / gC ind.-1)",
      main= "NBSS, All Atlantic, All data, n = 2840, C units")


for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
  
  points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
          pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
}

# Naive linear model (all data): ------------------
# matrix to to vector (transpose), TSWA -----------------------  
library (MASS)

x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
y1 <- c ( t(Data.All.Atl.NBSSmatrix))

df1 <- data.frame(x1, y1, ylog = log10(y1))
names(df1)

# -INF to NA!!!

# df1$ylog <- replace(df1$ylog, "-Inf", NA) 

df1[df1=="-Inf"]<-NA

df2<- na.omit(df1)
dim(df1)
dim(df2)

summary(df1)
summary(df2)

mean(df2$y1)

summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.954
summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
# slope = -0.955


abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
       lwd = 2.5, lty = 2, col = "darkorange" )





# define cutoff groups (considering taxon and size)  ------------------

Data.All.Atl$target_organisms_UVP_NET_CUTOFF <- Data.All.Atl$target_organisms_UVP_NET

Data.All.Atl$target.org.colour_CUTOFF <- Data.All.Atl$target.org.colour




# II. recalculate slopes (important for zooplankton and nekton)-------------


##############
#############
# # II.1 Phytoplankton ------------
# II.1 check slopes for phytoplankton and recalculate with Robust regression ---------

phytopl.All.Atl$new_slopes <- phytopl.All.Atl$Comment
hist(phytopl.All.Atl$NBSS_slope)


phytopl.All.Atl.NBSSmatrix <- as.matrix(phytopl.All.Atl[,37:91])
dim(phytopl.All.Atl.NBSSmatrix)
# 727 columns , 55 rows of carbon NBSS values (g C m-3 g-1C) 

sum(!is.na(phytopl.All.Atl.NBSSmatrix)) # 4847 bins, new: 7196 bins


# Calculate percent NA's per bin --------

# Convert to matrix (already done)
phytopl.Carbon.All.Atl.ok9997.NBSSmatrix <- as.matrix(phytopl.All.Atl.p.rsq.slope_ratio.ok[,37:91])
dim(phytopl.Carbon.All.Atl.ok9997.NBSSmatrix) # 311 rows

# Create index
index <- 1:ncol(phytopl.Carbon.All.Atl.ok9997.NBSSmatrix)


# Median per bin (column)
median_val <- apply(phytopl.Carbon.All.Atl.ok9997.NBSSmatrix, 2, median, na.rm = TRUE)

# Percent NA per bin (column)
percent_na <- apply(phytopl.Carbon.All.Atl.ok9997.NBSSmatrix, 2, function(x) {
  mean(is.na(x)) * 100
})

# Create bin index (columns now)
bin_index <- 1:ncol(phytopl.Carbon.All.Atl.ok9997.NBSSmatrix)

# Build table
result_table_bins311 <- data.frame(
  bin = bin_index,
  median_NBSS = median_val,
  percent_NA = percent_na
)

# Check
head(result_table_bins311)
dim(result_table_bins311)


#View(result_table_bins311) 

# bins [7] to [16]  are OK, with les thn 50% NAs! OK.

result_table_bins311

log10(X_vector_gCind[7:16]) # from -11.4389 to  -8.7305 log10(g C cell-1), n = 311 OK
length(log10(X_vector_gCind[7:16])) # 10 selected bins from -11.4 to  -8.7 log10(g C cell-1), n = 311,  OK

X_vector_gCind # 55 size bins g C ind.^-1, all size bins
X_vector_gCind_UVP <-  X_vector_gCind [X_vector_gCind > 10^-4.9] # 1.5300e-05 to 1.0240e+03
X_vector_gCind_UVP <-  X_vector_gCind_UVP [X_vector_gCind_UVP <  10^-2] # 1.5300e-05 to 1.0240e+03
X_vector_gCind_UVP # 10 classes for UVP, 0.0000153  to 0.0078130 g C ind.^-1
log10(0.0000153)
log10( 0.0078130)
X_vector_gCind[29:38] # # 10 classes for UVP, 0.0000153  to 0.0078130 g C ind.^-1

# Rationale 
# For phytoplankton NBSS, linear models were adjusted using all nano- and microphytoplankton 
# size bins in the range from -12 to -9 log10(gC cell-1), thus encompassing 10 size bins
# from -11.74 to -9.03  log10(gC cell-1) ). 
#   This most data-rich size range was selected to exclude picoplankton (which was not analysed consistently in all 
#   cruises and regions) and large-sized microplankton (which were rare in subsampling 
#   procedures and thus often had zero counts, i.e., empty bins) 
#   prior to linear model fitting. 

logCvector<- log10(X_vector_gCind)

logCvector[6] # first selcted bin,   phytoplankton
logCvector[15] # last selcted bin, phytoplankton
logCvector[6:15] #  10 selected bins, phytoplankton

# [bookm]

length(37:91)
length(X_vector_gCind)


dim(zoopl.All.Atl)
# 1579 samples, 100 variables, variables 37:91 are the NBSS data
# write.csv(zoopl.All.Atl, file = "zoopl_All_Atl_n1580_b.csv") 

names(zoopl.All.Atl [  c( 37:(37+5), (37+16):91 )] )   # 45 NBSS classes not used for phytopl NBSS (n = 10 size bins)


zoopl.All.Atl.NBSS.data <-  zoopl.All.Atl[, 37:91]

# CLEANUP - delete Phytopl outside size range -12 to -9 log10 (g C ind.^-1)
# clean all Phytopl data outside the size range

dim(phytopl.All.Atl)


phytopl.All.Atl.unfiltr.nolims727 <- phytopl.All.Atl

df <- phytopl.All.Atl


# [bookm]

# Selection of data-rich size range (nano- to microphytoplankton)


# Selected data-rich Phytopl. size range (nano- to microphytoplankton), for CLEANUP
indices_for_Phyto.cleanup <-  c( 37:(37+5), (37+16):91 ) # old and OK, gives n = 311 top linear samples, 
             # these indices give, after filtering,  a top-linear dataset of  311 samples 
# indices_for_Phyto.cleanup <- c(37:40, 52:91) # -12 to -9 logC
# indices_for_Phyto.cleanup <- c(37:38, 62:91) # test, (not used) gives > 370 top linear samples

log10(X_vector_gCind)[7]
# -11.4389 log10(gC cell-1): old and OK (gives n = 311) lower limit data-rich range (bin mid)
log10(X_vector_gCind)[16]
# -8.7305 log10(gC cell-1): old (gives n = 311) upper limit of data-rich range (bin mid)
7:17 # 11 size bins, selected data-rich range (nano- to microphytoplankton)
# selected size range:  -11.4 to -8.7 log10(gC cell-1)


for (i in indices_for_Phyto.cleanup) { 
  
  
  df[,i]    <- NA
  
  
}   

# View(df)

phytopl.All.Atl <- df

# Phytopl data finished, cleanup ok





phytopl.All.Atl$new_slopes <- phytopl.All.Atl$Comment
hist(phytopl.All.Atl$NBSS_slope)


phytopl.All.Atl.NBSSmatrix <- as.matrix(phytopl.All.Atl[,37:91])
dim(phytopl.All.Atl.NBSSmatrix)
# 727 columns , 55 rows of carbon NBSS values (g C m-3 g-1C) 






#first NBSS plots, phytopl. --------------


lnumb <- sample ( 1:300, 1)
plot (log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind) ,
      main = paste( phytopl.All.Atl$target_organisms[lnumb]  , " , ",  phytopl.All.Atl$CruiseID[lnumb] )     )

summary(mod1<-  lm(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind)))
mod1$coefficients[2]
mod1_slope <- as.numeric(mod1$coefficients[2])
abline(mod1)
phytopl.All.Atl$NBSS_slope[lnumb]
ratio_ralfvsChristina_slopes <- mod1_slope / phytopl.All.Atl$NBSS_slope[lnumb]
ratio_ralfvsChristina_slopes
# ratio = 1, slopes calculated by Cristina = slopes calclate by Ralf (lm, OLSR)
# 

# robust regression  (MASS) ---------

library(MASS)

modrr <- MASS::rlm(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind))
summary(modrr)

slope_rr <- as.numeric(modrr$coefficients[2])
intercept_rr <- as.numeric(modrr$coefficients[1])




# as function 1 (lm, OLSR, gives slope) ---------

fun.1slope.lm.olsr <- function(yvalues , xvalues) {
  mod1<-  lm( yvalues ~ xvalues)
  mod1$coefficients[2]
  mod1_slope <- as.numeric(mod1$coefficients[2])
    ; mod1_slope}


fun.1slope.lm.rr <- function(yvalues , xvalues) {
  mod1<-  rlm( yvalues ~ xvalues)
  mod1$coefficients[2]
  mod1_slope <- as.numeric(mod1$coefficients[2])
  ; mod1_slope}



fun.1bslope.lm.olsr.w.X_vec <- function(yvalues ) {
  mod1<-  lm( yvalues ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13, 1.82e-12, 3.64e-12, 7.28e-12, 1.46e-11, 2.91e-11, 5.82e-11, 1.16e-10, 2.33e-10, 4.66e-10, 9.31e-10, 1.86e-09, 3.73e-09, 7.45e-09, 1.49e-08, 2.98e-08, 5.96e-08, 1.19e-07, 2.38e-07, 4.77e-07, 9.54e-07, 1.91e-06, 3.81e-06, 7.63e-06, 1.53e-05, 3.05e-05, 6.1e-05, 0.000122, 0.000244, 0.000488, 0.000977, 0.001953, 0.003906, 0.007813, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)))
  mod1$coefficients[2]
  mod1_slope <- as.numeric(mod1$coefficients[2])
  ; mod1_slope}


fun.1bslope.lm.rr.w.X_vec <- function(yvalues ) {
  mod1<-  rlm( yvalues ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13, 1.82e-12, 3.64e-12, 7.28e-12, 1.46e-11, 2.91e-11, 5.82e-11, 1.16e-10, 2.33e-10, 4.66e-10, 9.31e-10, 1.86e-09, 3.73e-09, 7.45e-09, 1.49e-08, 2.98e-08, 5.96e-08, 1.19e-07, 2.38e-07, 4.77e-07, 9.54e-07, 1.91e-06, 3.81e-06, 7.63e-06, 1.53e-05, 3.05e-05, 6.1e-05, 0.000122, 0.000244, 0.000488, 0.000977, 0.001953, 0.003906, 0.007813, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)))
  mod1$coefficients[2]
  mod1_slope <- as.numeric(mod1$coefficients[2])
  ; mod1_slope}





lnumb <- 300
fun.1slope.lm.olsr(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))
fun.1bslope.lm.olsr.w.X_vec(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) )


fun.1slope.lm.rr(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))
fun.1bslope.lm.rr.w.X_vec(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) )


# as function 2 (gives a table with OLSR slope, OLSR intercept, OLSR p value, ----------------
                          # R-squared, AIC,  robust regression slope and intercept) ---------

fun.2.lm.olsr.rob.reg.complete.table <- function(yvalues , xvalues) {
  mod1<-  lm( yvalues ~ xvalues)
  
  
  slope_olsr <- as.numeric(mod1$coefficients[2])
  intercept_olsr <- as.numeric(mod1$coefficients[1])
  
  p_value <-  summary(mod1)$coefficients[2,4] 
  R_squared <-  summary(mod1)$r.squared
  AIC_model = AIC(mod1)
  modrr <- MASS::rlm(yvalues ~ xvalues)
  

  slope_robust_regr <- as.numeric(modrr$coefficients[2])
  intercept_robust_regr <- as.numeric(modrr$coefficients[1])
  
  
    
  mod1_table <- data.frame(   slope_olsr=slope_olsr, intercept_olsr=intercept_olsr,  p_value=p_value, 
                              R_squared = R_squared, AIC_model = AIC_model,  slope_robust_regr =slope_robust_regr , intercept_robust_regr =intercept_robust_regr     )
  
  
  ; mod1_table}



lnumb <- 300
fun.2.lm.olsr.rob.reg.complete.table(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))

res<- fun.2.lm.olsr.rob.reg.complete.table(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))

res[1] # slope of OLSR regr.

res[6] # slope of robust regr.


# apply function 1 to phytoplankton data (OLSR slopes recalculated, for cross-checking)------

out.slope.vec <- c(rep (NA, length(phytopl.All.Atl$Longitude)) )

for (i in 1: length(phytopl.All.Atl$Longitude)){

  tryCatch({
  
lnumb <- i
out.slope.vec[lnumb] <- fun.1slope.lm.olsr(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
  
}

hist(out.slope.vec)

phytopl.All.Atl$new_slopes <-out.slope.vec

phytopl.All.Atl$new_slopes / phytopl.All.Atl$NBSS_slope # all OK!cross-checking OK!
# Conclusion: cross-checking phytopl slopes OK! n = 565 useful slope data 

#  Compare OLSR slopes with robust regression and check "p" values --------------
# apply function 2 ----------------

out.p.vec <- c(rep (NA, length(phytopl.All.Atl$Longitude)) )

out.slope.rr.vec <- c(rep (NA, length(phytopl.All.Atl$Longitude)) )


for (i in 1: length(phytopl.All.Atl$Longitude)){
  
  
  tryCatch({
  
  lnumb <- i
  
  out.table <- fun.2.lm.olsr.rob.reg.complete.table(log10( phytopl.All.Atl.NBSSmatrix[lnumb,]) , log10(X_vector_gCind))

  out.p.vec[lnumb] <- out.table$p_value [1] 
  out.slope.rr.vec[lnumb] <-  out.table$slope_robust_regr[1]
  
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  

  
  }

hist(out.p.vec)
hist(out.slope.rr.vec)
hist (out.slope.vec )

out.slope.vec

summary(out.slope.vec)
summary(out.slope.rr.vec)



phytopl.All.Atl$rob_reg__slopes <- out.slope.rr.vec

phytopl.All.Atl$p_values_olsr_slopes <- out.p.vec


ratio_slopes_rr_olsr <- phytopl.All.Atl$rob_reg__slopes / phytopl.All.Atl$NBSS_slope # all OK!cross-checking OK!

hist(ratio_slopes_rr_olsr)

dim(phytopl.All.Atl)
summary(phytopl.All.Atl)
# pre-filter: 727 datasets, 163 could not be fitted any NBSS slope
727 - 163 # 564 with any kinf od NBSS slopes (pre-filter)
# USE ONKLY P < 0.05 models

phytopl.All.Atl.p.ok <- subset  (phytopl.All.Atl, p_values_olsr_slopes  < 0.05)

dim(phytopl.All.Atl.p.ok) # 311 phytopl. data sets with  p < 0.05 (depending on selected size range, this may be 372 samples)

ratio_slopes_rr_olsr_p_ok <- phytopl.All.Atl.p.ok$rob_reg__slopes / phytopl.All.Atl.p.ok$new_slopes # all OK!cross-checking OK!

hist(ratio_slopes_rr_olsr_p_ok)

ratio_slopes_rr_olsr_p_ok

summary(ratio_slopes_rr_olsr_p_ok)

hist(ratio_slopes_rr_olsr_p_ok)

which.max(ratio_slopes_rr_olsr_p_ok) # line 79, with 1.4 times higher rr slope (outlier?)
which.min(ratio_slopes_rr_olsr_p_ok) # line 143, with 0.85 times higher rr slope (outlier?) , 1.2 times higher olsr slope, OK
# acceptable bias of OLSR vs Rob_reg = 30%   (ratio < 0.67 or ratio > 1.5) 

+1/1.5

# plot

phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(phytopl.All.Atl[,37:91])
dim(phytopl.All.Atl.NBSSmatrix.1)

phytopl.All.Atl.NBSSmatrix.p.ok <- as.matrix(phytopl.All.Atl.p.ok[,37:91])
dim(phytopl.All.Atl.NBSSmatrix.p.ok)



lnumb <- 79
plot (log10( phytopl.All.Atl.NBSSmatrix.p.ok[lnumb,]) ~ log10(X_vector_gCind) ,
      main = paste( phytopl.All.Atl.p.ok$target_organisms[lnumb]  , 
                    " , st = ",  phytopl.All.Atl.p.ok$StationID[lnumb],
                    " , ",  phytopl.All.Atl.p.ok$CruiseID[lnumb] ,
                    " , line = "  , lnumb)     )

# action: delete one outlier (delete first NBSS data point of line 79)

phytopl.All.Atl.NBSSmatrix.p.ok2 <- phytopl.All.Atl.NBSSmatrix.p.ok

phytopl.All.Atl.NBSSmatrix.p.ok2[79,6] <- NA

phytopl.All.Atl.p.ok[79,6+36] <- NA

phytopl.All.Atl.p.ok[79, 37:100]
phytopl.All.Atl.NBSSmatrix.p.ok2 [79,]



lnumb <- 79
plot (log10( phytopl.All.Atl.NBSSmatrix.p.ok2[lnumb,]) ~ log10(X_vector_gCind) ,
      main = paste( phytopl.All.Atl.p.ok$target_organisms[lnumb]  , 
                    " , st = ",  phytopl.All.Atl.p.ok$StationID[lnumb],
                    " , ",  phytopl.All.Atl.p.ok$CruiseID[lnumb] ,
                    " , line = "  , lnumb)     )


# apply function 1 again to phytoplankton data (OLSR slopes recalculated, for cross-checking)------

out.slope.vec <- c(rep (NA, length(phytopl.All.Atl.p.ok$Longitude)) )

for (i in 1: length(phytopl.All.Atl.p.ok$Longitude)){
  
  lnumb <- i
  out.slope.vec[lnumb] <- fun.1slope.lm.olsr(log10( phytopl.All.Atl.NBSSmatrix.p.ok2[lnumb,]) , log10(X_vector_gCind))
}

hist(out.slope.vec)

phytopl.All.Atl.p.ok$new_slopes <-out.slope.vec

phytopl.All.Atl.p.ok$new_slopes / phytopl.All.Atl.p.ok$NBSS_slope # all OK!cross-checking OK!
# Conclusion: cross-checking phytopl slopes OK! n = 565 useful slope data 

#  Compare OLSR slopes with robust regression and check "p" values --------------
# apply function 2 ----------------

out.p.vec <- c(rep (NA, length(phytopl.All.Atl.p.ok$Longitude)) )

out.slope.rr.vec <- c(rep (NA, length(phytopl.All.Atl.p.ok$Longitude)) )

out.slope.r.sq.vec <- c(rep (NA, length(phytopl.All.Atl.p.ok$Longitude)) )



for (i in 1: length(phytopl.All.Atl.p.ok$Longitude)){
  
  lnumb <- i
  
  out.table <- fun.2.lm.olsr.rob.reg.complete.table(log10( phytopl.All.Atl.NBSSmatrix.p.ok2[lnumb,]) , log10(X_vector_gCind))
  
  out.p.vec[lnumb] <- out.table$p_value [1] 
  out.slope.rr.vec[lnumb] <-  out.table$slope_robust_regr[1]
  out.slope.r.sq.vec[lnumb] <-  out.table$R_squared [1]
  
}

hist(out.p.vec)
hist(out.slope.rr.vec)
hist (out.slope.vec )
hist (out.slope.r.sq.vec)

phytopl.All.Atl.p.ok$rob_reg__slopes <- out.slope.rr.vec

phytopl.All.Atl.p.ok$p_values_olsr_slopes <- out.p.vec

phytopl.All.Atl.p.ok$r.square.olsr <- out.slope.r.sq.vec



phytopl.All.Atl.p.ok$ratio_slopes_rr_olsr <- phytopl.All.Atl.p.ok$rob_reg__slopes / phytopl.All.Atl.p.ok$new_slopes # all OK!cross-checking OK!
ratio_slopes_rr_olsr <- phytopl.All.Atl.p.ok$ratio_slopes_rr_olsr

hist(ratio_slopes_rr_olsr)

summary(ratio_slopes_rr_olsr)

summary(phytopl.All.Atl.p.ok$p_values_olsr_slopes)
summary(phytopl.All.Atl.p.ok$r.square.olsr)


# removal of all rows where r_squared is lower than 50% -----------
dim(phytopl.All.Atl.p.ok) # 311  rows ((depending on selected size range, this may be 372 samples))
summary (phytopl.All.Atl.p.ok$ratio_slopes_rr_olsr)

# # 
  phytopl.All.Atl.p.rsqok <-  phytopl.All.Atl.p.ok 
#  phytopl.All.Atl.p.rsqok <-   subset (phytopl.All.Atl.p.ok, r.square.olsr > 0.5)
#  dim(phytopl.All.Atl.p.rsqok) # 459  rows
# 
# summary( phytopl.All.Atl.p.rsqok$ratio_slopes_rr_olsr)
 
 # removal of all rows where outlier bias is > 30%  ----------
# (slope ratio olsr/rob. reg.   has to be between  0.67  and 1.5)
phytopl.All.Atl.p.rsq.slope_ratio.ok <-  phytopl.All.Atl.p.rsqok
 phytopl.All.Atl.p.rsq.slope_ratio.ok <-   subset (phytopl.All.Atl.p.rsq.slope_ratio.ok,
                                                   ratio_slopes_rr_olsr < 1.5)

 
 
  hist(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
 hist (phytopl.All.Atl.p.rsq.slope_ratio.ok$r.square.olsr )
 hist (phytopl.All.Atl.p.rsq.slope_ratio.ok$ratio_slopes_rr_olsr)
 
 
 summary( phytopl.All.Atl.p.rsq.slope_ratio.ok$ratio_slopes_rr_olsr)
 
  dim(phytopl.All.Atl.p.rsq.slope_ratio.ok) # 311  rows 
 # 311 105
 311/545# 57% percentage of total were kept after filtering
 
 summary( phytopl.All.Atl.p.rsq.slope_ratio.ok$r.square.olsr)
 summary (phytopl.All.Atl.p.rsq.slope_ratio.ok$ratio_slopes_rr_olsr)
 summary (phytopl.All.Atl.p.rsq.slope_ratio.ok$p_values_olsr_slopes)# p: 0 to 0.037 
  summary( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # -2.2017 -1.4095 -1.2562 -1.2390 -1.0658 -0.4416 
   
  
  
#  
 dim(phytopl.All.Atl.p.rsq.slope_ratio.ok) 
 # new: 311 samples
 # OLD: 458  rows with R2 >0.5, rr/olsr < 30%, p < 0.05
 # 494  datsets with ALLL R2 , rr/olsr < 50%, ALL,  p < 0.05, 
 # in practice, xx were excluded due to bad  "p"

 
# Phytoplankton: 311 top-quality NBSS size spectra
# We omitted all linear models (delete rows), where p > 0.05, and
# # all  these 311 sampes were OK regrding outlier bias (slope ratio olsr/rob. reg.   between  0.67  and 1.5)  was above 50%
# up to one outlier per station  (obvious outliers, see supp mat), at up to three station were removed.
# Phytoplankton dataset: one obvious outlier was removed, for  ACEx/Simteco cruise station station 9 (see suppl mat, fig xx)

 
 
 phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix <- as.matrix(phytopl.All.Atl.p.rsq.slope_ratio.ok[,37:91])
 dim(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix)
 # 365 columns , 55 rows of carbon NBSS values (g C m-3 g-1C) 
 
 sum(!is.na(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix)) # 3564 bins
 

 # Create index
 index <- 1:nrow(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix)
 
 # Median per sample (row)
 median_val <- apply(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix, 1, median, na.rm = TRUE)
 
 # Percent NA per sample
 percent_na <- apply(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix, 1, function(x) {
   mean(is.na(x)) * 100
 })
 
 # Build table
 result_table365 <- data.frame(
   index = index,
   median_NBSS = median_val,
   percent_NA = percent_na
 )
 
 # View result
 head(result_table365)
 
 # Confirm size
 dim(result_table365)
 
 #View(result_table365) 
 
 
 
 
# analyse top-quality dataset, 458 top-quality NBSS size spectra

summary(phytopl.All.Atl.p.rsq.slope_ratio.ok) # n = 311 top-quality NBSS size spectra
summary(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) # n = 311 top-quality NBSS size spectra
 # n = 311 top-quality NBSS size spectra
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.1797 -1.2168 -1.0252 -1.0522 -0.8740 -0.3707 

# These linear models ( RR NBSS) yielded 311 RR NBSS slope values  that were all negative 
# and had a median that was virtually identical  to -1
# (median: -1.03, mean:  -1.05, range:    -2.18 to -0.37 ,  IQR: -1.21 to   -0.87 n = 311). 


# 
#          new_slopes      rob_reg__slopes   p_values_olsr_slopes
#       Min.   :-2.1962   Min.   :-2.2017   Min.   :0.000e+00   
#       1st Qu.:-1.3995   1st Qu.:-1.4095   1st Qu.:2.000e-08   
#       Median :-1.2558   Median :-1.2562   Median :4.034e-05   
#       Mean   :-1.2317   Mean   :-1.2390   Mean   :1.873e-03   
#       3rd Qu.:-1.0281   3rd Qu.:-1.0658   3rd Qu.:7.134e-04   
#       Max.   :-0.3708   Max.   :-0.4416   Max.   :3.722e-02

#  r.square.olsr    ratio_slopes_rr_olsr
#  Min.   :0.5010   Min.   :0.8649      
#  1st Qu.:0.7200   1st Qu.:0.9880      
#  Median :0.8210   Median :1.0015      
#  Mean   :0.8038   Mean   :1.0082      
#  3rd Qu.:0.9009   3rd Qu.:1.0261      
#  Max.   :0.9895   Max.   :1.2574 


# New summary functions with convenient outputs -----------

fun.summary_minmaxmedianmeanIQR <- function (x){
 
  first_quartile <- round (quantile(x, 0.25), 2)
  third_quartile <- round( quantile(x, 0.75), 2)
  inter.quart.text <- paste ("IQR:",  first_quartile  ,"to", third_quartile )
  
  range.text <- paste ("range:",  round(min(x),2)  ,"to", round(max(x),2) )
  
  
   head<- c( "min",   "max", "median", "mean"  )
  res <- c(  round (min(x),2) ,   round(max(x),2) , round(median    (x),2) , round( mean (x),2)  )
  
; print  (list( df <- data.frame(  head = head, res = res), 
                paste (inter.quart.text),
                paste (range.text) ))

  }

 

  summary( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
  # NEW carbon (n = 311, narrow sie range)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # -2.1797 -1.2168 -1.0252 -1.0522 -0.8740 -0.3707 
  
  # OLD carbon (all sizes...)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.2017 -1.4095 -1.2562 -1.2390 -1.0658 -0.4416 

  fun.summary_minmaxmedianmeanIQR( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
 
  # head   res
  # 1    min -2.20
  # 2    max -0.44
  # 3 median -1.26
  # 4   mean -1.24
  # 
  # [[2]]
  # [1] "IQR: -1.41 to -1.07"
  # 
  # [[3]]
  # [1] "range: -2.2 to -0.44"
  #  
  length( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope [phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope < -1])
  
  # 90  flatter than -1
  # 368 steeper than-1     # 368 with slope steer than -1 # new, OK: 372 steeper than-1(new, ok)
  # 458 top-quality NBSS size spectra
     90/ 458# 90  flatter than -1 (20%)
     368/ 458 # (80%) slopes steeper than -1
  
  
# Conclusion: cross-checking phytopl slopes OK! OLD n = 565 useful slope data 
#     OLD   458 top-quality NBSS size spectra (p <0.05, rsquared > 50%, outlier bias < 30%, only one outlier removed)
# NEW 494 RR NBSS slpes OK




# check top-quality phytoplankton NBSS (random plots) ---------

phytopl.NBSS.matrix.1 <- as.matrix(phytopl.All.Atl[,37:91])
dim(phytopl.NBSS.matrix.1)

phytopl.NBSS.matrix.ok3 <- as.matrix(phytopl.All.Atl.p.rsq.slope_ratio.ok[,37:91])
dim(phytopl.NBSS.matrix.ok3)



lnumb <- round( runif( n = 1 ,min = 1, max = 311), 0)
plot (log10(phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind),
      main = paste( phytopl.All.Atl.p.rsq.slope_ratio.ok$target_organisms[lnumb]  , ",",  
                    phytopl.All.Atl.p.rsq.slope_ratio.ok$CruiseID[lnumb], ", SST =",
                    round(phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu[lnumb],1 ) ),  
                     ylim = c(-20, 15))
      
                    # 
                    # " , st = ",  phytopl.All.Atl.p.ok$StationID[lnumb],
                    # ", line = " , lnumb  )

abline(lm(log10(phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind)), col = "red")
abline(MASS::rlm(log10(phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind)), col = "darkgreen")

# Checked, OK!

# Checked, OK!




# save complete File (phytoplankton only) with slopes and indices ---------------

 write.csv(phytopl.All.Atl.p.rsq.slope_ratio.ok, file = "outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_okv55.csv" )
# outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_ok.csv

write.csv (phytopl.All.Atl.unfiltr.nolims727 , file = "phytopl.All.Atl.unfiltr.nolims727.csv")


#phytopl.Carbon.All.Atl.ok3 <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_ok.csv")

phytopl.Carbon.All.Atl.ok3 <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_okv55.csv")

phytopl.Carbon.All.Atl.OLD.filt_458nolims <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_ok.csv")

phytopl.All.Atl.unfiltr.nolims727 <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/phytopl.All.Atl.unfiltr.nolims727.csv")


dim(phytopl.All.Atl.unfiltr.nolims727) # unfiltr., not cleaned at limits (n = 727)

phytopl.Carbon.All.Atl.OLD.filt_458nolims # old filter, not cleaned at limits (n = 458)

dim(phytopl.Carbon.All.Atl.ok3) # n = 311 filtered, cleaned at limits ( -12 to -9.4 log10(gC cell-1))



dim (phytopl.Carbon.All.Atl.ok3)
# get MLD and other paramters from Biovolume table --------------------


phytopl.Biovol.All.Atl.ok3 <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Biovol_phytopl_All_Atl_p_rsq_slope_ratio_ok_V6c.csv")

 phytopl.Carbon.All.Atl.ok3_OLD <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_ok.csv")

phytopl.Carbon.All.Atl.ok3 <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/outp_Carbon_phytopl_All_Atl_p_rsq_slope_ratio_okv55.csv")

phytopl.Biovol.All.Atl.ok3$ConsecutiveNumber <- as.character( phytopl.Biovol.All.Atl.ok3$ConsecutiveNumber)
phytopl.Carbon.All.Atl.ok3$ConsecutiveNumber <- as.character( phytopl.Carbon.All.Atl.ok3$ConsecutiveNumber)



summary(phytopl.Biovol.All.Atl.ok3)

summary(phytopl.Carbon.All.Atl.ok3)


phytopl.Carbon.All.Atl.ok3.NBSSmatrix <- as.matrix(phytopl.Carbon.All.Atl.ok3[,38:92])
dim(phytopl.Carbon.All.Atl.ok3.NBSSmatrix)
# 365 rows , 55 columns of carbon NBSS values (g C m-3 g-1C) 

length(38:92) # 55 columns

sum(!is.na(phytopl.All.Atl.p.rsq.slope_ratio.ok.NBSSmatrix)) # 3564 bins

# Check dimensions
dim(phytopl.Carbon.All.Atl.ok3.NBSSmatrix)
# 365 rows , 55 columns

# Percent non-empty per bin (column)
percent_non_empty <- colMeans(!is.na(phytopl.Carbon.All.Atl.ok3.NBSSmatrix)) * 100

# Percent empty
percent_empty <- 100 - percent_non_empty

# Index for bins
bin_index <- 1:ncol(phytopl.Carbon.All.Atl.ok3.NBSSmatrix)

# Flag column
flag <- ifelse(percent_empty >= 50, 0, 1)

# Build table
result_bins <- data.frame(
  bin_index = bin_index,
  percent_non_empty = percent_non_empty,
  percent_empty = percent_empty,
  keep_bin = flag
)

# Add bin mass vector
result_bins$X_vector_gCind_PHYTOPL <- X_vector_gCind_PHYTOPL

# Add log10 of bin mass
result_bins$log10_X_vector_gCind_PHYTOPL <- log10(X_vector_gCind_PHYTOPL)


# Check
head(result_bins)
dim(result_bins)
# View(result_bins)


log10(X_vector_gCind)
X_vector_gCind_PHYTOPL
log10(X_vector_gCind_PHYTOPL)



library(tidyverse) # merge the two tables (Carbon and Biovolume_with_MLD)

#phytopl.Carbon_andBiovol.All.Atl.ok3b <- full_join(phytopl.Biovol.All.Atl.ok3, phytopl.Carbon.All.Atl.ok3, 
  #        by = c("Longitude", "Latitude"))

#write.csv(phytopl.Carbon_andBiovol.All.Atl.ok3b, file = "phytopl_Carbon_andBiovol_All_Atl_ok3b.csv" )
#


# maps (phytoplankton) ----------


# maps --------------

plot (Data.All.Atl$Longitude ~ Data.All.Atl$Latitude)


#USING MAPS package (simple), top-quality NBSS slopes only
 library(maps)
#    maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), xlim=c(-95, 30), mar=c(0,0,0,0))
# points( phytopl.All.Atl.p.rsq.slope_ratio.ok$Longitude,
#         phytopl.All.Atl.p.rsq.slope_ratio.ok$Latitude,   col = alpha("darkgreen", 0.2), pch=16)
# 
# points( phytopl.All.Atl.unfiltr.nolims727 $Longitude,
#         phytopl.All.Atl.unfiltr.nolims727$Latitude,   col = alpha("darkgrey", 0.2), pch=16)
# 

#
# #
# #  #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
# #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
# symbols(phytopl.All.Atl.p.rsq.slope_ratio.ok$Latitude~ phytopl.All.Atl.p.rsq.slope_ratio.ok$Longitude,
#         circles = (10^(2-(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes))), fg= alpha("red", 0.6), inches=0.3,
#         add = TRUE)
# #


# N of Phytoplankton data  -----------------
#phytopl.All.Atl 

# phytopl.All.Atl.NBSSmatrix.1

dim(phytopl.All.Atl.NBSSmatrix.1)

phytopl.All.Atl.NBSSmatrix.1_no_na <-  phytopl.All.Atl.NBSSmatrix.1[rowSums(is.na(phytopl.All.Atl.NBSSmatrix.1)) != ncol(phytopl.All.Atl.NBSSmatrix.1), ]

dim(phytopl.All.Atl.NBSSmatrix.1_no_na)
# n = 545 Phytoplankton samples with data, unfiltered


# by regions (Phytoplankton, ALL data, not filteed) --------
# PHYTOPL. by region------------------

# DEFINE 4 regions, V5
# 4 key regions (where there is an exceptionally dense wealth of data)
# version v5:

# TSWA

#  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
# -14°S to -2.5°S
# -40°W to -26°W

# CCUS
# 12°N to 33°N
# -27.5°W to -9°W

#  rect( ybottom = 12 ,  ytop = 33,     xleft = -27.5 , xright =   -9 )


# EQU
#new3:
#  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
# -1°S to 6°N
# -32°W to -12°W


# BUS
# -17.5°S to -35°S
# 20°W to 10.2°W
#new 3:
#  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)



# TSWA ------------
#  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)

phytopl.unfiltr.ALL.TSWA.5 <- phytopl.All.Atl
phytopl.unfiltr.ALL.TSWA.5 <-  phytopl.All.Atl[ (phytopl.unfiltr.ALL.TSWA.5$Latitude > -14) & (phytopl.unfiltr.ALL.TSWA.5$Latitude < -2.5) &
                                (phytopl.unfiltr.ALL.TSWA.5$Longitude > -40) & (phytopl.unfiltr.ALL.TSWA.5$Longitude < -26),  ]

# symbols(phytopl.unfiltr.ALL.TSWA.5$Latitude~ phytopl.unfiltr.ALL.TSWA.5$Longitude, 
#         circles = (10^(-1* phytopl.unfiltr.ALL.TSWA.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
#         add = TRUE)

dim(phytopl.unfiltr.ALL.TSWA.5) # 110  samples!
 

# CCUS -------------

#  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )

phytopl.unfiltr.ALL.CCUS.5 <- phytopl.All.Atl
phytopl.unfiltr.ALL.CCUS.5 <-  phytopl.All.Atl[ (phytopl.unfiltr.ALL.CCUS.5$Latitude > 12) & (phytopl.unfiltr.ALL.CCUS.5$Latitude < 33) &
                                (phytopl.unfiltr.ALL.CCUS.5$Longitude > -27.5) & (phytopl.unfiltr.ALL.CCUS.5$Longitude < -9),  ]

# symbols(phytopl.unfiltr.ALL.CCUS.5$Latitude~ phytopl.unfiltr.ALL.CCUS.5$Longitude,
#         circles = (10^(-1* phytopl.unfiltr.ALL.CCUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
#         add = TRUE)
# 

dim (phytopl.unfiltr.ALL.CCUS.5) # 46 samples!

# EQU -----------
#  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )

phytopl.unfiltr.ALL.EQU.5 <- phytopl.All.Atl
phytopl.unfiltr.ALL.EQU.5 <-  phytopl.All.Atl[ (phytopl.unfiltr.ALL.EQU.5$Latitude > -1) & (phytopl.unfiltr.ALL.EQU.5$Latitude < 6) &
                               (phytopl.unfiltr.ALL.EQU.5$Longitude > -32) & (phytopl.unfiltr.ALL.EQU.5$Longitude < -12),  ]

# symbols(phytopl.unfiltr.ALL.EQU.5$Latitude~ phytopl.unfiltr.ALL.EQU.5$Longitude, 
#         circles = (10^(-1* phytopl.unfiltr.ALL.EQU.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
#         add = TRUE)

dim(phytopl.unfiltr.ALL.EQU.5)
#  29 samples


# BUS ---------------
#  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)


phytopl.unfiltr.ALL.BUS.5 <- phytopl.All.Atl
phytopl.unfiltr.ALL.BUS.5 <-  phytopl.All.Atl[ (phytopl.unfiltr.ALL.BUS.5$Latitude > -35) & (phytopl.unfiltr.ALL.BUS.5$Latitude < -17.5) &
                               (phytopl.unfiltr.ALL.BUS.5$Longitude < 20) & (phytopl.unfiltr.ALL.BUS.5$Longitude > 10.2),  ]

symbols(phytopl.unfiltr.ALL.BUS.5$Latitude~ phytopl.unfiltr.ALL.BUS.5$Longitude, 
        circles = (10^(-1* phytopl.unfiltr.ALL.BUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
        add = TRUE)

dim(phytopl.unfiltr.ALL.BUS.5)
# 59 samples



# T2a Plots All data Carbon vs Biovolume ---------------

# PLot ALL NBSS Data, Carbon  -------

nrow(phytopl.NBSS.matrix.ok3) # 311 high-quality datasets, OK

phytopl.All.Atl.unfiltr.nolims727 

# plot all points and rlm moldes (blue lines), and red median line  ----------

lnumb <- 34
plot (log10( phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind), 
      main = "ATLANTIC, n = 311",
      xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
      ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
      xlim = c(-15, -5), ylim = c(0, 14) ,
      pch = 16, col = "white")


phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727[,38:92])
dim(phytopl.NBSS.matrix.unfiltr727)


for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
{    
  points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("lightblue", 0.5))
}



for (lnumb in 1 : nrow(phytopl.NBSS.matrix.ok3) ) 
{ 
  abline( rlm(log10( phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind)), 
          col =alpha ("cornsilk4", 0.1))
}
  for (lnumb in 1 : nrow(phytopl.NBSS.matrix.ok3) ) 
  {    
  points (log10( phytopl.NBSS.matrix.ok3[lnumb,]) ~ log10(X_vector_gCind) ,
             pch = 16, col = alpha("forestgreen", 0.4))
}


# Carbon  (458 high-quality datasets, OK)
(median(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
# median RR slope : -1.256
(min(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)) 
(max(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
# min - max  RR slope : -0.44155
dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
summary(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)

dim(phytopl.NBSS.matrix.ok3)

# matrix to to vector (transpose), TSWA -----------------------  
x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.NBSS.matrix.ok3) )
y1 <- c ( t(phytopl.NBSS.matrix.ok3))

# regression line (Rob. regr.)
#plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
#points(  log10(y1) ~ x1 )
# abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
       lwd = 2.5, lty = 2, col = "darkorange" )

text(x = -12, y = 2, "n= 311, RR,  slope = -1.06", col = "darkorange"  )
text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )


#  Green dots: 311 samples and size bins filtered for the most data-rich size range (-12 to -9.4 log10 gC cell-1)
abline(v = -12)
abline(v = -9)

abline(v = -11.73, col = "red")
abline(v = -9.03, col = "red")

abline(v = log10(X_vector_gCind)[7], col = "darkgreen") # -11.4389 log(gC)
abline(v = log10(X_vector_gCind)[15], col = "darkgreen") # -9.03 log(gC)

abline(v = log10(X_vector_gCind)[7], col = "purple", lwd = 3) # -9.03 log(gC)
abline(v = log10(X_vector_gCind)[16], col = "purple",  lwd = 3) # -9.03 log(gC)

log10(X_vector_gCind)[7]
# -11.4389 log1-gC cell-1: old (gives n = 311) lower limit data-rich range (bin mid)
log10(X_vector_gCind)[16]
# -8.7305 log1-gC cell-1: old (gives n = 311) upper limit of data-rich range (bin mid)
7:17 # 11 size bins, secleted data-rich range (nano- to microphytoplankton)



# [bookm]

# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

mean_   <- function(...) mean(..., na.rm=T)
median_ <- function(...) median(..., na.rm=T)

mean_ (c(4,NA, 8, 34)) #15.33, OK
median_ (c(4,NA, 8,34)) #8, OK


# calculate median (or mean) only if less than 50% are NA ------------

vec1 <- c(4, 5,NA,NA, NA,  34)

  sum(is.na(vec1))
percNA <- sum(is.na(vec1))/ length(vec1)
  
  
mean_NA50   <- function(...)  {

  percNA <- sum(is.na(...))/ length(...)
 
  if(percNA <= 0.5) {
  mean(..., na.rm=T) }
  
  else{ ;NA }
}

median_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    median(..., na.rm=T) }
  
  else{ ;NA }
}

vec2 <- c(4, 5,NA,NA, NA,NA,  34)
vec3 <- c(4, 5,NA,34)

mean_NA50(vec2)
mean_NA50(vec3)


# apply median by columns

median_vec_ALLBiovol <- apply( phytopl.NBSS.matrix.ok3, 2, median_) 
mean_vec_ALLBiovol <- apply(phytopl.NBSS.matrix.ok3, 2, mean_) 




# lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)


# apply mean by columns (only if less thn 50% NAs)

# mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 

median_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, median_NA50) 

fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
  
max_vec_ALL <- apply(phytopl.NBSS.matrix.ok3, 2, fun95_quant ) 

length(median_vec_ALLBiovol2)

class(median_vec_ALLBiovol2)

median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))

#mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA



 median_vec_PHYTOP <-  median_vec_ALLBiovol2 

 X_vector_gCind_PHYTOPL <- X_vector_gCind


df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))

dfok<- df1
#dfok <- df1[complete.cases(df1$y1),]
#plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)

dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
        col = "purple" , lwd = 4.5)


abline ( v =  -11.15) # nanoplankton peak
# rlm based on medians ony

# abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
#               lwd = 2.5, lty = 2, col = "darkgreen")
#         
# summary(rlm_med)
# # slope:  -0.9629 (based on medians < 50% NAs )

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)

rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1))
summary(rlm_med)

  


      # col = "darkmagenta" , lwd = 2.5)


### Calculate median with ALL data (not only the ones filtered for linearity)


phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(phytopl.All.Atl[,38:92])
dim(phytopl.All.Atl.NBSSmatrix.1)

# PLot ALL NBSS Data, Carbon  -------

nrow(phytopl.NBSS.matrix.ok3) # 458 high-quality datasets, OK

# plot all points and rlm moldes (blue lines), and red median line  ----------

lnumb <- 34
plot (log10( phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind), 
      main = "ATLANTIC, n = 727",
      xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
      ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
      xlim = c(-15, -5), ylim = c(0, 14) ,
      pch = 16, col = "white")

for (lnumb in 1 : nrow(phytopl.NBSS.matrix.1) ) 
{    
  points (log10( phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("darkgrey", 0.2))
}

# Carbon  (458 high-quality datasets, OK)
(median(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
# median RR slope : -1.256
(min(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)) 
(max(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
# min - max  RR slope : -0.44155

dim(phytopl.NBSS.matrix.1)

# matrix to to vector (transpose), TSWA -----------------------  
x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.NBSS.matrix.1) )
y1 <- c ( t(phytopl.NBSS.matrix.1))

# regression line (Rob. regr.)
#plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
#points(  log10(y1) ~ x1 )
# abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
       lwd = 2.5, lty = 2, col = "darkorange" )
summary(rlm2)
# slope = -1.146 , using all green points, rlm
# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

# apply median by columns

median_vec_ALLBiovol <- apply( phytopl.NBSS.matrix.1, 2, median_) 
mean_vec_ALLBiovol <- apply(phytopl.NBSS.matrix.1, 2, mean_) 

# lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)


# apply mean by columns (only if less thn 50% NAs)

# mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 

median_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.1, 2, median_NA50) 

dim(median_vec_ALLBiovol2)

class(median_vec_ALLBiovol2)

median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))

#mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA

df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))

dfok<- df1
#dfok <- df1[complete.cases(df1$y1),]
#plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)


# rlm based on medians only

abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
        lwd = 2.5, lty = 2, col = "darkgreen")

summary(rlm_med)
# OLD slope: -1.1573 (based on medians < 50% NAs )
# slope: -0.7802, based on medians (-12 to -9)
lines(  log10(dfok$y1) ~ log10(dfok$x1) )

   text(-9, 13.5, "green line: rlm based on medians only, n = 727 samples")   
   text(-9, 13.5, "green line: rlm based on all 727 samples")   
   
         
#                col = "red" , lwd = 4.5)

# col = "darkmagenta" , lwd = 2.5)

   
   
#### Temperature and chla vs phytoplankton NBSS slope, plots and models ----------------

   
#  [bookmark]
   
### For PAPER ---------- for paper ------------ interaction plots -------

#1a.   THE BEST Best univarate model,log10(1+ Chl a in situ ) ----------
lmlogChla_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                         log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
# CARBON (new):  R-squared: 0.0285,  n = 311, p-value: 0.0029
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       log10(1+ phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu), 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")

dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)


cor.test (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ,
            log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu))
# Pearson, p-value = 0.002865


cor.test (method = "spearman", phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ,
          log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu))
# spearman, p-value = 0.01624

lmlogChla_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                         log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
rlmlogChla_insitu <- rlm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                           log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )

summary(rlmlogChla_insitu) #
summary(lmlogChla_insitu) #
#R-squared   lm:  0.0285

library(permuco)
 aovperm(rlmlogChla_insitu ,np = 2000)
#aovperm(rlmlogChla_insitu ,np = 20000)
# p = 0.0026

 
 df <- data.frame (x = log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) , 
                   y = phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) 
 df <- na.omit(df) # delete rows with NA
 
 # Calculate R-sq (rlm)
 y <- df$y; x <- df$x
 y_pred.rlm1 <- predict(rlm1 <- rlm(y ~ x))
 rss <- sum((y - y_pred.rlm1)^2)
 tss <- sum((y - mean(y))^2)
 (R_squared.rlm1 <- 1 - rss / tss)
 # #R-squared   rlm:  0.02515
 

# For paper

library(ggplot2)
df2 <- data.frame( NBSS_slope = phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes , 
                   Chla = log10( phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu))

ggplot(df2, aes(Chla, NBSS_slope ) ) +
  scale_x_continuous(name ="Chl a (log10 (mg m-3)), in situ", label = comma,
                     limits=c(-1.6, 1.6))+
  theme(axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10) ,
        axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10) )+ 
  geom_point( alpha = 0.3, col = "navy", size = 2.5) +
  geom_smooth(method = "lm")


log10(1)   # "0"  represnnts 1 mg m-1 on a log-scale


summary(m109 <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope ~ 
       log10( phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu )))
# p = 0.00481, R² = 0.0255,  only 2.5% of the variabilty is explained!
  
summary(rlm109 <- rlm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope ~ 
                     log10( phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu )))


cor.test( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope , 
          log10(phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu), method = "spearman")
#  sign for Chla , p-value = 0.01624, Spearman correlation

aovperm(rlm109)  
#  p = 0.005, R² = 0.0256 ,only 2.6% of the variabilty is explained!


cor.test( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope , 
          phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu, method = "spearman")
#  sign for Chla , p-value = 0.01624, Spearman correlation

#  cor.test( phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope , 
#            phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu, method = "spearman")
# # not sign for SST
#  
#  # p-value = 0.02388, weak but significant relationship
# # slighly flatter (more large-sized organisms) in colder waters


# summary( lm(zoopl_All_Atl4UVP$rob_reg__slope ~ zoopl_All_Atl4UVP$SST_C_insitu))




# 1b. phytopl slopes vs  Chl a (model)  not sign.   ---------- 

summary(lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
             log10(phytopl.All.Atl.p.rsq.slope_ratio.ok$chla)))
        
cor.test(method = "spearman", phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ,
           log10(phytopl.All.Atl.p.rsq.slope_ratio.ok$chla))    
# p-value = 0.0657, not sign.



#1b.  SST  in situ  not sign.   ----------

summary( phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    5.50   15.60   22.34   21.26   27.30   34.62      13 
   
lmSST_C_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                    phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu)
summary(lmSST_C_insitu) #OLD (biovolume , in situ SST data) R-squared:  0.349, 460 degrees of freedom, p-value: < 2.2e-16
# Carbon (new, insitu SST Data): R-squared:  0.338, 443 degrees of freedom, p-value: < 2.2e-16


plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmSST_C_insitu, lwd = 2.5, lty = 2, col = "darkorange")



#1c.  second best univarate model, SST in situ----------
lmSST_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                     phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu)
summary(lmSST_insitu) ##OLD (biovolume , model SST data)  R-squared:  0.347, 421 degrees of freedom, p-value: < 2.2e-16
# Carbon (new, insitu SST Data): R-squared:  0.338, 443 degrees of freedom, p-value: < 2.2e-16


plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmSST_C_insitu, lwd = 2.5, lty = 2, col = "darkorange")





#1d.  3rd univariate model, Temp  in situ at depth,  Not Sign.  ----------
lmTemp_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                      phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu)
summary(lmTemp_insitu) # R-squared:  0.26, 421 degrees of freedom, p-value: < 2.2e-16
# new (carbon): Multiple R-squared:  0.3056,  417 degrees of freedom, , p-value: < 2.2e-16
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmTemp_insitu, lwd = 2.5, lty = 2, col = "darkorange")




#2.   THE 3rd  worst univarate model,salinity) ----------
lm.sal <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
               (phytopl.All.Atl.p.rsq.slope_ratio.ok$salinity) )
summary(lm.sal) # R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
#  R-squared:  0.155, 456 DF,  p-value: < 2.2e-16


#3.   THE second worst univarate model NITRATE situ ) ----------

#univarate model NITRATE in situ
lm.nitr <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                ((phytopl.All.Atl.p.rsq.slope_ratio.ok$Nitrate__uM_L_1_insitu)) )
summary(lm.nitr) # (OLD, biovolume): R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
#NEW, Carbon:  R-squared:  0.02465,  298 DF,  p-value: 0.00643

#univarate model NITRATE Model
lm.nitr2 <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                (phytopl.All.Atl.p.rsq.slope_ratio.ok$nitrate) )
summary(lm.nitr2) #
#Carbon, new,   R-squared:  0.0475, 456 DF,  p-value: 2.507e-06



#4.   THE worst (!) univarate model,MLD ----------
# lm.mld <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
#                (phytopl.All.Atl.p.rsq.slope_ratio.ok$MLD_Mathilde) )
#summary(lm.mld) #OLD (biovolume): R-squared:  0.03, 349 degrees of freedom, p-value: 0.0008



### PLOTs For PAPER ----------

# Chlorophyll__mg_m_3_insitu
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       log10(1+ phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu), 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")



# Temperature__C_insitu
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmTemp_insitu, lwd = 2.5, lty = 2, col = "darkorange")



# robust regression (RR) univariate models -------
# Short name for dataframe: "PP.ATL.df4" (all phytopl. slopes, n = 458 linear NBSS datasets)
PP.ATL.df4 <- phytopl.All.Atl.p.rsq.slope_ratio.ok
dim (PP.ATL.df4)# 458  linear-shaped phytoplankton NBSS datasets(p < 0.05, R-sq > 0.5, rr/olsr < 1.3 )

library(MASS) # rlm # works OK, bu no "p" values, so I  had to do "a posteriori" permutation test
library(robustbase) # lmrob gives very low "p" values, based on M estimator... not reliable!
library(permuco)
#permuco::lmperm()
#permuco::aovperm()

# Chlorophyll__mg_m_3_insitu
plot(PP.ATL.df4$rob_reg__slopes ~ 
       log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu), 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
abline(rlm(PP.ATL.df4$rob_reg__slopes ~ 
             log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu)) 
  , col = "red")
abline(lmrob(PP.ATL.df4$rob_reg__slopes ~ 
             log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu)) 
       , col = "green")
abline(lmperm(PP.ATL.df4$rob_reg__slopes ~ 
                log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu)) 
       , col = "darkmagenta")


# all four regression lines are 100% identical! (lm, rlm, lmrob, lmperm)
# for a large sample, OK

# ?permuco::aovperm()
aovperm(np = 10000, rlm(PP.ATL.df4$rob_reg__slopes ~ 
              log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu)))



# Temperature__C_insitu
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmTemp_insitu, lwd = 2.5, lty = 2, col = "darkorange")
abline(rlm(PP.ATL.df4$rob_reg__slopes ~ 
             (PP.ATL.df4$Temperature__C_insitu)) 
       , col = "red")
abline(lmrob(PP.ATL.df4$rob_reg__slopes ~ 
               ( PP.ATL.df4$Temperature__C_insitu)) 
       , col = "green")
# all three regression lines are identical! (lm, rlm, lmrob)

aovperm(np = 10000, rlm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                          phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu))



# small test dataset ----------
# Temperature__C_insitu
plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes[3:25] ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu[3:25], 
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lmTemp_insitu,  col = "red")
abline(rlm(PP.ATL.df4$rob_reg__slopes[3:25] ~ 
             (PP.ATL.df4$Temperature__C_insitu[3:25])) 
       , lty = 2, lwd = 2, col = "darkorange"  )
abline(lmrob(PP.ATL.df4$rob_reg__slopes[3:25] ~ 
               ( PP.ATL.df4$Temperature__C_insitu[3:25])) 
       , col = "green")
abline(lmperm(PP.ATL.df4$rob_reg__slopes[3:25] ~ 
               ( PP.ATL.df4$Temperature__C_insitu[3:25])) 
       , col = "darkmagenta")

# 23 datapoints only, 21 degrees of freedom, very small dataset for testing 


# test the significance of the rlm model with permuco::aovperm !!!
# ?aovperm # {permuco}
aovperm(rlm(PP.ATL.df4$rob_reg__slopes[3:25] ~ 
              log10(1+ PP.ATL.df4$Temperature__C_insitu[3:25])),
        np = 50000) # p = 0.03, oK
summary(lm(PP.ATL.df4$rob_reg__slopes[3:25] ~ 
         log10(1+ PP.ATL.df4$Temperature__C_insitu[3:25]))) # p = 0.04, ok



# The  four regression lines are NOT identical! (lm, rlm, lmrob, lmperm)
# orange: lm
# red: rlm
# green: lmrob
# magenta: lmperm
# darkblue:  lmPerm:: lmp)

#  bootcoefs {complmrob} ----------------
# Bootstrap the regression coefficients for a robust linear regression model ---------
# based on lmrob (robustbase) outputs ---------
# library(complmrob)
# ? bootcoefs



# Multivarariate models (Phytoplankton) ----------------

# Best Multivariate models ----------------

# NBSS slope vs SST_C_insitu and  Chlorophyll__mg_m_3_insitu -----------

# Short name for dataframe: "PP.ATL.df4" (all phytopl. slopes, n = 458 linear NBSS datasets)
PP.ATL.df4 <- phytopl.All.Atl.p.rsq.slope_ratio.ok
dim (PP.ATL.df4)# 458  linear-shaped phytoplankton NBSS datasets(p < 0.05, R-sq > 0.5, rr/olsr < 1.3 )


lm4a.SSTno.int <-  lm(PP.ATL.df4$rob_reg__slopes ~ 
                     PP.ATL.df4$SST_C_insitu +  
                    log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu))
summary(lm4a.SSTno.int)
# Interaction is Not .Sign.! ----------------
# Chla is  sign. ! p <0.0001 ------------------------
#  "SST in situ"  is signif! p <0.0001 -----------------
# NO INTERACTION!
# BOTH VARIABLES are highly signif.!!! p <0.0001

aovperm( np = 10000,
rlm(PP.ATL.df4$rob_reg__slopes ~ 
     PP.ATL.df4$SST_C_insitu +  
      log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu))) #OK
# BOTH VARIABLES are highly signif.!!! p <0.0001, OK
# NO INTERACTION!


#  test for interaction terms -----------

lm4a.SSTw.int <-  lm(PP.ATL.df4$rob_reg__slopes ~ 
             PP.ATL.df4$SST_C_insitu *  
               log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu))
summary(lm4a.SSTw.int)
# Chla is highly sign. ! ------------------------
# Only "SST in situ"  is signif! -----------------
# There is NO INTERACTION!


aovperm( 
  rlm(PP.ATL.df4$rob_reg__slopes ~ 
        PP.ATL.df4$SST_C_insitu *  
        log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu))) #OK
# BOTH VARIABLES are highly signif.!!! p <0.0001, OK
# NO INTERACTION!, OK

# Test for Collinearity (reject if r > 0.8) ----------------------

cor.test( PP.ATL.df4$SST_C_insitu , log10(1+ PP.ATL.df4$Chlorophyll__mg_m_3_insitu) )
# -0.5153801 

# relative importance within TOP Model ---------
library(relaimpo)

relaimpo::calc.relimp(lm4a.SSTno.int) # two variable TOP model (SST, chl a )
# Phytoplankton 



# NBSS slope vs Temperature__C_insitu and  Chlorophyll__mg_m_3_insitu -----------

PP.ATL.df4 <- phytopl.All.Atl.p.rsq.slope_ratio.ok

# model with no inteaction 
lm4b.no.int <-  lm(PP.ATL.df4$rob_reg__slopes ~ 
                    PP.ATL.df4$Temperature__C_insitu +  
                    PP.ATL.df4$Chlorophyll__mg_m_3_insitu)
summary(lm4b.no.int)
# Interaction is Not .Sign.! ----------------
# BOTH VARIABLES are highly signif.!!! p <0.0001, OK
# NO INTERACTION!


aovperm( 
  rlm(PP.ATL.df4$rob_reg__slopes ~ 
        PP.ATL.df4$Temperature__C_insitu +  
        PP.ATL.df4$Chlorophyll__mg_m_3_insitu)) #OK
# BOTH VARIABLES are highly signif.!!! p <0.0001, OK
# NO INTERACTION!, OK


lm4b.w.int <-  lm(PP.ATL.df4$rob_reg__slopes ~ 
             PP.ATL.df4$Temperature__C_insitu *  
             PP.ATL.df4$Chlorophyll__mg_m_3_insitu)
summary(lm4b.w.int)
# Interaction is Not .Sign.! ----------------
# Chla is not sign. ! ------------------------
# Only "Temp in situ" is signif! -----------------
# interaction term is not signif. 

# Phytoplankton:
# BIOVOLUME:INTERACTION SST vs Chla IS SIGNIFICANT for BIOVOLUME NSS , but not for carbon NBSS!!!
# CARBON: INTERACTION SST vs Chla IS NOT SIGNIFICANT for Carbon NBSS 

aovperm( 
  rlm(PP.ATL.df4$rob_reg__slopes ~ 
        PP.ATL.df4$Temperature__C_insitu *  
        PP.ATL.df4$Chlorophyll__mg_m_3_insitu)) #OK
# Interaction is Not .Sign.! ----------------
# Chla is not sign. ! ------------------------
# Only "Temp in situ" is signif! -----------------

# 3d plots --------------------

library(car)
scatter3d(x = PP.ATL.df4$SST_C_insitu, y =  PP.ATL.df4$rob_reg__slopes, 
          z = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu), 
          grid = FALSE, fit = "smooth", point.col="darkorange") # bg = "black")

# 3d surface is not twisted! NO visible interactions for Carbon NBSS


library(car)
scatter3d(x = PP.ATL.df4$Temperature__C_insitu, y =  PP.ATL.df4$rob_reg__slopes, 
          z = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu), 
          grid = FALSE, fit = "smooth", point.col="darkorange") # bg = "black")

# 3d surface is not twisted! NO visible interactions for Carbon NBSS



# 3D plot with coluors

# library(scatterplot3d)
library(plot3D) # colured 3D plots

x = PP.ATL.df4$Temperature__C_insitu
z =  PP.ATL.df4$rob_reg__slopes
y = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu)
        
  ramp.col

# plot3D::scatter3D(x, y, z,  
#           col = ramp.col(c("darkgreen", "khaki", "darkred")), bty="b2", 
#           theta = 595, phi = 35, ticktype = "detailed", d=200,
#           xlab = "", ylab = "", zlab = "",  
#           surf = list(x = x.pred, y = y.pred, z = z.pred,  
#                       facets=T, border="black"), main = "Chl_a in situ, Temp. in situ")
# 
# points3D(x, y, z, pch = 16, color="black", alpha = 0.8, add=TRUE)
# 



### T-S Digram, t-s


# T-S diagram for Paper ! ---------------------------
# T-S diagram with circle size = Chla  ------------------
# with GGPLOT!

library(ggplot2)

# with chla in situ (ok, but only 313 good datasets)
df0insitu <- data.frame ( Chla_in_situ =   Data.All.Atl$Chlorophyll__mg_m_3_insitu,
                          Salinity = Data.All.Atl$Salinity_insitu,
                          SST =  Data.All.Atl$SST_C_insitu  ) 
dim(df0insitu)# N = 2840   datasets (with NAs)
df0insitu <- df0insitu[complete.cases(df0insitu),]
dim(df0insitu)
# N = 313   datasets (without NAs)

summary(df0insitu)

# with chla from model (bad...), with # 575 complete datsets
df0 <- data.frame (Chla =   Data.All.Atl$chla,
                   Salinity = Data.All.Atl$Salinity_insitu,
                   SST =  Data.All.Atl$SST_C_insitu  ) 
dim(df0)# N = 2840   datasets (with NAs)
summary(df0)
df0 <- df0[complete.cases(df0),]
dim(df0) # 575 complete datsets
summary(df0)

summary(Data.All.Atl$chla)

summary(as.numeric(Data.All.Atl$Chlorophyll__mg_m_3_insitu))


#cor.test (Chlorophyll__mg_m_3_insitu, chla_model) 

#plot(Chlorophyll__mg_m_3_insitu ~ chla_model) 

#abline(rlm(Chlorophyll__mg_m_3_insitu ~ chla_model ) , col = "darkgreen")



# FOR PAPER --------------
library(ggplot2)
library(scales)
library(car)
scatter3d(x = PP.ATL.df4$Temperature__C_insitu, y =  PP.ATL.df4$rob_reg__slopes, 
          z = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu), 
          grid = FALSE, fit = "smooth", point.col="darkorange") # bg = "black")

# 3d surface is not twisted! NO visible interactions for Carbon NBSS

# colour plot with ggplot2 ----------------

summary(PP.ATL.df4$Salinity_insitu)
summary(Chlorophyll__mg_m_3_insitu)
summary(SST_C_insitu)


df4insitu <- data.frame ( SST = PP.ATL.df4$SST_C_insitu,
                          Sal =  PP.ATL.df4$Salinity_insitu,
                          Temp = PP.ATL.df4$Temperature__C_insitu,
                          NBSS_slopes =  PP.ATL.df4$rob_reg__slopes, 
                         logChl_a = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu))

dim(df4insitu)# N = 458   datasets (with NAs)
summary(df4insitu)
df4insitu <- df4insitu[complete.cases(df4insitu),]
dim(df4insitu) # 205 complete datsets (no NAs in any in situ paramers)
summary(df4insitu)
# N = 205   datasets (without NAs)
df4insitu <- sort_by(df4insitu, df4insitu$SST )


b <- ggplot(df4insitu, aes(x = logChl_a, y = SST)) #+ xlim(31.5, 38) + ylim( 5,30 )
# Basic scatter plot
b+ geom_point() 

b <- ggplot(df4insitu, aes(x = logChl_a, y = SST)) #+ xlim(31.5, 38) + ylim( 5,30 )
# Basic scatter plot
b+ geom_point(aes(colour = df4insitu$NBSS_slopes))

# T-S diagram, points with colours, FOR PAPER ! ----------------
b + geom_point(aes(color = (NBSS_slopes)) , size = 3) +
  scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme(legend.position = "right")

library(scales)
colors3 = c( (alpha("#FC4E07", alpha = 0.8)) ,
             (alpha("#E7B800", alpha = 0.6)),
              (alpha("#00AFBB", alpha = 0.4)),
             (alpha("forestgreen", alpha = 0.4)), 
             (alpha("darkgreen", alpha = 0.4)) )
  
colors4 = c(  (alpha("darkorange", alpha = 0.6)),
             (alpha("#E7B800", alpha = 0.6)),
              (alpha("yellowgreen", alpha = 0.6)),
             (alpha("forestgreen", alpha = 0.4)), 
             (alpha("darkgreen", alpha = 0.4)),
             (alpha("#00AFBB", alpha = 0.4)),
             (alpha("blue", alpha = 0.4)),
             (alpha("darkblue", alpha = 0.4)),
             (alpha("purple", alpha = 0.4)))
             

colors4b = c(  (alpha("darkorange", alpha = 0.9)),
               (alpha("#E7B800", alpha = 0.9)),
               (alpha("yellowgreen", alpha = 0.9)),
               (alpha("forestgreen", alpha = 0.9)), 
               (alpha("darkgreen", alpha = 0.9)),
               (alpha("#00AFBB", alpha = 0.9)),
               (alpha("blue", alpha = 0.9)),
               (alpha("darkblue", alpha = 0.9)),
               (alpha("purple", alpha = 0.9)))


colors4c = c(  (alpha("darkred", alpha = 0.9)),
              (alpha("darkorange", alpha = 0.9)),
              (alpha("#E7B800", alpha = 0.9)),
              (alpha("yellowgreen", alpha = 0.9)),
              (alpha("forestgreen", alpha = 0.9)), 
              (alpha("darkgreen", alpha = 0.9)),
              (alpha("#00AFBB", alpha = 0.9)),
              (alpha("blue", alpha = 0.9)),
              (alpha("darkblue", alpha = 0.9)),
              (alpha("purple", alpha = 0.9)))



colors5 = c( (alpha("darkred", alpha = 0.8)),
              (alpha("darkorange", alpha = 0.7)),
             (alpha("orange2", alpha = 0.6)),
               (alpha("orange", alpha = 0.6)),
                     (alpha("lightsalmon", alpha = 0.6)),
             (alpha("yellow2", alpha = 0.6)),
              (alpha("yellowgreen", alpha = 0.6)),
              (alpha("forestgreen", alpha = 0.6)), 
               (alpha("darkgreen", alpha = 0.5)),
              (alpha("navyblue", alpha = 0.5)),
                 (alpha("navyblue", alpha = 0.5)))
              



# nice plots for paper! ---------              
# T-S diagram, points with colours, FOR PAPER ! ----------------
b <- ggplot(df4insitu, aes(x = logChl_a, y = SST)) #+ xlim(31.5, 38) + ylim( 5,30 )
b + geom_point(aes(color = (NBSS_slopes)) , size = 2) +
  scale_color_gradientn(colors = colors3) +
  theme(legend.position = "right")

##
b + geom_point(aes(color = (NBSS_slopes)) , size = 2) +
  scale_color_gradientn(colors = colors4) +
  theme(legend.position = "right")

#############
## T-S diagram --- points with colours, FINAL PLOT FOR PAPER !!!! ----------------
library(RColorBrewer)
library(wesanderson)


names(df4insitu)

c <- ggplot(df4insitu, aes(x = Sal, y = SST))+ #+ xlim(31.5, 38) + ylim( 5,30 )
 geom_point(aes(color = logChl_a) , size = 3) +
  scale_color_gradientn(colors = colors4b) +
  theme(legend.position = "right")




##
#  plots for paper OK, FOR PAPER---------              
## NBSS_slopes
c1 <- ggplot(df4insitu, aes(x = logChl_a, y = SST)) #+ xlim(31.5, 38) + ylim( 5,30 )

c1 + geom_point(aes(color = NBSS_slopes) , size = 3) +
  scale_color_gradientn(colors = colors5) +
  theme(legend.position = "right")

#[bookm]
### FINAL PLOT FOR PAPER ---------------------------------
#  plot for paper OK, FOR PAPER ---------              
## NBSS_slopes
c2 <- ggplot(df4insitu, aes(x = logChl_a, y = Temp)) #+ xlim(31.5, 38) + ylim( 5,30 )


c2 + geom_point(aes(color = NBSS_slopes) , size = 3) +
  scale_color_gradientn(colors = colors5) +
  theme(legend.position = "right")

# add rectangles---------------
# temperatures for classes: 3 to 20, 20 to 27.8, < 27.8 (27.8 to 31)  ----------

#?geom_rect
 # c2 + 
 #   geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),  alpha=0.1, inherit.aes = FALSE) +
 #   geom_rect(data=t2.rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes = FALSE) +
 #   geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.1, inherit.aes = FALSE)
 # 


##############

# symbol plots (for rectangles) ----------


symbols(Data.All.Atl$SST_C_insitu ~ Data.All.Atl$Salinity_insitu, main = "circles = Chla, in situ", 
        circles = (Data.All.Atl$Chlorophyll__mg_m_3_insitu), fg= alpha("navyblue", 0.6), inches=0.2,
        xlim = c(31,38),
        ylim = c(0, 32))
# rect(31.5, 10, 38, 27)
# rect(36.9, 27, 38, 32)
# rect(31.5, 3, 38, 10)


# limit at 20 degress Celsius
rect(36, 27.8, 38, 32) # TW, oligotrophic tropical water (high salinity high temp, low nutr, low chl a)
rect(34.5, 27.8, 36, 32) # Warm Coastal waters, warm, nutrient-rich waters (e.g., productive coastal waters in summer") 
rect(34.5, 20, 38, 27.8)  # mixed warm waters 
rect(31.5, 3, 34, 20) # upwelling, extremely productive, cold waters, temp < 20 deg. celsius, salinity < 34, upwelling
rect(34, 3, 36, 20) # productive, cold waters, temp < 20 deg. celsius, salinity 34 to 36, upwelling
rect(36, 3, 37, 20) # less productive cold waters, salinity > 36

# temperatures for classes: 3 to 20, 20 to 27.8, < 27.8 (27.8 to 31)  ----------


#### By Chla. classes -------------------

#############
# three Chla groups FOR PAPER -----------------------
# classified by  Chla  in situ !

# Separate OLigo-, meso-, and  eutrophic samples --------
# boundaries: 0.2 and 1 mg Chla m-3 in situ

# find boundaries -------------------
length(phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 727 datasets
# 458 total datasets



oligotr.phytopl.All.Atl <- subset(phytopl.All.Atl.p.ok, Chlorophyll__mg_m_3_insitu < 0.2 )
summary(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
length(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 88 datasets


mesotr.phytopl.All.Atl <- subset(phytopl.All.Atl.p.ok, Chlorophyll__mg_m_3_insitu > 0.2 )
mesotr.phytopl.All.Atl <- subset(mesotr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 2 )#>2 mg m3 Chla 
summary(mesotr.phytopl.All.Atl$SST_C_insitu) # 
length(mesotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
# 301 datasets


eutr.phytopl.All.Atl <-subset(phytopl.All.Atl.p.ok, Chlorophyll__mg_m_3_insitu > 2 )
eutr.phytopl.All.Atl <-subset(eutr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 31 )
summary(eutr.phytopl.All.Atl$SST_C_insitu) # < 31 mg 
length(eutr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 88 datasets
# 55 datasets












# Backwards selection of variables  
# full models ----------------


# Finding the best linear model
# "Selecting a subset of predictor variables from a larger set 
# (e.g., stepwise selection) is a controversial topic. 
# You can perform stepwise selection (forward, backward, both) 
# using the stepAIC( ) function from the MASS package.
# stepAIC( ) performs stepwise model selection by exact AIC. 
# "forward" method does not work here (retains all variables)
# same results are obtained for "backward" and "both" methods", OK. 


Data.lin.mod.insitu <-   cbind( PP.ATL.df4[,c(30, 31, 97, 102)], PP.ATL.df4$Nitrate__uM_L_1_insitu)
#names(PP.ATL.df4) #six in_situ variables and rob_reg__slopes"
names(Data.lin.mod.insitu)
names(Data.lin.mod.insitu) <- c( "SST_C_insitu"   ,     "Temperature__C_insitu" ,           
                                   "log10_1X_chla_insitu" ,   "rob_reg__slopes"  ,                
                                  "Nitrate_insitu" )


# TEST FOR COLLINEARITY ----------------

names(Data.lin.mod.insitu)
cor.test(Data.lin.mod.insitu$Temperature__C_insitu , Data.lin.mod.insitu$SST_C_insitu)
# r = 0.92 , collineaity, choose one variable

cor.test(PP.ATL.df4$nitrate , PP.ATL.df4$Nitrate__uM_L_1_insitu)
# r = 0.61 ,  no collineaity (! strange!), but however,  choose one variable (in situ)

cor.test(PP.ATL.df4$SST_C_insitu , PP.ATL.df4$Temperature__C_insitu)
# r = 0.92 ,  with collineaity  choose one variable (Temp in situ) !!!

cor.test(PP.ATL.df4$Chlorophyll__mg_m_3_insitu , PP.ATL.df4$Temperature__C_insitu)
# r = -0.33 , no collinearity , OK!

cor.test(PP.ATL.df4$Chlorophyll__mg_m_3_insitu , PP.ATL.df4$SST_C_insitu)
# r = -0.39 , no collinearity , OK!

cor.test(log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu) , PP.ATL.df4$SST_C_insitu)
# r = -0.51 ,  no collinearity , OK!

cor.test(log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu) , PP.ATL.df4$Temperature__C_insitu)
# r = -0.46 ,  no collinearity , OK!


plot(log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu) , PP.ATL.df4$Temperature__C_insitu, 
     col = alpha ("darkgreen", 0.3), pch = 16,
     ylab= " Temperature in situ (°C) ",
     xlab = "log10( 1+Chlorophyll in situ, mg m-3 )")
abline(lmSST_C_insitu, lwd = 2.5, lty = 2, col = "darkorange")





names( Data.lin.mod.insitu)
dim(Data.lin.mod.insitu) 
Data.lin.mod.na.omit <- na.omit(Data.lin.mod.insitu) # remove lines with NAs
dim(Data.lin.mod.na.omit) 
# n = 404 complete NBSS datasets (with model nitrate , no in situ nitrate  .. , )
# n = 287    complete NBSS datasets (with in situ nitrate and  model nitrate, model silicate  )
# n = 287    complete NBSS datasets (with in situ nitrate only, no model data , no model silicate )



library(MASS)
#?stepAIC

full.lm1 <- lm((rob_reg__slopes) ~ ., data = Data.lin.mod.na.omit)
summary(full.lm1) # all variables are signif. except silicate
full.lm2 <- stepAIC(full.lm1, trace = FALSE)
full.lm2$anova
# all variables are retained except SALIINITY (n = 287 with in situ nitrate) !!!




# BEST MODEL (R-square = 0.4), 399 DF, p-value: < 2.2e-16 --------
summary(full.lm2)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.809567   0.063543 -12.741  < 2e-16 ***
#   nitrate               -0.022834   0.005178  -4.410 1.33e-05 ***
#   SST_C_insitu          -0.010395   0.004151  -2.504  0.01268 *  
#   Temperature__C_insitu -0.013277   0.004258  -3.118  0.00195 ** 
#   log10_1X_chla_insitu   0.383606   0.059855   6.409 4.14e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1921 on 399 degrees of freedom
# Multiple R-squared:  0.4039,	Adjusted R-squared:  0.3979 
# F-statistic: 67.58 on 4 and 399 DF,  p-value: < 2.2e-16


# #NEW:
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.904492   0.071889 -12.582  < 2e-16 ***
#   Temperature__C_insitu -0.021783   0.002781  -7.833 9.67e-14 ***
#   log10_1X_chla_insitu   0.361951   0.090168   4.014 7.64e-05 ***
#   Nitrate_insitu        -0.005392   0.002627  -2.052   0.0411 *  
# 
  AIC( full.lm2) # -85.25768
 # OLD: -179.324, Negative AIC indicates less information loss than a positive AIC and therefore a better model

 # BEST MODEL (R-square = 0.4), 399 DF, p-value: < 2.2e-16 --------
 summary(full.lm2)
 # Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)           -0.809567   0.063543 -12.741  < 2e-16 ***
 #   nitrate               -0.022834   0.005178  -4.410 1.33e-05 ***
 #   SST_C_insitu          -0.010395   0.004151  -2.504  0.01268 *  
 #   Temperature__C_insitu -0.013277   0.004258  -3.118  0.00195 ** 
 #   log10_1X_chla_insitu   0.383606   0.059855   6.409 4.14e-10 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Residual standard error: 0.1921 on 399 degrees of freedom
 # Multiple R-squared:  0.4039,	Adjusted R-squared:  0.3979 
 # F-statistic: 67.58 on 4 and 399 DF,  p-value: < 2.2e-16
 
 # VERY, Very Best Best model , without collinearity, only three explanatory variables,  ----
 # SST in situ, nitrate , log10(x+1 Chl a in situ) 
 full.lm3 <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
                  Data.lin.mod.na.omit$Nitrate_insitu+
                Data.lin.mod.na.omit$SST_C_insitu+
                  Data.lin.mod.na.omit$log10_1X_chla_insitu,
                                   data = Data.lin.mod.na.omit)

 summary(full.lm3)# Adjusted R-squared:  0.3848 , slightly worse Adjusted R-squared: than with  SST and Temp at depth
  
 AIC( full.lm3) # -179.324, Negative AIC indicates less information loss than a positive AIC and therefore a better model
 # -171.5964,  slightly worse AIC th an with  SST and Temp at depth
 
 ##################
 # TOP MODEL (full.lm4TOP) !!!! -----------
 # VERY, Very, very  Best Best  TOP model, without collinearity, only three explanatory variables,  ----
 # Temperature__C_insitu, nitrate , log10(x+1 Chl a in situ) 
 full.lm4TOP <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
                  Data.lin.mod.na.omit$Nitrate_insitu+
                  Data.lin.mod.na.omit$Temperature__C_insitu+
                  Data.lin.mod.na.omit$log10_1X_chla_insitu,
                data = Data.lin.mod.na.omit)
 
 summary(full.lm4TOP)
 # Multiple R-squared:  0.3945,	Adjusted R-squared:   0.39 ,
 # best R-squared for parsimonious model without collinearity! 
 # Explains about 39 % of the variability
 #  slightly worse Adjusted R-squared: than with  SST and Temp at depth
 # 400 degrees of freedom
 
 AIC( full.lm4TOP) # -179.324, Negative AIC indicates less information loss than a positive AIC and therefore a better model
 # -171.5964,  slightly worse AIC th an with  SST in situ  and Temp at depth in situ
 # -175.0244, slightly better  AIC than with  SST in situ 
 library(MASS)
 library(permuco)
 
 full.rlm4TOP <- rlm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
                     Data.lin.mod.na.omit$Nitrate_insitu+
                     Data.lin.mod.na.omit$Temperature__C_insitu+
                     Data.lin.mod.na.omit$log10_1X_chla_insitu,
                   data = Data.lin.mod.na.omit)
 summary(full.rlm4TOP)
 # 
 # Data.lin.mod.na.omit$nitrate                -0.0166   
 # Data.lin.mod.na.omit$Temperature__C_insitu  -0.0196   
 # Data.lin.mod.na.omit$log10_1X_chla_insitu    0.4335  
 
 # ratio abs(coef(log10_1X_chla_insitu)) / abs( coef (nitrate)        )
 abs(0.4335 )/  abs(-0.0166 )
 26.1 # the effect of Chla is 26 times stonger than for nitrate!
 
aovperm(full.rlm4TOP)# OK, all three variables are signif, with p < 0.001 
 
# relative importance within TOP Model ---------
library(relaimpo)

relaimpo::calc.relimp(full.lm4TOP) # three variable TOP model (nitrate,temp, chl a )
# Phytoplankton 
# lmg
# Data.lin.mod.na.omit$nitrate               0.019 # nitrate explains only 2% of variablity in NBSS Slope
# Data.lin.mod.na.omit$Temperature__C_insitu 0.22# temp in situ explains  22% of variability  in NBSS Slope
# Data.lin.mod.na.omit$log10_1X_chla_insitu  0.12 # chl a explains  12% of variability  in NBSS Slope
# Analysis based on 287 observations  (with in situ nitrate data, no filter for NBSS)

 cor.test ( Data.lin.mod.na.omit$Temperature__C_insitu , Data.lin.mod.na.omit$log10_1X_chla_insitu )
# -0.448031
 cor.test ( Data.lin.mod.na.omit$Nitrate_insitu , Data.lin.mod.na.omit$log10_1X_chla_insitu )
# -0.009446719
 cor.test ( Data.lin.mod.na.omit$Nitrate_insitu , Data.lin.mod.na.omit$Temperature__C_insitu )
 # -0.009446719
 
 



# # Using new Mass ratio productivity indices
# 
# #index based on effects
# 
# 
# Data.lin.mod.na.omit$prod_index_26_logX_1__mass_ratio4 <- 
#   (26*Data.lin.mod.na.omit$log10_1X_chla_insitu)  + ( Data.lin.mod.na.omit$nitrate)
# 
# Data.lin.mod.na.omit$prod_index_logX_1_logX_1_mass_ratio5 <- 
#   Data.lin.mod.na.omit$log10_1X_chla_insitu  + log10(1+ Data.lin.mod.na.omit$nitrate)
# 
# Data.lin.mod.na.omit$prod_index_2logX_1_logX_1_mass_ratio6 <- 
#   Data.lin.mod.na.omit$log10_1X_chla_insitu  + log10(1+ Data.lin.mod.na.omit$nitrate)
# 

# 
# 
# full.lm4 <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
#                    Data.lin.mod.na.omit$prod_index_26_logX_1__mass_ratio4,
#                                           data = Data.lin.mod.na.omit)
# summary(full.lm4)
#chla only: Multiple R-squared:  0.2535,	Adjusted R-squared:  0.2516 
# Multiple R-squared:  0.247,	Adjusted R-squared:  0.2451 , worse than Chla only!


# 
# full.lm5 <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
#                  Data.lin.mod.na.omit$prod_index_logX_1_logX_1_mass_ratio5,
#                data = Data.lin.mod.na.omit)
# summary(full.lm5)
# #chla only: Multiple R-squared:  0.2535,	Adjusted R-squared:  0.2516 
# # Multiple R-squared:  0.1578,	Adjusted R-squared:  0.155, much worse than Chla only!
# 
# 
# full.lm6 <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
#                  Data.lin.mod.na.omit$prod_index_2logX_1_logX_1_mass_ratio6,
#                data = Data.lin.mod.na.omit)
# summary(full.lm6)
#chla only: Multiple R-squared:  0.2535,	Adjusted R-squared:  0.2516 
# Multiple R-squared:  0.1578,	Adjusted R-squared:  0.1557 , much worse than Chla only!

lm.chla <- lm((Data.lin.mod.na.omit$rob_reg__slopes) ~ 
                 Data.lin.mod.na.omit$log10_1X_chla_insitu,
               data = Data.lin.mod.na.omit)
summary(lm.chla)
#chla only: Multiple R-squared:  0.2535,	Adjusted R-squared:  0.2516 


##############
# Carbon content of typical pico- and nanoplankton -----------

# Picoplancton Size and Carbon
#phytoplankton with sizes between 0,2 µm and 2 µm.--> Synechococcus 1  µm = 0.4 pgC/cell
#picogram (plural picograms) (metrology) An SI unit of mass equal to 10−12 grams. Symbol: pg.
#Picoplankton
0.4e-12 #g C
4e-13
log10(4e-13)
# -12.39 log10 g C ind. for Synechococcus

#Nanoplankton:
  
#Nanoplankton: phytoplankton with sizes between 2 µm and 20 µm. --> E.g.,
#coccolithophore: Emiliania huxleyi 8 µm  = 25 pgC/cell; 
#Coccolithus braarudii 20 µm  = 200 pgC/cell


#Emiliania huxleyi, 8 µm  = 25 pgC/cell; 
25e-12 #g C
2.5e-11 #g C
log10(25e-12) 
# -10.60206 log10 g C ind. for Emiliania huxleyi


#Coccolithus braarudii 20 µm  = 200 pgC/cell
200e-12 #g C
2e-10 #g C
log10(2e-10) # -9.7 log10 g C ind. for 

# Picopl. approx.  -13 log10gC cell.
# Nanopl. ca. -12 to -10 log10gC cell.
# Microplankton > -10 log10gC cell.


  #############
# three Chla groups FOR PAPER -----------------------
# classified by  Chla  in situ !

# Separate OLigo-, meso-, and  eutrophic samples --------

# find boundaries -------------------
length(phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) # 88 datasets
# 458 total datasets

plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
       phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu)

abline (v = 0.2)
abline (v = 2)

oligotr.phytopl.All.Atl <- subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, Chlorophyll__mg_m_3_insitu < 0.2 )
summary(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
length(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 88 datasets

mesotr.phytopl.All.Atl <- subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, Chlorophyll__mg_m_3_insitu > 0.2 )
mesotr.phytopl.All.Atl <- subset(mesotr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 2 )
summary(mesotr.phytopl.All.Atl$SST_C_insitu) # 20.15 to 27.75 degr Celsius
length(mesotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
# 301 datasets


eutr.phytopl.All.Atl <-subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, Chlorophyll__mg_m_3_insitu > 2 )
eutr.phytopl.All.Atl <-subset(eutr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 31 )
summary(eutr.phytopl.All.Atl$SST_C_insitu) # 27.97 to 29.28 degr Celsius
length(eutr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 88 datasets
# 55 datasets


# 1. oligotrophic waters --------

# the hotter, the steeper, with decreasing temperature, more large cells.
# with increasing temperature, less large-sized cells (in carbon)
summary(oligotr.phytopl.All.Atl$SST_C_insitu) # 11.04 to 34.6 degr Celsius, oligotr.

lm1.olig <-  lm(oligotr.phytopl.All.Atl$slope ~ (oligotr.phytopl.All.Atl$SST_C_insitu ))
summary(lm1.olig)

# (Intercept)                          -1.047371   0.106384  -9.845 9.38e-16 ***
#   oligotr.phytopl.All.Atl$SST_C_insitu -0.011589   0.004191  -2.765  0.00696 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1794 on 86 degrees of freedom
# Multiple R-squared:  0.08165,	Adjusted R-squared:  0.07097 
# F-statistic: 7.646 on 1 and 86 DF,  p-value: 0.006961

# library(broom)
# broom::tidy(lm1.olig)



plot(oligotr.phytopl.All.Atl$slope ~ (oligotr.phytopl.All.Atl$SST_C_insitu), 
     xlim = c(3, 32), main = "oligotrophic (0 to 0.19 mg m-3 Chla)",
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lm1.olig, lwd = 2.5, lty = 2, col = "darkorange")

# Flatter with increasing chl a ... !!! 
# relatively more large-sized  cells  with increasing Chla ! 
# very steep!


# 2. mesotrophic waters --------

# the hotter, the steeper, with decreasing temperature, more large cells.
# with increasing temperature, less large-sized cells (in carbon)
summary(mesotr.phytopl.All.Atl$SST_C_insitu) # 11.04 to 34.6 degr Celsius, oligotr.

lm1.mesotr <-  lm(mesotr.phytopl.All.Atl$slope ~ (mesotr.phytopl.All.Atl$SST_C_insitu ))
summary(lm1.mesotr)

# #                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.781346   0.042508  -18.38   <2e-16 ***
#   mesotr.phytopl.All.Atl$SST_C_insitu -0.022551   0.001938  -11.64   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2123 on 287 degrees of freedom
# (12 observations deleted due to missingness)
# Multiple R-squared:  0.3206,	Adjusted R-squared:  0.3182 
# F-statistic: 135.4 on 1 and 287 DF,  p-value: < 2.2e-16

#broom::tidy(lm1.mesotr)


plot(mesotr.phytopl.All.Atl$slope ~ (mesotr.phytopl.All.Atl$SST_C_insitu), 
     xlim = c(3, 32), main = "mesotrophic (2 to 1.9 mg m-3 Chla)",
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lm1.mesotr, lwd = 2.5, lty = 2, col = "darkorange")

# very steep!


# 3. eutrophic waters --------

# All sampes have very flat NBSS slopes
# no relationship with temperature (mostl cold samples, below 20ºC)

summary(eutr.phytopl.All.Atl$SST_C_insitu) # 11.04 to 34.6 degr Celsius, oligotr.

lm1.eutroph <-  lm(eutr.phytopl.All.Atl$slope ~ (eutr.phytopl.All.Atl$SST_C_insitu ))
summary(lm1.eutroph)
# NOT SIGNIFICANT!!!
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       -0.791507   0.135900  -5.824 3.44e-07 ***
#   eutr.phytopl.All.Atl$SST_C_insitu -0.004557   0.008536  -0.534    0.596    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1575 on 53 degrees of freedom
# Multiple R-squared:  0.005349,	Adjusted R-squared:  -0.01342 
# F-statistic: 0.285 on 1 and 53 DF,  p-value: 0.5957

# broom::tidy(lm1.eutroph)

summary(eutr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu)

plot(eutr.phytopl.All.Atl$slope ~ (eutr.phytopl.All.Atl$SST_C_insitu), 
     xlim = c(3, 32), main = "eutrophic (2 to 15 mg m-3 Chla)",
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
#abline(lm1.eutroph, lwd = 2.5, lty = 2, col = "darkorange")

# very steep!


# output lm summary as table, with broom---------------
# library(broom)
# model <- lm(mpg~cyl, mtcars)
# broom::tidy(model)



#############
# three Temperature groups FOR PAPER -----------------------
# classified by  SST in situ !


### slope vs Chla, in three SST groups --------------------

###  slope vs Chla, by Temp_class --------------------

length(phytopl.All.Atl.p.rsq.slope_ratio.ok $Date) # 462 phytoplankton
summary(phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu) # 13 NAs
# range from 5.50  to  29.28 ºC.

# Cold: 3 to 20 degr. Celsius
# Suptropical: 20 to 27.8 degr. Celsius 
# Tropical: > 27.8 (27.8 to 31) degr. Celsius

Cold.phytopl.All.Atl <- subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, SST_C_insitu < 20 )
summary(Cold.phytopl.All.Atl$SST_C_insitu) # 5.5 to 19.8 degr Celsius
summary(Cold.phytopl.All.Atl$rob_reg__slopes) #  
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7789 -1.2098 -1.0498 -1.0661 -0.8879 -0.4416 

# median Cold RR slope: -1.05

Cold.phytopl.All.Atl.unfiltr <- subset(phytopl.All.Atl, SST_C_insitu < 20 )
dim(Cold.phytopl.All.Atl)# 196 Cold water datasets, top linear 
dim(Cold.phytopl.All.Atl.unfiltr)# 371 Cold water datasets , ALL (linear and messy)
+(371-196)/371 # 47.2 % percent non-linear datasets  (Cold waters only)


SubTrop.phytopl.All.Atl <- subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, SST_C_insitu > 20 )
SubTrop.phytopl.All.Atl <- subset(SubTrop.phytopl.All.Atl, SST_C_insitu < 27.8 )
summary(SubTrop.phytopl.All.Atl$SST_C_insitu) # 20.15 to 27.75 degr Celsius
summary(SubTrop.phytopl.All.Atl$rob_reg__slopes) #  
dim(SubTrop.phytopl.All.Atl) # 152 datasets, top-linear only, 220 ALL datasets
# median -1.28

SubTrop.phytopl.All.Atl.unfiltr <-subset(phytopl.All.Atl, SST_C_insitu > 20 )
SubTrop.phytopl.All.Atl.unfiltr <-subset(SubTrop.phytopl.All.Atl.unfiltr, SST_C_insitu < 27.8 )
dim(SubTrop.phytopl.All.Atl.unfiltr)# 152 tropical top-linear datasets, 220 ALL datasets, 

+(220-152)/220 # 30.9 % percent non-linear datasets  (Subtropical)


Tropic.phytopl.All.Atl <-subset(phytopl.All.Atl.p.rsq.slope_ratio.ok, SST_C_insitu > 27.8 )
Tropic.phytopl.All.Atl <-subset(Tropic.phytopl.All.Atl, SST_C_insitu < 31 )
summary(Tropic.phytopl.All.Atl$SST_C_insitu) # 27.97 to 29.28 degr Celsius
summary(Tropic.phytopl.All.Atl$rob_reg__slopes) #
dim(Tropic.phytopl.All.Atl)# 95 tropical top-linear datasets
# median NBSS slope (RR): -1.40

Tropic.phytopl.All.Atl.unfiltr <-subset(phytopl.All.Atl, SST_C_insitu > 27.8 )
Tropic.phytopl.All.Atl.unfiltr <-subset(Tropic.phytopl.All.Atl.unfiltr, SST_C_insitu < 31 )
dim(Tropic.phytopl.All.Atl.unfiltr)# 95 tropical top-linear datasets, 110 ALL datasets, 

+(110-95)/110 # 13.6 % percent non-linear datasets  (Tropical)



# Cold waters (temperate & Upwelling): 3 to 20 degr. Celsius -----------------------
#Cold.phytopl.All.Atl <- subset(phytopl.All.Atl, SST_C_insitu < 20 )
summary(Cold.phytopl.All.Atl$SST_C_insitu) # 5. to 19.8 degr Celsius

lm1.cold <-  lm(Cold.phytopl.All.Atl$slope ~ log10(Cold.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu))
summary(lm1.cold)
# old (Biovolume)  R-squared:  0.3586,  p-value: < 2.2e-16, 193 degrees of freedom
# NEW, C units:  R-squared:  0.314,  194 DF,  p-value: < 2.2e-16
plot(Cold.phytopl.All.Atl$slope ~ log10(Cold.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu), 
     xlim = c(-2, 1.5),
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lm1.cold, lwd = 2.5, lty = 2, col = "darkorange")

summary(Cold.phytopl.All.Atl$rob_reg__slopes) # 
# Old, biovolume:
#    Min.   1st Qu. Median    Mean 3rd Qu.    Max. 
# -1.4554 -1.0457 -0.9208 -0.9240 -0.7797 -0.3081  

# new, Carbon:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.7789 -1.2098 -1.0498 -1.0661 -0.8879 -0.4416 

# Flatter with increasing chl a ... !!! 
# relatively more large-sized  cells  with increasing Chla ! 
# very steep!


# 
# Subtropical  waters: 20 to  27.8 -----------------------
summary(SubTrop.phytopl.All.Atl$SST_C_insitu) # 20 to 27.8 degr Celsius

lm1.subtrop <-  lm(SubTrop.phytopl.All.Atl$rob_reg__slopes ~ log10(SubTrop.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu))
summary(lm1.subtrop)
# Old, (Biovolume): R-squared:  0.04376, p = 0.008, 155 degrees of freedom
# New, carbon: R-squared:  0.05918, 149 DF,  p-value: 0.002613
plot(SubTrop.phytopl.All.Atl$rob_reg__slopes ~ 
       log10(SubTrop.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu),
     xlim = c(-2,  1.5),
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
abline(lm1.subtrop, lwd = 2.5, lty = 2, col = "darkorange")

summary(SubTrop.phytopl.All.Atl$rob_reg__slopes) # 
# OLD biovolume:
#     Min.     1st Qu.  Median   Mean    3rd Qu      Max. 
#   -1.4092   -1.1593   -1.1079 -1.0992 -1.0376 -0.8184 

# NEW (Carbon)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.6464 -1.3593 -1.2828 -1.2840 -1.2045 -0.8323 

# Flatter with increasing chl a ... !!! 
  # relatively more large-sized  cells  with increasing Chla ! 
  


# Tropical  waters: > 27.8 -----------------------
summary(Tropic.phytopl.All.Atl$SST_C_insitu) # > 27.8 degr Celsius
summary(Tropic.phytopl.All.Atl$rob_reg__slopes) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.8343 -1.4394 -1.3019 -1.3125 -1.1853 -0.8082 
# 

lm1.trop <-  lm(Tropic.phytopl.All.Atl$rob_reg__slopes ~ log10(1+Tropic.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu))
summary(lm1.trop) # ns.

rlm1.trop <-  rlm(Tropic.phytopl.All.Atl$rob_reg__slopes ~ log10(1+Tropic.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu))
aovperm(rlm1.trop, np = 5000) #  p = 0.043


plot(Tropic.phytopl.All.Atl$rob_reg__slopes ~
       log10(1+Tropic.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu),
     xlim = c(-2,  1.5),
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
#abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")


plot(Tropic.phytopl.All.Atl$rob_reg__slopes ~
       log10(1+Tropic.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu),
    # xlim = c(-2,  1.5),
     ylim = c(-2.2, 0.1),
     col = alpha ("darkgreen", 0.3), pch = 16)
   abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")

 # nearly flat shape!  But... steeper with more chl a ... !!! 
   # relatively more small-sized  cells  with increasing Chla ! 
   
 
   
   # Divide into 9 groups(7 actually exist with real samples) ----------
  #plot with 9 sections (cold, eutrophic etc.)
   
   plot(phytopl.All.Atl$Chlorophyll__mg_m_3_insitu ,
        phytopl.All.Atl$SST_C_insitu,
        #               xlim = c(-2,  1.5),
        ylim = c(3, 30),
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
   abline(h = c(20, 27.5))
   abline(v = c(2, 0.2))
   
   
   
   plot(log10(1+phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) ,
        phytopl.All.Atl$SST_C_insitu,
        #               xlim = c(-2,  1.5),
        ylim = c(3, 30),
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
   abline(h = c(20, 27.5), col = alpha("darkgrey", 0.6), lwd = 2)
   abline(v = c(  0.0792,  0.477), col = alpha("darkgrey", 0.6), lwd = 2)
   log10(1+2) # 0.477
   log10(1+0.2) # 0.0792
   
   
   
   plot(log10(1+phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) ,
        phytopl.All.Atl$SST_C_insitu,
        #               xlim = c(-2,  1.5),
        ylim = c(3, 30),
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
   abline(h = c(20, 27.5), col = alpha("darkgrey", 0.6), lwd = 2)
   abline(v = c(  0.0792,  0.477), col = alpha("darkgrey", 0.6), lwd = 2)
   log10(1+2) # 0.477
   log10(1+0.2) # 0.0792
   
   
   # K. Calculate % linear a datasets within each temp and chla range --------
   
   # K.1 Calculate % linear a datasets within  Chla ranges --------
   
   # Define the breakpoints (ranges)
   breaks.log.Chla <- c(0, 0.25, 0.5, 0.75, 1.25)  # 
   diff.breaks.log.Chla <-    # 
     
   
   
   #  Chl a vector (in situ)
   x <- log10(1+phytopl.All.Atl$Chlorophyll__mg_m_3_insitu)
   length(x)# all data, n = 727
 
     x2.unfiltr <- log10(1+phytopl.All.Atl.unfiltr.nolims727$Chlorophyll__mg_m_3_insitu)
   length(x2.unfiltr)# all data, n = 727, linear and non-linear
   length(x2.unfiltr[!is.na(x2.unfiltr)])
   # all data,with Chla in situ data n = 562, linear and non-linear
   x2.unfiltr.wChla <- x2.unfiltr[!is.na(x2.unfiltr)] # n = 562
   summary(x2.unfiltr.wChla)
   length(x2.unfiltr.wChla)
   
    x3.linear <- log10(1+phytopl.Carbon.All.Atl.ok3$Chlorophyll__mg_m_3_insitu)
   length(x3.linear)# all linear data, n = 311
   length(x3.linear[!is.na(x3.linear)])
   # all data,with Chla in situ data n = 310, linear top-quality oly
   x3.linear.wChla <- x3.linear[!is.na(x3.linear)] # # all linear data, n = 311
     summary(x3.linear.wChla)
   length(x3.linear.wChla)
   
   # Use cut() to categorize the data into intervals
   x_cut2 <- cut(x2.unfiltr.wChla, breaks = breaks.log.Chla, right = FALSE)
   
   # Count the number of observations in each range
   range_counts2 <- table(x_cut2)
   sum(range_counts2)
   
   # Print the result
   range_counts2
   

   # Use cut() to categorize the data into intervals
   x_cut3 <- cut(x3.linear.wChla, breaks = breaks.log.Chla, right = FALSE)
   
   # Count the number of observations in each range
   range_counts3 <- table(x_cut3)
   sum(range_counts3)
   
   # Print the result
   range_counts3
   
   
   # Calculate percentage in each Chla range
   
   (x3.2.ratio  <- range_counts3/range_counts2 )
   (perc.linear  <- range_counts3/range_counts2 * 100)
   
   ## FOR PAPER - BARPLOT for suppl mat ---------- 
   
    barplot(x3.2.ratio, ylab = "% linear-shaped NBSS", 
            xlab = "log10 ((1+Chla in situ), mg m-3)", col = "lightblue")

     
   # Compute midpoints of each bin
   midpoints.log.Chla <- (head(breaks.log.Chla, -1) + tail(breaks.log.Chla, -1)) / 2
   
   # Combine into a data frame
   res.Chla <- data.frame(
       midpoint = midpoints.log.Chla,
     perc.linear = as.numeric(perc.linear)
   )
   
   print(res.Chla)
 
   permuco::aovperm(lm( as.numeric(perc.linear) ~ midpoints.log.Chla ))  
   # 0.04167,  OK

   # cor.test (as.numeric(perc.linear) , midpoints.log.Chla, 
   #           method = "spearman")   

   
   
   
   # K.2 Calculate % linear a datasets within  Temp. ranges --------
   
   # Define the breakpoints (ranges)
   breaks.log.Temp <- c(5, 15, 20, 25, 37)  # 
     
     
     #  Temp  as vector (SST in situ)
     Temp <- (phytopl.All.Atl$SST_C_insitu)
   length(Temp)# all data, n = 727
 
   summary(Temp) # 13 NA's
     
   Temp2.unfiltr <- (phytopl.All.Atl.unfiltr.nolims727$SST_C_insitu)
   length(Temp2.unfiltr)# all data, n = 727, linear and non-linear
   length(Temp2.unfiltr[!is.na(Temp2.unfiltr)])
   # all data,with Temp in situ data n = 714, linear and non-linear
   Temp2.unfiltr.wTemp <- Temp2.unfiltr[!is.na(Temp2.unfiltr)] # n = 562
   summary(Temp2.unfiltr.wTemp)
   length(Temp2.unfiltr.wTemp)
   
   Temp3.linear <- (phytopl.Carbon.All.Atl.ok3$SST_C_insitu)
   length(Temp3.linear)# all linear data, n = 311
   length(Temp3.linear[!is.na(Temp3.linear)])
   # all data,with Temp in situ data n = 306, linear top-quality only
   Temp3.linear.wTemp <- Temp3.linear[!is.na(Temp3.linear)] # # all linear data, n = 311
   summary(Temp3.linear.wTemp)
   length(Temp3.linear.wTemp)
   
   # Use cut() to categorize the data into intervals
   Temp_cut2 <- cut(Temp2.unfiltr.wTemp, breaks = breaks.log.Temp, right = FALSE)
   
   # Count the number of observations in each range
   range_counts2 <- table(Temp_cut2)
   sum(range_counts2)
   
   # Print the result
   range_counts2
   
   
   # Use cut() to categorize the data into intervals
   Temp_cut3 <- cut(Temp3.linear.wTemp, breaks = breaks.log.Temp, right = FALSE)
   
   # Count the number of observations in each range
   range_counts3 <- table(Temp_cut3)
   sum(range_counts3)
   
   # Print the result
   range_counts3
   
   
   # Calculate percentage in each Temp range
   
   (Temp3.2.ratio  <- range_counts3/range_counts2 )
   (perc.linear  <- range_counts3/range_counts2 * 100)
   
   ##  BARPLOT for suppl mat ---------- 
   
   barplot(Temp3.2.ratio, ylab = "% linear-shaped NBSS", 
           xlab = "(SST in situ), mg m-3", col = "lightblue")
   
   
   # Compute midpoints of each bin
   midpoints.log.Temp <- (head(breaks.log.Temp, -1) + tail(breaks.log.Temp, -1)) / 2
   
   # Combine into a data frame
   res.Temp <- data.frame(
     midpoint = midpoints.log.Temp,
     perc.linear = as.numeric(perc.linear)
   )
   
   print(res.Temp)
   
   permuco::aovperm(lm( as.numeric(perc.linear) ~ midpoints.log.Temp ))  
   # 0.04167,  OK
   
   # cor.test (as.numeric(perc.linear) , midpoints.log.Chla, 
   #           method = "spearman")   
   
   
   # L. Picoplankton: analyze samples with picoplankton ------------
   
   # L.1 select samples with picoplankton ----------------
   
   
   phytopl.All.Atl.unfiltr.nolims727b <- phytopl.All.Atl.unfiltr.nolims727
   
   # calculate mean picoplankton biomass for each sample
   
   
  dim( phytopl.NBSS.matrix.unfiltr727)
  
   names(phytopl.All.Atl.unfiltr.nolims727 [  c( 38:(38+5), (38+16):91 )] )  
   # 45 NBSS classes not used for phytopl NBSS (n = 10 size bins)
   
   names(phytopl.All.Atl.unfiltr.nolims727 [  c( 38:(38+4))] )  
   # first 5 NBSS classes  used for Picoplancton counts (n = 5 size bins)
   
   logCvector[1:5] #  5 selected bins, picoplaknton

      # Count non-NA values in columns 1 to 5 
   non_na_counts <- rowSums(!is.na(phytopl.NBSS.matrix.unfiltr727[, 1:5]))
   
   # number of  datsets with picoplankton (Count non-zero values)
   non_zero_count <- sum(non_na_counts != 0)
   non_zero_count
   # 510 datasets with picoplankton (at least one value)
  
   at.least3_count <- sum(non_na_counts > 2)
   at.least3_count
   # only 309 datasets with 3 or more picoplankton non-empty bins
   
   at.least4_count <- sum(non_na_counts > 3)
   at.least4_count
   # only 279 datasets (49%) with 4 or more picoplankton non-empty bins
   +279/565
    
   # Select  rows with at least 4 picoplankton data points (non-empty bins)
   
  ( Pico.na_counts <- apply(phytopl.NBSS.matrix.unfiltr727[, 1:5], 1, function(x) sum(!is.na(x))) )
  (  Pico.NOT.NA_counts <- apply(phytopl.NBSS.matrix.unfiltr727[, 1:5], 1, function(x) sum(!is.na(x))) )
   non_zero_count <- sum(Pico.NOT.NA_counts > 0)
   non_zero_count # 510 datasets with picoplankton, OK
   
   # median and mean PICO NBSS values by sample
   
   ( Pico.row_medians <- apply(phytopl.NBSS.matrix.unfiltr727[, 1:5], 1, function(x) median(x, na.rm = TRUE))  )
   ( Pico.row_means <- apply(phytopl.NBSS.matrix.unfiltr727[, 1:5], 1, function(x) mean(x, na.rm = TRUE))  )
   
   

   phytopl.All.Atl.unfiltr.nolims727b$Pico.NOT.NA_counts <- Pico.NOT.NA_counts
   phytopl.All.Atl.unfiltr.nolims727b$Pico.row_medians <- Pico.row_medians
   phytopl.All.Atl.unfiltr.nolims727b$Pico.row_means <- Pico.row_means
   
   
   df <- phytopl.All.Atl.unfiltr.nolims727b
   
   # select samples with at least 4 Picopl. data
   phytopl.All.Atl.atleast.4.PICO <- df[df$Pico.NOT.NA_counts > 3, ]
    dim(phytopl.All.Atl.atleast.4.PICO)  
   # 279 datasets with at least 4 Picopl. data, OK
   names(phytopl.All.Atl.atleast.4.PICO)
    
   # L.2 analyse  Picoplankton slopes vs temp and chl a  ----------------
   
   
    # Calculate NBSS slopes with Picopl only ------------
     
   # Define functions that calculate  picoplankton NBSS SLOPES and p values
   #   Use first 5 bins only (PICOPL.)
   # Fit functions for OLSR and RR models
   
    fun.1cslope.lm.olsr.w.X_vec.PICO <- function(yvalues ) {
      mod1<-  lm( log10(yvalues) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
      mod1$coefficients[2]
      mod1_slope <- as.numeric(mod1$coefficients[2])
      ; mod1_slope}
    
    
    fun.1cslope.lm.rr.w.X_vec.PICO <- function(yvalues ) {
      mod1<-  rlm( log10(yvalues) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
      mod1$coefficients[2]
      mod1_slope <- as.numeric(mod1$coefficients[2])
      ; mod1_slope}
   
    fun.1cINTERC.lm.rr.w.X_vec.PICO <- function(yvalues ) {
      mod1<-  rlm( log10(yvalues) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
      mod1$coefficients[2]
      mod1_interc <- as.numeric(mod1$coefficients[1])
      ; mod1_interc}
    
    
    
    fun.1d.P_Val_lm.olsr.w.X_vec.PICO <- function(yvalues ) {
      mod1<-  lm( log10(yvalues) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
      mod1$coefficients[2]
      #mod1_slope <- as.numeric(mod1$coefficients[2])
      p_value_olsr <-  summary(mod1)$coefficients[2,4] 
      ; p_value_olsr}
    
    fun.1d.P_Val_lm.rr.w.X_vec.PICO <- function(yvalues ) {
      mod1<-  rlm( log10(yvalues) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)),
                                   na.action = "na.omit")
      mod1$coefficients[2]
      #mod1_slope <- as.numeric(mod1$coefficients[2])
      p_value_rr <-  permuco::aovperm(mod1, nperm = 2000) 
      ; p_value_rr}
    
    
    
    ( Pico.olsr.slopes <- apply(phytopl.All.Atl.atleast.4.PICO[, 38:(38+4)], 1, function(x) fun.1cslope.lm.olsr.w.X_vec.PICO(x))  )
    non_zero_count <- sum(Pico.olsr.slopes != 0)
    non_zero_count # 510 datasets with picoplankton olsr SLOPES, OK
    summary(Pico.olsr.slopes)
    
    ( Pico.rr.slopes <- apply(phytopl.All.Atl.atleast.4.PICO[, 38:(38+4)], 1, function(x) fun.1cslope.lm.rr.w.X_vec.PICO(x))  )
    non_zero_count <- sum(Pico.olsr.slopes != 0)
    non_zero_count # 510 datasets with picoplankton olsr SLOPES, OK
    summary(Pico.rr.slopes)
  
      ( Pico.rr.interc <- apply(phytopl.All.Atl.atleast.4.PICO[, 38:(38+4)], 1, function(x) fun.1cINTERC.lm.rr.w.X_vec.PICO(x))  )
    non_zero_count <- sum(Pico.rr.interc != 0)
    non_zero_count # 510 datasets with picoplankton olsr SLOPES, OK
    summary(Pico.rr.interc)
    
    
    
    ( Pico.olsr.p_val <- apply(phytopl.All.Atl.atleast.4.PICO[, 38:(38+4)], 1, function(x) fun.1d.P_Val_lm.olsr.w.X_vec.PICO(x))  )
    non_zero_count <- sum(Pico.olsr.p_val != 0)
    non_zero_count # 510 datasets with picoplankton olsr SLOPES, OK
    summary(Pico.olsr.p_val)
    
     
    
    length(Pico.olsr.p_val[Pico.olsr.p_val < 0.05 ] ) # only 119 signif. linear models....
    
    # lnumb = 4
    # fun.1d.P_Val_lm.rr.w.X_vec.PICO(as.numeric( phytopl.All.Atl.atleast.4.PICO[lnumb, 38:(38+4)]))
    # 
    # ( Pico.rr.p_val <- apply(phytopl.All.Atl.atleast.4.PICO[, 38:(38+4)], 1, function(x) fun.1d.P_Val_lm.rr.w.X_vec.PICO(x))  )
    # non_zero_count <- sum(Pico.rr.p_val != 0)
    # non_zero_count # 510 datasets with picoplankton olsr SLOPES, OK
    # summary(Pico.rr.p_val)
    
      phytopl.All.Atl.atleast.4.PICO$Pico.olsr.slopes <- Pico.olsr.slopes
      phytopl.All.Atl.atleast.4.PICO$Pico.rr.slopes <- Pico.rr.slopes
      phytopl.All.Atl.atleast.4.PICO$Pico.olsr.p_val <- Pico.olsr.p_val
      phytopl.All.Atl.atleast.4.PICO$Pico.rr.interc <- Pico.rr.interc
      
      
      
    plot(Pico.rr.slopes ~ phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)
    plot(Pico.rr.slopes ~ phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)
    
    plot(Pico.rr.slopes ~ 
           log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
    abline (h = -1)
    
    
      plot(log10(phytopl.All.Atl.atleast.4.PICO$Pico.row_medians) ~ 
           log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
abline( lm (log10(phytopl.All.Atl.atleast.4.PICO$Pico.row_medians) ~ 
              log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) )
      
plot((Pico.olsr.slopes) ~ 
       log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
abline( lm (log10(phytopl.All.Atl.atleast.4.PICO$Pico.row_medians) ~ 
              log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) )

  
  
    summary(lm(Pico.rr.slopes ~ 
                 log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)))
    # p-value: p-value: < 2.2e-16
    
  cor.test  (Pico.rr.slopes , 
        log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu),
  method = "spearman")
  # p-value < 2.2e-16
        
    
    plot(log10( phytopl.All.Atl.atleast.4.PICO$Pico.row_medians) ~ 
           log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
      
       
    # Picoplancton NBSS slopes vs Temp. is SIGNIFFICANT!
    plot(Pico.rr.slopes ~ 
                 (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
    abline (h = -1)
    
     summary(m101Temp <- lm(Pico.rr.slopes ~ 
                 (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
    # p-value:  < 2.2e-16
    #  R-squared:  0.560
     
     abline (m101Temp)     
     
     
      summary(m101Chl <- lm(as.numeric(Pico.rr.slopes) ~ 
          log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)))
     # p-value: < 2.2e-16
     # R-squared: only      0.45
      plot(as.numeric(Pico.rr.slopes) ~ 
        log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
  
abline (m101Chl)     
      
     
     aovperm(lm(as.numeric(Pico.rr.slopes) ~ 
                  (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
     # p-value:  2e-04
  
     aovperm(lm(as.numeric(Pico.rr.slopes) ~ 
                  log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)))
     # p-value:  0.0014
     
     
     
    cor.test  (Pico.rr.slopes , 
               (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu),
               method = "spearman")
    # p-value < 2.2e-16
    
  
    # Picoplancton NBSS median biomass  vs Temp. is SIGNIFFICANT!
    plot(phytopl.All.Atl.atleast.4.PICO$Pico.row_medians ~ 
      (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
    
     summary(lm(phytopl.All.Atl.atleast.4.PICO$Pico.row_medians ~ 
                 (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
    # p-value: ..
    
    cor.test  (phytopl.All.Atl.atleast.4.PICO$Pico.row_medians , 
               (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu),
               method = "spearman")
    # p-value = 0.0001685
    
    
    # L.3 multivariate Picopl. models -------------

    # no interaction (Adjusted R-squared:  0.06)
    m55 <- (lm(as.numeric(Pico.rr.slopes) ~ 
                 log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)+
                 phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
  summary(m55)
   # Multiple R-squared:  0.6307,	Adjusted R-squared:  0.6279 
  
    aovperm(m55)
   # SST is signifficant  < 2e-16 ***
    # Chla is  also signifficant, 3.94e-12 ***
    # Overall p-value: p-value: < 2.2e-16
    
    
    # Including another variable improved the adjusted R² and the AIC of the model.
    
    relaimpo::calc.relimp(m55)
    # lmg
    # log10(1 + Chlorophyll__mg_m_3_insitu) 0.26 (26% is explained)
    # SST_C_insitu                         0.37 (37% is explained)
    # 
    
    
    
    # linear model with interaction   ----------
    
    m57 <- (lm(as.numeric(Pico.rr.slopes) ~ 
                 log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu) *
                 phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
    summary(m57)
   # p-value: < 2.2e-16
    # Multiple R-squared:  0.6318,	Adjusted R-squared:  0.6276 
    aovperm(m57, nperm = 2000)
    
    # Including another variable with interaction 
    # did NOT improve the adjusted R² and the AIC of the model.
    
    
    # SST is highly  sign. (!)
    # Chl a is NOT signif.,
    # Interaction Chl a x SST is NOT signiff. (
    
    # prefer the model without interaction!!!, same adjsted R-squ (!)
      
    
    gam58b <- (mgcv::gam(as.numeric(Pico.rr.slopes) ~ 
                          s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu) ))) 
                            
    summary(gam58b)
    
    #R-sq.(adj) =  0.0731   Deviance explained = 8.09%, p = 1e-04 ***
    plot(gam58b)
   
     gam58c <- (mgcv::gam(as.numeric(Pico.rr.slopes) ~ 
                           s((phytopl.All.Atl.atleast.4.PICO$SST_C_insitu) ))) 
     summary(gam58c)
     # R-sq.(adj) =  0.625   Deviance explained = 63.7%
     # p = 3.04e-05 **
    plot(gam58c)
    
    
   
    gam59 <- (mgcv::gam(as.numeric(Pico.rr.slopes) ~ 
                 s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu) ,
                 phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
    #plot(gam59)
    summary(gam59)
    # R-sq.(adj) =  0.689   Deviance explained = 71.2% - BEST univariate model !!!!
    
    #?mgcv::gam
    
    gam61.noint <- (mgcv::gam(as.numeric(Pico.rr.slopes) ~ 
                          s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) +
                            s(phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
    #plot(gam59)
    summary(gam61.noint)
    # R-sq.(adj) =  0.654   Deviance explained = 66.6%
    # SST and  Chl are signif
    # R sq is worse thn univariate SST model   !
   
    
    gam61.w.int <- (mgcv::gam(as.numeric(Pico.rr.slopes) ~ 
                                s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) +
                                s(phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)+
                                +ti(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu),phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
    #plot(gam59)
    summary(gam61.w.int)
    # R-sq.(adj) =  0.657   Deviance explained =   67%
    # SST and  Chl are signif
    # interactioonis not  signif.
    # bad model.. 
    
    # the BEST picopl. model is the univariate gam SST model (Deviance explained = 71.2%)
    
  
    
    # * L.4 now only with  significant PICOPL. NBSS models (p <005) * -------------
    
    phytopl.All.Atl.PICO.p_OK <- subset (phytopl.All.Atl.atleast.4.PICO , Pico.olsr.p_val < 0.05)
dim(phytopl.All.Atl.PICO.p_OK) 
# 119 samples with at least 4 Picopl. bins and p < 0.05
    
summary(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.851  -3.330  -2.708  -2.487  -2.171   1.068 
# > 


plot(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
       log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) )
    # , col = phytopl.All.Atl.PICO.p_OK$SST_C_insitu)
symbols(log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu),
        phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes,
        circles  = (phytopl.All.Atl.PICO.p_OK$SST_C_insitu/3000), inches = 0.2)
abline(h = -1)

m62Chla <- lm(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
                (log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) ))
summary(m62Chla)
# R-squared:  0.3986,	 p-value: 1.251e-05, Chla  explains 40% variabilty , n = 40

plot(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
       (phytopl.All.Atl.PICO.p_OK$SST_C_insitu) )
symbols((phytopl.All.Atl.PICO.p_OK$SST_C_insitu),
        phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes,
        circles  = log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu/3000), inches = 0.3)

abline(h = -1)

m62Temp <- lm(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
                 (phytopl.All.Atl.PICO.p_OK$SST_C_insitu) )
summary(m62Temp)
# R-squared:  0.4831, p-value: 9.123e-07
#   SST explains 48% of the variability , n = 40 top PICO lm models

# linear model with interaction, n = 40

m63.w.int <- (lm(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
             log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) *
             phytopl.All.Atl.PICO.p_OK$SST_C_insitu))
summary(m63.w.int)
#Multiple R-squared:  0.5587,	Adjusted R-squared:  0.5208 
# p-value: 2.234e-06 
# only SST is signifficant... Not a good model , 
# Chl a is not siognifficant, intaraction is not signif
# Not a good model , 

aovperm(m63.w.int, nperm = 2000)

# linear model without interaction, n = 119
m63.no.int <- (lm(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
                   log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) +
                   phytopl.All.Atl.PICO.p_OK$SST_C_insitu))
summary(m63.no.int)
# Multiple R-squared:  0.5569,	Adjusted R-squared:  0.5322 , explains 55% of the variabilty!
# p-value: 4.343e-07
#  SST is highly signifficant !!!   0.000767 ***
# Chl  is also siognifficant! 0.0193
# interaction is not signif !!!
# A very good model 
# explains 55% of the variability!
# Adding Chla improves the model! Both variable are important
# N = 40 top picoplankton NBSS

aovperm(m63.no.int, nperm = 2000)


# non-linear models (GAM), picopl. only (n = 40) ----------

# non-linear model with interaction, n = 40

gam65.w.int <- (mgcv::gam(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ) ~ 
                            s(log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu)) +
                            s(phytopl.All.Atl.PICO.p_OK$SST_C_insitu)+
                            +ti(log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu),phytopl.All.Atl.PICO.p_OK$SST_C_insitu)))
#plot(gam59)
summary(gam65.w.int)
# R-sq.(adj) =  0.917   Deviance explained =   96%
# SST is not sign. (!)
# interaction SST x Chl a  is not sign. (!)
# R-sq.(adj) =  0.917   Deviance explained =   96%


gam65.no.int <- (mgcv::gam(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ) ~ 
                            s(log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu)) +
                            s(phytopl.All.Atl.PICO.p_OK$SST_C_insitu)) )
                            
#plot(gam59)
summary(gam65.no.int)
# R-sq.(adj) =  0.917   Deviance explained =   96%
# SST is not sign. (!),  p = 0.487
# CHla is highly signif, p < 2e-16 ***
# 
#  R-sq.(adj) =  0.885   Deviance explained = 91.5%


# gam  model with Chl:
gam65.Chla <- (mgcv::gam(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ) ~ 
                          s(log10( 1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu))) )
summary(gam65.Chla)
# R-sq.(adj) =    0.4   Deviance explained = 42.5%, 4.85e-05 ***


# Best gam model (SST): Deviance explained = 91.3%
gam65.SST <- (mgcv::gam(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ) ~ 
                             s(phytopl.All.Atl.PICO.p_OK$SST_C_insitu)) )
summary(gam65.SST)

# 
# R-sq.(adj) =  0.887   Deviance explained = 91.3%
# p<2e-16 ***
#  N = 40 top picoplankton NBSS

    
    
  # * M.    Pico+Nano+Micropl. models * -------------
    # n = 279 datasets (49%) with 4 or more picoplankton non-empty bins 
    # redo all multivariate models with all size ranges, 
    # but only for 279 samples where PICOPL is present (at least 4 picopl data)
    
    
# M 1. Univariate linear models, n = 279 datasets, all Sizes 
#   only for 279 samples where PICOPL is present in at least 4 Picopl. bins 
# all Phytoplanktn size classes

plot(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
       log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu))
  
 lm71.Chl <-  lm (phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
              log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) 
summary(lm71.Chl)
 # Multiple R-squared:  0.4683,	Adjusted R-squared:  0.4664 
 # 274 DF,  p-value: < 2.2e-16

cor.test( phytopl.All.Atl.atleast.4.PICO$slope.clean ,
             log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu),
          method = "spearman")
# p-value: p-value < 2.2e-16

# Phytoplankton NBSS slopes vs Chla is SIGNIFFICANT!
# Chla explains 47% of the variablity
#  n = 279 samples where PICOPL is present (at least 4 picopl data)


plot(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
       (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))

lm71.Temp <-  lm (phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
                   (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)) 
abline(lm71.Temp)
summary(lm71.Temp)
#Multiple R-squared:  0.5314,	Adjusted R-squared:  0.5296 
# F-statistic:   305 on 1 and 269 DF,  p-value: < 2.2e-16
# Temp explains 53% of the variablity


cor.test( phytopl.All.Atl.atleast.4.PICO$slope.clean ,
          (phytopl.All.Atl.atleast.4.PICO$SST_C_insitu),
          method = "spearman")
# p-value: p-value < 2.2e-16

# Phytoplankton NBSS slopes vs Temp (SST)  is SIGNIFFICANT!
# Temp explains 53% of the variablity in Phytopl slope
#  n = 279 samples where PICOPL is present (at least 4 picopl data)


aovperm(lm71.Chl)
# p-value:  2e-04

aovperm(lm71.Temp)
# p-value:  2e-04


# M.2 multivariate Phytopl. models -------------
#   only for 279 samples where PICOPL is present in at least 4 Picopl. bins 
# all Phytoplanktn size classes

# BEST LINEAR MODEL:  no interaction (Adjusted R-squared:  0.62)
#   only for 279 samples where PICOPL is present in at least 4 Picopl. bins 

m75.no.int <- (lm(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
             log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)+
               phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
summary(m75.no.int)
# Multiple R-squared:  0.6226,	Adjusted R-squared:  0.6197 
# F-statistic: 218.6 on 2 and 265 DF,  p-value: < 2.2e-16
# explains 62 % of the variablity in   phytopl. slope !

aovperm(m75.no.int)
# SST is signifficant  p < 2e-16 ***
# Chla is also signifficant p = 3.63e-14 ***
# (no interction model) 
# Overall p-value:   < 2.2e-16

# Including another variable improved the adjusted R² and the AIC of the model.

relaimpo::calc.relimp(m75.no.int)
# lmg
# log10(1 + Chlorophyll__mg_m_3_insitu) 0.2787 (28 % is explained by Chl a)
# SST_C_insitu                          0.343 (34 % is explained by SST)
# 



# linear model with interaction   ----------
#   only for 279 samples where PICOPL is present in at least 4 Picopl. bins 
# , all Phytoplanktn size classes

# w interaction (Adjusted R-squared:  0.06)
m75.w.int <- (lm(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
                    log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)*
                    phytopl.All.Atl.atleast.4.PICO$SST_C_insitu))
summary(m75.w.int)
# Multiple R-squared:  0.6231,	Adjusted R-squared:  0.6188 
# F-statistic: 145.5 on 3 and 264 DF,  p-value: < 2.2e-16
# explains 62 % of the variablity in   phytopl. slope !

aovperm(m75.w.int)
# SST is signifficant  p < 2e-16 ***
# Chla is also signifficant p =0.022
# interaction  is  NOT sig. (prefer the no. int. model !) 
# Overall p-value:   < 2.2e-16

# Including interaction did not improve the adjusted R² and the AIC of the model.
# interaction  is  NOT sign. .. Conclusion: prefer the no. int. model !)

relaimpo::calc.relimp(m75.w.int)
# lmg
# log10(1 + Chlorophyll__mg_m_3_insitu) 0.2787 (28 % is explained by Chl a)
# SST_C_insitu                          0.343 (34 % is explained by SST)
# interaction explains  only 0.04 % of the variability ! 
# Conclusion: prefer the model without  interaction !


# Including another variable with interaction improved the adjusted R² and the AIC of the model.

# SST is not sign. (!)
# Chl a is signif., p = 0.0036
# Interaction Chl a x SST i signiff. ( p = 0.046)

relaimpo::calc.relimp(m57)
# lmg
# log10(1 +Chlorophyll__mg_m_3_insitu)                  3.4%      0.03412297
# SST_C_insitu                                          3.8%    0.03878442
# log10(1 + Chlorophyll__mg_m_3_insitu) :SST_C_insitu   1.4% 0.01436110
# non-linear models with interction  


# M.3 non-linear models (GAM), all Phytoplanktn size classes ----------
#   only for 279 samples where PICOPL is present in at least 4 Picopl. bins 

gam78bChla <- (mgcv::gam(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
                       s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu) ))) 

summary(gam78bChla)
#R-sq.(adj) =  0.543   Deviance explained = 55.3%
# p = <2e-16 ***

plot(gam78bChla)


gam78bSST <- (mgcv::gam(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
                       s((phytopl.All.Atl.atleast.4.PICO$SST_C_insitu) ))) 
summary(gam78bSST)
# R-sq.(adj) =   0.65   Deviance explained = 66.1%
# p = <2e-16 ***
plot(gam78bSST)


# GAM, no interaction
# only for 279 samples where PICOPL is present in at least 4 Picopl. bins 

gam79.no.int <- (mgcv::gam(phytopl.All.Atl.atleast.4.PICO$slope.clean ~ 
                             s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) +
                             s(phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
                             
summary(gam79.no.int)
#R-sq.(adj) =  0.724   Deviance explained = 73.8%
AIC(gam79.no.int)
# -223.9437

# GAM, with interaction
# only for 279 samples where PICOPL is present in at least 4 Picopl. bins 

gam79.with.int <- (mgcv::gam(phytopl.All.Atl.atleast.4.PICO$slope.clean ~      
                          s(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu)) +
                           s(phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)+
                            +ti(log10(1+phytopl.All.Atl.atleast.4.PICO$Chlorophyll__mg_m_3_insitu),phytopl.All.Atl.atleast.4.PICO$SST_C_insitu)))
                        
summary(gam79.with.int)
# p-value
# s(log10(1 + Chlorophyll__mg_m_3_insitu))                     p =   0.0468
# s(SST_C_insitu)                                               p   <2e-16
# interaction                                                   0.2122 (not sign.)



# R-sq.(adj) =  0.742   Deviance explained = 76.1%
# R squared is not kch better  with interaction than without.
# Conclusion: the BEST lineaer model is  without interaction , 
# the  best GAM model is with interaction r
AIC(gam79.with.int)

# Model gam79.with.int (with AIC = -236.6408) is better,
# because it has a lower AIC than Model gam79.no.int.
# Also the GAM model with interaction has a better R2
# However, in the GAM model the interaction was also NOT signif. !!!

#  A rule of thumb: A difference in AIC > 2 is often considered meaningful. 
# In this case, the difference is about 12.7 


# N. Picopl. plots FOR PAPER !!! ----------------- 

summary(phytopl.All.Atl.PICO.p_OK)
dim(phytopl.All.Atl.PICO.p_OK)

# plot WITH CIRCLES 
plot(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
       log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) )
# , col = phytopl.All.Atl.PICO.p_OK$SST_C_insitu)
symbols((phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu),
        phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes,
        circles  = (phytopl.All.Atl.PICO.p_OK$SST_C_insitu/3000), inches = 0.1)
abline(h = -1)

symbols(log10(phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu),
        phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes,
        circles  = (phytopl.All.Atl.PICO.p_OK$SST_C_insitu/3000), inches = 0.1)
abline(h = -1)



## PLOT WITH circles! 
plot(as.numeric(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes) ~ 
       (phytopl.All.Atl.PICO.p_OK$SST_C_insitu) )
symbols((phytopl.All.Atl.PICO.p_OK$SST_C_insitu),
        phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes,
        circles  = log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu/3000), inches = 0.3)

abline(h = -1)
#abline(h = -1.4)
abline(v = 18.4)


# nbss picopl plot
# select Picopl spectra warm and cood waters (< 20 degrees Celsius >)

phytopl.All.Atl.PICO.p_OK.warm <- subset (phytopl.All.Atl.PICO.p_OK,phytopl.All.Atl.PICO.p_OK$SST_C_insitu   > 18.5 )
dim(phytopl.All.Atl.PICO.p_OK.warm)
# 91 samples

phytopl.All.Atl.PICO.p_OK.cold <- subset (phytopl.All.Atl.PICO.p_OK,phytopl.All.Atl.PICO.p_OK$SST_C_insitu   < 18.5 )
dim(phytopl.All.Atl.PICO.p_OK.cold)
# 20 samples

phytopl.All.Atl.PICO.p_OK.oligotr  <- subset (phytopl.All.Atl.PICO.p_OK, phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu < (1) )
dim(phytopl.All.Atl.PICO.p_OK.oligotr)
# 100 samples

phytopl.All.Atl.PICO.p_OK.mesotr  <- subset (phytopl.All.Atl.PICO.p_OK, phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu > (1) )
dim(phytopl.All.Atl.PICO.p_OK.mesotr)
# 17 samples




summary(phytopl.All.Atl.PICO.p_OK)
dim(phytopl.All.Atl.PICO.p_OK)
# 119 samples are OK

summary(phytopl.All.Atl.PICO.p_OK$SST_C_insitu)[]  # 8 NA's

summary(phytopl.All.Atl.PICO.p_OK.warm)
dim(phytopl.All.Atl.PICO.p_OK.warm)
# 92 WARM WATER PICOPL OK SAMPLES

summary(phytopl.All.Atl.PICO.p_OK.cold)
dim(phytopl.All.Atl.PICO.p_OK.cold)
# 20 COLD WATER PICOPL OK SAMPLES

summary(phytopl.All.Atl.PICO.p_OK)
dim(phytopl.All.Atl.PICO.p_OK)
# 119 PICOPL OK SAMPLES

P.DATA <- phytopl.All.Atl.PICO.p_OK # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED
  
phytopl.All.Atl.PICO.p_OK.MATRIX <- as.matrix(P.DATA[,38:92])
dim(phytopl.All.Atl.PICO.p_OK.MATRIX)
names(phytopl.All.Atl.PICO.p_OK[,38:92])


# 3D plots 


# 3d plots --------------------

library(car)
scatter3d(x = phytopl.All.Atl.PICO.p_OK$SST_C_insitu, y =  phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes, 
          z = log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu), 
          grid = FALSE, fit = "smooth", point.col="darkorange") # bg = "black")

# test for collnearity 
cor.test (log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu), phytopl.All.Atl.PICO.p_OK$SST_C_insitu) 

# 3d surface is not twisted! NO visible interactions for Carbon NBSS


# check the effect of chla variabllity  in cold waters....

m103cold_chl.effect <- lm (   phytopl.All.Atl.PICO.p_OK.cold$Pico.rr.slopes ~     log10(1+phytopl.All.Atl.PICO.p_OK.cold$Chlorophyll__mg_m_3_insitu)) 
summary(m103cold_chl.effect)
# NOT sign. !! CHla has no effect, in cold waters  !!!

phytopl.All.Atl.PICO.p_OK.oligotr  <- subset (phytopl.All.Atl.PICO.p_OK, phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu < (1) )
dim(phytopl.All.Atl.PICO.p_OK.oligotr)
# 100 samples


phytopl.All.Atl.PICO.p_OK.mesotr  <- subset (phytopl.All.Atl.PICO.p_OK, phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu > 1 )
dim(phytopl.All.Atl.PICO.p_OK.mesotr)
# 17 samples


fun.inverse_log10_1plusx <- function(y) {  10^y - 1 }

fun.log10_1plusx <- function(x) {  log10 (1+x) }

fun.inverse_log10_1plusx(0) 
# 0.3  -> 1 log10(1+x) 
#      log10(1+1)  = 0.3
#       log10( 1 ) = 0 

log10(2) # 0.3
log10(1) # 0



dim(phytopl.All.Atl.PICO.p_OK.oligotr)

m103ologtophicsamples_SST.effect <- lm (   phytopl.All.Atl.PICO.p_OK.oligotr$Pico.rr.slopes ~ 
                                             (phytopl.All.Atl.PICO.p_OK.oligotr$SST_C_insitu)) 

summary(m103ologtophicsamples_SST.effect)
# NSST is sign. !! SST  has an  effect, even for  oligtrophic waters !!!
# p-value: 2.458e-14


m103mesotophicsamples_SST.effect <- lm (   phytopl.All.Atl.PICO.p_OK.mesotr$Pico.rr.slopes ~ 
                                             (phytopl.All.Atl.PICO.p_OK.mesotr$SST_C_insitu)) 

summary(m103mesotophicsamples_SST.effect)
# SST is not sign, low N..., only 17 samples,  and all COLD!, all below 18.6 derees Cellsius
summary(phytopl.All.Atl.PICO.p_OK.mesotr$SST_C_insitu)



dim(phytopl.All.Atl.PICO.p_OK.cold)

m103low_SST.effect <- lm (   phytopl.All.Atl.PICO.p_OK.cold$Pico.rr.slopes ~     log10(1+phytopl.All.Atl.PICO.p_OK.cold$Chlorophyll__mg_m_3_insitu)) 
summary(m103cold_chl.effect)
# NOT sign. !! CHla has no effect, in cold waters  !!!
# in cold water, picopl.  NBSS was  flat, regrdless of Chla !!!, in spite of large vaROVLTY IN CHLa
summary(phytopl.All.Atl.PICO.p_OK.cold$Chlorophyll__mg_m_3_insitu)

summary(phytopl.All.Atl.PICO.p_OK.cold$Pico.rr.slopes)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.9645 -1.0561 -0.8453 -0.5293  0.3207  1.0685 
# For Results text:   median:  -0.845
# (-1.96 to 1.07, median: b = -0.845)



  
dim(phytopl.All.Atl.PICO.p_OK.warm) # 91 samples

summary(phytopl.All.Atl.PICO.p_OK.warm$Pico.rr.slopes)
# Min.    1st Qu.  Median    Mean  3rd Qu.    Max. 
# -3.851  -3.416  -2.892  -2.897  -2.410    -1.483 
# For Results text:   median:  -2.89
# (-3.85 to -1.48, median: b = -2.89)
# varying from (-3.85 to -1.48, with an extremely steep median NBSS slope of b = -2.89


m103high_SST.effect <- lm (   phytopl.All.Atl.PICO.p_OK.warm$Pico.rr.slopes ~   
                                log10(1+phytopl.All.Atl.PICO.p_OK.warm$Chlorophyll__mg_m_3_insitu)) 
summary(m103high_SST.effect)
# NOT sign. !! CHla has little effect, in warm waters  !!! expains only 15% variablity
plot( (   phytopl.All.Atl.PICO.p_OK.warm$Pico.rr.slopes ~   
            log10(1+phytopl.All.Atl.PICO.p_OK.warm$Chlorophyll__mg_m_3_insitu)) )





# PLot ALL picop  NBSS Data, Carbon  -------

nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) # 40 high-quality datasets, OK

phytopl.All.Atl.PICO.p_OK$CruiseID 
phytopl.All.Atl.PICO.p_OK$StationID 

# MAP the 40 PICOPL samples
 
 # library(maps)
 #     maps::map(main = "all PHYTOPL. SAMPLES, N = 565" ,"world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), xlim=c(-95, 30), mar=c(0,0,0,0))
 #  points( phytopl.All.Atl.p.rsq.slope_ratio.ok$Longitude,
 #         phytopl.All.Atl.p.rsq.slope_ratio.ok$Latitude,   col = alpha("darkgreen", 0.2), pch=16)

   # dev.off()

  # library(maps)
  # maps::map(main = "all picopl. ok SAMPLES, N = " ,
  #           "world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), xlim=c(-95, 30), mar=c(0,0,0,0))
  # points( phytopl.All.Atl.PICO.p_OK$Longitude,
  #         phytopl.All.Atl.PICO.p_OK$Latitude,   col = alpha("red", 0.5), pch=19, cex = 2)

    # dev.off()

  # library(maps)
  # maps::map(main = "N = 91 picopl. WARM SAMPLES" ,
  #           "world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), xlim=c(-95, 30), mar=c(0,0,0,0))
  # points( phytopl.All.Atl.PICO.p_OK.warm$Longitude,
  #         phytopl.All.Atl.PICO.p_OK.warm$Latitude,   col = alpha("red", 0.5), pch=19, cex = 2)

  # dev.off()

# 
#   library(maps)
#   maps::map(main = "N = 20 picopl. COLD SAMPLES" ,
#             "world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), xlim=c(-95, 30), mar=c(0,0,0,0))
#   points( phytopl.All.Atl.PICO.p_OK.cold$Longitude,
#           phytopl.All.Atl.PICO.p_OK.cold$Latitude,   col = alpha("red", 0.5), pch=19, cex = 2)
#  
  # dev.off()

  
  
  
############ NBSS PLOT ########
### COMPLEX NBSS PLOT - ALL PICOPLANKTON  --------- 
# plot all points and rlm moldes (blue lines), and red median line  ----------


P.DATA <- phytopl.All.Atl.PICO.p_OK # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED


lnumb <- 3
 
 mod1<-  lm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
abline(mod1) 
abline(v = -12) 

 
plot (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind), 
      main = "ATLANTIC, PICOPLANKTON, n = 40",
      xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
      ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
      xlim = c(-15, -5), ylim = c(0, 14) ,
      pch = 16, col = "white")

# Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 

abline(a = -3.2915, b = -0.9548 ,  
       lwd = 2.5, lty = 2, col = "grey44" )


# phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727b[,38:92])
# dim(phytopl.NBSS.matrix.unfiltr727)
# 
# 
# for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
# {    
#   points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
#           pch = 16, col = alpha("lightblue", 0.5))
# }



for (lnumb in 1 : nrow(P.DATA) ) 
{ 
  # abline( rlm(log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind)), 
  #         col =alpha ("cornsilk4", 0.1))
  
  mod1<-  rlm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~
                 log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)), na.action = "na.omit")
  abline( mod1, col =alpha ("cornsilk4", 0.1))
  
}

for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
{    
  points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("lightgrey", 0.4))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10(phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("lightblue", 0.5))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,1:5]) ~ log10(X_vector_gCind[1:5]) ,
          pch = 16, col = alpha("forestgreen", 0.7))
}

abline (v = -12)


# Carbon  (458 high-quality datasets, OK)
(median( P.DATA$Pico.rr.slopes) )
dim(P.DATA)
summary(P.DATA$Pico.rr.slopes)

dim(phytopl.All.Atl.PICO.p_OK.MATRIX)

# matrix to to vector (transpose),  -----------------------  
x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) )
y1 <- c ( t(phytopl.All.Atl.PICO.p_OK.MATRIX))

# regression line PICOPL, Rob. regr.

med.slope.rr.PICO <- median(P.DATA$Pico.rr.slopes)
med.interc.rr.PICO <- median(P.DATA$Pico.rr.interc)

abline( a = med.interc.rr.PICO, b = med.slope.rr.PICO)

# regression line (Rob. regr.)
#plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
#points(  log10(y1) ~ x1 )
# abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
       lwd = 2.5, lty = 2, col = "darkorange" )


text(x = -12, y = 2, "n= 40, RR,  slope = -1.06", col = "darkorange"  )
#text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )


# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

mean_   <- function(...) mean(..., na.rm=T)
median_ <- function(...) median(..., na.rm=T)

mean_ (c(4,NA, 8, 34)) #15.33, OK
median_ (c(4,NA, 8,34)) #8, OK


# calculate median (or mean) only if less than 50% are NA ------------

vec1 <- c(4, 5,NA,NA, NA,  34)

sum(is.na(vec1))
percNA <- sum(is.na(vec1))/ length(vec1)


mean_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    mean(..., na.rm=T) }
  
  else{ ;NA }
}

median_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    median(..., na.rm=T) }
  
  else{ ;NA }
}

vec2 <- c(4, 5,NA,NA, NA,NA,  34)
vec3 <- c(4, 5,NA,34)

mean_NA50(vec2)
mean_NA50(vec3)


# apply median by columns

median_vec_ALLBiovol <- apply( phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_) 
mean_vec_ALLBiovol <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, mean_) 


# lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)


# apply mean by columns (only if less thn 50% NAs)

# mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 

median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_NA50) 

fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }

max_vec_ALL <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, fun95_quant ) 

length(median_vec_ALLBiovol2)

class(median_vec_ALLBiovol2)

median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))

#mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA



median_vec_PHYTOP <-  median_vec_ALLBiovol2 

X_vector_gCind_PHYTOPL <- X_vector_gCind


df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))

dfok<- df1
#dfok <- df1[complete.cases(df1$y1),]
#plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)

dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
        col = "purple" , lwd = 4.5)


############ END OF NBSS PLOT ########




############ NBSS PLOT ########
# NBSS FIG. PICOPL. OK (p <0.05) WARM Waters ( SST in situ < 18.4 degrees calsius) only, n = 15 samples --------------

P.DATA <- phytopl.All.Atl.PICO.p_OK.warm # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED

dim(phytopl.All.Atl.PICO.p_OK.warm) # 91 samples
  # plot all points and rlm moldes (blue lines), and red median line  ----------



phytopl.All.Atl.PICO.p_OK.MATRIX <- as.matrix(P.DATA[,38:92])

lnumb <- 3

mod1<-  lm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
abline(mod1) 
abline(v = -12) 


plot (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind), 
      main = "ATLANTIC, PICOPLANKTON, n = 91, warm wates",
      xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
      ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
      xlim = c(-15, -5), ylim = c(0, 14) ,
      pch = 16, col = "white")

# Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 

abline(a = -3.2915, b = -0.9548 ,  
       lwd = 2.5, lty = 2, col = "grey44" )


# phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727b[,38:92])
# dim(phytopl.NBSS.matrix.unfiltr727)
# 
# 
# for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
# {    
#   points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
#           pch = 16, col = alpha("lightblue", 0.5))
# }



for (lnumb in 1 : nrow(P.DATA) ) 
{ 
  # abline( rlm(log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind)), 
  #         col =alpha ("cornsilk4", 0.1))
  
  mod1<-  rlm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~
                 log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)), na.action = "na.omit")
  abline( mod1, col =alpha ("cornsilk4", 0.4))
  
}

for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
{    
  points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("lightgrey", 0.2))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10(phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("pink", 0.8))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,1:5]) ~ log10(X_vector_gCind[1:5]) ,
          pch = 16, col = alpha("forestgreen", 0.7))
}

abline (v = -12)


# Carbon  (458 high-quality datasets, OK)
(median( P.DATA$Pico.rr.slopes) )
dim(P.DATA)
summary(P.DATA$Pico.rr.slopes)

dim(phytopl.All.Atl.PICO.p_OK.MATRIX)

# matrix to to vector (transpose),  -----------------------  
x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) )
y1 <- c ( t(phytopl.All.Atl.PICO.p_OK.MATRIX))

# regression line PICOPL, Rob. regr.

med.slope.rr.PICO <- median(P.DATA$Pico.rr.slopes)
med.interc.rr.PICO <- median(P.DATA$Pico.rr.interc)

abline( a = med.interc.rr.PICO, b = med.slope.rr.PICO, lwd = 2.5, lty = 2, col = "forestgreen")

# regression line (Rob. regr.)
#plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
#points(  log10(y1) ~ x1 )
# abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
       lwd = 2.5, lty = 2, col = "darkorange" )


text(x = -12, y = 2, "n= 40, RR,  slope = -1.06", col = "darkorange"  )
#text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )


# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

mean_   <- function(...) mean(..., na.rm=T)
median_ <- function(...) median(..., na.rm=T)

mean_ (c(4,NA, 8, 34)) #15.33, OK
median_ (c(4,NA, 8,34)) #8, OK


# calculate median (or mean) only if less than 50% are NA ------------

vec1 <- c(4, 5,NA,NA, NA,  34)

sum(is.na(vec1))
percNA <- sum(is.na(vec1))/ length(vec1)


mean_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    mean(..., na.rm=T) }
  
  else{ ;NA }
}

median_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    median(..., na.rm=T) }
  
  else{ ;NA }
}

vec2 <- c(4, 5,NA,NA, NA,NA,  34)
vec3 <- c(4, 5,NA,34)

mean_NA50(vec2)
mean_NA50(vec3)


# apply median by columns

median_vec_ALLBiovol <- apply( phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_) 
mean_vec_ALLBiovol <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, mean_) 


# lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)


# apply mean by columns (only if less thn 50% NAs)

# mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 

median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_NA50) 

fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }

max_vec_ALL <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, fun95_quant ) 

length(median_vec_ALLBiovol2)

class(median_vec_ALLBiovol2)

median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))

#mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA



median_vec_PHYTOP <-  median_vec_ALLBiovol2 

X_vector_gCind_PHYTOPL <- X_vector_gCind


df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))

dfok<- df1
#dfok <- df1[complete.cases(df1$y1),]
#plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)

dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
        col = "purple" , lwd = 4.5)


############ END OF NBSS PLOT ########



############ NBSS PLOT ########
# NBSS FIG. PICOPL. OK (p <0.05) COLD Waters ( SST in situ < 18.4 degrees calsius) only, n = 15 samples --------------

P.DATA <- phytopl.All.Atl.PICO.p_OK.cold # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED

dim(phytopl.All.Atl.PICO.p_OK.cold) # 20 samples
# plot all points and rlm moldes (blue lines), and red median line  ----------



phytopl.All.Atl.PICO.p_OK.MATRIX <- as.matrix(P.DATA[,38:92])

lnumb <- 3

mod1<-  lm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~ log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)))
abline(mod1) 
abline(v = -12) 


plot (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind), 
      main = "ATLANTIC, PICOPLANKTON, n = 20, cold waters",
      xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
      ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
      xlim = c(-15, -5), ylim = c(0, 14) ,
      pch = 16, col = "white")

# Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 

abline(a = -3.2915, b = -0.9548 ,  
       lwd = 2.5, lty = 2, col = "grey44" )


# phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727b[,38:92])
# dim(phytopl.NBSS.matrix.unfiltr727)
# 
# 
# for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
# {    
#   points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
#           pch = 16, col = alpha("lightblue", 0.5))
# }



for (lnumb in 1 : nrow(P.DATA) ) 
{ 
  # abline( rlm(log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind)), 
  #         col =alpha ("cornsilk4", 0.1))
  
  mod1<-  rlm( log10(as.numeric(P.DATA[lnumb, 38:(38+4)])) ~
                 log10(c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13)), na.action = "na.omit")
  abline( mod1, col =alpha ("cornsilk4", 0.6))
  
}

for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
{    
  points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("lightgrey", 0.2))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10(phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind) ,
          pch = 16, col = alpha("pink", 0.7))
}


for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
{    
  points (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,1:5]) ~ log10(X_vector_gCind[1:5]) ,
          pch = 16, col = alpha("forestgreen", 0.7))
}

abline (v = -12)


# Carbon  (458 high-quality datasets, OK)
(median( P.DATA$Pico.rr.slopes) )
dim(P.DATA)
summary(P.DATA$Pico.rr.slopes)

dim(phytopl.All.Atl.PICO.p_OK.MATRIX)

# matrix to to vector (transpose),  -----------------------  
x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) )
y1 <- c ( t(phytopl.All.Atl.PICO.p_OK.MATRIX))

# regression line PICOPL, Rob. regr.

med.slope.rr.PICO <- median(P.DATA$Pico.rr.slopes)
med.interc.rr.PICO <- median(P.DATA$Pico.rr.interc)

abline( a = med.interc.rr.PICO, b = med.slope.rr.PICO, lwd = 2.5, lty = 2, col = "forestgreen")

# regression line (Rob. regr.)
#plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
#points(  log10(y1) ~ x1 )
# abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
       lwd = 2.5, lty = 2, col = "darkorange" )


text(x = -12, y = 2, "n= 40, RR,  slope = -1.06", col = "darkorange"  )
#text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )


# means and medians without considering NAs  
mean (c(4,NA, 18,34)) # NA
median( c(4,NA, 18, 34)  ) #NA

mean_   <- function(...) mean(..., na.rm=T)
median_ <- function(...) median(..., na.rm=T)

mean_ (c(4,NA, 8, 34)) #15.33, OK
median_ (c(4,NA, 8,34)) #8, OK


# calculate median (or mean) only if less than 50% are NA ------------

vec1 <- c(4, 5,NA,NA, NA,  34)

sum(is.na(vec1))
percNA <- sum(is.na(vec1))/ length(vec1)


mean_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    mean(..., na.rm=T) }
  
  else{ ;NA }
}

median_NA50   <- function(...)  {
  
  percNA <- sum(is.na(...))/ length(...)
  
  if(percNA <= 0.5) {
    median(..., na.rm=T) }
  
  else{ ;NA }
}

vec2 <- c(4, 5,NA,NA, NA,NA,  34)
vec3 <- c(4, 5,NA,34)

mean_NA50(vec2)
mean_NA50(vec3)


# apply median by columns

median_vec_ALLBiovol <- apply( phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_) 
mean_vec_ALLBiovol <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, mean_) 


# lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)


# apply mean by columns (only if less thn 50% NAs)

# mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 

median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_NA50) 

fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }

max_vec_ALL <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, fun95_quant ) 

length(median_vec_ALLBiovol2)

class(median_vec_ALLBiovol2)

median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))

#mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA



median_vec_PHYTOP <-  median_vec_ALLBiovol2 

X_vector_gCind_PHYTOPL <- X_vector_gCind


df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))

dfok<- df1
#dfok <- df1[complete.cases(df1$y1),]
#plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")

lines(  log10(dfok$y1) ~ log10(dfok$x1), 
        col = "red" , lwd = 4.5)

dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
        col = "purple" , lwd = 4.5)


############ END OF NBSS PLOT ########


# Sizerange of Picoplankton: -13.6 to -11.6 log₁₀ gC cell⁻¹ 
# Organism	         Carbon per cell	           log₁₀ gC cell⁻¹
# Synechococcus      	0.1 – 1.5 pg            	–13.00 to –11.82
# Prochlorococcus	    66 – 158 fg              	–13.18 to –12.80
# 
#  Picoeukaryotes (fixed)	1.978×10⁻¹²            	−11.70
# Picoeukaryotes	24–2590 fg	                 –13.62 to –11.59

# For paper

library(ggplot2)
df2 <- data.frame( NBSS_slope = phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes , 
                   Chla = (phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu),
                   log1plusChla = log10(1+ phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu))

ggplot(df2, aes(log1plusChla, NBSS_slope ) ) +
  scale_x_continuous(name ="Chl a (log10 (x+1, mg m-3)), in situ", label = comma,
                     limits=c(0, 1))+
  theme(axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10) ,
        axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10) )+ 
  geom_point( alpha = 0.3, col = "navy", size = 2.5) +
  geom_smooth(method = "lm")

ggplot(df2, aes(Chla, NBSS_slope ) ) +
  scale_x_continuous(name ="Chl a (mg m-3)), in situ", label = comma,
                     limits=c(0, 8))+
  theme(axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10) ,
        axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10) )+ 
  geom_point( alpha = 0.3, col = "navy", size = 2.5)  +
  geom_smooth() # (method = "lm")




# For paper

library(ggplot2)
df2 <- data.frame( NBSS_slope = phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes , 
                   SST = (phytopl.All.Atl.PICO.p_OK$SST_C_insitu))

ggplot(df2, aes(SST, NBSS_slope ) ) +
  scale_x_continuous(name ="SST, in situ", label = comma,
                     limits=c(15, 30))+
  theme(axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10) ,
        axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10) )+ 
  geom_point( alpha = 0.3, col = "navy", size = 2.5) +
  geom_smooth(method = "lm")









   # [bookmark]
   
   
# II.2 Zooplankton ------------
# II.2  recalculate slopes for Zooplankton and recalculate with Robust regression 
# Recalculate  for zooplankton (use data within pre-defined size ranges only) ----------------

# UVP vs Net zoopl ------------------------------------
   
   
#   5.1  Cut UVP data at limits -4.9 to -2 log10 (g C ind.^-1) -------------------
   
   10^-4.9 # 1.26 e-05 gC
   10^-2 # 0.01 gC
  
   # X_vector_gCind = Bin mids for NBSS, in g C ind.^-1
   
   X_vector_gCind # 55 size bins g C ind.^-1, all size bins
   X_vector_gCind_UVP <-  X_vector_gCind [X_vector_gCind > 10^-4.9] # 1.5300e-05 to 1.0240e+03
   X_vector_gCind_UVP <-  X_vector_gCind_UVP [X_vector_gCind_UVP <  10^-2] # 1.5300e-05 to 1.0240e+03
   X_vector_gCind_UVP # 10 classes for UVP, 0.0000153  to 0.0078130 g C ind.^-1
   log10(0.0000153)
   log10( 0.0078130)
   X_vector_gCind[29:38] # # 10 classes for UVP, 0.0000153  to 0.0078130 g C ind.^-1
   
   length(38:92)
   length(X_vector_gCind)
      
   dim(zoopl.All.Atl)
   # 1579 samples, 100 variables, variables 37:91 are the NBSS data
    # write.csv(zoopl.All.Atl, file = "zoopl_All_Atl_n1580_b.csv") 
   
   names(zoopl.All.Atl [37:91]) # all NBSS (n = 55 size bins)
   names(zoopl.All.Atl [(37+28):(37+37)]) # 10 classes UVP OK NBSS (n = 10 size bins)
   names(zoopl.All.Atl [(37):(37+27)]) # NBSS smaller than UVP NBSS (n = 10 size bins)
   names(zoopl.All.Atl [(37+38):91]) # NBSS larger than UVP NBSS (n = 10 size bins)
   names(zoopl.All.Atl [  c( 37:(37+27), (37+38):91 )] )   # 45 NBSS classes not used for UVP NBSS (n = 10 size bins)
   
   # zoopl.All.Atl
   
   
   zoopl.All.Atl.NBSS.data <-  zoopl.All.Atl[, 37:91]
  
       
       
   
   # CLEANUP - delete UVP outside size range -4.9 to -2 log10 (g C ind.^-1)
   # clean all UVP NBSS data outside the size range
   
   # LOOP (loops through the column and replaces NBSS data with NA)
   # for UVP data only
   # # 45 NBSS classes not used for UVP NBSS
   
   df <- zoopl.All.Atl
   
   indices_for_UVPcleanup <- c( 37:(37+27), (37+38):91 )
   
   for (i in indices_for_UVPcleanup) { 
     
    
     df[,i]   [df$Gear == 'UVP'] <- NA
     
     
   }   
   
   # View(df)
   
   zoopl.All.Atl <- df
   
  # UVP data finished, cleanup ok

       
   # 5.2 Calculate Zooplankton NBSS slopes-------------------------------------------
   
   
   
   # Calculate new slopes from maximum to first non-empty bin -----------------
   
   # Define Functions for Zooplankton NBSS slopes (OLSR and Rob. Regr.) ------------
   #  Define  functions that select from maximum to last non_empty bin -----------
   # and fit a linear model
   # Inputs: two vectors of the same length, "bin.means", "raw.NBSS"
   # Output:  results  with "NBSS.slope", "NBSS.intercept", 
   # "NBSS.p.value.of.Intercept",  "NBSS.p.value.of.Slope", "r-squared", 
   # "first.bin.used". "last.bin.used" , "N.of.bins.used" , etc.
   
   
   # define function NBSS.select.w.lm.MS.A , with OLSR lm only----------
   
   NBSS.select.w.lm.MS.A   <- function (bin.means = bin.means, raw.NBSS = raw.NBSS) {
     
     # 
     # bin.means = c(2^1, 2^2, 2^3, 2^4, 2^5 , 2^6, 2^7, 2^8, 2^9)
     # raw.NBSS <- c (NA, 10,   28,  3,  0.4,  0.04 ,  NA,  3, NA)  
     
     
     #  mean individual biomass (g.C ind.-1) in each size class (X-Axis)
     mean.ind.biomass  <- bin.means
     # Comment: a geometric series whose terms are the successive powers of two 
     
     
     #  Normalized Biomass (NBSS)  in each size class (y-Axis) unit: "gC * m-3 * gC-1 * ind."   or simply "ind. m-3"
     NBSS  <- raw.NBSS
     # Comment: each value represents the normalized biomass in each size class
     
     d1_all <- data.frame(mean.ind.biomass, NBSS)
     
     n.NBSS.withdata <-  length (which(!is.na(  ( zoopl.NBSSmatrix_no_na[lnumb,]))))
     # number of non-NA values
     
     
     #plot (log10(mean.ind.biomass), log10(NBSS))  
     # 3. Automatic selection (from the maximum to the last non-empty bin before empty bins occur) ------------------------------
     
     # 3.0 preparations -------------------------------------------------------------
     
     # 3.0.a create log10 data, for convenience --------------------------------------
     z <- log10(NBSS)
     z[z== -Inf] <- NA   # replace "-Inf" with NA
     
     # 3.0.b create an index vector (1,2,3,..) with the positions of the data,  for convenience ----------
     
     d1_index <- 1:(length(z))
     d1_index
     #  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14
     
     # 3.0.b build a data.frame (a table)
     
     d1_all <- data.frame(d1_index, mean.ind.biomass, NBSS, z )
     
     # 3.1 find the maximum 
     
     d_1MAXindex <- which.max(z) # position of the maximum
     
     d_1MAXindex # 3 (the maximum is located at the position 3)
     
     # 3.2 select only the data from the maximum onwards
     
     d1_all.2 <- d1_all[ d1_all$d1_index[d1_all$d1_index > (d_1MAXindex-1) ] ,]
     
     # 3.3 find the last non-empty bin before empty bins occur
     
     NAindex <- which(is.na(d1_all.2$z)) # locate ALL NA values in vector z
     firstNA <- min(NAindex)# locate the first NA value in vector z 
     firstNA # position of first NA value in vector z
     LastNonNA <-firstNA - 1 # position of last non_NA value in vector z
     LastNonNA  # position last non-empty bin before empty bins occur, index as for d1_all.2
     
     # 3.4. select the data (based on automatically selected indices)
     
     d1_all.3 <- d1_all.2[ 1:LastNonNA ,]
     
     # position last non-empty bin before empty bins occur, index as for d1_all.2
     #ALL.LastNonNA <- max(d1_all.3$d1_index) # maximum index 
     # position last non-empty bin before empty bins occur, index as for d1_all.2
     
     
     lm.sel <- lm(log10(d1_all.3$NBSS) ~ log10(d1_all.3$mean.ind.biomass))
     summary(lm.sel)
     coef(lm.sel)
     
     # abline(lm.sel)
     
     NBSS.Intercept<- as.numeric(coef(lm.sel)[1]) # NBSS Intercept
     NBSS.slope <- as.numeric(coef(lm.sel)[2] )# NBSS Slope
     
     sum.table <- as.list (summary(lm.sel))
     
     p.values <- summary(lm.sel)$coefficients[,4]  # p-values 
     r.squared <- summary(lm.sel)$r.squared # R-squared
     
     N.bins.used <- length(d1_all.2[ 1:LastNonNA ,]$z)
    
     
     # res <- c(NBSS.slope, NBSS.Intercept, p.values, r.squared, d_1MAXindex, 
     #          ALL.LastNonNA, N.bins.used)
     
     
     results <- list(slope = NBSS.slope, intercept = NBSS.Intercept, p.values = p.values, 
                     r.squared = r.squared, index.of.maximum = d_1MAXindex, 
                     x.value.of.maximum = d1_all$mean.ind.biomass[d_1MAXindex],
                     y.value.of.maximum =  d1_all$NBSS[d_1MAXindex],
                     index.of.start.bin <- d_1MAXindex, 
                     x.value.of.start.bin = d1_all$mean.ind.biomass[d_1MAXindex],
                     #endpoint = ALL.LastNonNA, 
                     # x.value.of.end.bin = d1_all$mean.ind.biomass[ALL.LastNonNA],
                     N.bins.used = N.bins.used,
                     n.NBSS.withdata = n.NBSS.withdata)
     
     
     names(results)
     
     res <- list(results = results, sum.table = sum.table, selected.data = d1_all.3 , all.data = d1_all)
     
     
     
     return(res)
     
     
   }
   
   #  Function NBSS.select.w.lm.MS.B.RR , OLSR and Robust Regression------------------
   # with robust regression, output includes Rsquared, "p" values
   
   NBSS.select.w.lm.MS.B.RR   <- function (bin.means = bin.means, raw.NBSS = raw.NBSS) {
     
     # 
     # bin.means = c(2^1, 2^2, 2^3, 2^4, 2^5 , 2^6, 2^7, 2^8, 2^9)
     # raw.NBSS <- c (NA, 10,   28,  3,  0.4,  0.04 ,  NA,  3, NA)  
     
     
     #  mean individual biomass (g.C ind.-1) in each size class (X-Axis)
     mean.ind.biomass  <- bin.means
     # Comment: a geometric series whose terms are the successive powers of two 
     
     
     #  Normalized Biomass (NBSS)  in each size class (y-Axis) unit: "gC * m-3 * gC-1 * ind."   or simply "ind. m-3"
     NBSS  <- raw.NBSS
     # Comment: each value represents the normalized biomass in each size class
     
     d1_all <- data.frame(mean.ind.biomass, NBSS)
     
     n.NBSS.withdata <-  length (which(!is.na(  ( zoopl.NBSSmatrix_no_na[lnumb,]))))
            # number of non-NA values
     
     
     #plot (log10(mean.ind.biomass), log10(NBSS))  
     # 3. Automatic selection (from the maximum to the last non-empty bin before empty bins occur) ------------------------------
     
     # 3.0 preparations -------------------------------------------------------------
     
     # 3.0.a create log10 data, for convenience --------------------------------------
     z <- log10(NBSS)
     z[z== -Inf] <- NA   # replace "-Inf" with NA
     
     # 3.0.b create an index vector (1,2,3,..) with the positions of the data,  for convenience ----------
     
     d1_index <- 1:(length(z))
     d1_index
     #  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14
     
     # 3.0.b build a data.frame (a table)
     
     d1_all <- data.frame(d1_index, mean.ind.biomass, NBSS, z )
     
     # 3.1 find the maximum 
     
     d_1MAXindex <- which.max(z) # position of the maximum (index)
     
     maxim_Xvalue.meanbiom <- mean.ind.biomass[d_1MAXindex] # position of the NBSS maximum (X value, in mean.ind.biomass units)
      
    maxim_Yvalue.log10NBSS <-  z[d_1MAXindex] # position of the NBSS maximum (y value, in log10 norm. biomass units)
       
       
     #d_1MAXindex # 3 (the maximum is located at the position 3)
     
     # 3.2 select only the data from the maximum onwards
     
     d1_all.2 <- d1_all[ d1_all$d1_index[d1_all$d1_index > (d_1MAXindex-1) ] ,]
     
     # 3.3 find the last non-empty bin before empty bins occur
     
     NAindex <- which(is.na(d1_all.2$z)) # locate ALL NA values in vector z
     firstNA <- min(NAindex)# locate the first NA value in vector z 
     firstNA # position of first NA value in vector z
     LastNonNA <-firstNA - 1 # position of last non_NA value in vector z
     LastNonNA  # position last non-empty bin before empty bins occur, index as for d1_all.2
   
     
     
   
     Last_Xvalue_meanbiom <- d1_all.2$mean.ind.biomass[LastNonNA] # position of the last used X value (in mean.ind.biomass units)
     
     Last_Yvalue_log10NBSS <-  d1_all.2$z[LastNonNA] # position of the last used Y value (in log10 norm. biomass units)
     
         
     # 3.4. select the data (based on automatically selected indices)
     
     d1_all.3 <- d1_all.2[ 1:LastNonNA ,]
     
     # position last non-empty bin before empty bins occur, index as for d1_all.2
     #ALL.LastNonNA <- max(d1_all.3$d1_index) # maximum index 
     # position last non-empty bin before empty bins occur, index as for d1_all.2
     
     
     lm.sel <- lm(log10(d1_all.3$NBSS) ~ log10(d1_all.3$mean.ind.biomass))
     summary(lm.sel)
     coef(lm.sel)
     
     rlm.sel <- MASS::rlm(log10(d1_all.3$NBSS) ~ log10(d1_all.3$mean.ind.biomass))
     summary(rlm.sel)
     coef(rlm.sel)
     
     
     
     # abline(lm.sel)
     
     
     NBSS.Intercept<- as.numeric(coef(lm.sel)[1]) # NBSS Intercept
     NBSS.slope <- as.numeric(coef(lm.sel)[2] )# NBSS Slope
     
     rrNBSS.Intercept<- as.numeric(coef(rlm.sel)[1]) # NBSS Intercept, Robust regression
     rrNBSS.slope <- as.numeric(coef(rlm.sel)[2] )# NBSS Slope, Robust regression
     
     
     
     sum.table <- as.list (summary(lm.sel))
     
     p.values <- summary(lm.sel)$coefficients[,4]  # p-values 
     r.squared <- summary(lm.sel)$r.squared # R-squared
     
     #N.bins.used <- length(d_1MAXindex:ALL.LastNonNA)
     
     
     
     # res <- c(NBSS.slope, NBSS.Intercept, p.values, r.squared, d_1MAXindex, 
     #          ALL.LastNonNA, N.bins.used)
     
     
     
     
     results <- list(slope_OLSR = NBSS.slope, intercept_OLSR = NBSS.Intercept, 
                     slope_RR = rrNBSS.slope,
                     intercept_RR = rrNBSS.Intercept,
                     p.values = p.values, 
                     r.squared = r.squared, index.of.maximum = d_1MAXindex, 
                     x.value.of.maximum = d1_all$mean.ind.biomass[d_1MAXindex], #X value of maximum, in mean.ind.biomass units
                     y.value.of.maximum =  d1_all$NBSS[d_1MAXindex],#y value of maximum, in log10NBSS units
                     index.of.start.bin <- d_1MAXindex, 
                     x.value.of.start.bin = d1_all$mean.ind.biomass[d_1MAXindex],
                     Last_Xvalue_meanbiom = d1_all.2$mean.ind.biomass[LastNonNA], # last used X value (in mean.ind.biomass units)
                     Last_Yvalue_log10NBSS =  d1_all.2$z[LastNonNA] ,# last used Y value (in log10 norm. biomass units)
                     n.NBSS.withdata = n.NBSS.withdata
                      #endpoint = ALL.LastNonNA, 
                     # x.value.of.end.bin = d1_all$mean.ind.biomass[ALL.LastNonNA],
                     #N.bins.used = N.bins.used)
     )
     
     names(results)
     
     res <- list(results = results, sum.table = sum.table, selected.data = d1_all.3 , all.data = d1_all)
     
     
     
     return(res)
     
     
   }
   
  
    
   # create matrix to fit slopes ---------
   
   # zoopl.NBSSmatrix ---------
   
   #View(zoopl.All.Atl[35:97])
   
   dim(zoopl.All.Atl)
   # 1580 samples
   # write.csv(zoopl.All.Atl, file = "zoopl_All_Atl_n1580.csv") 
   
   dim(zoopl.All.Atl)
   
   names(zoopl.All.Atl[37:92])
   length(names(zoopl.All.Atl[37:92])) # 56 NBSS data columns
   
   zoopl.NBSSmatrix <- as.matrix(zoopl.All.Atl[,37:91])
   
   
   zoopl.NBSSmatrix_no_na <-  zoopl.NBSSmatrix[rowSums(is.na(zoopl.NBSSmatrix)) != ncol(zoopl.NBSSmatrix), ]
   
   
   bin.means <-  exp( X_vector_gCind)
   
   X_vector_gCind
   
   # first plots
   
   plot(zoopl.NBSSmatrix_no_na [22, ]  ~   logCvector)
   
   length(zoopl.NBSSmatrix_no_na [22, ])  
   length(logCvector)  
   
    
   lnumb = 30
   #lnumb = sample( 1: nrow(zoopl.NBSSmatrix_no_na) , 1)
   
   plot(log10(zoopl.NBSSmatrix_no_na [lnumb, ])  ~   logCvector,
        main = paste (lnumb))
   
   abline (lm (  log10(zoopl.NBSSmatrix_no_na [lnumb, ])  ~   logCvector ) )
   
   
   # apply function NBSS.select.w.lm.MS.A-------------------------
   # test
   results <- NBSS.select.w.lm.MS.A ( 10^logCvector , 
                                      zoopl.NBSSmatrix_no_na[lnumb,] )
   
   plot(log10(zoopl.NBSSmatrix_no_na [lnumb , ])  ~    (logCvector),
        main = paste (lnumb)    )
   
   b =  results$results$slope
   #[1]  -0.6459868
   
   a = results$results$intercept
   #[1] -0.7477856
   
   nbins <-results$results$N.bins.used
   nbins
   
   abline( a= a , b= b, col = "red")
   
   
   #   raw.NBSSi <-  zoopl.NBSSmatrix[4,]
   #   raw.NBSSi[is.na(raw.NBSSi)] <- 0
   #   
   #   NBSS.select.w.lm.MS.A(as.vector(bin.means), as.vector(raw.NBSSi)) 
   
   
   # apply function NBSS.select.w.lm.MS.B.RR (with Rob. Regr. and OLSR-------------------------
   # test
   library(MASS)
   
   lnumb = 78
   # lnumb = sample( 1: nrow(zoopl.NBSSmatrix_no_na) , 1)
   
   results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                         zoopl.NBSSmatrix_no_na[lnumb,] )
 
   
   length (which(!is.na(  ( zoopl.NBSSmatrix_no_na[lnumb,]))))
   
   results$results$n.NBSS.withdata
     

      plot(log10(zoopl.NBSSmatrix_no_na [lnumb , ])  ~    (logCvector),
        main = paste ("line =", lnumb , ", ", "Gear = ", zoopl.All.Atl$Gear[lnumb] ,  
                      ", ", zoopl.All.Atl$CruiseID[lnumb])  ,
        xlim = c(-10, 7), ylim = c(-2, 4))
   
   b =  results$results$slope_OLSR
   #[1]  -0.6459868
   
   a = results$results$intercept_OLSR
   #[1] -0.7477856
   
   abline( a= a , b= b, col = "red")
   
   brr =  results$results$slope_RR
   #[1]  -0.6459868
   
   arr = results$results$intercept_OLSR
   #[1] -0.7477856
   
   abline( a = arr , b= brr, col = "darkgreen")
   
   
   p.val <- as.numeric( results$results$p.values[2])
   p.val
   
   r.sq <- as.numeric( results$results$r.squared)
   r.sq
   
   olsr_rr_slope_ratio <-  b / brr
   olsr_rr_slope_ratio
   
   text(x = -6, y =  3 , paste ( "OLSR slope =", round(b, 3))  )
   text(x = -6, y =  2.5 , paste ( "Rob. Regr. slope =", round(brr, 3))  )
   text(x = -6, y =  2 , paste ( "R² =", round(r.sq, 8))  )
   text(x = -6, y =  1.5 , paste ( " p =", round(p.val , 3))  )
   text(x = -6, y =  1 , paste ( "slope OLSR/RR =", round(olsr_rr_slope_ratio, 3))  )
   
   
   
   
   
   
   ####
   #  Four Tasks (Zooplankton NBSS) ----------------------------------------
   #  Task 1 --------------------
   # 1a Recalculate  all zoopl slopes  with point selection and Robust Regression 
   # 1b Apply a function with robust regression, output includes Rsquared, "p" values, and rr/olsr ratio  ...
   # 1c Create an output table, by applying the function NBSS.select.w.lm.MS.B 
   # 1d  Calcute  rr/olsr slope ratio  ...
   # 1e Filter for high-quality models only
   # 1f  Describe (median, etc.) high-quality slolpes for UVP zoopl. and UVP zoopl 
   
   #  Task 2 --------------------
   # plots (all Atlantic), carbon and  biovolume
   # compare carbon vs   biovolume slopes
   # plots (all Atlantic), NBSS in carbon and in biovolume, 
   # maps and NBSS plots for  net-caught  zoopl. vs UVP zoopl. 
   # compare carbon vs   biovolume slopes
   # compare  net-caught  zoopl. vs  UVP zoopl. 
   
   
   #  Task 3 --------------------
   # analyse driving factors for Zoopl. NBSS
   # net-caught  zoopl., UVP zoopl. and UVP zoopl.+detr.
   
   #  Task 4 --------------------
   # Zoopl. NBSS by 4 regions
   
   
   
   
   #####
   ###
   # ' # TASK 1 (Zooplankton)
   #  Task 1 --------------------
   
   
   #  Task 1 --------------------
   # 1a Recalculate  all zoopl slopes  with point selection and Robust Regression 
   # 1b Apply a function with robust regression, output includes Rsquared, "p" values, and rr/olsr ratio  ...
   # 1c Create an output table, by applying the function NBSS.select.w.lm.MS.B 
   # 1d  Calcute  rr/olsr slope ratio  ...
   
   
   # 1a Recalculate  all zoopl. slopes  with point selection and Robust Regression ---------------
   # 1b Apply a function with robust regression, output includes Rsquared, "p" values, and rr/olsr ratio  ...-----------
   # 1c Create an output table, by applying the function NBSS.select.w.lm.MS.B ------------------
   # 1d  Calcute  rr/olsr slope ratio  ...-----------------
   
   
   #  Compare OLSR slopes with robust regression and check "p" values --------------
   # apply function "NBSS.select.w.lm.MS.B.RR()"  ----------------
   
   # first quality selection (p and R-squared)   -------------------- 
   
   # first step> calculate p vales and R-squared (loop 1) ----
   
   out.slope.OLSR.vec <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )
   out.interc.OLSR.vec <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )
   out.p.vec <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )
   out.r.sq.vec <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )
   out.N.bins.used <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )
  
   out.x.value.of.maximum <- c(rep (NA, length(zoopl.All.Atl$Longitude)) ) #X value of maximum, in mean.ind.biomass units
   out.y.value.of.maximum<- c(rep (NA, length(zoopl.All.Atl$Longitude)) )#y value of maximum, in log10NBSS units
   out.Last_Xvalue_meanbiom <- c(rep (NA, length(zoopl.All.Atl$Longitude)) ) # last used X value (in mean.ind.biomass units)
   out.Last_Yvalue_log10NBSS <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )# last used Y value (in log10 norm. biomass units)
   out.n.NBSS.withdata <- c(rep (NA, length(zoopl.All.Atl$Longitude)) )# last used Y value (in log10 norm. biomass units)
   
 ##  lnumb = 78
   lnumb = 97
   
   
    length(zoopl.All.Atl$slope.clean)
   
   for (i in 1: length(zoopl.All.Atl$slope.clean)){
     
     tryCatch({
       
       lnumb <- i
       
       results <- NBSS.select.w.lm.MS.A ( 10^logCvector , 
                                             zoopl.NBSSmatrix[lnumb,] )
       
       
         
       b <-  results$results$slope
       
       a <- results$results$intercept
       
       p.val <- as.numeric( results$results$p.values[2])
       
       r.sq <- as.numeric( results$results$r.squared)
       
       N.bins.used <- as.numeric( results$results$N.bins.used)
       
       
       out.slope.OLSR.vec[lnumb] <- b
       out.interc.OLSR.vec[lnumb] <- a 
       out.p.vec[lnumb] <- p.val
       out.r.sq.vec[lnumb] <- r.sq
       out.N.bins.used[lnumb] <- N.bins.used
       
       out.x.value.of.maximum[lnumb] <- results$results$x.value.of.maximum    #X value of maximum, in mean.ind.biomass units
       # out.y.value.of.maximum[lnumb] <- results$results$y.value.of.maximum #y value of maximum, in log10NBSS units
       #out.Last_Xvalue_meanbiom[lnumb] <- results$results$Last_Xvalue_meanbiom  # last used X value (in mean.ind.biomass units)
       # out.Last_Yvalue_log10NBSS[lnumb] <- results$results$Last_Yvalue_log10NBSS   # last used Y value (in log10 norm. biomass units)
       out.n.NBSS.withdata[lnumb] <- results$results$n.NBSS.withdata  # N of bins in NBSS (not NA)
       
       
     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
     
   }
   
   out.slope.OLSR.vec 
   out.n.NBSS.withdata
   
   summary(out.slope.OLSR.vec) #  1358  NA's
   length(out.slope.OLSR.vec) # 1579 data with NA's
   1580 - 1358 # 222 OLSR slope data without NA's
   
   hist(out.slope.OLSR.vec)  
   hist(out.interc.OLSR.vec)
   hist (out.p.vec )
   hist (out.r.sq.vec)
   hist (out.N.bins.used)
   
   hist (log10(out.x.value.of.maximum ))  #X value of maximum, in mean.ind.biomass units
   # hist (out.y.value.of.maximum)  #y value of maximum, in log10NBSS units
   # hist (out.Last_Xvalue_meanbiom)   # last used X value (in mean.ind.biomass units)
   # hist (out.Last_Yvalue_log10NBSS) # last used Y value (in log10 norm. biomass units)
   # 
   
     
   
   zoopl.All.Atl$slope_OSLR_max.to.NA <- out.slope.OLSR.vec
   
   zoopl.All.Atl$intercept_OSLR_max.to.NA <- out.interc.OLSR.vec
   
   zoopl.All.Atl$p_value_olsr <- out.p.vec
   
   zoopl.All.Atl$r.square.olsr <- out.r.sq.vec
  
   zoopl.All.Atl$N.bins.used <- out.N.bins.used
   
    
      zoopl.All.Atl$x.value.of.maximum  <-out.x.value.of.maximum
      # zoopl.All.Atl$y.value.of.maximum <-out.y.value.of.maximum
      # zoopl.All.Atl$Last_Xvalue_meanbiom <-out.Last_Xvalue_meanbiom
      # zoopl.All.Atl$Last_Yvalue_log10NBSS <-out.Last_Yvalue_log10NBSS
      # 
   
      zoopl.All.Atl$n.NBSS.withdata  <- out.n.NBSS.withdata
      
      
      # View(zoopl.All.Atl)
      length  (   zoopl.All.Atl$slope_OSLR_max.to.NA  [zoopl.All.Atl$Gear == "Bongo" ])
     # 35  ABRACOS Bongo samples, 21 NA, total = 14 samples with carbon NBSS
      zoopl.All.Atl$CruiseID  [zoopl.All.Atl$Gear == "Bongo" ]
      zoopl.All.Atl$StationID  [zoopl.All.Atl$Gear == "Bongo" ]
      zoopl.All.Atl$slope_OSLR_max.to.NA  [zoopl.All.Atl$Gear == "Bongo" ]
      zoopl.All.Atl$N.bins.used  [zoopl.All.Atl$Gear == "Bongo" ]
      zoopl.All.Atl$p_value_olsr  [zoopl.All.Atl$Gear == "Bongo" ] 
      # all p  values are OK!
      # total = 14 ABRACOS samples with carbon NBSS, all OLSR slopes are OK
      
       summary(   zoopl.All.Atl$slope_OSLR_max.to.NA  [zoopl.All.Atl$Gear == "Bongo" ])
   # [1] -1.0544088         NA -0.9131190         NA -0.9591326         NA -0.8167237         NA -0.7489733
   # [10]         NA -1.1517037         NA -0.8220523         NA -1.1749962 -1.3245572         NA         NA
   # [19]         NA -1.6570372         NA         NA         NA         NA         NA         NA -1.2366240
   # [28]         NA -1.0022773         NA -1.2859452         NA         NA         NA -1.8491355  
  # 14 slope data from ABRACOS    
 
   
     dim (zoopl.All.Atl) # 1579  datasets, 107 variables
   
     dim(zoopl.NBSSmatrix_no_na)
      
      # N = 904 datasets, zoopl. (UVP+NET, not fitltered yet for "p" , Rsq, and linearity)   
     
      
      # NET zoopl -------------------
      zoopl_All_Atl4.unfiltr_NET <- subset (zoopl.All.Atl, Gear != "UVP")
      dim( zoopl_All_Atl4.unfiltr_NET)
      # Carbon: 838 NET  unfiltered samples, many are NA
      # # 169 zoopl. NET samples with slope data, 669 with NA's only
      # Carbon: 134 Net-caught high-quality linear filtered samples, OK
      #Biovolume: 111 Net  high-quality linear filtered samples, OK
      
     summary(as.factor( zoopl_All_Atl4.unfiltr_NET$CruiseID))
      
      summary( zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA) # 669 NA's
      +169+669 # 169 zoopl. NET samples with slope data, 669 with NA's only
      
      hist(zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA)
       
     #    write.csv(zoopl_All_Atl4.unfiltr_NET, file = "zoopl_All_Atl4.unfiltr_NETv8b.csv", sep =";")
      
      names( zoopl_All_Atl4.unfiltr_NET[37:91])
       
      zoopl_All_Atl4.unfiltr_NET.MATRIX <- as.matrix(zoopl_All_Atl4.unfiltr_NET[37:91])
      
      zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na <-  zoopl_All_Atl4.unfiltr_NET.MATRIX[rowSums(is.na(zoopl_All_Atl4.unfiltr_NET.MATRIX)) != ncol(zoopl_All_Atl4.unfiltr_NET.MATRIX), ]
      rowSums(is.na(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na))
      rowSums(!is.na(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na))
      
      
      dim( zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na)  # ne zoopl N = 169 Net  unfiltered samples, no NA's OK
      # Carbon: 169 Net  unfiltered samples, no NA's OK
      # Carbon: 134 Net  high-quality linear filtered samples, OK
      # Biovolume: 111 Net  high-quality linear filtered samples, OK
    
      #phytopl.All.Atl.MATRIX.unfiltr.nolims727<-  phytopl.All.Atl.unfiltr.nolims727
      phytopl.All.Atl.MATRIX.unfiltr.nolims727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727[,38:92])
   
         dim(phytopl.All.Atl.MATRIX.unfiltr.nolims727) # 727
      
       zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na727 <-  phytopl.All.Atl.MATRIX.unfiltr.nolims727 [rowSums(is.na(phytopl.All.Atl.MATRIX.unfiltr.nolims727)) != ncol(phytopl.All.Atl.MATRIX.unfiltr.nolims727), ]
       dim(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na727)
       
      
    
       
    # UVP---------  
      
      zoopl_All_Atl.unfiltr_UVP <- subset (zoopl.All.Atl, Gear == "UVP")
      dim( zoopl_All_Atl.unfiltr_UVP) 
    
      rowSums(!is.na(zoopl_All_Atl.unfiltr_UVP)) # OK
      
      
      # Carbon: TOTAL = 741  unfiltered UVP samples (linear and nonlinear),OK
      # Carbon: 651 zoopl. UVP samples with OLSR slope data, 90 with NA's only
      # Carbon: 271  UVP  high-quality linear filtered samples,OK
      # Biovolume: 484 UVP high-quality linear filtered samples
      
      
      hist( zoopl_All_Atl.unfiltr_UVP$n.NBSS.withdata)
      summary( zoopl_All_Atl.unfiltr_UVP$n.NBSS.withdata) # 90 NA's in Slope OLSR from max
      741-90 
      # 651 UVP profiles could be fitted with a NBSS slope (-4.9 to -2)
      # (all kinds of models.. low "n", bad "p", etc)
      
      
        hist( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
            summary( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA) # 90 NA's in Slope OLSR from max
      741-90 
      # 651 UVP profiles could be fitted with a NBSS slope (all kinds of models.. low "n", bad "p", etc)
      
      
      # summary( zoopl_All_Atl.unfiltr_UVP$slope.clean) # 46 NA's in OLD slope (not from maximum )
      # 741-46  
       
     hist(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
      
     # t.test(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA, 
     #        zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA) 
   
     zoopl_All_Atl4NET <-  read.csv("outp_zoopl_All_Atl4f_Net_CARB.csv", row.names = NULL)
     dim( zoopl_All_Atl4NET)  
     # 111  Net  samples,OK
     
      
     # COMPARE PHYTOPL vs UVP ZOOPL  slopes
      wilcox.test(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes, 
                  zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA) 
     #  p-value = 3.497e-12 (NEW, nano-to-microphyt.  vs UVP zoopl  )
      #   311  nano-to-microphyt.  vs 157  UVP zoopl  )
      
     
      
      # COMPARE PHYTOPL vs NET ZOOPL slopes
      wilcox.test (zoopl_All_Atl4NET$rob_reg__slope,
                           phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
      # p-value = 1.51e-08 (NEW, nano-to-microphyt.  vs NET zoopl  )
      #   311  nano-to-microphyt.  vs 117  NET zoopl  )
      1.51e-08
      0.000000015
       
      # COMPARE UVP  vs NET ZOOPL slopes (unfiltr. ALLL, OLSR)
      wilcox.test (zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA,
                   zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
     # p-value = 9.772e-05, OK
          
      # COMPARE UVP  vs NET ZOOPL slopes (filtr. top models, RR)
     #  wilcox.test (zoopl_All_Atl4NET$rob_reg__slope,
     #               zoopl_All_AtlUVP $rob_reg__slope)
     # #  p-value = 0.6836, n.s.
      
      
      
      options(scipen = 999) 
      9.772e-05 #  9.772e-05 = 0.00009772 
      options(scipen = 0)
     # The OLSR slope of UVP (OLSR, unfiltr.) was significant flatter than the 
      # slope of the NET Zoopl.
      # (Mann-Whitney U test,p < 0.0001)
      
      summary(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
      summary(zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA)
      
    
    #  UVP slopes are flatter than Zoopl. net samples, p < p-value = 5.66e-10
     # N UVP sampes = 651 zoopl. UVP samples with OLSR slope (max to non-empty)
     # N NET-caught samples = 169 NET-caught samples with OLSR slope (max to non-empty)
     
     summary( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA) # 90 NA's
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -4.1642 -1.0572 -0.7640 -0.8470 -0.5381  0.1482      90 
    
    length(( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA [zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA < -1]))
    length(( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA [zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA > -1]))
   
     
     summary( zoopl_All_Atl.unfiltr_UVP$intercept_OSLR_max.to.NA) # 90 NA's
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -15.2226  -2.3674  -1.3224  -1.4349  -0.2563   2.9971        6    
     
     
     summary( zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA) # 669 NA's
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -3.0297 -1.0881 -0.9156 -0.9845 -0.8159 -0.2018     669  

     # MULTIPLE COMPARISONS  ------------
     # K-W ANOVA-----------
     # Permanova--------------
     
            # K-W ANOVA
      # Permanova
    
      
      boxplot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes,
               zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA,
              zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA) 
      
      abline (h = -1, col = "green")    
      abline (h = -0.9)    
          abline (h = -0.8)    
  abline (h = -0.7)  
  abline (h = -0.6)  
  
  # Creating variable
  PP= phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes
  NETZ= zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA
  UVPZ= zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA
  
  # Finding maximum length
  max_ln1 <- max(c(length(PP), length(NETZ)))
  max_ln2 <- max(c(length(UVPZ)))
  max_ln<-max(max_ln2,max_ln1)
  gfg_data<- data.frame(PP = c(PP,rep(NA, max_ln - length(PP))),
                        NETZ = c(NETZ,rep(NA, max_ln - length(NETZ))),
                        UVPZ = c(UVPZ,rep(NA, max_ln - length(UVPZ))))
  
  gfg_data
  is.data.frame((gfg_data))
   summary(gfg_data)    
  
 boxplot(gfg_data)
 
 library(PMCMRplus)
   
  kruskalTest(gfg_data) # p-value < 2.2e-16
  kruskal.test(gfg_data) # p-value < 2.2e-16
  
  PMCMRplus::kwAllPairsNemenyiTest(gfg_data)
  
  #           1       2   
  # 2      3.9e-08   -   
  # 3       1.1e-14   0.15
  # 
    
     boxplot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes,
             zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA,
             zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA) 
     
     abline (h = -1, col = "green")    
     abline (h = -0.9)    
     abline (h = -0.8)    
     abline (h = -0.7)  
     abline (h = -0.6)  
     
     
     
     
     
     boxplot (zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA, 
              zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA) 
     
     
      zoopl_All_Atl4.unfiltr_UVP.MATRIX <- as.matrix(zoopl_All_Atl.unfiltr_UVP[37:91])
      
      zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na <-  zoopl_All_Atl4.unfiltr_UVP.MATRIX[rowSums(is.na(zoopl_All_Atl4.unfiltr_UVP.MATRIX)) != ncol(zoopl_All_Atl4.unfiltr_UVP.MATRIX), ]
      rowSums(is.na(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na))
      
      
      # Categories of UVP Data - what made so many  be filtered out?
      dim( zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na)  
      # Carbon: 735  No-NA, unfiltered UVP samples (linear and nonlinear),OK
        # OLD Carbon: 271  UVP  high-quality linear filtered samples,OK
      # OLD Biovolume: 484 UVP high-quality linear filtered samples
      
      # among these, 46 had less than 4 bins with data

      # 202 high-quality
     summary( zoopl_All_Atl.unfiltr_UVP)
   
      write.csv(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na ,  file = "zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_naV49j.csv")
          
      # HISTOGRAM of maxima, unflitered -----------
      
      dim(zoopl_All_Atl.unfiltr_UVP) # 741 unfiltered samples
      logposit.Xof_maximaUVP_unfiltr <-  log10(zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum)
      zoopl_All_Atl.unfiltr_UVP$logposit.Xof_maximaUVP_unfiltr <-  log10(zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum)
      
     hist( logposit.Xof_maximaUVP_unfiltr,
           breaks = seq( -6 , -1 ,by = 0.1)  )
   abline(v = -4.9)
     
     #  # new variable ,  max_ini vs max_bump, cutoff is -5.5 logC ind.-1
     # zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini <- zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum
     # 
     # zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini <- zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum
     # 
     # zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini[zoopl_All_Atl.unfiltr_UVP$logposit.Xof_maximaUVP_unfiltr > -5.5] <- "Bump"
     # zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini[zoopl_All_Atl.unfiltr_UVP$logposit.Xof_maximaUVP_unfiltr < -5.5] <- "Max_ini"
     # zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini <- as.factor(zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini)
     # 
     # summary(zoopl_All_Atl.unfiltr_UVP$logposit.Xof_maximaUVP_unfiltr )
     #  summary(zoopl_All_Atl.unfiltr_UVP$bump_vs_max_ini)
    
     
      
      # unfiltered UVP samples------------------  
      # Bump: 394 samples with max in the middle of the spectrum at biomass larger -5.5. logC ind.-1 )
      # Dome: xxx max inthe middle , subsequent  linear decline = dome-shaped : xxx samples
      # Max_ini: 345 samples with max in the middle of the spectrum
      # Max_ini linear: xxx samples with max in the middle of the spectrum and liear shape
       #   Bump Max_ini    NA's 
    #    394     345       2 
       
      
     dim( zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA)
      dim( zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
      
       plot(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA~
              zoopl_All_Atl.unfiltr_UVP$SST_C_insitu)
      
     mod1 <-   lm(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA~
              zoopl_All_Atl.unfiltr_UVP$SST_C_insitu)
summary(mod1) #    p-value: 0.3052 , not signif. 
  #abline(mod1)     


plot(zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum~ log10(
       zoopl_All_Atl.unfiltr_UVP$SST_C_insitu))      
mod2 <-   lm(zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA~
               zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum)
summary(mod2) #    p-value: 0.3052 , not signif. 
#abline(mod1)     


      
       # FILTER, CLEANUP ---------------
   # Cleanup--------
   # Apply filter (select ) for "p", "R-squared" and "nbins used" ----------
   
   zoopl.All.Atl_p_rsq_OK <-  zoopl.All.Atl[   (zoopl.All.Atl$p_value_olsr < 0.05) &
                                               (zoopl.All.Atl$N.bins.used > 4),]
   
   dim (zoopl.All.Atl)
   # BIOVOLUME: 1579 zooplankton samples, total
   # CARBON: 1579 zooplankton samples, total
   
   dim (zoopl.All.Atl_p_rsq_OK) 
   # # CARBON:1467   zoopl. samples, R-sq > 0.5,  p < 0.05, nbins > 4
   
   # # BIOVOLUME: 1078 samples, R-sq > 0.5,  p < 0.05, nbins > 4
   # # BIOVOLUME: 1078 zopl. samples, filtered for R-sq > 0.5,  p < 0.05, nbins > 4
   
   # # CARBON:1194   zoopl. samples, R-sq > 0.5,  p < 0.05, nbins > 4
   # # CARBON:1194    zopl. samples, filtered for R-sq > 0.5,  p < 0.05, nbins > 4
   
   summary(zoopl.All.Atl_p_rsq_OK)
  
      # 1467 samples (NET+UVP) 
   #  1344     NA rows
   1467-1344 # 123 OK sampes with rows
   
   # intercept_OSLR_max.to.NA  p_value_olsr       r.square.olsr     N.bins.used    
   # Min.   :-2.753           Min.   :0.0000000   Min.   :0.5787   Min.   : 5.000  
   # 1st Qu.: 1.291           1st Qu.:0.0000837   1st Qu.:0.8621   1st Qu.: 6.000  
   # Median : 1.739           Median :0.0006373   Median :0.9189   Median : 7.000  
   # Mean   : 1.641           Mean   :0.0044103   Mean   :0.8994   Mean   : 7.256  
   # 3rd Qu.: 2.104           3rd Qu.:0.0044965   3rd Qu.:0.9557   3rd Qu.: 8.000  
   # Max.   : 3.291           Max.   :0.0483672   Max.   :0.9987   Max.   :15.000  
   # NA's   :2                NA's   :2           NA's   :2        NA's   :2 
   
   
   # as matrix (NBBS only), # 669 samples, R-sq > 0.5and p < 0.05,  al least 5 bins, OK
   
   
   zoopl.All.Atl_p_rsq_OK3 <- zoopl.All.Atl_p_rsq_OK
   
   dim(zoopl.All.Atl_p_rsq_OK3) # 1467 samples, ok
   # 
   zoopl.NBSSmatrix_no_na.3.ok <- as.matrix(zoopl.All.Atl_p_rsq_OK3[,37:92])
   
   
   # loop2 (with Robust regression, requires good data, with sufficient data points)
   
   out.slope.rr.vec <- c(rep (NA, length(zoopl.All.Atl_p_rsq_OK3$Longitude)) )
   out.interc.rr.vec <- c(rep (NA, length(zoopl.All.Atl_p_rsq_OK3$Longitude)) )
   out.slope.slope_ratio.vec <- c(rep (NA, length(zoopl.All.Atl_p_rsq_OK3$Longitude)) )
   
   
   for (i in 1: length(zoopl.All.Atl_p_rsq_OK3$slope.clean)){
     
     tryCatch({
       
       lnumb <- i
       
       results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                             zoopl.NBSSmatrix[lnumb,] )
       
       b <-  results$results$slope_OLSR
       
       brr <-  results$results$slope_RR
       
       arr <- results$results$intercept_RR
       
       olsr_rr_slope_ratio <-  b / brr
       
       
       out.slope.rr.vec[lnumb]  <- brr
       out.interc.rr.vec[lnumb] <- arr
       out.slope.slope_ratio.vec[lnumb] <- olsr_rr_slope_ratio
       
     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   }
   
   out.slope.rr.vec
   
   length (out.slope.rr.vec  )# 1467 zooplankton data 
   
   summary(out.slope.rr.vec) # 1246 NA's
   
   1467-1246
     length( out.slope.rr.vec[!is.na(out.slope.rr.vec)] )
   #carbon: 221 useful Rob. Reg.  data
   # carbon: 450 useful Rob. reg. slopes
   
   hist(out.slope.rr.vec)
   hist(out.interc.rr.vec)
   hist(out.slope.slope_ratio.vec)
   
   
   
   summary(out.interc.rr.vec)
   summary(out.interc.OLSR.vec)
   
   zoopl.All.Atl_p_rsq_OK3$rob_reg__slope <- out.slope.rr.vec
   zoopl.All.Atl_p_rsq_OK3$rob_reg_intercept <- out.interc.rr.vec
   zoopl.All.Atl_p_rsq_OK3$slope_ratio_OLSR_RR <- out.slope.slope_ratio.vec
   
   summary(zoopl.All.Atl_p_rsq_OK3)
   
   # Task 1e Filter for high-quality models only---------------------
   
   
   ### Quality indicators used for filter (nbins, Rsquared, "p" , RR/OLSR ratio) ----------
   
   # histograms before filter (zoopl, UVP+Net) --------
   hist (out.p.vec )
   hist (out.r.sq.vec)
   hist (out.N.bins.used)
   hist(out.slope.slope_ratio.vec)
   
   
   
   # summary before filter (zoopl. UVP+Net) --------
   summary (out.p.vec )
   summary (out.r.sq.vec)
   summary (out.N.bins.used)
   
   length(which (!is.na(zoopl.All.Atl_p_rsq_OK3$rob_reg__slope)))
   # 1579 zoopl. (UVP+Net) dataset before filter
   # 221 with rob_reg__slope
   
     length(which (!is.na(zoopl.All.Atl_p_rsq_OK3$slope_ratio_OLSR_RR)))
      # 1579 zoopl. (UVP+Net) dataset before filter
   # 221 with slope_ratio_OLSR_RR
   
  hist(zoopl.All.Atl_p_rsq_OK3$slope_ratio_OLSR_RR)
   
   dim(zoopl.All.Atl_p_rsq_OK3) # 597 samples
   # 1194 samples CARBON
   # 1326 samples BIOVOLUME, R-sq > 0.5and p < 0.05,  
   #      at least 5 bins, OK
   # includes, 744 NA's
   
      summary(zoopl.All.Atl_p_rsq_OK3$rob_reg__slope) #  744 NA's
      summary(zoopl.All.Atl_p_rsq_OK3$slope_ratio_OLSR_RR) #  744 NA's
     
       dim(zoopl.All.Atl_p_rsq_OK3)
      
      
      # delete 744 rows with NA's in "rob_reg__slope"
      
      #df[!is.na(df$B), ] # delete rows where there is NA in column "B" 
      
      zoopl.All.Atl_p_rsq_OK3b <-  zoopl.All.Atl_p_rsq_OK3[!is.na(zoopl.All.Atl_p_rsq_OK3$rob_reg__slope), ]
  
      dim(zoopl.All.Atl_p_rsq_OK3b)
      # 221 samples with rob_reg__slope, no NA's in slope
    
     
         
   # 
   # FILTER for slope_ratio_OLSR_RR ( 0.667 to 1.5, = <50% bias )
   
   zoopl.All.Atl_p_rsq_RR_OK4 <-  zoopl.All.Atl_p_rsq_OK3b[ (zoopl.All.Atl_p_rsq_OK3b$slope_ratio_OLSR_RR < 1.5) ,]
   zoopl.All.Atl_p_rsq_RR_OK4 <-  zoopl.All.Atl_p_rsq_RR_OK4[ (zoopl.All.Atl_p_rsq_RR_OK4$slope_ratio_OLSR_RR > 0.67) ,]
   
   dim(zoopl.All.Atl_p_rsq_RR_OK4)
   # 220 high-quality samples with rob_reg__slope, no NA's in slope
   # OLD 357 high-quality samples with rob_reg__slope, no NA's in slope
   # OlD 350 top-quality samples with rob_reg__slope, no NA's in slope, slope_ratio_OLSR_RR ( 0.77 to 1.3, = <30% bias )
   #  (17 samples were removed, with extreme slope_ratio_OLSR_RR (0.77 to 1.3, = <30% bias )  )
   
   
   dim (zoopl.All.Atl_p_rsq_RR_OK4) 
   
    
   hist (log10(zoopl.All.Atl_p_rsq_RR_OK4$x.value.of.maximum ))
   hist ((zoopl.All.Atl_p_rsq_RR_OK4$N.bins.used ))
   median_(zoopl.All.Atl_p_rsq_RR_OK4$N.bins.used) # median = 9 bins used, UVP & NET zoopl.
    
     
    #      R-sq > 0.5, p < 0.05,  
   #      at least 5 bins, 
   # slope ratio OLSR/RR   0.77<x< 1.3, OK, 4 smaples delted due to bad OLSR/RR ratio (above 1.3, strong outlier(s) effects)
   
   # write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "outp_Biovol_zoopl_All_Atl_p_rsq_slope_ratio_ok_n595.csv")
   
   # write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "outp_Biovol_zoopl_All_Atl4.csv")
   
   # write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "zoopl.All.Atl_p_rsq_RR_OK4b_n1322")
   
   # write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "zoopl.All.Atl_p_rsq_RR_OK4b_n1322.csv")
   
  #  write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "zoopl.All.Atl_p_rsq_RR_OK4c_n1183_CARBON.csv")
 
   #    write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "zoopl.All.Atl_p_rsq_RR_OK4g_n220_CARBON.csv")
   
       write.csv(zoopl.All.Atl_p_rsq_RR_OK4, file = "zoopl.All.Atl_p_rsq_RR_OK4h_n220_CARBON.csv")
       
   #zoopl_All_Atl4 <- zoopl.All.Atl_p_rsq_RR_OK4
   
    #zoopl_All_Atl4 <-  read.csv("zoopl.All.Atl_p_rsq_RR_OK4e_n433_CARBON.csv")

   
    getwd()
  
    
    # reanlyze the old, non-filred dataset
     
        zoopl_All_Atl4 <-  read.csv("zoopl.All.Atl_p_rsq_RR_OK4e_n433_CARBON.csv")
   
    
   
   length(zoopl_All_Atl4$Longitude) 
   # 350 non-NA, high-quality  filtered zooplankton NBSS data, OK, CARBON units 
   summary(zoopl_All_Atl4)
   
   # CARBON:
 #   target_organisms_UVP_NET  slope.clean      slope_OSLR_max.to.NA intercept_OSLR_max.to.NA
 #   Length:433               Min.   :-3.2862   Min.   :-2.1133      Min.   :-11.4336        
 #   Class :character         1st Qu.:-1.1856   1st Qu.:-0.9207      1st Qu.: -3.1900        
 #   Mode  :character         Median :-0.9347   Median :-0.7283      Median : -1.6915        
 #   Mean   :-0.9194   Mean   :-0.7468      Mean   : -2.0425        
 #   3rd Qu.:-0.6495   3rd Qu.:-0.5315      3rd Qu.: -0.8235        
 #   Max.   : 1.1509   Max.   :-0.1814      Max.   :  0.7160        
 #   NA's   :40        NA's   :28           NA's   :28              
 #  p_value_olsr      r.square.olsr     N.bins.used     rob_reg__slope    rob_reg_intercept     
 # Min.   :0.000000   Min.   :0.5026   Min.   : 5.000   Min.   :-4.1642   Min.   :-22.0150  
 # 1st Qu.:0.000036   1st Qu.:0.7388   1st Qu.: 7.000   1st Qu.:-0.9801   1st Qu.: -2.8203  
 # Median :0.000647   Median :0.8536   Median : 9.000   Median :-0.6239   Median : -1.1696  
 # Mean   :0.004517   Mean   :0.8213   Mean   : 9.286   Mean   :-0.7624   Mean   : -1.8786  
 # 3rd Qu.:0.004684   3rd Qu.:0.9258   3rd Qu.:11.000   3rd Qu.:-0.3845   3rd Qu.: -0.2669  
 # Max.   :0.048619   Max.   :0.9944   Max.   :18.000   Max.   : 0.1470   Max.   :  1.8706  
 # NA's   :28         NA's   :28       NA's   :28                                           
 #   slope_ratio_OLSR_RR
 #   Min.   :0.7954     
 #   1st Qu.:1.0000     
 #   Median :1.0000     
 #   Mean   :1.0087     
 #   3rd Qu.:1.0056     
 #   Max.   :1.2984   
   
   
   
   
   # BIOVOLUME:
   #   slope_OSLR_max.to.NA intercept_OSLR_max.to.NA  p_value_olsr    r.square.olsr     N.bins.used   
   #   Min.   :-2.5205      Min.   :-2.753           Min.   :0.0000   Min.   :0.5787   Min.   : 5.00  
   #   1st Qu.:-1.1206      1st Qu.: 1.133           1st Qu.:0.0000   1st Qu.:0.8616   1st Qu.: 6.00  
   #   Median :-0.8832      Median : 1.630           Median :0.0006   Median :0.9190   Median : 7.00  
   #   Mean   :-0.9580      Mean   : 1.473           Mean   :0.0042   Mean   :0.8988   Mean   : 7.54  
   #   3rd Qu.:-0.7262      3rd Qu.: 2.029           3rd Qu.:0.0043   3rd Qu.:0.9545   3rd Qu.: 9.00  
   #   Max.   :-0.3458      Max.   : 3.291           Max.   :0.0484   Max.   :0.9981   Max.   :15.00  
   #   NA's   :798          NA's   :798              NA's   :798      NA's   :798      NA's   :798    
   # rob_reg__slope    rob_reg_intercept     slope_ratio_OLSR_RR
   # Min.   :-4.0338   Min.   :-0.0042   Min.   :0.3368     
   # 1st Qu.:-1.1107   1st Qu.: 1.3057   1st Qu.:1.0000     
   # Median :-0.8085   Median : 1.6708   Median :1.0000     
   # Mean   :-0.9051   Mean   : 1.7123   Mean   :1.0007     
   # 3rd Qu.:-0.6155   3rd Qu.: 2.0934   3rd Qu.:1.0000     
   # Max.   :-0.0562   Max.   : 4.8950   Max.   :1.2689     
   # NA's   :796       NA's   :796       NA's   :796     
   
   
   # Task 1f  Describe (median, etc.) high-quality slopes for net zoopl., UVP zoopl and UVP zoopl+detr.-------------
   
hist(zoopl_All_Atl4$rob_reg__slope)   
   summary(zoopl_All_Atl4$rob_reg__slope)
   length(zoopl_All_Atl4$rob_reg__slope)
   # 350 non-NA, high-quality  filtered zooplankton NBSS data, OK, CARBON units 
   
   # summary(zoopl_All_Atl4$rob_reg__slope)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -4.1642 -0.9801 -0.6239 -0.7624 -0.3845  0.1470 
   
   # n = 350 high-qualty zoopl. data 
   
   
   # Task 1g  separate UVP ws Net-Caught samples
   
   
   zoopl_All_Atl4_UVP <- subset (zoopl_All_Atl4, Gear == "UVP")
   dim( zoopl_All_Atl4_UVP) 
   # Carbon: 202  UVP samples,OK
   # Biovolume: 484(?)  UVP samples
   
   zoopl_All_Atl4_Net <- subset (zoopl_All_Atl4, Gear != "UVP")
   dim( zoopl_All_Atl4_Net)  
   # Carbon: 117 Net  samples, OK
   #Biovolume: 111 Net  samples, OK
   
   
   
  # write.csv(zoopl_All_Atl4_UVP, file = "outp_zoopl_All_Atl4f_UVP_CARB.csv")
   
   zoopl_All_Atl4UVP <-  read.csv("outp_zoopl_All_Atl4f_UVP_CARB.csv", row.names = NULL)
   dim( zoopl_All_Atl4UVP)  
   # 271  UVP  samples,OK
   
   write.csv(zoopl_All_Atl4_Net, file = "outp_zoopl_All_Atl4f_Net_CARB.csv")
   
   zoopl_All_Atl4NET <-  read.csv("outp_zoopl_All_Atl4f_Net_CARB.csv", row.names = NULL)
   dim( zoopl_All_Atl4NET)  
   # 111  Net  samples,OK
   
   
   # Task 1h  plot zooplankton NBSS (All data) --------------------
   
   # for Paper  -----------
   # Figure FOR PAPER , ALL Atlantic Zooplankton NBSS -----------------
   # PLot ALL NBSS Data, Biovolume -------
   # with red median line, darkorange dashed line, lots of blue lines, darkgreen points -----------
   # BIOVOLUME units --------------
   
   
   # NET Zooplankton ------------- 
   # plot all points and rlm moldes (blue lines), and red median line  ----------
   
   # Net-caught zopl., 111 samples
   # zoopl_All_Atl4NET
   
   summary(zoopl_All_Atl4NET)
   length(zoopl_All_Atl4NET$Longitude) # n = 117, OK
   
   summary(zoopl_All_Atl4NET$rob_reg__slope)
   # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   # -4.16424 -1.08440 -0.84045 -0.86508 -0.46466 -0.02725 
   
   fun.summary_minmaxmedianmeanIQR(zoopl_All_Atl4NET$rob_reg__slope)
   
   # head   res
   # 1    min -4.16
   # 2    max -0.03
   # 3 median -0.84
   # 4   mean -0.87
   # 
   # [[2]]
   # [1] "IQR: -1.08 to -0.46"
   # 
   # [[3]]
   # [1] "range: -4.16 to -0.03"
    
   # COMPARE PHYTOPL vs ZOOPL slopes
   wilcox.test (zoopl_All_Atl4NET$rob_reg__slope,
                phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)
   
   
length( zoopl_All_Atl4NET$rob_reg__slope [zoopl_All_Atl4NET$rob_reg__slope > -1])
   # 80 with slope flaatter than -1
80/ 117   # 68% of the slopes are flatter than -1


   summary(zoopl_All_Atl4NET)
   length(zoopl_All_Atl4NET$Longitude) # n = 117, OK
   
   names (zoopl_All_Atl4NET[,39:93])
   
   zoopl.NBSSmatrix_4.NET <- as.matrix(zoopl_All_Atl4NET[,39:93])
   dim(zoopl.NBSSmatrix_4.NET)
   # 111 samples
   
   
   # plot a few randomly selected NET zoopl NBSS plots --------------
   
   lnumb = sample( 1: nrow(zoopl.NBSSmatrix_4.NET) , 1)
   
   plot(log10(zoopl.NBSSmatrix_4.NET [lnumb, ])  ~   logCvector,
        xlim = c(-7, 1), 
        ylim = c(-3, 5) ,
        main = paste (lnumb))
   
   results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                         zoopl.NBSSmatrix_4.NET[lnumb,] )
   
   b =  results$results$slope_RR
   #[1]  -0.6459868
   
   a = results$results$intercept_RR
   #[1] -0.7477856
   
   abline( a= a , b= b, col =alpha ("grey65", 0.3))
   
   
   
   
   
   
   lnumb <- 34
   plot (log10(zoopl.NBSSmatrix_4.NET[lnumb, ]) ~ logCvector, 
         main = "ATLANTIC, Net-caught Zooplankton",
         xlab = "log10( C Biomass (gC ind.-1)) ",
         ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
         xlim = c(-8.5, 1), 
         ylim = c(-3, 5) ,
         pch = 16, col = "white")
   for (lnumb in 1 : nrow(zoopl.NBSSmatrix_4.NET) ) 
   {    points (log10(zoopl.NBSSmatrix_4.NET[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("navy", 0.1))

     results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                           zoopl.NBSSmatrix_4.NET[lnumb,] )
   
     b =  results$results$slope_OLSR
     #[1]  -0.6459868
     
     a = results$results$intercept_OLSR
     #[1] -0.7477856
     
     abline( a= a , b= b, col =alpha ("grey65", 0.3))
     
     
          
     # abline( rlm(log10( zoopl.NBSSmatrix_4.NET[lnumb,]) ~ logCvector),
     #         col =alpha ("grey65", 0.3))
   }
   
   names(zoopl_All_Atl4NET)
   median(zoopl_All_Atl4NET$intercept_OSLR_max.to.NA)
   # medioan lm intercept :-3.464819  carbon
   median(zoopl_All_Atl4NET$rob_reg_intercept)
   
   
   # Net zoopl., n = 111
   # Biovolume  (462 useful Data)
   (median(zoopl_All_Atl4NET$rob_reg__slope) )
   # median RR slope : -0.914 (bivol)
   # median RR slope : -0.7288432 (CARBON)
   
   (min(zoopl_All_Atl4NET$rob_reg__slope)) 
   (max(zoopl_All_Atl4NET$rob_reg__slope) )
   # min to max  RR slope 
   # -2.1248 to -0.6658 
   
   summary(zoopl_All_Atl4NET$rob_reg__slope)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -2.1248 -1.0195 -0.9147 -0.9821 -0.8328 -0.6658 
  
   # 
   # abline( a = median(zoopl_All_Atl4NET$rob_reg_intercept) ,
   #         b =  median(zoopl_All_Atl4NET$rob_reg__slope)   ,
   #   lwd = 2.5, lty = 2, col = "darkorange" )
   # 
   
   
  # Analyze N bins and max position (X value of NBSS Maximum)  ----------  
  # Net zoopl (Carbon):
    hist (log10(zoopl_All_Atl4NET$x.value.of.maximum ))
   hist (log10(zoopl_All_Atl4NET$x.value.of.maximum ),
         breaks = seq(-8, -3, by =0.1) )
  # Two peaks, two maxima at -5.7  and -4.8 log10(gC ind.-1)
  # separate into NetZoopeakmesozoo and NetZoopeakmacrozoo?
   # diff. between regionsbetween NetZoopeakmesozoo and NetZoopeakmacrozoo ?
   # diff. in mesh size or gear  between NetZoopeakmesozoo and NetZoopeakmacrozoo?? 
       # diff. in NBSS slope between NetZoopeakmesozoo and NetZoopeakmacrozoo?? 
   
summary(log10(zoopl_All_Atl4NET$x.value.of.maximum ))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.924  -5.719  -5.419  -5.227  -4.815  -4.215 
      hist ((zoopl_All_Atl4NET$N.bins.used ))
   median_(zoopl_All_Atl4NET$N.bins.used) # median = 9 bins used,NET zoopl.
   max(zoopl_All_Atl4NET$N.bins.used)#  max 18 bins
   
    
   # matrix to to vector (transpose), TSWA -----------------------  
   x1 <- rep( logCvector, nrow(zoopl.NBSSmatrix_4.NET) )
   y1 <- c ( t(zoopl.NBSSmatrix_4.NET))
   
   # regression line (Rob. regr.)
   #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
   #points(  log10(y1) ~ x1 )
   # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
  
# ALL Data, no selection (bad)-----------   
  #  abline(rlm (log10(y1) ~ x1, , 
  #        lwd = 2.5, lty = 2, col = "pink" )

   
   # means and medians without considering NAs  
   mean (c(4,NA, 18,34)) # NA
   median( c(4,NA, 18, 34)  ) #NA
   
   mean_   <- function(...) mean(..., na.rm=T)
   median_ <- function(...) median(..., na.rm=T)
   
   mean_ (c(4,NA, 8, 34)) #15.33, OK
   median_ (c(4,NA, 8,34)) #8, OK
   
   # apply median by columns
   
   # median_vec_ALLBiovol <- apply( zoopl.NBSSmatrix_4.NET, 2, median_) 
   # mean_vec_ALLBiovol <- apply(zoopl.NBSSmatrix_4.NET, 2, mean_) 
   
    
  # lines(  log10(median_vec_ALLBiovol) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate median (or mean) only if less than 50% are NA ------------
   
   vec1 <- c(4, 5,NA,NA, NA,  34)
   
   sum(is.na(vec1))
   percNA <- sum(is.na(vec1))/ length(vec1)
   
   
   mean_NA50   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.5) {
       mean(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   median_NA50   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.5) {
       median(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   
  # median_vec_ALLBiovol_NA <- apply( zoopl.NBSSmatrix_4.NET, 2, median_NA50) 
  # mean_vec_ALLBiovol_NA <- apply(zoopl.NBSSmatrix_4.NET, 2, median_NA50) 
   
   median_vec_ALLNET.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na, 2, median_NA50) 
   
   
   
   #lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate rlm and lm from medians only (for comparison) -----------------
   # 
   # 
   # lm8_from_medians <- lm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # abline(lm8_from_medians, , 
   #        lwd = 2.5, lty = 2, col = "darkred")
   # 
   # # rlm from Net zoop medians-------
   # rlm7_from_medians <- rlm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # 
   # abline(rlm7_from_medians, 
   # lwd = 2.5, lty = 2, col = "darkgreen")
   # 
   # 
   # lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,
   #         pch = 16, col = "red" , lwd = 3.5)
    
   
   
   # medians (curved line) ------
   logmedians5 <- log10(median_vec_ALLNET.NA_unfiltr)
   length(logmedians5) # 55 data points with medians
   (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
   logmedians5_large <- logmedians5[30:55]
   logmedians5_small <- logmedians5[1:30]
   logCvector_large <- logCvector[30:55]
   logCvector_small <- logCvector[1:30]
   
   logCvector_largeNET <- logCvector_large
   logCvector_smallNET <- logCvector_small
   
   logmedians5_largeNET <- logmedians5_large
   logmedians5_smallNET <- logmedians5_small
   
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   # lm from Net zoop medians-------
   # lm8b_from_medians <- lm(  logmedians5_large ~ logCvector_large)
   # abline(lm8b_from_medians,  
   #        lwd = 2.5, lty = 2, col = "darkblue")
   # 
   
   # rlm from Net zoop medians-------
   rlm7b_from_medians <- rlm(  logmedians5_large ~ logCvector_large )
   abline(rlm7b_from_medians, 
          lwd = 2.5, lty = 2, col = "darkorange")
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   #Regression equation (dashed line)
    
   summary(rlm7b_from_medians)
   # Value    Std. Error t value 
   # (Intercept)       -3.6171   0.1245   -29.0457
   # logCvector_large  -0.9938   0.0353   -28.1711
   aovperm(rlm(logmedians5_large[1:8] ~ logCvector_large[1:8]
   )     )
    # p < 0.0001, OK  
   
   
   
   #####################
   # LINEAR MODELS -----------------
   ##################
   # NET Zooplankton ---------------
   # Linear models for ZOOplankton  -------------
   # Factors affecting the zooplankton NBSS slopes ------------------
   
   
   #### Temperature and chla vs zooplankton NBSS slope, plots and models ----------------


  # zoopl.All.Atl_p_rsq_RR_OK4
   summary( zoopl_All_Atl4NET$rob_reg__slope)
   summary( zoopl_All_Atl4NET$slope_ratio_OLSR_RR)
   dim (zoopl_All_Atl4NET$rob_reg__slope  )
   # 134 top-quality net zooplankton samples (carbon)
   
   
   summary( zoopl_All_Atl4NET$chla)# OK
   # (bad data, however N seems OK, 100%)
   summary( zoopl_All_Atl4NET$Chlorophyll__mg_m_3_insitu)# OK
   # (nice data in situ, however only 28 data...)
   summary( zoopl_All_Atl4NET$SST)# OK
   # (bad data, however N seems OK, 100%)
   plot(zoopl_All_Atl4NET$SST ~zoopl_All_Atl4NET$SST_C_insitu)
   plot(Data.All.Atl$SST ~Data.All.Atl$SST_C_insitu)
   # (bad data, however N seems OK, 100%)
   summary( zoopl_All_Atl4NET$primprod)# OK
   # (bad data, however N seems OK, 100%)
   
   #10.    univarate model with Latitude ----------
  
    lmLati <- lm(zoopl_All_Atl4NET$rob_reg__slope ~ 
                            (zoopl_All_Atl4NET$Latitude) )
    summary(lmLati) 
    # Multiple R-squared:  0.05699,	Adjusted R-squared:  0.04984 
    # F-statistic: 7.977 on 1 and 132 DF,  p-value: 0.005474
     # Carbon: Signifficant!!!p-value: 0.005474
   
      plot(zoopl_All_Atl4NET$rob_reg__slope ~ 
          (zoopl_All_Atl4NET$Latitude), 
        col = alpha ("darkgreen", 0.3), pch = 16)
   abline(lmLati, lwd = 2.5, lty = 2, col = "darkorange")
   #HORRIBLE PLOT!
   # ridiculous low R-squared, R-squared:  0.056
   # not useful!
   
   
   #1a.    univarate model with chl a ,log10(1+ Chl a in situ ) ----------
   lmlogChla_insitu <- lm(zoopl_All_Atl4NET$rob_reg__slope ~ 
                            log10(1+zoopl_All_Atl4NET$Chlorophyll__mg_m_3_insitu) )
   summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
   # CARBON (new): 26 DF,  p-value: 0.9049
   plot(zoopl_All_Atl4NET$rob_reg__slope ~ 
          log10(1+ zoopl_All_Atl4NET$Chlorophyll__mg_m_3_insitu), 
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
   # FEW data (only 28 data with chl a in situ!)
   
   
   
   
   #1b.    univarate model with chl a ,log10(1+ Chl a MODEL ) ----------
   lmlogChla_insitu <- lm(zoopl_All_Atl4NET$rob_reg__slope ~ 
                            (1+zoopl_All_Atl4NET$chla) )
   summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
   # CARBON (new): 26 DF,  p-value: 0.9049
   plot(zoopl_All_Atl4NET$rob_reg__slope ~ 
          log10(1+ zoopl_All_Atl4NET$chla), 
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
   # lots of  data (only 134 data with chl a model!)
   # but... awkward satellite/model data, not significant!
   # NOT SIGNIF!
   
   
   #1c.  Best univariate model, SST  model ----------
   
   summary( zoopl_All_Atl4NET$SST)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   #    5.50   15.60   22.34   21.26   27.30   34.62      13 
 
   dim(zoopl_All_Atl4NET) 
#   132 DF , 132 top-quelty filtered net zoopl. datasets
    
   
   
   lmSST_C_MODEL <- lm(zoopl_All_Atl4NET$rob_reg__slope ~ 
                          zoopl_All_Atl4NET$SST)
   summary(lmSST_C_MODEL) #OLD (biovolume , model SST data) R-squared:  0.349, 460 degrees of freedom, p-value: < 2.2e-16
   # #  100% data with SST from MODEL!)
   #  model data !!!
   # SIGNIF. , p = value: 0.007857
   
   rlmSST_C_MODEL <- rlm(zoopl_All_Atl4NET$rob_reg__slope ~ 
                         zoopl_All_Atl4NET$SST)
   summary(rlmSST_C_MODEL) #OLD (biovolume , model SST data) R-squared:  0.349, 460 degrees of freedom, p-value: < 2.2e-16
  aovperm(rlmSST_C_MODEL, nperm = 50000)

 # View(zoopl_All_Atl4NET)
     
    # #  100% data with SST from MODEL!)
   #  model data !!!
   # SIGNIF. , p = value: 0.008
   
   plot(zoopl_All_Atl4NET$rob_reg__slope ~ 
          zoopl_All_Atl4NET$SST, 
        col = alpha ("darkgreen", 0.3), pch = 16)
   abline(rlmSST_C_MODEL, lwd = 2.5, lty = 2, col = "darkorange")
   # SIGNIF. , p = value: 0.008
   # weird result! warmer = flatter shape! 
   # Model with negilible explantion power!
   # R-squared:  0.05!!
   # explains omnly 5% of all variabilty in zoopl. NBSS!
   
   
   # look at the intercept ---------

   summary( lm ( zoopl_All_Atl4NET$rob_reg_intercept ~ 
                   zoopl_All_Atl4NET$SST ))
      plot(zoopl_All_Atl4NET$rob_reg_intercept ~ 
          zoopl_All_Atl4NET$SST, 
        col = alpha ("darkgreen", 0.3), pch = 16)
   abline(rlmSST_C_MODEL, lwd = 2.5, lty = 2, col = "darkorange")
   # weird result! warmer = more biomass! 
   # Model with negligible explanatory power!
   # R-squared:  0.109!!
   # explains only 11% of all variability in zoopl. NBSS!
   
   plot(log10(zoopl_All_Atl4_Net$x.value.of.maximum ) ~ zoopl_All_Atl4_Net$SST,
        col = alpha ("navy", 0.3), pch = 16,
        xlab = "SST (°C, in situ)"  ,
        ylab = "ind.biomass at NBSS maximum (log10gC ind-1)"   )
   abline (v = c(20, 27.8), lwd = 2.5, lty = 2, col = "darkgrey")
   cor.test(log10(zoopl_All_Atl4_Net$x.value.of.maximum ) , zoopl_All_Atl4_Net$SST_C_insitu,  
            method = "spearman") #  not signif.
   length(zoopl_All_Atl4_Net$x.value.of.maximum) # 
   # NET zoopl: 134 samples (filtered for log-linear sections) 
   # UVP: (271 UVP profiles (filtered for log-linear sections)  )
   
   
   
   #MAP ---------------------
   
   ##################################
   ### MAPS ######################
   #MAP NET zoopl  data  --------
   # # with and without bumps ----------
   # # # NET zoopl --------
   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   
   points( zoopl_All_Atl4_UVP$Longitude,zoopl_All_Atl4_UVP$Latitude,   col="navy", pch=16)

   points( zoopl_All_Atl4_UVP$Longitude,zoopl_All_Atl4_UVP$Latitude,   col="magenta", pch=16)
   

   # #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   # #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
   # symbols(zoopl_All_Atl4NET$Latitude ~ zoopl_All_Atl4NET$Longitude,
   #         circles = (10^(-1* zoopl_All_Atl4NET$slope_OSLR_max.to.NA)), fg= alpha("navy", 0.6), inches=0.3,
   #         add = TRUE)
   # 
   
   
   
   # DEFINE 4 regions, V5
   # 4 key regions (where there is an exceptionally dense wealth of data)
   # version v4:
   # TSWA
   #new3
   #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
   # -2°S to -12°S
   # -38°W to -25°W
   
   # CCUS
   #old2  #  rect( ybottom = 16.9 ,  ytop = 30,           xleft = -26 , xright =   -18 )
   # 16.9°N to 30°N
   # -26°W to -18°W
   #new4:
   #  rect( ybottom = 12 ,  ytop = 33,     xleft = -27.5 , xright =   -9)
   
   
   # EQU
   #new3:
   #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
   # -1°S to 8°N
   # -33°W to -20°W
   
   
   # BUS
   # (old) #  rect( ybottom = -17.5,  ytop = -23.1,           xleft =  10.2, xright =   12.8 )
   
   #old 2:  #  rect( ybottom = -17.5,  ytop = -35,           xleft =  17, xright =   10.2)
   # -17.5°S to -28.5°S
   # 10.2°W to 17°W
   #new 3:
   #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
   
   
   
   
 #  [bookmark UVP]  ------------------
   
   
    #####################  
   ################
   # UVP
   #######
   # UVP zoopl --------
   
   # plot all points and rlm models (blue lines), and red median line  ----------
   
   # UVP Zoopl. , 484 samples
   # zoopl_All_Atl4NET
  
   
   dim(zoopl_All_Atl4UVP)
    # 202  UVP  samples,OK
   
   
   summary(zoopl_All_Atl4UVP)
   length(zoopl_All_Atl4UVP$Longitude) 
   # Median RR slope :-0.8591
   
   h <- hist(zoopl_All_Atl4UVP$rob_reg__slope)
   h$counts
   h$mids
   ( df  <- data.frame( c = h$counts , m   =  h$mids))
   +90+53 # 143 UVP samples with RR slope > -1
   143/202
   sum(df$c)
   summary(zoopl_All_Atl4UVP$rob_reg__slope)
   
   summary(zoopl_All_Atl4UVP$rob_reg__slope)
   
    
   
   length( zoopl_All_Atl4UVP$rob_reg__slope [zoopl_All_Atl4UVP$rob_reg__slope > -1])
   # 143 with slope flatter than -1
   143/ 202   # 71% of the slopes are flatter than -1
   
#   mean RR NBSS slope of -0.88 
   # (median:  -0.72, 
   # range: -3.34 to -0.07,
   # IQR: -1.13 to  -0.49, n = 202, ind. Biomass range analyzed = -4.9 to -2.2 log10 (gC ind-1)).
   
   fun.summary_minmaxmedianmeanIQR(zoopl_All_Atl4UVP$rob_reg__slope)
   
   
   summary(zoopl_All_Atl4UVP$slope_ratio_OLSR_RR)
   summary(zoopl_All_Atl4UVP$N.bins.used)
   
   
   
   
   
   zoopl.NBSSmatrix_4.UVP <- as.matrix(zoopl_All_Atl4UVP[,39:93])
   dim(zoopl.NBSSmatrix_4.UVP)
   # 202 sampes CARBON (271 OLD Carbon)
   # 111 samples biovolume
   
   
    
   
   
   # plot a few randomly selected UVP NBSS plots --------------
   
   lnumb = sample( 1: nrow(zoopl.NBSSmatrix_4.UVP) , 1)
   
   plot(log10(zoopl.NBSSmatrix_4.UVP [lnumb, ])  ~   logCvector,
        xlim = c(-8.5, 1), 
        ylim = c(-3, 5) ,
         main = paste (c("lnumb = ", lnumb, "Temp =  ",
                         round (zoopl_All_Atl4UVP$SST_C_insitu[lnumb], 1)  ))  )
   
   results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                         zoopl.NBSSmatrix_4.UVP[lnumb,] )
   
   b =  results$results$slope_OLSR
   
   a = results$results$intercept_OLSR
  
   b =  results$results$slope_OLSR
   #[1]  -0.6459868
   
   a.RR = results$results$intercept_RR
   b.RR = results$results$slope_RR
   
     
   abline( a= a , b= b, col =alpha ("darkorange", 0.3))
   abline( a = a.RR , b= b.RR, col = alpha ("darkgreen"))
   
      
   text( x = -3,y = 4.3, 
         paste("OLSR slope = ", round(results$results$slope_OLSR, 2))) 
   text( x = -3,y = 3.5, 
         paste("RR slope = ", round(results$results$slope_RR, 2))) 
   text(cex = 0.7, x = -5,y = -2, 
         paste("p = ", round(results$results$p.values[2],6) ) )
   
   results$results$p.values[2]
   
   
   # Plot all UVP Data ----------------
   
   lnumb <- 34
   plot (log10(zoopl.NBSSmatrix_4.UVP[lnumb, ]) ~ logCvector, 
         main = "ATLANTIC, UVP Zooplankton",
         xlab = "log10( C Biomass (gC ind.-1)) ",
         ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
         xlim = c(-8.5, 1), 
         ylim = c(-3, 5) ,
         pch = 16, col = "white")
   
   for (lnumb in 1 : nrow(zoopl.NBSSmatrix_4.UVP) ) 
   {    points (log10(zoopl.NBSSmatrix_4.UVP[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("navy", 0.3))
     
     results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                           zoopl.NBSSmatrix_4.UVP[lnumb,] )
     
     b =  results$results$slope_OLSR
     #[1]  -0.6459868
     
     a = results$results$intercept_OLSR
     #[1] -0.7477856
     
     abline( a= a , b= b, col =alpha ("grey65", 0.2))
     
   }
   
   # plot with depth or temperature as colour code ----------
   
   library(colourvalues)
   zoopl_All_Atl4_UVP$z_colours_fromTEMPinsitu   <-  colour_values(zoopl_All_Atl4_UVP$Temperature__C_insitu, palette = "viridis")
   zoopl_All_Atl4_UVP$z_colours_fromSSTinsitu   <-  colour_values(2+zoopl_All_Atl4_UVP$SST_C_insitu^8, palette = "viridis")
   
   plot(zoopl_All_Atl4_UVP$Temperature__C_insitu, col = alpha(  zoopl_All_Atl4_UVP$z_colours_fromTEMPinsitu , 0.7))
   
   plot(zoopl_All_Atl4_UVP$Catching_Depth_max_m, col = alpha(  zoopl_All_Atl4_UVP$z_colours_fromTEMPinsitu , 0.7))
   
   plot(zoopl_All_Atl4_UVP$Catching_Depth_min, col = alpha(  zoopl_All_Atl4_UVP$z_colours_fromTEMPinsitu , 0.7))
   
   # SUBSETS by depth?? -  select by depth strata, only samples with ??
   
   summary(zoopl_All_Atl4_UVP$Catching_Depth_max_m)
   hist(zoopl_All_Atl4_UVP$Catching_Depth_max_m)
   
   # Big question: 
   # deep sea samples! How can they be more biomass per m3, on average, than in  net samples????  
   
   
   lnumb <- 34
   plot (log10(zoopl.NBSSmatrix_4.UVP[lnumb, ]) ~ logCvector, 
         main = "ATLANTIC, UVP Zooplankton",
         xlab = "log10( C Biomass (gC ind.-1)) ",
         ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
         xlim = c(-8.5, 1), 
         ylim = c(-3, 5) ,
         pch = 16, col = "white")
   
   for (lnumb in 1 : nrow(zoopl.NBSSmatrix_4.UVP) ) 
   {    points (log10(zoopl.NBSSmatrix_4.UVP[lnumb,]) ~ logCvector ,
                pch = 16, 
            col = alpha(  zoopl_All_Atl4_UVP$z_colours_fromSSTinsitu , 0.04))
     
      }
   
   
   
   
   # UVP zoopl.,
   # Biovolume  
   (median(zoopl_All_Atl4UVP$rob_reg__slope) )
   (median(zoopl_All_Atl4UVP$rob_reg_intercept) )

   abline(a = median(zoopl_All_Atl4UVP$rob_reg_intercept), 
          b = median(zoopl_All_Atl4UVP$rob_reg__slope)  , 
   lwd = 2.5, lty = 2, col = "darkorange" )
   
    # median RR slope : -0.914 Net zoopl,  -0.859 UVP zoopl
   (min(zoopl_All_Atl4UVP$rob_reg__slope)) 
   (max(zoopl_All_Atl4UVP$rob_reg__slope) )
   
   summary(zoopl_All_Atl4UVP$rob_reg__slope)
   # UVP (carbon) - very flat !!!
      # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -4.1642 -0.9448 -0.5315 -0.7273 -0.3130  0.1470 
   
   # UVP (biovolume)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -2.5205 -1.1599 -0.8591 -0.9591 -0.6836 -0.3019 
   
   # NET
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -2.1248 -1.0195 -0.9147 -0.9821 -0.8328 -0.6658 
   
   ####################################
   # Analyze N bins and max position (X value of NBSS Maximum)  ----------  
   # UVP zoopl (Carbon):
   hist(log10(zoopl_All_Atl4UVP$x.value.of.maximum ))
  
   
    hist(log10(zoopl_All_Atl4UVP$x.value.of.maximum ),
         breaks = seq(-8, -2, by = 0.1 ))
   abline (v = -5.5, col = "red")
   abline (v = -4.9, col = "green")
   abline (v = -2, col = "green")
   
   hist(log10(zoopl_All_Atl4UVP$x.value.of.maximum ),
        breaks = seq(-8, -2, by = 0.5 ))
   abline (v = -5.5, col = "red")
  
   # Categories of UVP profiles ------------
   
dim(   zoopl_All_Atl.unfiltr_UVP) # 741 UVP samples
   summary(zoopl_All_Atl.unfiltr_UVP)
   
   hist(zoopl_All_Atl4_UVP$N.bins.used) # 5 to 10 bins used
   
   hist(zoopl_All_Atl.unfiltr_UVP$N.bins.used)
   
length(zoopl_All_Atl.unfiltr_UVP$N.bins.used [zoopl_All_Atl.unfiltr_UVP$N.bins.used > 5  ])   
# only 249 had 5 or more  useful bins (max to first non-empty bin)  


   # Cat A high-quality linear samples  n = 202
  dim( zoopl_All_Atl4UVP) # 202 high-quality linear UVP samples (-4.9 to -2)
  
  # As (A small-sized NBSS maximum)
  length(   log10(zoopl_All_Atl4UVP$x.value.of.maximum) [log10(zoopl_All_Atl4UVP$x.value.of.maximum) <  -4.4 ])
 # 90 small-sided maximum samples , all with SST equal o warmer 15
  # mean SST : 25.4 C
  
  length(   log10 (zoopl_All_Atl4UVP$x.value.of.maximum) [log10(zoopl_All_Atl4UVP$x.value.of.maximum) <  -4.4 ])
  # 90 small-sized maximum samples, <  -4.4, =  -4.4 to -2,   (tested range: -4.9 to -2)
  # among the 202 high-qualoiitylnear samples, all with SST > 15 degr 
  # (84 have data , 6 don't have SST data)
   summary( x <- (zoopl_All_Atl4UVP$SST_C_insitu [log10(zoopl_All_Atl4UVP$x.value.of.maximum) <  -4.4 ] ))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  # 14.99   25.33   25.93   25.37   26.86   28.59       6 
  
  # AL (A Large-sized NBSS maximum)
  
  length(   log10(zoopl_All_Atl4UVP$x.value.of.maximum) [log10(zoopl_All_Atl4UVP$x.value.of.maximum) >  -4.4 ])
  # 112 small-sided maximum samples, <  -4.4, =  -4.4 to -2,   (tested range: -4.9 to -2)
  # among the 202 high-qualoiitylnear samples, all at them > 15 degr
  # ALL have SST data , includes cold waters < 15 C
   summary( x <- (zoopl_All_Atl4UVP$SST_C_insitu [log10(zoopl_All_Atl4UVP$x.value.of.maximum) >  -4.4 ] ))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 13.03   20.85   25.78   23.42   26.90   28.61 


   plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$rob_reg__slope )
   plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$rob_reg_intercept )
   
   plot(log10(zoopl_All_Atl4_Net$x.value.of.maximum ) ~ zoopl_All_Atl4_Net$rob_reg__slope )
   plot(log10(zoopl_All_Atl4_Net$x.value.of.maximum ) ~ zoopl_All_Atl4_Net$rob_reg_intercept )
  
   plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$Temperature__C_insitu )
   
   plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu )
   abline (v = c(20, 27.8))
   
   summary(mod1 <- lm(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu ))
   abline (v = c(20, 27.8))
   abline (mod1)
   
   # 202 high-quallity UVP profiles (filtered for log-linear sections), 6 without SST data
   
   # FOR PAPER!!---------------
   # UVP PLOT!!! --------------
   # NBSS maximum vs STT ! -----------------
   
   plot((zoopl_All_Atl4UVP$rob_reg__slope ) ~ zoopl_All_Atl4UVP$SST_C_insitu)
   summary(lm(zoopl_All_Atl4UVP$rob_reg__slope  ~ zoopl_All_Atl4UVP$SST_C_insitu))
   # n.s.
  cor.test(zoopl_All_Atl4UVP$rob_reg__slope  , 
           zoopl_All_Atl4UVP$SST_C_insitu, method = "spearman")
  # p-value = 0.02388. cold is flatter!
  library(mgcv)
bt <- mgcv::gam( zoopl_All_Atl4UVP$rob_reg__slope  ~ s(zoopl_All_Atl4UVP$SST_C_insitu))  
plot(bt)  
anova(bt)
   
     plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu,
        col = alpha ("navy", 0.3), pch = 16,
        xlab = "SST (°C, in situ)"  ,
        ylab = "ind.biomass at NBSS maximum (log10gC ind-1)"   )
   abline (v = c(20, 27.8), lwd = 2.5, lty = 2, col = "darkgrey")
cor.test(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) , zoopl_All_Atl4UVP$SST_C_insitu,  
      method = "spearman") #  p-value = 0.009701
   length(zoopl_All_Atl4UVP$x.value.of.maximum) # 202 UVP profiles (filtered for log-linear sections), 6 without SST data
   zoopl_All_Atl4UVP$SST_C_insitu
   hist(zoopl_All_Atl4UVP$x.value.of.maximum)
   
  summary(( lm(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu)))    
   
   
   library(ggplot2)
   df1 <- data.frame( UVP_posMAX = log10(zoopl_All_Atl4UVP$x.value.of.maximum) , SST = zoopl_All_Atl4UVP$SST_C_insitu)
  
    ggplot(df1, aes(SST, UVP_posMAX ) ) +
      geom_point( alpha = 0.3, col = "navy", size = 2.5) +
      geom_smooth(span = 0.7)
    
      ggplot(df1, aes(SST, UVP_posMAX ) ) +
     theme(axis.text.y = element_text(face="bold", color="#993333", 
                                        size=10), 
           axis.text.x = element_text(face="bold", color="#993333", 
                                      size=10))     +
     scale_x_continuous(name ="SST (°C)", label = comma,
                        limits=c(12, 30))+
     scale_y_continuous( name ="ind. biomass at NBSS maximum (log10 gC ind-1) ", 
                         label = comma)+
          geom_point( alpha = 0.3, col = "navy", size = 2.5) +
     geom_smooth(span = 0.7)
   
   ### SST vs UVP NBSS Slope, # p-value = 0.02388, spearman rank cor
      
      # For paper
   
      library(ggplot2)
      df2 <- data.frame( NBSS_slope = zoopl_All_Atl4UVP$rob_reg__slope , SST = zoopl_All_Atl4UVP$SST_C_insitu)
      
      ggplot(df2, aes(SST, NBSS_slope ) ) +
        scale_x_continuous(name ="SST (°C)", label = comma,
                           limits=c(12, 30))+
  theme(axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10) ,
        axis.text.x = element_text(face="bold", color="#993333", 
                                        size=10) )+ 
        geom_point( alpha = 0.3, col = "navy", size = 2.5) +
        geom_smooth(method = "lm")
      
      
      
      cor.test( zoopl_All_Atl4UVP$rob_reg__slope , 
                zoopl_All_Atl4UVP$SST_C_insitu, method = "spearman")
        # p-value = 0.02388, weak but significant relationship
         # slighly flaater (more large-sized organisms) in colder waters
      
      
     summary( lm(zoopl_All_Atl4UVP$rob_reg__slope ~ zoopl_All_Atl4UVP$SST_C_insitu))
      
   
      # two UVP groups ("linear decline" and "macrozoopl. bump"), ----------
   # peaks at -6.1 and -4.3 log10(gC ind.-1) ------------
   # 1.) samples peaks at -6.1 and"log-linear decline"   from lower detection range, as expected
   # 2.) and samples with "macrozoopl. bump" at -4.3 ("macrozoopl. bump") log10(gC ind.-1) ------------
  
    # in Cold subtropical waters (SST from 20 to 27.8 ºC): both types of samples
   # in Tropical waters (SST >27.5°C): mostly  samples without "macrozoopl. bump" at -4.3 !!!
   # in cold waters (SST <20°C): mostly  samples with "macrozoopl. bump" at -4.3 !!!
   
   #    cut at -5.5 
   logCvector
   # separate into UVP_peak meso and UVP_peak macro based on Max_NBSS?
  # diff. between regions between NetZoopeakmesozoo and NetZoopeakmacrozoo ?
   # diff. in mesh size or gear  between NetZoopeakmesozoo and NetZoopeakmacrozoo?? 
   # diff. in NBSS slope between NetZoopeakmesozoo and NetZoopeakmacrozoo?? 
  
   
   summary(log10(zoopl_All_Atl4UVP$x.value.of.maximum ))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # -6.924  -5.719  -5.419  -5.227  -4.815  -4.215 
   hist ((zoopl_All_Atl4UVP$N.bins.used ))
   median_(zoopl_All_Atl4UVP$N.bins.used) # median = 9 bins used,NET zoopl.
   max(zoopl_All_Atl4UVP$N.bins.used)#  max 18 bins
   
  ## [bookmark]
   
   
   #########################
   # Analyze N bins and max position (X value of NBSS Maximum)  ----------  
       # # Net zoopl (Carbon):
   # hist (log10(zoopl_All_Atl4NET$x.value.of.maximum ))
   # summary(log10(zoopl_All_Atl4NET$x.value.of.maximum ))
   # # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # # -6.924  -5.719  -5.419  -5.227  -4.815  -4.215 
   # hist ((zoopl_All_Atl4NET$N.bins.used ))
   # median_(zoopl_All_Atl4NET$N.bins.used) # median = 9 bins used,NET zoopl.
   # max(zoopl_All_Atl4NET$N.bins.used)#  max 16 bins
   # 
   
   
   # matrix to to vector (transpose), TSWA -----------------------  
   x1 <- rep( logCvector, nrow(zoopl.NBSSmatrix_4.UVP) )
   y1 <- c ( t(zoopl.NBSSmatrix_4.UVP))
   
   # regression line (Rob. regr.)
   #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
   #points(  log10(y1) ~ x1 )
   # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
   # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
   #        lwd = 2.5, lty = 2, col = "salmon" )
   # 
   # means and medians without considering NAs  
   mean (c(4,NA, 18,34)) # NA
   median( c(4,NA, 18, 34)  ) #NA
   
   mean_   <- function(...) mean(..., na.rm=T)
   median_ <- function(...) median(..., na.rm=T)
   
   mean_ (c(4,NA, 8, 34)) #15.33, OK
   median_ (c(4,NA, 8,34)) #8, OK
   
   # apply median by columns
   
   median_vec_UVPcarb <- apply( zoopl.NBSSmatrix_4.UVP, 2, median_) 
   mean_vec_UVPcarb <- apply(zoopl.NBSSmatrix_4.UVP, 2, mean_) 
   
   median_vec_UVPcarb_NA <- apply( zoopl.NBSSmatrix_4.UVP, 2, median_NA50) 
   mean_vec_UVPcarb_NA <- apply(zoopl.NBSSmatrix_4.UVP, 2, mean_NA50) 
   
   
   #lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   # medians (curved line) ------
   logmedians5 <- log10(median_vec_UVPcarb_NA)
   length(logmedians5) # 55 data points with medians
   (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
   logmedians5_large <- logmedians5[30:55]
   logmedians5_small <- logmedians5[1:30]
   logCvector_large <- logCvector[30:55]
   logCvector_small <- logCvector[1:30]
  
    logCvector_largeUVP <- logCvector_large
   logCvector_smallUVP <- logCvector_small
   
   logmedians5_largeUVP <- logmedians5[42:55]
   logmedians5_smallUVP <- logmedians5[1:42]
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   
   
   # naive regresson acrss all data (yellow)
   abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
          lwd = 2.5, lty = 2, col = "yellow" )
   
   
   
# Since the beginning,  when reading and analyzing the zooplanktn data, 
# we selected only top-quality and selected range, as in Net Zoopl.
# after fitting rlm (and lm) models to selected sections, 
#  we divided into UVP and NET zoopl!
# So all rlm and lm slopes are from max - fporst non_NA sections, also for UVP!
   
   
   #####################
   # LINEAR MODELS -----------------
   ##################
   # UVP Zooplankton ---------------
   # Linear models for UVP  ZOOplankton  -------------
   # Factors affecting the UVP zooplankton NBSS slopes ------------------
   
   
   #### Temperature and chla vs zooplankton NBSS slope, plots and models ----------------
   
   
   # zoopl.All.Atl_p_rsq_RR_OK4
   summary( zoopl_All_Atl4_UVP$rob_reg__slope)
   summary( zoopl_All_Atl4_UVP$slope_ratio_OLSR_RR)
   dim (zoopl_All_Atl4_UVP )
   # 271 top-quality UVP zooplankton samples (carbon)
   
   
   summary( zoopl_All_Atl4_UVP$Chlorophyll__mg_m_3_insitu)# OK
   # (bad data, 100% NA)
   summary( zoopl_All_Atl4_UVP$chla)# OK
   # (bad  data, and N only about 60 chla model data)
   
   summary( zoopl_All_Atl4_UVP$SST)# OK
   # (nice data in situ, however only about 60  data...)
   summary( zoopl_All_Atl4_UVP$SST_C_insitu)# OK
   # NICE DATA! About 250 useful SST data!
   
   summary( zoopl_All_Atl4_UVP$Temperature__C_insitu)# OK
   # NICE DATA! About 250 useful Temp in situ data!
 
   summary( zoopl_All_Atl4_UVP$Temperature__C_insitu)# OK
   # NICE DATA! About 250 useful Temp in situ data!
   
   
     
   plot( zoopl_All_Atl4_UVP$Temperature__C_insitu , zoopl_All_Atl4_UVP$SST_C_insitu)
  
   plot( zoopl_All_Atl4_UVP$Temperature__C_insitu , 
           zoopl_All_Atl4_UVP$Oxygen__umol_kg_insitu)
   
    
   
   #1a.    univarate model with chl a ,log10(1+ Chl a in situ ) ----------
   # lmlogChla_insitu <- lm(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
   #                          log10(1+zoopl_All_Atl4_UVP$Chlorophyll__mg_m_3_insitu) )
   # summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
   # # CARBON (new): 26 DF,  p-value: 0.9049
   # plot(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
   #        log10(1+ zoopl_All_Atl4_UVP$Chlorophyll__mg_m_3_insitu), 
   #      col = alpha ("darkgreen", 0.3), pch = 16)
   # #abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
   # FEW data (only 28 data with chl a in situ!)
   
   
   #1b.    univarate model with chl a ,log10(1+ Chl a MODEL ) ----------
   lmlogChla_insitu <- lm(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
                            (1+zoopl_All_Atl4_UVP$chla) )
   summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
   # CARBON (new): 26 DF,  p-value: 0.9049
   plot(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
          log10(1+ zoopl_All_Atl4_UVP$chla), 
        col = alpha ("darkgreen", 0.3), pch = 16)
   #abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
   # lots of  data (only 134 data with chl a model!)
   # but... awkward satellite/model data, not significant!
   # NOT SIGNIF!
   
   
   #1c.   univarate model, SST  model ----------
   
   summary( zoopl_All_Atl4_UVP$SST)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   #    5.50   15.60   22.34   21.26   27.30   34.62      13 
   
   lmSST_C_MODEL <- lm(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
                         zoopl_All_Atl4_UVP$SST)
   summary(lmSST_C_MODEL) #OLD (biovolume , model SST data) R-squared:  0.349, 460 degrees of freedom, p-value: < 2.2e-16
   # #  100% data with SST from MODEL!)
   #  model data !!!
   # NOT SIGNIF. , p-value: 0.6403

   
   rlmSST_C_MODEL <- rlm(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
                           zoopl_All_Atl4_UVP$SST)
   summary(rlmSST_C_MODEL) #OLD (biovolume , model SST data) R-squared:  0.349, 460 degrees of freedom, p-value: < 2.2e-16
   aovperm(rlmSST_C_MODEL, nperm = 50000)
   
   # #  100% data with SST from MODEL!)
   #  model data !!!
   # SIGNIF. , p = value: 0.008
   
   plot(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
          zoopl_All_Atl4_UVP$SST, 
        col = alpha ("darkgreen", 0.3), pch = 16)
   abline(rlmSST_C_MODEL, lwd = 2.5, lty = 2, col = "darkorange")
   # SIGNIF. , p = value: 0.008
   # Model with negilible explantion power!
   # R-squared:  0.05!!
   # explains omnly 5% of all variabilty in zoopl. NBSS!
 
     plot(zoopl_All_Atl4_UVP$rob_reg__slope ~ 
          zoopl_All_Atl4_UVP$SST, 
        col = alpha ("darkgreen", 0.3), pch = 16,
        ylim = c(-2,0.5) )
   abline(rlmSST_C_MODEL, lwd = 2.5, lty = 2, col = "darkorange")
   
   
   # look at the intercept---------
   
   summary( lm ( zoopl_All_Atl4_UVP$rob_reg_intercept ~ 
                   zoopl_All_Atl4_UVP$SST ))
   plot(zoopl_All_Atl4_UVP$rob_reg_intercept ~ 
          zoopl_All_Atl4_UVP$SST, 
        col = alpha ("darkgreen", 0.3), pch = 16)
   abline(rlmSST_C_MODEL, lwd = 2.5, lty = 2, col = "darkorange")
   # weird result! warmer = more bomass! 
   # Model with negilible explantion power!
   # R-squared:  0.109!!
   # explains omnly 11% of all variabilty in zoopl. NBSS!
   
   
   
   
   # Use only UVP with max at < -5.5  log10mgC ind-1 ------------
   # separate (subsests ) by positiin of maximum -------------
   # Use only loglinearly declining NBSS , declining from mesozooplankton size ("no bump")
   
   
   zoopl_All_Atl4_UVP_no_bump <-  subset (zoopl_All_Atl4_UVP, 
                                  log10(zoopl_All_Atl4_UVP$x.value.of.maximum) < -5.5)
   
   dim(zoopl_All_Atl4_UVP_no_bump)
   
   # FOR PAPER ---------------
   ### FOR PAPER ----------------
   ### UVP bump-free  profiles (n = 157) , RR slope vs SST in situ ------------- 
 
  # summary(lm55 <- lm(zoopl_All_Atl4_UVP_no_bump$rob_reg__slope ~zoopl_All_Atl4_UVP_no_bump$SST_C_insitu)   )
  
   
   
   # Estimate Std. Error t value Pr(>|t|)  
   # (Intercept)                              0.91674    0.75194   1.219   0.2250  
   # zoopl_All_Atl4_UVP_no_bump$SST_C_insitu -0.05928    0.02785  -2.129   0.0351 *
   #   ---
   #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
   # 
   # Residual standard error: 0.5746 on 132 degrees of freedom
   # (23 observations deleted due to missingness)
   # Multiple R-squared:  0.03319,	Adjusted R-squared:  0.02587 
   # F-statistic: 4.532 on 1 and 132 DF,  p-value: 0.03513
   # 
   
   
  #    plot(  zoopl_All_Atl4_UVP_no_bump$rob_reg__slope ~zoopl_All_Atl4_UVP_no_bump$SST_C_insitu)
  #  abline(v = c(20, 27.5))
  #  abline (h = -0.4)
  # 
  #  plot(  zoopl_All_Atl4_UVP_no_bump$rob_reg__slope ~zoopl_All_Atl4_UVP_no_bump$SST_C_insitu,
  #     col = alpha ("navy", 0.3), pch = 16,
  #  xlab = "SST (°C, in situ)"  )
  #  #abline (h = -0.35, lwd = 2.5, lty = 2, col = alpha ("darkgrey", 0.6))
  #  abline ( v = c( 28.1), lwd = 2.5, lty = 2, col = alpha ("darkgrey", 0.9))
  # 
  #  y = zoopl_All_Atl4_UVP_no_bump$rob_reg__slope
  #  x =   zoopl_All_Atl4_UVP_no_bump$SST_C_insitu
  # dfxy <- data.frame( y = y , x = x)
  #   qplot(x,y, geom='smooth', span =0.5)

   # In ggplot2 you can do smooths in a number of ways, for example:
   #   
   #   library(ggplot2)
   # ggplot(dfxy, aes(x, y)) + geom_point() +
   #   geom_smooth(method = "gam", formula = y ~ poly(x, 2)) 
   # ggplot(dfxy, aes(x, y)) + geom_point() +
   #   geom_smooth(method = "loess", span = 0.3, se = FALSE)  
   # 
   # 
   # 
   #   cor.test(zoopl_All_Atl4_UVP_no_bump$rob_reg__slope ,
   #            zoopl_All_Atl4_UVP_no_bump$SST_C_insitu,
   #            method =  "spearman") # p-value = 0.007593
   #   aovperm(np= 20000,  rlm(zoopl_All_Atl4_UVP_no_bump$rob_reg__slope ~zoopl_All_Atl4_UVP_no_bump$SST_C_insitu))
   #   # p = 0.033
   # # 132 degrees of freedom, 134 samples with SST
   # NBSS slopes are significant steeper in warmer  waters (p = 0.007 Spearman rank test). 
     # VERY FLAT UVP NBSS slope around -0.3 in waters colder than  24°C  ("no bump" UVP samples only)
     #  Much STEEPER UVP NBSS with slopes around -2  to -0.4 in waters warmer  than 28°C  ("no bump" UVP samples only)
     # n = 134 "no-bump" completely linear UVP profiles
     
     
     ##################################
    ### MAPS ######################
      #MAP UVP data  --------
     # with and without bumps ----------
     # UVP  zoopl --------
     #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
    # points( zoopl_All_Atl4_UVP$Longitude,zoopl_All_Atl4_UVP$Latitude,   col="navy", pch=16)


     #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
     #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
     # symbols(zoopl_All_Atl4_UVP$Latitude ~ zoopl_All_Atl4_UVP$Longitude,
     #         circles = (10^(-1* zoopl_All_Atl4_UVP$slope_OSLR_max.to.NA)), fg= alpha("navy", 0.6), inches=0.3,
     #         add = TRUE)

     
     # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
     # points( zoopl_All_Atl4NET$Longitude,zoopl_All_Atl4NET$Latitude,   col="navy", pch=16)
     # 
     # 
     # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
     # #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
     # symbols(zoopl_All_Atl4NET$Latitude ~ zoopl_All_Atl4NET$Longitude, 
     #         circles = (10^(-1* zoopl_All_Atl4NET$slope_OSLR_max.to.NA)), fg= alpha("navy", 0.6), inches=0.3,
     #         add = TRUE)
     
     
     
     
     
     ## FINAL NBSS plots ---------  
     # For Paper ---------------
   ########
   #PLOT WITH ALL DATA and Zoopl (net zoopl. and UVP) regression 
   ######
   # plot zoopl (Net Zoopl.) regression on ALL data (incl phytopl.) ------------
   
   
   
   library(scales)
   
   lnumb.all <- length(Data.All.Atl$Longitude)
   lnumb.all# 2840 data (lines) in C format
   
   lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
   
   plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
         ylim = c(-20, 13),
         xlab = "log10(gC ind.-1)",
         ylab = "log10(gC m-3 / gC ind.-1)",
         main= "NBSS, All Atlantic, All data, n = 2840, C units")
   
   
   
   for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
     
     points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("darkgreen", 0.1))
   }
   
   # Naive linear model (all data): ------------------
   # matrix to to vector (transpose), TSWA -----------------------  
   library (MASS)
   
   x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
   y1 <- c ( t(Data.All.Atl.NBSSmatrix))
   
   df1 <- data.frame(x1, y1, ylog = log10(y1))
   names(df1)
   
   # -INF to NA!!!
   
   # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
   
   df1[df1=="-Inf"]<-NA
   
   df2<- na.omit(df1)
   dim(df1)
   dim(df2)
   
   summary(df1)
   summary(df2)
   
   mean(df2$y1)
   
   summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
   # slope = -0.954
   summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
   # slope = -0.955
   
   # naive regresson acrss all data (yellow)
   abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
          lwd = 2.5, lty = 2, col = "yellow" )
   
   
   # medians (curved line) ------
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   # lm from Net zoop medians-------
   lm8b_from_medians <- lm(  logmedians5_large ~ logCvector_large)
   abline(lm8b_from_medians,  
          lwd = 2.5, lty = 2, col = "darkblue")
   
   
   # rlm from Net zoop medians-------
   rlm7b_from_medians <- rlm(  (logmedians5_large) ~ logCvector_large )
   abline(rlm7b_from_medians, 
          lwd = 2.5, lty = 2, col = "darkorange")
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, col = "darkmagenta" , lwd = 2.5)
   
   
   
   
   ####
   # MAPS (Zoopl, filtered for linear slopes, UVP vs NET) ---------------
   
   
   #USING MAPS package (simple) --------
   library(maps)
   
   # NET zoopl --------
   # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   # points( zoopl_All_Atl4NET$Longitude,zoopl_All_Atl4NET$Latitude,   col="navy", pch=16)
   # 
   # 
   # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   # #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
   # symbols(zoopl_All_Atl4NET$Latitude ~ zoopl_All_Atl4NET$Longitude, 
   #         circles = (10^(-1* zoopl_All_Atl4NET$slope_OSLR_max.to.NA)), fg= alpha("navy", 0.6), inches=0.3,
   #         add = TRUE)
   
   summary(zoopl_All_Atl4NET$Latitude)
   
   
   
   
   
   
   
   
   
   
   #par(opar)
   
   # 4 key regions (where there is an exceptionally dense wealth of data)
   
   
   
   # subsets by regions  (Net zoopl)--------------------
   
   zoo_NET.5 <- zoopl_All_Atl4NET
   
   # TSWA ------------
   #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
   # -2°S to -12°S
   # -38°W to -25°W
   
   zoo_NET.TSWA.5 <- zoo_NET.5
   zoo_NET.TSWA.5 <-  zoo_NET.5[ (zoo_NET.TSWA.5$Latitude > -14) & (zoo_NET.TSWA.5$Latitude < -2.5) &
                                   (zoo_NET.TSWA.5$Longitude > -40) & (zoo_NET.TSWA.5$Longitude < -26),  ]
   
   symbols(zoo_NET.TSWA.5$Latitude~ zoo_NET.TSWA.5$Longitude, 
           circles = (10^(-1* zoo_NET.TSWA.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   dim(zoo_NET.TSWA.5) # only 5 samples!
   
   # CCUS -------------
   
   #  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )
   # 16.9°N to 30°N
   # -26°W to -18°W
   
   zoo_NET.CCUS.5 <- zoo_NET.5
   zoo_NET.CCUS.5 <-  zoo_NET.5[ (zoo_NET.CCUS.5$Latitude > 12) & (zoo_NET.CCUS.5$Latitude < 33) &
                                   (zoo_NET.CCUS.5$Longitude > -27.5) & (zoo_NET.CCUS.5$Longitude < -9),  ]
   
   # symbols(zoo_NET.CCUS.5$Latitude~ zoo_NET.CCUS.5$Longitude, 
   #         circles = (10^(-1* zoo_NET.CCUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
   #         add = TRUE)
   # 
   
   dim (zoo_NET.CCUS.5) # 28 samples!
   
   # EQU -----------
   #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
   # -1°S to 8°N
   # -33°W to -20°W
   
   zoo_NET.EQU.5 <- zoo_NET.5
   zoo_NET.EQU.5 <-  zoo_NET.5[ (zoo_NET.EQU.5$Latitude > -1) & (zoo_NET.EQU.5$Latitude < 6) &
                                  (zoo_NET.EQU.5$Longitude > -32) & (zoo_NET.EQU.5$Longitude < -12),  ]
   
   symbols(zoo_NET.EQU.5$Latitude~ zoo_NET.EQU.5$Longitude, 
           circles = (10^(-1* zoo_NET.EQU.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   dim(zoo_NET.EQU.5)
   #  56 samples
   
   
   
   # BUS ---------------
   #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
   # -17.5°S to -28.5°S
   # 10.2°W to 17°W
   
   
   zoo_NET.BUS.5 <- zoo_NET.5
   zoo_NET.BUS.5 <-  zoo_NET.5[ (zoo_NET.BUS.5$Latitude > -35) & (zoo_NET.BUS.5$Latitude < -17.5) &
                                  (zoo_NET.BUS.5$Longitude < 20) & (zoo_NET.BUS.5$Longitude > 10.2),  ]
   
   symbols(zoo_NET.BUS.5$Latitude~ zoo_NET.BUS.5$Longitude, 
           circles = (10^(-1* zoo_NET.BUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   dim(zoo_NET.BUS.5)
   # only 6 samples...
   
   
   # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   # #points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="red", pch=16)
   # symbols(zoopl_All_Atl4UVP$Latitude ~ zoopl_All_Atl4UVP$Longitude, 
   #         circles = (10^(-1* zoopl_All_Atl4UVP$rob_reg__slope)), fg= alpha("navy", 0.6), inches=0.3,
   #         add = TRUE)
   # map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   # points( zoopl_All_Atl4UVP$Longitude,zoopl_All_Atl4UVP$Latitude,   col="navy", pch=16)
   # 
   # 
   summary(zoopl_All_Atl4NET$Latitude)
   
   
   #par(opar)
   
   # 4 key regions (where there is an exceptionally dense wealth of data)
   
     
   # subsets by regions  (UVP zoopl)--------------------
   
   zoo_UVP.5 <- zoopl_All_Atl4UVP
   
   # TSWA ------------
   #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
   # -2°S to -12°S
   # -38°W to -25°W
   
   zoo_UVP.TSWA.5 <- zoo_UVP.5
   zoo_UVP.TSWA.5 <-  zoo_UVP.5[ (zoo_UVP.TSWA.5$Latitude > -14) & (zoo_UVP.TSWA.5$Latitude < -2.5) &
                                   (zoo_UVP.TSWA.5$Longitude > -40) & (zoo_UVP.TSWA.5$Longitude < -26),  ]
   
   symbols(zoo_UVP.TSWA.5$Latitude~ zoo_UVP.TSWA.5$Longitude, 
           circles = (10^(-1* zoo_UVP.TSWA.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   dim(zoo_UVP.TSWA.5)
   # # carbon:  34 samples
   # biovol:51 samples
   
   # CCUS -------------
   
   #  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )
   # 16.9°N to 30°N
   # -26°W to -18°W
   
   zoo_UVP.CCUS.5 <- zoo_UVP.5
   zoo_UVP.CCUS.5 <-  zoo_UVP.5[ (zoo_UVP.CCUS.5$Latitude > 12) & (zoo_UVP.CCUS.5$Latitude < 33) &
                                   (zoo_UVP.CCUS.5$Longitude > -26.5) & (zoo_UVP.CCUS.5$Longitude < -9),  ]
   
   symbols(zoo_UVP.CCUS.5$Latitude~ zoo_UVP.CCUS.5$Longitude, 
           circles = (10^(-1* zoo_UVP.CCUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   
   dim(zoo_UVP.CCUS.5)
   # # carbon: 38 samples
   # biovol: 45 samples
   
   
   # EQU -----------
   #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
   # -1°S to 8°N
   # -33°W to -20°W
   
   zoo_UVP.EQU.5 <- zoo_UVP.5
   zoo_UVP.EQU.5 <-  zoo_UVP.5[ (zoo_UVP.EQU.5$Latitude > -1) & (zoo_UVP.EQU.5$Latitude < 6) &
                                  (zoo_UVP.EQU.5$Longitude > -32) & (zoo_UVP.EQU.5$Longitude < -12),  ]
   
   symbols(zoo_UVP.EQU.5$Latitude~ zoo_UVP.EQU.5$Longitude, 
           circles = (10^(-1* zoo_UVP.EQU.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   
   dim(zoo_UVP.EQU.5)
   ## carbon: 72 samples
   # biovol: 122 samples
   
   
   # BUS ---------------
   #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
   # -17.5°S to -28.5°S
   # 10.2°W to 17°W
   
   
   zoo_UVP.BUS.5 <- zoo_UVP.5
   zoo_UVP.BUS.5 <-  zoo_UVP.5[ (zoo_UVP.BUS.5$Latitude > -35) & (zoo_UVP.BUS.5$Latitude < -17.5) &
                                  (zoo_UVP.BUS.5$Longitude < 20) & (zoo_UVP.BUS.5$Longitude > 10.2),  ]
   
   symbols(zoo_UVP.BUS.5$Latitude~ zoo_UVP.BUS.5$Longitude, 
           circles = (10^(-1* zoo_UVP.BUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
           add = TRUE)
   
   
   dim(zoo_UVP.BUS.5)
   # carbon: 17 samples
   # biovol: 31 samples
   
   
   
   
   
   
   
   #####
   ###
   # ' # TASK 3 (Zooplankton)
   #  Task 3 --------------------
   # analyse driving factors for Zoopl. NBSS
   
   
   
   # Zoopl. T-S diagram with circle size = NBSS slope -------------------
   
   

 # [bookmark]
   
   
    # II.3    
   # Mesopelagics () ----------------------
    # mesop_fish.All.Atl
   # micronekton.All.Atl
   # MESOPELAGICS ------------------
   
   
   
   ##################################
   ### MAPS ######################
   #meso_fish  --------
   #  ----------
   # meso_fish --------
   #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
   #points( mesop_fish.All.Atl$Longitude,mesop_fish.All.Atl$Latitude,   col="navy", pch=16)
   
   
   # micronekton --------
   #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-70, 60), xlim=c(-95, 30), mar=c(0,0,0,0))
  # points( micronekton.All.Atl$Longitude,micronekton.All.Atl$Latitude,   col="navy", pch=16)
   
   
   
   # Mesopelagic fish ----------
   
   dim (   mesop_fish.All.Atl) # 390 samples, NA's   :360 are NA, only 30  samples with slopes
   dim ( micronekton.All.Atl) # 140 samples, NA's   :93, only 47  samples with slopes     
   
summary (   mesop_fish.All.Atl)
   summary( micronekton.All.Atl)
   
    
   names(mesop_fish.All.Atl[37:91])
   length(names(mesop_fish.All.Atl[37:91])) # 56 NBSS data columns
   names(mesop_fish.All.Atl[37:91])
   
   
   meso_fi.NBSSmatrix <- as.matrix(mesop_fish.All.Atl[,37:91])
   dim(meso_fi.NBSSmatrix)
 rowSums(is.na(meso_fi.NBSSmatrix))
     
   meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
   rowSums(is.na(meso_fi.NBSSmatrix_no_na))
   
   
   
   dim(meso_fi.NBSSmatrix_no_na)# 86  samples with NBSS data

   
   # CLEANUP - delete  outside size range -1.3 to 0.3 log10 (g C ind.^-1)
   # clean all UVP NBSS data outside the size range
   
   # LOOP (loops through the column and replaces NBSS data with NA)
   # for UVP data only
   # # 45 NBSS classes not used for UVP NBSS
   
   df <- mesop_fish.All.Atl
  df2 <-  df[,c(41:46)]
  
   
   indices_for_mesop_fishcleanup <- c( (37:(39+37)), ((46+37):91) )
  
  # indices_for_mesop_fishcleanup <- 37:69
  
   for (i in indices_for_mesop_fishcleanup) { 
     
     
     df[,i]    <- NA
     
     
   }   
   
   # View(df)
   # View(df2)
   
   dim(df)
    
   mesop_fish.CLEAN.All.Atl <- df
   meso_fi.CLEAN.NBSSmatrix <- as.matrix(mesop_fish.CLEAN.All.Atl[,37:91])
   
  
   yvec <- as.vector(t(meso_fi.CLEAN.NBSSmatrix))
  
   df4 <- data.frame( y = yvec, x = logCvector)
   dim(df4)
   
   # linear modell based on ALL BINS (not by sample)
   
plot(log10( df4$y) ~df4$x)   

summary (mod88 <- lm(log10( df4$y) ~df4$x) )# p-value: < 2.2e-16  
rmod88 <- rlm(log10( df4$y) ~df4$x) 
aovperm(rmod88)# 2e-04,  p = 0.0002, resampled P

# p-value: < 2.2e-16  


# Bootstrap confidence interval for slope:

# confidence interval for slope:

# Bootstrap function ------------
# creating vectors to store bootstrap regression coefficients
# As function 
fun.bootstr.rr.slope <- function (x,y, nruns) {

# B = 2500 # n of runs
 B <- nruns # n of runs

df <- data.frame (x = x , y = y) 
df <- na.omit(df) # delete rows with NA

boot.out <- rep(NA, B)
vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
### LOOP ###
for(i in 1:B){ #starting loop
  
  ##creating samples of observation IDs with replacement, of same size    as original sample
  boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
  
  #matching response and explanatory variable values to bootstrap sample   observation IDs
   ysam <- df$y[boot.id]
   xsam <- df$x[boot.id]
  
  #generating bootstrap SLR model for each bootstrap sample
  boot.mod <-rlm( log10( ysam) ~  xsam )       
  
  #storing regression coefficient values for each bootstrap SLR
  boot.out[i] <- coef(boot.mod)[2]
}
  #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
  (medianboot <-median(boot.out))
  # -0.883927
  ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
  # 2.5%      97.5% 
  #   -1.0279427 -0.7463216 
  

;  boot.CI 
}

fun.bootstr.rr.slope(df4$x, df4$y, 200)

summary(lm( log10(df4$y) ~df4$x))

# fun.bootstr.rr.slope(df4$x, df4$y, 20000)


# Boostrap 95% CI (lm):
# 2.5%        50%      97.5% 
# -0.98 -0.84 -0.71 

# OLSR CI 95% (summary, mean +_ (se*1.96)) ,
mean.slope = -0.84749 ; se.slope =  0.06547
lower.lim = mean.slope - (1.96 * se.slope);   upper.lim = mean.slope + (1.96 * se.slope);   
(lower.lim); (upper.lim)
# -0.98 to -0.72, mean = -0.85

dim(df4)

# RLM:
# Boostrap 95% CI (rlm, nruns = 20000):
# 2.5%          50%         97.5% 
# -1.03       -0.89       -0.74 



   
      ##########################
   ### PLOT meso_fi (for Paper)--------------

   
   lnumb <- 34
   plot (log10(meso_fi.NBSSmatrix_no_na[lnumb, ]) ~ logCvector, 
         main = "ATLANTIC, Mesopelagic fish",
         xlab = "log10( C Biomass (gC ind.-1)) ",
         ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
         xlim = c(-4, 2), 
         ylim = c(-7, 2) ,
         pch = 16, col = "white")
   for (lnumb in 1 : nrow(meso_fi.NBSSmatrix_no_na) ) 
   {    points (log10(meso_fi.NBSSmatrix_no_na[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("navy", 0.1))
   
   }
   
   for (lnumb in 1 : nrow(meso_fi.CLEAN.NBSSmatrix) ) 
   {    points (log10(meso_fi.CLEAN.NBSSmatrix[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("darkgreen", 0.1))
     
   }
   
   
  # 
   
   
   # matrix to to vector (transpose), TSWA -----------------------  
   x1 <- rep( logCvector, nrow(meso_fi.NBSSmatrix_no_na) )
   y1 <- c ( t(meso_fi.NBSSmatrix_no_na))
   
   # regression line (Rob. regr.)
   #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
   #points(  log10(y1) ~ x1 )
   # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
   
   # ALL Data, no selection (bad)-----------   
   #  abline(rlm (log10(y1) ~ x1, , 
   #        lwd = 2.5, lty = 2, col = "pink" )
   
   
   # means and medians without considering NAs  
   mean (c(4,NA, 18,34)) # NA
   median( c(4,NA, 18, 34)  ) #NA
   
   mean_   <- function(...) mean(..., na.rm=T)
   median_ <- function(...) median(..., na.rm=T)
   
   mean_ (c(4,NA, 8, 34)) #15.33, OK
   median_ (c(4,NA, 8,34)) #8, OK
   
   # apply median by columns
   
   median_vec_ALLBiovol <- apply( meso_fi.NBSSmatrix_no_na, 2, median_) 
   mean_vec_ALLBiovol <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_) 
   
   # lines(  log10(median_vec_ALLBiovol) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate median (or mean) only if less than 50% are NA ------------
   
   vec1 <- c(4, 5,NA,NA, NA,  34)
   
   sum(is.na(vec1))
   percNA <- sum(is.na(vec1))/ length(vec1)
   
   
   mean_NA10   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.1) {
       mean(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   median_NA20   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.2) {
       median(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   
   median_vec_ALLBiovol <- apply( meso_fi.NBSSmatrix_no_na, 2, median_) 
   mean_vec_ALLBiovol <- apply(meso_fi.NBSSmatrix_no_na, 2, median_) 
   
      
   median_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, median_NA50) 
   mean_vec_ALLBiovol_NA50 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA50) 
   mean_vec_ALLBiovol_NA10 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA10) 
   
   
   
   
   length(meso_fi.NBSSmatrix_no_na[,2])
dim(meso_fi.NBSSmatrix_no_na)   
   
      
   #lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate rlm and lm from medians only (for comparison) -----------------
   # 
   # 
   # lm8_from_medians <- lm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # abline(lm8_from_medians, , 
   #        lwd = 2.5, lty = 2, col = "darkred")
   # 
   # # rlm from Net zoop medians-------
   # rlm7_from_medians <- rlm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # 
   # abline(rlm7_from_medians, 
   # lwd = 2.5, lty = 2, col = "darkgreen")
   # 
   # 
     lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,
             pch = 16, col = "red" , lwd = 3.5)
   # # 
   
   
   # medians (curved line) ------
   logmedians5 <- log10(median_vec_ALLBiovol_NA)
   length(logmedians5) # 55 data points with medians
   (cut_point_Xval <- logCvector[41]) #  cutoff at -1.2 log10 C ind.-1
   logmedians5_large <- logmedians5[41:55]
   logmedians5_small <- logmedians5[1:41]
   logCvector_large <- logCvector[41:55]
   logCvector_small <- logCvector[1:41]
   
   logCvector_largeMESO_FI <- logCvector_large
   logCvector_smallMESO_FI <- logCvector_small
   
   logmedians5_largeMESO_FI <- logmedians5_large
   logmedians5_smallMESO_FI <- logmedians5_small
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   # lm from Net zoop medians-------
   # lm8b_from_medians <- lm(  logmedians5_large ~ logCvector_large)
   # abline(lm8b_from_medians,  
   #        lwd = 2.5, lty = 2, col = "darkblue")
   # 
   
   # rlm from Net zoop medians-------
   rlm7b_from_medians <- rlm(  logmedians5_large ~ logCvector_large )
   abline(rlm7b_from_medians, 
          lwd = 2.5, lty = 2, col = "darkorange")
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   #Regression equation (dashed line)
   
   summary(rlm7b_from_medians)
   # Value    Std. Error t value 
   # (Intercept)       -3.6171   0.1245   -29.0457
   # logCvector_large  -0.9938   0.0353   -28.1711
   aovperm(rlm(logmedians5_large[1:8] ~ logCvector_large[1:8]
   )     )
   # p = 0.0017, parametric  
   # p = 0.0028, OK  
   
   rlm7b_from_medians <- rlm(  logmedians5_large ~ logCvector_large )
   aovperm(rlm7b_from_medians)
  
    lm7b_from_medians <- lm(  logmedians5_large ~ logCvector_large )
   summary(lm7b_from_medians)
  # p-value: 0.0017
   
   
   # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
   
   abline(a = -3.2915, b = -0.9548 ,  
          lwd = 2.5, lty = 2, col = "grey44" )
   
   
   
   # naive regresson acrss all data (yellow)
   # abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
   #        lwd = 2.5, lty = 2, col = "yellow" )
   # 
   
   # Linear Regression using aall data (not medians) in useful size range
   # cutoff at -1.2 log10 C ind.-1
   abline ( v = -1.2)
   abline ( v = 0.3)
   logCvector_large
   
   
  # points(log10( df4$y) ~df4$x   )
  
   
   
   # Maximum carrying capacity spectrum (98% quantiles)
   
   fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
   fun98_quant <- function(x) {quantile(x, probs = 0.98, na.rm =  TRUE) }
    
   
 
   quant95_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, fun95_quant) 
   quant98_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, fun98_quant) 
   
   lines(   log10(quant95_vec_ALLBiovol_NA) ~ logCvector,
            pch = 16, col = "red" , lwd = 5.5)
   
   lines(   log10(quant98_vec_ALLBiovol_NA) ~ logCvector,
            pch = 16, col = "purple" , lwd = 5.5)
   
   
   length(quant98_vec_ALLBiovol_NA)
   
     
   #################
   ###############
   #### Micronekton ---------------------------------------
   ### Micronekton (fish and invertebrates) ---------
   
   # micronekton.All.Atl
   
   dim (   mesop_fish.All.Atl) # 390 samples
   dim ( micronekton.All.Atl) # 140 samples
    
   
   names(micronekton.All.Atl[37:91])
   length(names(micronekton.All.Atl[37:91])) # 56 NBSS data columns
   names(micronekton.All.Atl[37:91])
   
   
   micronekt.NBSSmatrix <- as.matrix(micronekton.All.Atl[,37:91])
   dim(micronekt.NBSSmatrix)
   rowSums(is.na(micronekt.NBSSmatrix))
   
   micronekt.NBSSmatrix_no_na <-  micronekt.NBSSmatrix[rowSums(is.na(micronekt.NBSSmatrix)) != ncol(micronekt.NBSSmatrix), ]
   rowSums(is.na(micronekt.NBSSmatrix_no_na))
   
   dim(micronekt.NBSSmatrix_no_na)# 85  micronekton samples with NBSS data
   
   #   dim(meso_fi.NBSSmatrix_no_na)# 86  samples with NBSS data
   
   
   # CLEANUP - delete  outside size range -1.3 to 0.3 log10 (g C ind.^-1)
   # clean all UVP NBSS data outside the size range
   
   # LOOP (loops through the column and replaces NBSS data with NA)
   # for UVP data only
   # # 45 NBSS classes not used for UVP NBSS
   
   df <- micronekton.All.Atl
   df2 <-  df[,c(41:46)]
   
   
   indices_for_micronekton.cleanup <- c( (37:(39+38)), ((46+36):91) )
   
   # indices_for_mesop_fishcleanup <- 37:69
   
   for (i in indices_for_micronekton.cleanup) { 
     
     
     df[,i]    <- NA
     
     
   }   
   
   # View(df)
   # View(df2)
   
   
   micronektCLEAN.All.Atl <- df
   micronekt.CLEAN.NBSSmatrix <- as.matrix(micronektCLEAN.All.Atl[,37:91])
   
   
   yvec <- as.vector(t(micronekt.CLEAN.NBSSmatrix))
   
   df4 <- data.frame( y = yvec, x = logCvector)
   dim(df4)
   
   # linear modell based on ALL BINS (not by sample)
   
   plot(log10( df4$y) ~df4$x)   
   
   
   summary (mod89 <- lm(log10( df4$y) ~df4$x) )# p-value: < 2.2e-16  
   rmod89 <- rlm(log10( df4$y) ~df4$x) 
   permuco::aovperm(rmod89) #, nperm = 2000) # 2e-04,  p = 0.0002, resampled P2e-04
   #permuco::aovperm()
   # p-value: < 2.2e-16  
   
   # confidence interval for slope:
   
   # creating vectors to store bootstrap regression coefficients
   
  #  B = 2500 # n of runs
   B = 200 # n of runs
   
   boot.out <- rep(NA, B)
   vector.id <- 1:length(df4$x)   # vector of observation IDs (indices)
   ### LOOP ###
   for(i in 1:B){ #starting loop
     
     ##creating samples of observation IDs with replacement, of same size    as original sample
     boot.id <- sample(vector.id, length(df4$y), replace=TRUE) 
     
     #matching response and explanatory variable values to bootstrap sample   observation IDs
     ysam <- df4$y[boot.id]
     xsam <- df4$x[boot.id]
     
     #generating bootstrap SLR model for each bootstrap sample
     boot.mod <-rlm( log10( ysam) ~  xsam )       
     
     #storing regression coefficient values for each bootstrap SLR
     boot.out[i] <- coef(boot.mod)[2]
   }
   
   #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
   (medianboot <-median(boot.out))
   # -0.9044162
   ( boot.CI <- quantile(boot.out, c(0.025, 0.5, 0.975), type = 2) )
   # 2.5%      97.5% 
  #  -1.1175059 -0.6847944 
   
   
   
   
   
   
   ##########################
   ### PLOT micronekt (for Paper)--------------
   
   
   lnumb <- 4
   plot (log10(micronekt.NBSSmatrix_no_na[lnumb, ]) ~ logCvector, 
         main = "ATLANTIC, Micronekton",
         xlab = "log10(C Biomass (gC ind.-1)) ",
         ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
         xlim = c(-4, 2),  # xlim = c(-6, 2), 
         ylim = c(-7, 2),   # ylim = c(-8, 5) ,
               pch = 16, col = "white")
   for (lnumb in 1 : nrow(micronekt.NBSSmatrix_no_na) ) 
   {    points (log10(micronekt.NBSSmatrix_no_na[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("navy", 0.1))
     
   }
   
   
   
   for (lnumb in 1 : nrow(micronekt.CLEAN.NBSSmatrix) ) 
   {    points (log10(micronekt.CLEAN.NBSSmatrix[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("darkgreen", 0.2))
     
   }
   
   
   # 
   
   
   # matrix to to vector (transpose), TSWA -----------------------  
   x1 <- rep( logCvector, nrow(micronekt.NBSSmatrix_no_na) )
   y1 <- c ( t(micronekt.NBSSmatrix_no_na))
   
   # regression line (Rob. regr.)
   #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
   #points(  log10(y1) ~ x1 )
   # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
   
   # ALL Data, no selection (bad)-----------   
   #  abline(rlm (log10(y1) ~ x1, , 
   #        lwd = 2.5, lty = 2, col = "pink" )
   
   
   # means and medians without considering NAs  
   mean (c(4,NA, 18,34)) # NA
   median( c(4,NA, 18, 34)  ) #NA
   
   mean_   <- function(...) mean(..., na.rm=T)
   median_ <- function(...) median(..., na.rm=T)
   
   mean_ (c(4,NA, 8, 34)) #15.33, OK
   median_ (c(4,NA, 8,34)) #8, OK
   
   # apply median by columns
   
   median_vec_ALLBiovol <- apply( micronekt.NBSSmatrix_no_na, 2, median_) 
   mean_vec_ALLBiovol <- apply(micronekt.NBSSmatrix_no_na, 2, mean_) 
   
   # lines(  log10(median_vec_ALLBiovol) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate median (or mean) only if less than 50% are NA ------------
   
   vec1 <- c(4, 5,NA,NA, NA,  34)
   
   sum(is.na(vec1))
   percNA <- sum(is.na(vec1))/ length(vec1)
   
   
   mean_NA10   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.1) {
       mean(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   median_NA20   <- function(...)  {
     
     percNA <- sum(is.na(...))/ length(...)
     
     if(percNA <= 0.2) {
       median(..., na.rm=T) }
     
     else{ ;NA }
   }
   
   
   median_vec_ALLBiovol <- apply( micronekt.NBSSmatrix_no_na, 2, median_) 
   mean_vec_ALLBiovol <- apply(micronekt.NBSSmatrix_no_na, 2, median_) 
   
   
   median_vec_ALLBiovol_NA <- apply( micronekt.NBSSmatrix_no_na, 2, median_NA50) 
   mean_vec_ALLBiovol_NA50 <- apply(micronekt.NBSSmatrix_no_na, 2, mean_NA50) 
   mean_vec_ALLBiovol_NA10 <- apply(micronekt.NBSSmatrix_no_na, 2, mean_NA10) 
   
   
   length(micronekt.NBSSmatrix_no_na[,2])
   dim(micronekt.NBSSmatrix_no_na)   
   
   
   #lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,pch = 16, col = "red" , lwd = 2.5)
   
   
   # calculate rlm and lm from medians only (for comparison) -----------------
   # 
   # 
   # lm8_from_medians <- lm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # abline(lm8_from_medians, , 
   #        lwd = 2.5, lty = 2, col = "darkred")
   # 
   # # rlm from Net zoop medians-------
   # rlm7_from_medians <- rlm(  log10(median_vec_ALLBiovol_NA) ~ logCvector)
   # 
   # abline(rlm7_from_medians, 
   # lwd = 2.5, lty = 2, col = "darkgreen")
   # 
   # 
   lines(  log10(median_vec_ALLBiovol_NA) ~ logCvector,
           pch = 16, col = "red" , lwd = 3.5)
   # # 
   logCvector
   
   # medians (curved line) ------
   logmedians5 <- log10(median_vec_ALLBiovol_NA)
   length(logmedians5) # 55 data points with medians
   (cut_point_Xval <- logCvector[42]) #  cutoff at -0.9 log10 C ind.-1
   logmedians5_large <- logmedians5[42:55]
   logmedians5_small <- logmedians5[1:42]
   logCvector_large <- logCvector[42:55]
   logCvector_small <- logCvector[1:42]
   
   logCvector_largeMicronekt <- logCvector_large
   logCvector_smallMicronekt <- logCvector_small
   
   logmedians5_largeMicronekt <- logmedians5[42:55]
   logmedians5_smallMicronekt <- logmedians5[1:42]
   
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   
   
   # lm from Net zoop medians-------
   # lm8b_from_medians <- lm(  logmedians5_large ~ logCvector_large)
   # abline(lm8b_from_medians,  
   #        lwd = 2.5, lty = 2, col = "darkblue")
   # 
   
   # rlm from Net zoop medians-------
   rlm7b_from_medians <- rlm(  logmedians5_large ~ logCvector_large )
   abline(rlm7b_from_medians, 
          lwd = 2.5, lty = 2, col = "darkorange")
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(  logmedians5_large ~ logCvector_large,
           pch = 16, col = "red" , lwd = 3.5)
   
   # medians with selection cutoff point  
   lines(  logmedians5_small ~ logCvector_small,
           pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   

      # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
   
   abline(a = -3.2915, b = -0.9548 ,  
          lwd = 2.5, lty = 2, col = "grey44" )
   
   
   # Add legend -----------
   
   leg.names <- c(   "Micronekton", "All Taxa"  ) 
   col.taxa  <- c(  "darkorange" , "grey44"  ) 
   
   legend("topright", leg.names, lty = 2,lwd = 2.5,
          col = col.taxa)
   
   
   
   
     #Regression equation (dashed line)
   
   summary(rlm7b_from_medians)
   # Value    Std. Error t value 
   # (Intercept)       -3.6171   0.1245   -29.0457
   # logCvector_large  -0.9938   0.0353   -28.1711
   aovperm(rlm(logmedians5_large[1:8] ~ logCvector_large[1:8]
   )     )
   # p < 0.0001, OK  
  
   #abline ( v = -1.2)
   abline ( v =  -0.903090)
   abline (v = 0) 
   # abline ( v = 0.3)
   # 
   # # naive regresson acrss all data (yellow)
   # abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
   #        lwd = 2.5, lty = 2, col = "darkgreen" )
   # 
   
   # Maximum carrying capacity spectrum (98% quantiles)
  
   fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
  # max95_vec_ALL.micronekt.all <- apply(micronekt.9.MATRIX, 2, fun95_quant ) 
   
  #
      # lines(   log10(max95_vec_ALL.micronekt.all) ~ log10(X_vector_gCind),
      #       pch = 16, col = "red" , lwd = 5.5)
      # 
      #   lines(   log10(max98_vec_ALL.micronekt.all) ~ log10(X_vector_gCind),
      #       pch = 16, col = "purple" , lwd = 5.5)
      # 
    
        
        # quant95_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun95_quant) 
        # quant98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
        # 
        # lines(   log10(quant95_vec_ALLBiovol_NA) ~ logCvector,
        #          pch = 16, col = "red" , lwd = 5.5)
        # 
        # lines(   log10(quant98_vec_ALLBiovol_NA) ~ logCvector,
        #          pch = 16, col = "purple" , lwd = 5.5)
        # 
        # 
   
   
   
 # [bookm]
   
    # [bookmark ECOSYSTEM  NBSS] -------------------
   
   ### ECOSYSTEM  NBSS, ALL ATLANTIC --------------- 
   #### ALL TAXA ------------------
   # For Paper ALL TAxa ------------
   # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
   # FOR PAPER -  grey plot wit coloured medians 
   ### Plot All data (grey) with medians and regression lines ------------------- 

   # plot ALL NBSS Data (loop) --------------
   
   library(scales)
   
   lnumb <- length(Data.All.Atl$Longitude)
   lnumb# 2840 data (lines) in C format
   
   lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
   
   plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
         ylim = c(-7, 13),col = "white",
         xlim = c(-14, 3),
         xlab = "log10(gC ind.-1)",
         ylab = "log10(gC m-3 / gC ind.-1)",
         main= "NBSS, All Atlantic, All data, n = 2840, C units")
   
   for (i in 1: (nrow(Data.All.Atl.NBSSmatrix)) ){
     
     points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
 #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
             pch = 16, col = alpha("darkgrey", 0.05))
 
  }
#    [bookmark]
   
   # Naive linear model (all data): ------------------
   # matrix to to vector (transpose), TSWA -----------------------  
   library (MASS)
   
   x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
   y1 <- c ( t(Data.All.Atl.NBSSmatrix))
   
   df1 <- data.frame(x1, y1, ylog = log10(y1))
   names(df1)
   
   # -INF to NA!!!
   
   # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
   
   df1[df1=="-Inf"]<-NA
   
   df2<- na.omit(df1)
   dim(df1)
   dim(df2)
   
   summary(df1)
   summary(df2)
   
   mean(df2$y1)
   
   summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
   # slope = -0.954
   summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
   # slope = -0.955
   
   rlm1.ALL.ATL.ALLTaxa.naive <- rlm1
   
   # 95% CI (bootstrap, converting y to "log y")
   fun.bootstr.rr.slope(df2$x1, df2$y1, 100)
      #   fun.bootstr.rr.slope(df2$x1, df2$y1, 800)
   #   fun.bootstr.rr.slope(df2$x1, df2$y1, 50000)
   # 
   # 2.5%        50%      97.5% 
   # -0.9602939 -0.9548072 -0.9488146 
   # 
 #  permuco::aovperm(rlm1.ALL.ATL.ALLTaxa.naive)
   
   
   
   abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
          lwd = 2.5, lty = 2, col = "darkorange" )
   
   abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
          lwd = 2.5, lty = 2, col = "grey44" )
   
   summary(rlm1.ALL.ATL.ALLTaxa.naive)
   
   # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
   
   abline(a = -3.2915, b = -0.9548 ,  
          lwd = 2.5, lty = 2, col = "grey44" )
   
   
   # bootstrap the Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
   #    [bookmark]
   
   
   # medians (curved line) ------
   
   # Phtopl - ALL Areas - unfiltered OK , only medians with < 50% NA
   median_vec_PHYTOP <-  median_vec_ALLBiovol2 
   
   logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
   
   lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
           pch = 16, col = "forestgreen" , lwd = 5.5)
   
   
   
   # Net Zoopl. , medians,  - ALL Areas -------
   
 
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
  
  # Filtered NET zoopl. samples only-------------
    # lines(   logmedians5_largeNET ~ logCvector_largeNET,
    #         pch = 16, col = "dodgerblue" , lwd = 5.5)
    
   
    # Unfiltered NET zoopl. samples -------
     
   median_vec_ALLNET.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na, 2, median_NA50) 
   
   
     logmedians5 <- log10(median_vec_ALLNET.NA_unfiltr)
     length(logmedians5) # 55 data points with medians
     (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
     logmedians5_large <- logmedians5[30:55]
     logmedians5_small <- logmedians5[1:30]
     logCvector_large <- logCvector[30:55]
     logCvector_small <- logCvector[1:30]
     
     logCvector_largeNET <- logCvector_large
     logCvector_smallNET <- logCvector_small
     
     logmedians5_largeNET <- logmedians5_large
     logmedians5_smallNET <- logmedians5_small
   
     
    
     
     # Unfiltered NET zoopl. samples, ALL --------
     
      lines(   logmedians5_largeNET ~ logCvector_largeNET,
              pch = 16, col = "dodgerblue" , lwd = 5.5)
     
     
     # UVP ----------
     
   # UVP  Zoopl. , medians,  - ALL Areas -------
   
     # Unfiltered UVP zoopl. samples -------
     
     median_vec_ALLUVP.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na, 2, median_NA50) 
     
     
     logmedians5 <- log10(median_vec_ALLUVP.NA_unfiltr)
     
   
     lines(  log10(median_vec_ALLUVP.NA_unfiltr) ~ logCvector,
             pch = 16, col = "darkmagenta" , lwd = 5.5)
     
     
        
    # Mesop. fish, medians,  - ALL Areas -------
     
  
       # all medians (no cutoff selection), but olnly < 50% NA
   lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
           pch = 16, col = "darkorange" , lwd = 5.5)

   
   # Micronekt, medians, - ALL Areas -------
   
   
   # all medians (no cutoff selection), but olnly < 50% NA
   lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
            pch = 16, col = "darkred" , lwd = 5.5)
   
      
   # # medians with selection cutoff point  
   # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
   #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   # 
   
   
   # Add legend -----------
   
   leg.names <- c(  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
   col.taxa  <- c(  "forestgreen" , "dodgerblue" ,  "darkmagenta", "darkorange" , "darkred"  ) 
   
   legend("topright", leg.names, lty = 1,
                    lwd = 5.5, 
          col = col.taxa)

   ### end of grey plot with coloured medians (ALL Atlantic) ###
      
   
   
   
   
   # lm from Net zoop medians-------
   # lm8b_from_medians <- lm(  logmedians5_large ~ logCvector_large)
   # abline(lm8b_from_medians,  
   #        lwd = 2.5, lty = 2, col = "darkblue")
   # 
   
   # rlm from Net zoop medians-------
   # rlm7b_from_medians <- rlm(  (logmedians5_large) ~ logCvector_large )
   # abline(rlm7b_from_medians, 
   #        lwd = 2.5, lty = 2, col = "darkorange")
   # 
   # # all medians (no cutoff selection), but olnly < 50% NA
   # lines(  logmedians5_large ~ logCvector_large,
   #         pch = 16, col = "red" , lwd = 3.5)
   # 
   # # medians with selection cutoff point  
   # lines(  logmedians5_small ~ logCvector_small,
   #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
   # 
   
   
   ##################################
   ### MAPS ######################
   #MAP ALL DATA --------
   #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-65, 65), xlim=c(-85, 29), mar=c(0,0,0,0))
  
   # points( phytopl.All.Atl$Longitude,phytopl.All.Atl$Latitude,   col="darkgreen", pch=16)
   # 
   #  
   #  
   #  points( micronekton.All.Atl$Longitude,micronekton.All.Atl$Latitude,   col="darkred", pch=16)
   # 
   #  points( zoopl_All_Atl4UVP$Longitude,zoopl_All_Atl4UVP$Latitude,   col="darkmagenta", pch=16)
   #  
   #  points( mesop_fish.All.Atl$Longitude,mesop_fish.All.Atl$Latitude,   col="darkorange", pch=16)
   #  
   #  points( zoopl_All_Atl4_Net$Longitude,zoopl_All_Atl4_Net$Latitude,   col="navy", pch=16)
   #  
    
    # DEFINE 4 regions, V5
    # 4 key regions (where there is an exceptionally dense wealth of data)
    # version v5:
   
     # TSWA
   
    #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
   
    
    # CCUS
  
    #  rect( ybottom = 12 ,  ytop = 33,     xleft = -27.5 , xright =   -9 )
    
    
    # EQU
    
    #  rect( ybottom = -1,  ytop = 7,           xleft = -32 , xright =   -16 )
     
    # BUS
       #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
    
   
   
   ####GREY NBSS PLOTS by REGION -----------
    ### Regionwise comparisons and plots -------------------
    
    ##   Subsets by regions -----------
       # Mesopelagic fish  -------------
       
       # TSWA ------------
       #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
       
        
       Meso_fiALL.TSWA.5 <- mesop_fish.All.Atl
       Meso_fiALL.TSWA.5 <-  Meso_fiALL.TSWA.5[ (Meso_fiALL.TSWA.5$Latitude > -14) & (Meso_fiALL.TSWA.5$Latitude < -2.5) &
                                                     (Meso_fiALL.TSWA.5$Longitude > -40) & (Meso_fiALL.TSWA.5$Longitude < -26),  ]
       
       # symbols(Meso_fiALL.TSWA.5$Latitude~ Meso_fiALL.TSWA.5$Longitude, 
       #         circles = (10^(-1* Meso_fiALL.TSWA.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
       #         add = TRUE)
       
       dim(Meso_fiALL.TSWA.5) # only 12  samples!
       
       
       # CCUS -------------
       
       #  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )
       
       Meso_fiALL.CCUS.5 <- mesop_fish.All.Atl
       Meso_fiALL.CCUS.5 <-  Meso_fiALL.CCUS.5[ (Meso_fiALL.CCUS.5$Latitude > 12) & (Meso_fiALL.CCUS.5$Latitude < 33) &
                                                     (Meso_fiALL.CCUS.5$Longitude > -27.5) & (Meso_fiALL.CCUS.5$Longitude < -9),  ]
       
       # symbols(Meso_fiALL.CCUS.5$Latitude~ Meso_fiALL.CCUS.5$Longitude,
       #         circles = (10^(-1* Meso_fiALL.CCUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
       #         add = TRUE)
       # 
       
       dim (Meso_fiALL.CCUS.5) # 53 samples
       
       
       # EQU -----------
       #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
       
       Meso_fiALL.EQU.5 <- mesop_fish.All.Atl
       Meso_fiALL.EQU.5 <-  Meso_fiALL.EQU.5[ (Meso_fiALL.EQU.5$Latitude > -1) & (Meso_fiALL.EQU.5$Latitude < 6) &
                                                    (Meso_fiALL.EQU.5$Longitude > -32) & (Meso_fiALL.EQU.5$Longitude < -12),  ]
       
       # symbols(Meso_fiALL.EQU.5$Latitude~ Meso_fiALL.EQU.5$Longitude, 
       #         circles = (10^(-1* Meso_fiALL.EQU.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
       #         add = TRUE)
       
       dim(Meso_fiALL.EQU.5)
       #  26 samples
       
       
       # BUS ---------------
       #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
       
       
       Meso_fiALL.BUS.5 <- mesop_fish.All.Atl
       Meso_fiALL.BUS.5 <-  Meso_fiALL.BUS.5[ (Meso_fiALL.BUS.5$Latitude > -35) & (Meso_fiALL.BUS.5$Latitude < -17.5) &
                                                    (Meso_fiALL.BUS.5$Longitude < 20) & (Meso_fiALL.BUS.5$Longitude > 10.2),  ]
       
       # points(Meso_fiALL.BUS.5$Latitude~ Meso_fiALL.BUS.5$Longitude, 
       #          col = alpha("darkmagenta", 0.6), inches=0.3,
       #         add = TRUE)
       # 
       summary(Meso_fiALL.BUS.5) 
       dim(Meso_fiALL.BUS.5)
       
       # 14 samples
       
    
          
       ##   Subsets by regions -----------
                   # Micronekton  -----------------
       
       # TSWA ------------
       #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
       
       micronekton.All.Atl
       
       Micronekton.ALL.TSWA.5 <- micronekton.All.Atl
       Micronekton.ALL.TSWA.5 <-  Micronekton.ALL.TSWA.5[ (Micronekton.ALL.TSWA.5$Latitude > -14) & (Micronekton.ALL.TSWA.5$Latitude < -2.5) &
                                                         (Micronekton.ALL.TSWA.5$Longitude > -40) & (Micronekton.ALL.TSWA.5$Longitude < -26),  ]

       # points(Micronekton.ALL.TSWA.5$Latitude~ Micronekton.ALL.TSWA.5$Longitude,
       #         col = alpha("red", 0.6), inches=0.3,
       #         add = TRUE)

       dim(Micronekton.ALL.TSWA.5) # 83  samples!
       
       
       # CCUS -------------
       
       #  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )
       
       Micronekton.ALL.CCUS.5 <- micronekton.All.Atl
       Micronekton.ALL.CCUS.5 <-  Micronekton.ALL.CCUS.5[ (Micronekton.ALL.CCUS.5$Latitude > 12) & (Micronekton.ALL.CCUS.5$Latitude < 33) &
                                                         (Micronekton.ALL.CCUS.5$Longitude > -27.5) & (Micronekton.ALL.CCUS.5$Longitude < -9),  ]
       
       # symbols(Micronekton.ALL.CCUS.5$Latitude~ Micronekton.ALL.CCUS.5$Longitude,
       #         circles = (10^(-1* Micronekton.ALL.CCUS.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
       #         add = TRUE)
       # 
       
       dim (Micronekton.ALL.CCUS.5) # 115 samples!
       
       # EQU -----------
       #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
       
       Micronekton.ALL.EQU.5 <- micronekton.All.Atl
       Micronekton.ALL.EQU.5 <-  Micronekton.ALL.EQU.5[ (Micronekton.ALL.EQU.5$Latitude > -1) & (Micronekton.ALL.EQU.5$Latitude < 6) &
                                                        (Micronekton.ALL.EQU.5$Longitude > -32) & (Micronekton.ALL.EQU.5$Longitude < -12),  ]
       
       # symbols(Micronekton.ALL.EQU.5$Latitude~ Micronekton.ALL.EQU.5$Longitude, 
       #         circles = (10^(-1* Micronekton.ALL.EQU.5$rob_reg__slope)), fg= alpha("green", 0.6), inches=0.3,
       #         add = TRUE)
       
       dim(Micronekton.ALL.EQU.5)
       #  only 7 samples!
       
       
       # BUS ---------------
       #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
       
       
       Micronekton.ALL.BUS.5 <- micronekton.All.Atl
       Micronekton.ALL.BUS.5 <-  Micronekton.ALL.BUS.5[ (Micronekton.ALL.BUS.5$Latitude > -35) & (Micronekton.ALL.BUS.5$Latitude < -17.5) &
                                                        (Micronekton.ALL.BUS.5$Longitude < 20) & (Micronekton.ALL.BUS.5$Longitude > 10.2),  ]
       
       # points(Micronekton.ALL.BUS.5$Latitude~ Micronekton.ALL.BUS.5$Longitude, 
       #         col= alpha("red", 0.6), inches=0.3,
       #         add = TRUE)
       # 
       dim(Micronekton.ALL.BUS.5)
       # 24 samples

       
       # 
       # ### ALL Taxa --------------
       # # subsets by region, ALL Taxa , All samples unfiltered -----
       # 
       # # Data.All.Atl
       # 
       #  length(Data.All.Atl$ConsecutiveNumber) # 2840  datasets, ALL Data
       # 
       # 
       # #   maps::map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-65, 65), xlim=c(-85, 29), mar=c(0,0,0,0))
       # 
       # points( Data.All.Atl$Longitude,Data.All.Atl$Latitude,   col="forestgreen", pch=16)
       # 
       # 
       # # TSWA ------------
       # #  rect( ybottom = -14,  ytop = -2.5,           xleft =  -40, xright =   -26)
       # 
       # 
       Data.All.ALL.TSWA.5 <- Data.All.Atl
       Data.All.ALL.TSWA.5 <-  Data.All.ALL.TSWA.5[ (Data.All.ALL.TSWA.5$Latitude > -14) & (Data.All.ALL.TSWA.5$Latitude < -2.5) &
                                                            (Data.All.ALL.TSWA.5$Longitude > -40) & (Data.All.ALL.TSWA.5$Longitude < -26),  ]

       # points(Data.All.ALL.TSWA.5$Latitude~ Data.All.ALL.TSWA.5$Longitude,
       #        col="red", pch=16,
       #        add = TRUE)
       # 
       # dim(Data.All.ALL.TSWA.5) # 83  samples!
       # 
       # 
       # # CCUS -------------
       # 
       # #  rect( ybottom = 12 ,  ytop = 33,           xleft = -27.5 , xright =   -9 )
       # 
       Data.All.ALL.CCUS.5 <- Data.All.Atl
       Data.All.ALL.CCUS.5 <-  Data.All.ALL.CCUS.5[ (Data.All.ALL.CCUS.5$Latitude > 12) & (Data.All.ALL.CCUS.5$Latitude < 33) &
                                                            (Data.All.ALL.CCUS.5$Longitude > -27.5) & (Data.All.ALL.CCUS.5$Longitude < -9),  ]
       # 
       # points(Data.All.ALL.CCUS.5$Latitude~ Data.All.ALL.CCUS.5$Longitude,
       #        col="red", pch=16,
       #          add = TRUE)
       #  
       # 
       # dim (Data.All.ALL.CCUS.5) # 115 samples!
       # 
       # # EQU -----------
       # #  rect( ybottom = -1,  ytop = 6,           xleft = -32 , xright =   -12 )
       # 
       Data.All.ALL.EQU.5 <- Data.All.Atl
       Data.All.ALL.EQU.5 <-  Data.All.ALL.EQU.5[ (Data.All.ALL.EQU.5$Latitude > -1) & (Data.All.ALL.EQU.5$Latitude < 6) &
                                                          (Data.All.ALL.EQU.5$Longitude > -32) & (Data.All.ALL.EQU.5$Longitude < -12),  ]
       # 
       # points(Data.All.ALL.EQU.5$Latitude~ Data.All.ALL.EQU.5$Longitude,
       #        col="red", pch=16,
       #         add = TRUE)
       # 
       # dim(Data.All.ALL.EQU.5)
       # #  only 7 samples!
       # 
       # 
       # # BUS ---------------
       # #  rect( ybottom = -17.5,  ytop = -35,           xleft =  20, xright =   10.2)
       # 

       Data.All.ALL.BUS.5 <- Data.All.Atl
       Data.All.ALL.BUS.5 <-  Data.All.ALL.BUS.5[ (Data.All.ALL.BUS.5$Latitude > -35) & (Data.All.ALL.BUS.5$Latitude < -17.5) &
                                                          (Data.All.ALL.BUS.5$Longitude < 20) & (Data.All.ALL.BUS.5$Longitude > 10.2),  ]
       # 
       # points(Data.All.ALL.BUS.5$Latitude~ Data.All.ALL.BUS.5$Longitude, 
       #        col="red", pch=16,
       #        add = TRUE)
       # 
       # dim(Data.All.ALL.BUS.5)
       # # 24 samples
       # 
       # 
       
       
       #### PLOT NBSS by AREA, ALL groups,  ----------
       ####GREY AND COLOURED NBSS PLOTS by REGION -----------
       ### Regionwise comparisons and plots -------------------
       
       #### PLOT NBSS by AREA, ALL groups, and by group ----------
       # NBSS by areas , ALL data, colour codes --------------       
       
       
       # plot ALL NBSS Data (loop) --------------
       
       
       Data.All.Atl.NBSSmatrix <- as.matrix(Data.All.Atl[,37:91])
       
       
       library(scales)
       
       lnumb.all <- length(Data.All.Atl$Longitude)
       lnumb.all# 2840 data (lines) in C format
       
       lnumb <-  3
       
       plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-8, 13.5),
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "NBSS, All Atlantic, All data, n = 2840, C units")
       
       
       
       for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
         
         points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                 # pch = 16, col = alpha("darkgreen", 0.1))
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
       x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
       y1 <- c ( t(Data.All.Atl.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))

       #   [bookmark] ALL ATLANTIC LINEAR MODEL
     
            # ALL DATA, ALL taxa , RR linear model (rlm)
            # slope = -0.9548 +-   0.0026 (SE)
            #  Intercept)   -3.2915  +-   0.0191 (SE) 
            # R2: 0.86,  df: 19422
       
       #fun.bootstr.rr.summary.simple.complete.p.R2.95CI <- function (x, y , n.runsBOOT = 1000, nruns.aovperm = 2000) 
        # fun.bootstr.rr.summary.simple.complete.p.R2.95CI( y = log10(df2$y1), x = (df2$x1)   )
       
       
       abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "darkorange" )
    
       text (x = -2.5, y = 10, "RR slope, ALL =  -0.955" )
       
          
       # TSWA ----------
       
       
       # plot ALL TSWA NBSS Data (loop) --------------
       
       Data.All.ALL.TSWA.5.NBSSmatrix <- as.matrix(Data.All.ALL.TSWA.5 [,37:91])
       
       dim(Data.All.ALL.TSWA.5.NBSSmatrix) # 319 samples, all taxa,unfiltered
       
       library(scales)
       
       lnumb.all <- length(Data.All.ALL.TSWA.5$Longitude)
       lnumb.all# 2840 data (lines) in C format
       
       lnumb <-  3
       
       plot (log10(Data.All.ALL.TSWA.5.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-8, 13.5),
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "NBSS, TSWA, All data, n = 319, C units")
       
       
       
       for (i in 1:nrow(Data.All.ALL.TSWA.5.NBSSmatrix)){  
         
         points (log10(Data.All.ALL.TSWA.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 pch = 16, col = alpha(Data.All.ALL.TSWA.5$target.org.colour[i], 0.1))
         # pch = 16, col = alpha("darkgreen", 0.1))
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.TSWA.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.TSWA.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope (all Atlantic) = -0.955
       # slope  (TSWA) = -0.9085, flatter than average
      
      # (Intercept)   -3.6148    0.0566
       # slope    -0.9085    0.0072 
        #R-squared:  0.86
       # 2810 degrees of freedom
       
       abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "darkorange" )
       
       text (x = -2.5, y = 10, "RR slope, TSWA = -0.909 +- 0.007" )
       
     #  fun.bootstr.rr.summary.simple.complete.p.R2.95CI( y = log10(df2$y1), x = (df2$x1)   )
       
       
       
       
       
       # CCUS ----------
       
       
       # plot ALL NBSS Data (loop) --------------
       
       library(scales)
       
       Data.All.ALL.CCUS.5.NBSSmatrix <- as.matrix(Data.All.ALL.CCUS.5 [,37:91])
       
       dim(Data.All.ALL.CCUS.5.NBSSmatrix) # 463 samples, all taxa,unfiltered
      
       
        lnumb.all <- length(Data.All.ALL.CCUS.5$Longitude)
       lnumb.all# 463 data (lines) in C format
       
       lnumb <-  3
       
       plot (log10(Data.All.ALL.CCUS.5.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-8, 13.5),
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "NBSS, CCUS, All data, n = 463, C units")
       
       
       for (i in 1:nrow(Data.All.ALL.CCUS.5.NBSSmatrix)){  
         
         points (log10(Data.All.ALL.CCUS.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 pch = 16, col = alpha(Data.All.ALL.CCUS.5$target.org.colour[i], 0.1))
         # pch = 16, col = alpha("darkgreen", 0.1))
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.CCUS.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.CCUS.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # Intercept   -3.2412    0.0478 
       #     # slope = -0.9172    0.00875
     
       # Intercept)   -3.2463   +- 0.0495 
       # slope       -0.9408    +-0.0079   
       # R-squared:  0.8169
       # 2670 degrees of freedom
       
       abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "darkorange" )
       
       text (x = -2.5, y = 10, "RR slope, CCUS = -0.940 +- 0.008" )
       
       #  fun.bootstr.rr.summary.simple.complete.p.R2.95CI( y = log10(df2$y1), x = (df2$x1)   )
       
       
       fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
     #  max95_vec_ALL.micronekt.all <- apply(micronekt.9.MATRIX, 2, fun95_quant ) 
       
       #
       # lines(   log10(max95_vec_ALL.micronekt.all) ~ log10(X_vector_gCind),
       #       pch = 16, col = "red" , lwd = 5.5)
       # 
       #   lines(   log10(max98_vec_ALL.micronekt.all) ~ log10(X_vector_gCind),
       #       pch = 16, col = "purple" , lwd = 5.5)
       # 
       
       Micronekton.ALL.CCUS.5
       
       #quant95_vec_ALLBiovol_NA.CCUS <- apply( Micronekton.ALL.CCUS.5, 2, fun95_quant) 
     #  quant98_vec_ALLBiovol_NA.CCUS <- apply( Micronekton.ALL.CCUS.5, 2, fun98_quant) 
       
       # lines(   log10(quant95_vec_ALLBiovol_NA) ~ logCvector,
       #          pch = 16, col = "red" , lwd = 5.5)
       # 
       lines(   log10(quant98_vec_ALLBiovol_NA) ~ logCvector,
                pch = 16, col = "purple" , lwd = 5.5)
       
       
       
       
         
       # EQU
       
       
       # plot ALL NBSS Data (loop) --------------
       
       library(scales)
    
       Data.All.ALL.EQU.5.NBSSmatrix <- as.matrix(Data.All.ALL.EQU.5 [,37:91])
       
       dim(Data.All.ALL.EQU.5.NBSSmatrix) # 406 samples, all taxa,unfiltered
       
          
       lnumb.all <- length(Data.All.ALL.EQU.5$Longitude)
       lnumb.all# 406 data (lines) in C format
       
       lnumb <-  3
       
       plot (log10(Data.All.ALL.EQU.5.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-8, 13.5),
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "NBSS, EQU, All data, n = 406, C units")
       
       for (i in 1:nrow(Data.All.ALL.EQU.5.NBSSmatrix)){  
         
           
         points (log10(Data.All.ALL.EQU.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 pch = 16, col = alpha(Data.All.ALL.EQU.5$target.org.colour[i], 0.1))
         }
         
         # pch = 16, col = alpha("darkgreen", 0.1))
       
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.EQU.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.EQU.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1eq <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.921494   +- 0.009654 
       summary(rlm1eq <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # Intercept =    -3.2412    0.0478
         # slope = -0.9172 +-  0.0087
       
       
       abline(rlm1eq, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "darkorange" )
       
       
       text (x = -2.5, y = 10, "RR slope, EQU = -0.917  +- 0.009  " )
       
       #  fun.bootstr.rr.summary.simple.complete.p.R2.95CI( y = log10(df2$y1), x = (df2$x1)   )
       
       
       
       
       
       #BUS -----------------
       
       
       # plot ALL NBSS Data (loop) --------------
       
       library(scales)
       
       Data.All.ALL.BUS.5.NBSSmatrix <- as.matrix(Data.All.ALL.BUS.5 [,37:91])
       
       dim(Data.All.ALL.BUS.5.NBSSmatrix) # 202 samples, all taxa,unfiltered
       
       
       lnumb.all <- length(Data.All.ALL.BUS.5$Longitude)
       lnumb.all# 202 data (lines) in C format
       
       lnumb <-  3
       
       plot (log10(Data.All.ALL.BUS.5.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-8, 13.5),
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "NBSS, BUS, All data, n = 202, C units")
       
       
       for (i in 1:nrow(Data.All.ALL.BUS.5.NBSSmatrix)){
         
         points (log10(Data.All.ALL.BUS.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 pch = 16, col = alpha(Data.All.ALL.BUS.5$target.org.colour[i], 0.1))
         # pch = 16, col = alpha("darkgreen", 0.1))
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.BUS.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.BUS.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # sllpe all = -0.955
       # slope BUS =    -0.9997   
       
       # (Intercept)   -3.4922  +-  0.0595  
       # df2$x1        -0.9997  +-   0.0088  
       # 1982 DF
       # R-squared:  0.868
       
       
       abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "darkorange" )
       
       
       text (x = -2.5, y = 10, "RR slope, BUS = -1.000 +- 0.009  " )
       
       #  fun.bootstr.rr.summary.simple.complete.p.R2.95CI( y = log10(df2$y1), x = (df2$x1)   )
       
       
       #####################################
       #####################################
       # COMPARISONS AND TESTS - for paper #
       #####################################
       #####################################
       
       ### Compare slopes and intercepts between regions (ANCOVA & aovperm) ------------
       
      #  slopes and intercepts (ALL Data, ALL TAXA, by regions)
 
       df <- data.frame(
         Region = c("TSWA", "EQU", "CCUS", "BUS"),
         slope = c(-0.909, -0.917, -0.941, -1.000),
         slope_se = c(0.007, 0.009, 0.008, 0.009),
         intercept = c(-3.61, -3.24, -3.25, -3.5),
         intercept_se = c(0.06, 0.05, 0.05, 0.06),
         n = c(2812, 2877, 2672, 1986)
       )  
       
       # Calculate SDs
       df$slope_sd     <- df$slope_se * sqrt(df$n)
       df$intercept_sd <- df$intercept_se * sqrt(df$n)
       
       # View result
       print(df)
       
       # with t-test (preliminary)
       library(BSDA)
       # t.test, BUS slope vs TSWA slope
       tsum.test(mean.x = -0.909, mean.y =  -1,
                 s.x = 0.371198, s.y = 0.40108,
                 n.x = 2812 , n.y = 1986 )
       # t-test: p-value = 1.869e-15
       
       # t.test, EQU slope vs TSWA slope
       tsum.test(mean.x = -0.909, mean.y =   -0.941,
                 s.x = 0.371198, s.y = 0.482739,
                 n.x = 2812 , n.y = 2877 )
       # t-test:  p-value = 1.757e-15
   
       # t.test, EQU slope vs CCUS slope
       tsum.test(mean.x = -0.941, mean.y =  -1.000,
                 s.x = 0.482739 , s.y = -0.917,
                 n.x = 2812 , n.y = 2877 )
       # t-test: p-value = 8.754e-07
      
     region1 = "TSWA" ; region2 = "BUS"
          mean.1 = df$slope[df$Region == region1]
        mean.2 = df$slope[df$Region == region2]
        s.1 =  df$slope_sd [df$Region == region1]
        s.2 = df$slope_sd[df$Region == region2]
        n.1 =  df$n[df$Region == region1]
        n.2 = df$n[df$Region == region2]
        
       tsum.test(mean.x = mean.1, mean.y =  mean.2,
                 s.x = s.1 , s.y = s.2,
                 n.x = n.1 , n.y = n.2 )
     #  p-value = 1.869e-15, t-test compariing slopes
      
       region1 = "CCUS" ; region2 = "BUS"
       mean.1 = df$slope[df$Region == region1]
       mean.2 = df$slope[df$Region == region2]
       s.1 =  df$slope_sd [df$Region == region1]
       s.2 = df$slope_sd[df$Region == region2]
       n.1 =  df$n[df$Region == region1]
       n.2 = df$n[df$Region == region2]
       
       tsum.test(mean.x = mean.1, mean.y =  mean.2,
                 s.x = s.1 , s.y = s.2,
                 n.x = n.1 , n.y = n.2 )
       #  p-value = p-value = 9.949e-07 , t-test comparing slopes
       
       region1 = "CCUS" ; region2 = "EQU"
       mean.1 = df$slope[df$Region == region1]
       mean.2 = df$slope[df$Region == region2]
       s.1 =  df$slope_sd [df$Region == region1]
       s.2 = df$slope_sd[df$Region == region2]
       n.1 =  df$n[df$Region == region1]
       n.2 = df$n[df$Region == region2]
       
       tsum.test(mean.x = mean.1, mean.y =  mean.2,
                 s.x = s.1 , s.y = s.2,
                 n.x = n.1 , n.y = n.2 )
       #  p-value = p-value = 0.0463 , t-test comparing slopes
       # CCUS and EQU were not signiifocnt consideng a Bonferroni correction
        # 5 groups = 10 comparisons, p crit = 0.05 / 10 :  p crit = 0.005
       # p corr = 0.463 ,  n.s.
       
       # Pemutation(Aovperm))
       
       sam1 <- rnorm ( n = n.1 , mean =   mean.1, sd =  s.1)
        sam2 <- rnorm ( n = n.2 , mean =   mean.2, sd =  s.2)
        library(coin)
        
      wilcox.test(sam1, sam2) # p-value = 0.8865, n.s.
         ##independence_test( sam1, sam2)
       
      
      region1 = "TSWA" ; region2 = "CCUS"
      mean.1 = df$slope[df$Region == region1]
      mean.2 = df$slope[df$Region == region2]
      s.1 =  df$slope_sd [df$Region == region1]
      s.2 = df$slope_sd[df$Region == region2]
      n.1 =  df$n[df$Region == region1]
      n.2 = df$n[df$Region == region2]
      
      tsum.test(mean.x = mean.1, mean.y =  mean.2,
                s.x = s.1 , s.y = s.2,
                n.x = n.1 , n.y = n.2 )
      #  p-value = p-value = 0.002622, ok
      
      # Pemutation(Aovperm))
      
      sam1 <- rnorm ( n = n.1 , mean =   mean.1, sd =  s.1)
      sam2 <- rnorm ( n = n.2 , mean =   mean.2, sd =  s.2)
      library(coin)
      
      wilcox.test(sam1, sam2) # p-value = 0.8865, n.s.
      ##independence_test( sam1, sam2)
      
      
      region1 = "TSWA" ; region2 = "EQU"
      mean.1 = df$slope[df$Region == region1]
      mean.2 = df$slope[df$Region == region2]
      s.1 =  df$slope_sd [df$Region == region1]
      s.2 = df$slope_sd[df$Region == region2]
      n.1 =  df$n[df$Region == region1]
      n.2 = df$n[df$Region == region2]
      
      tsum.test(mean.x = mean.1, mean.y =  mean.2,
                s.x = s.1 , s.y = s.2,
                n.x = n.1 , n.y = n.2 )
      #  p-value = p-value = 0.48  n.s.
      
      # Pemutation(Aovperm))
      
      sam1 <- rnorm ( n = n.1 , mean =   mean.1, sd =  s.1)
      sam2 <- rnorm ( n = n.2 , mean =   mean.2, sd =  s.2)
      library(coin)
      
      wilcox.test(sam1, sam2) # p-value = 0.8865, n.s.
      ##independence_test( sam1, sam2)
      
      
      region1 = "TSWA" ; region2 = "BUS"
      mean.1 = df$slope[df$Region == region1]
      mean.2 = df$slope[df$Region == region2]
      s.1 =  df$slope_sd [df$Region == region1]
      s.2 = df$slope_sd[df$Region == region2]
      n.1 =  df$n[df$Region == region1]
      n.2 = df$n[df$Region == region2]
      
      tsum.test(mean.x = mean.1, mean.y =  mean.2,
                s.x = s.1 , s.y = s.2,
                n.x = n.1 , n.y = n.2 )
      #  p-value = p-value = p-value = 1.869e-15, ok
      
      
      region1 = "CCUS" ; region2 = "BUS"
      mean.1 = df$slope[df$Region == region1]
      mean.2 = df$slope[df$Region == region2]
      s.1 =  df$slope_sd [df$Region == region1]
      s.2 = df$slope_sd[df$Region == region2]
      n.1 =  df$n[df$Region == region1]
      n.2 = df$n[df$Region == region2]
      
      tsum.test(mean.x = mean.1, mean.y =  mean.2,
                s.x = s.1 , s.y = s.2,
                n.x = n.1 , n.y = n.2 )
      #  p-value = p-value = p-value = p-value = 9.949e-07, ok
      
      
      region1 = "EQU" ; region2 = "BUS"
      mean.1 = df$slope[df$Region == region1]
      mean.2 = df$slope[df$Region == region2]
      s.1 =  df$slope_sd [df$Region == region1]
      s.2 = df$slope_sd[df$Region == region2]
      n.1 =  df$n[df$Region == region1]
      n.2 = df$n[df$Region == region2]
      
      tsum.test(mean.x = mean.1, mean.y =  mean.2,
                s.x = s.1 , s.y = s.2,
                n.x = n.1 , n.y = n.2 )
      #  p-value = p-value = p-value = p-value = 9.949e-07, ok
     
      
      # Pemutation(Aovperm))
      
      sam1 <- rnorm ( n = n.1 , mean =   mean.1, sd =  s.1)
      sam2 <- rnorm ( n = n.2 , mean =   mean.2, sd =  s.2)
      library(coin)
      
      wilcox.test(sam1, sam2) # p-value = 0.8865, n.s.
     
         df99 <- data.frame(
          value = c(sam1, sam2),
          group = as.factor( c(rep("sam1", length(sam1)), rep("sam2", length(sam2))) )
        )
      
      # View the data frame
      print(df99)
      
      # permutation test
       independence_test( df99$value ~df99$group) # p-value = 4.094e-09
    
      # aovperm(  df99$value ~ df99$group,  np = 20000)   
      
        aovperm(  df99$value ~ df99$group,  np = 2000)  # p = 5e-04 
       
       
       #####################################
       #####################################
        
       # Plots For Paper ------------
       # TSWA
       # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
       # FOR PAPER -  grey plot wit coloured medians 
       ### Plot TSWA with medians and regression lines ------------------- 
       
       # plot TSWA  NBSS Data (loop) --------------
       
       library(scales)
       
       lnumb <- length(Data.All.Atl$Longitude)
       lnumb# 2840 data (lines) in C format
       
       lnumb <-3
       
       plot (log10(Data.All.ALL.TSWA.5.NBSSmatrix [lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-7, 13),col = "white",
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "TSWA")
       
       
       for (i in 1:nrow(Data.All.ALL.TSWA.5.NBSSmatrix)){  
         
         points (log10(Data.All.ALL.TSWA.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                 pch = 16, col = alpha("darkgrey", 0.3))
         
       }
      
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), TSWA -----------------------  
       library (MASS)
       
    
          x1 <- rep( logCvector, nrow(Data.All.ALL.TSWA.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.TSWA.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.955
       
       rlm1TSWA.All.naive <- rlm1
       
       
       abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("slateblue4", 0.4) )
       
       abline(rlm1TSWA.All.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("darkorange", 0.9) )
       
        
       # grey medians for ALL data (backround grey curved lines) ------
 
       phytopl.All.Atl.NBSSmatrix.1.TSWA <- as.matrix(phytopl.unfiltr.ALL.TSWA.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.TSWA)
      # 110 samples
       
       median_vec_ALLBiovol2.TSWA <- apply(phytopl.All.Atl.NBSSmatrix.1.TSWA, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.TSWA)
       
       median_vec_ALLBiovol2.TSWA <- as.vector(unlist(median_vec_ALLBiovol2.TSWA))
       
       median_vec_PHYTOP.TSWA <-  median_vec_ALLBiovol2.TSWA 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
    
      
          
       lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = alpha("grey27", 0.5), lwd = 4.5)
      
       
       # Darkgrey median lines as background ---------- 
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col =alpha("grey27", 0.5), lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
                pch = 16, col = alpha("grey27", 0.5) , lwd = 4.5)
     
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeNET ~ logCvector_largeNET,
                pch = 16, col = alpha("grey27", 0.6) , lwd = 4.7)
     
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, 
               col = alpha("grey27", 0.55), lwd = 4.6)
    
    
       # median lines for TSWA area ----------
       
       phytopl.All.Atl.NBSSmatrix.1.TSWA <- as.matrix(phytopl.unfiltr.ALL.TSWA.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.TSWA)
       # 110 samples
       
       median_vec_ALLBiovol2.TSWA <- apply(phytopl.All.Atl.NBSSmatrix.1.TSWA, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.TSWA)
       
       median_vec_ALLBiovol2.TSWA <- as.vector(unlist(median_vec_ALLBiovol2.TSWA))
       
       median_vec_PHYTOP.TSWA <-  median_vec_ALLBiovol2.TSWA 
           
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       lines(  log10(median_vec_PHYTOP.TSWA) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = "forestgreen" , lwd = 5.5) # OK
       
       
       # Net Zoopl. , medians,  TSWA , unfiltered-------
       
       # all medians 
       
       zoopl_All_Atl.unfiltrTSWA <-  subset( Data.All.ALL.TSWA.5, target_organisms == "zooplankton")  
       dim(zoopl_All_Atl.unfiltrTSWA) # 114
       
       zoopl_All_Atl.unfiltrTSWA_NET <- subset (zoopl_All_Atl.unfiltrTSWA, Gear != "UVP")
        dim(zoopl_All_Atl.unfiltrTSWA_NET)  # 35 samples in TSWA, unfiltered     
   

       zoo_NET.NBSSmatrixTSWA <- as.matrix(zoopl_All_Atl.unfiltrTSWA_NET[,37:91])
     names(as.data.frame (zoo_NET.NBSSmatrixTSWA))
         dim(zoo_NET.NBSSmatrixTSWA) # 83 samples
       rowSums(is.na(zoo_NET.NBSSmatrixTSWA))
       
       zoo_NET.NBSSmatrixTSWA_no_na <-  zoo_NET.NBSSmatrixTSWA[rowSums(is.na(zoo_NET.NBSSmatrixTSWA)) != ncol(zoo_NET.NBSSmatrixTSWA), ]
       rowSums(is.na(zoo_NET.NBSSmatrixTSWA_no_na))
       dim(zoo_NET.NBSSmatrixTSWA_no_na)
       # 14 samples with real data (no NA's), unfiltered
       
       median_vec_zoo_NET.NBSSmatrixTSWA_no_na <- apply( zoo_NET.NBSSmatrixTSWA_no_na, 2, median_NA50) 
       
       
       length(zoo_NET.NBSSmatrixTSWA_no_na[,2])
       dim(zoo_NET.NBSSmatrixTSWA_no_na)   
       # 14 samples with real data (no NA's), unfiltered
      
        
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_zoo_NET.NBSSmatrixTSWA_no_na)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[30:55]
       logmedians5_small <- logmedians5[1:30]
       logCvector_large <- logCvector[30:55]
       logCvector_small <- logCvector[1:30]
       
       logCvector_largezoo_NET <- logCvector_large
       logCvector_smallzoo_NET <- logCvector_small
       
       logmedians5_largezoo_NET_TSWA <- logmedians5_large
       
       
       lines(   logmedians5_largezoo_NET_TSWA ~ logCvector_largeNET,
                pch = 16, col = "dodgerblue" , lwd = 5.5)
       
      
      
        # UVP  Zoopl., medians,  TSWA, unfiltered -------
      
       zoopl_All_Atl.unfiltrTSWA <-  subset( Data.All.ALL.TSWA.5, target_organisms == "zooplankton")  
       dim(zoopl_All_Atl.unfiltrTSWA) # 114
       
       zoopl_All_Atl.unfiltrTSWA_UVP <- subset (zoopl_All_Atl.unfiltrTSWA, Gear == "UVP")
       dim(zoopl_All_Atl.unfiltrTSWA_UVP)  # 79 samples in TSWA, unfiltered     
       
       
       zoo_UVP.NBSSmatrixTSWA <- as.matrix(zoopl_All_Atl.unfiltrTSWA_UVP[,37:91])
       names(as.data.frame (zoo_UVP.NBSSmatrixTSWA))
       dim(zoo_UVP.NBSSmatrixTSWA) # 79 samples
       rowSums(is.na(zoo_UVP.NBSSmatrixTSWA))
       
       zoo_UVP.NBSSmatrixTSWA_no_na <-  zoo_UVP.NBSSmatrixTSWA[rowSums(is.na(zoo_UVP.NBSSmatrixTSWA)) != ncol(zoo_UVP.NBSSmatrixTSWA), ]
       rowSums(is.na(zoo_UVP.NBSSmatrixTSWA_no_na))
       dim(zoo_UVP.NBSSmatrixTSWA_no_na)
       # 79 samples with real data (no NA's), unfiltered
       
       median_vec_zoo_UVP.NBSSmatrixTSWA_no_na <- apply( zoo_UVP.NBSSmatrixTSWA_no_na, 2, median_NA50) 
       
       
       # length(zoo_UVP.NBSSmatrix_no_na[,2])
       # dim(zoo_UVP.NBSSmatrix_no_na)   
       # # 40 samples with data, fltered, vs 79 unfiltred
        
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_zoo_UVP.NBSSmatrixTSWA_no_na)
       
       
       median_vec_UVPcarb_NA_TSWA.unfiltered <- median_vec_zoo_UVP.NBSSmatrixTSWA_no_na
       
       lines(  log10(median_vec_UVPcarb_NA_TSWA.unfiltered) ~ logCvector,pch = 16, col = "darkmagenta" , lwd = 5.5)
       
       # overplotting
       lines(   logmedians5_largezoo_NET_TSWA ~ logCvector_largeNET,
                pch = 16, col = "dodgerblue" , lwd = 5.5)
       
       
       
       # Mesop. fish, medians,  TSWA -------
       
       
       
       meso_fi.NBSSmatrix <- as.matrix(Meso_fiALL.TSWA.5[,37:91])
       dim(meso_fi.NBSSmatrix)
       rowSums(is.na(meso_fi.NBSSmatrix))
       
       meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
       rowSums(is.na(meso_fi.NBSSmatrix_no_na))
       
       median_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, median_NA50) 
       mean_vec_ALLBiovol_NA50 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA50) 
       mean_vec_ALLBiovol_NA10 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA10) 
       
       
       length(meso_fi.NBSSmatrix_no_na[,2])
       dim(meso_fi.NBSSmatrix_no_na)   
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[41]) #  cutoff at -1.2 log10 C ind.-1
       logmedians5_large <- logmedians5[41:55]
       logmedians5_small <- logmedians5[1:41]
       logCvector_large <- logCvector[41:55]
       logCvector_small <- logCvector[1:41]
       
       logCvector_largeMESO_FI <- logCvector_large
       logCvector_smallMESO_FI <- logCvector_small
       
       logmedians5_largeMESO_FI <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col = "darkorange" , lwd = 5.5)
       
       
       # Micronekt, medians,  TSWA -------
       
       
       Micronekton.NBSSmatrix <- as.matrix(Micronekton.ALL.TSWA.5[,37:91])
       dim(Micronekton.NBSSmatrix) # 83 samples
       rowSums(is.na(Micronekton.NBSSmatrix))
       
       Micronekton.NBSSmatrix_no_na <-  Micronekton.NBSSmatrix[rowSums(is.na(Micronekton.NBSSmatrix)) != ncol(Micronekton.NBSSmatrix), ]
       rowSums(is.na(Micronekton.NBSSmatrix_no_na))
       dim(Micronekton.NBSSmatrix_no_na)
       # 40 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(Micronekton.NBSSmatrix_no_na[,2])
       dim(Micronekton.NBSSmatrix_no_na)   
       # 40 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[42]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[42:55]
       logmedians5_small <- logmedians5[1:42]
       logCvector_large <- logCvector[42:55]
       logCvector_small <- logCvector[1:42]
       
       logCvector_largeMicronekton <- logCvector_large
       logCvector_smallMicronekton <- logCvector_small
       
       logmedians5_largeMicronekton <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekton ~ logCvector_largeMicronekt,
                pch = 16, col = "darkred" , lwd = 5.5)
       
       
       # # medians with selection cutoff point  
       # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
       #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
       # 
 
      
       max98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
       lines(   log10(max98_vec_ALLBiovol_NA) ~ logCvector,
                pch = 16, col = "purple" , lwd = 5.5)
       
       
       # Add legend -----------
       
       leg.names <- c( "TSWA", "ALL Areas",  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
       col.taxa  <- c(  "darkorange" , alpha("slateblue4", 0.4),  "forestgreen" , "dodgerblue" ,  "darkmagenta", "darkorange" , "darkred"  ) 
       linetypes <- c(  2,2, 1,1,1,1,1   ) 
       
       legend("topright", leg.names, lty = linetypes,
              lwd = 2.5, 
              col = col.taxa)
      
       
       
        
       ### end of grey plot with coloured medians (TSWA) ###
       
       ###########
       ##########
       
       
       
       # For Paper ------------
       # CCUS
       # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
       # FOR PAPER -  grey plot wit coloured medians 
       ### Plot CCUS with medians and regression lines ------------------- 
       
       # plot CCUS  NBSS Data (loop) --------------
       
       library(scales)
       
       #lnumb <- length(Data.All.Atl$Longitude)
       #lnumb# 2840 data (lines) in C format
       
       lnumb <-3
       
       plot (log10(Data.All.ALL.CCUS.5.NBSSmatrix [lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-7, 13),col = "white",
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "CCUS")
       
       
       for (i in 1:nrow(Data.All.ALL.CCUS.5.NBSSmatrix)){  
         
         points (log10(Data.All.ALL.CCUS.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                 pch = 16, col = alpha("darkgrey", 0.3))
         
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), CCUS -----------------------  
       library (MASS)
       
       
       
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.CCUS.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.CCUS.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.955
       
       rlm1CCUS.All.naive <- rlm1
       
       
       abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("slateblue4", 0.4) )
       
       abline(rlm1CCUS.All.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("darkorange", 0.9) )
       
       
       # grey medians for ALL data (backround grey curved lines) ------
       
       phytopl.All.Atl.NBSSmatrix.1.CCUS <- as.matrix(phytopl.unfiltr.ALL.CCUS.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.CCUS)
       # 110 samples
       
       median_vec_ALLBiovol2.CCUS <- apply(phytopl.All.Atl.NBSSmatrix.1.CCUS, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.CCUS)
       
       median_vec_ALLBiovol2.CCUS <- as.vector(unlist(median_vec_ALLBiovol2.CCUS))
       
       median_vec_PHYTOP.CCUS <-  median_vec_ALLBiovol2.CCUS 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       
       
       lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = alpha("grey27", 0.5), lwd = 4.5)
       
       
       # Darkgrey median lines as background ---------- 
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col =alpha("grey27", 0.5), lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
                pch = 16, col = alpha("grey27", 0.5) , lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeNET ~ logCvector_largeNET,
                pch = 16, col = alpha("grey27", 0.6) , lwd = 4.7)
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, 
               col = alpha("grey27", 0.55), lwd = 4.6)
       
       
       # median lines for CCUS area ----------
       
       phytopl.All.Atl.NBSSmatrix.1.CCUS <- as.matrix(phytopl.unfiltr.ALL.CCUS.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.CCUS)
       # 110 samples
       
       median_vec_ALLBiovol2.CCUS <- apply(phytopl.All.Atl.NBSSmatrix.1.CCUS, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.CCUS)
       
       median_vec_ALLBiovol2.CCUS <- as.vector(unlist(median_vec_ALLBiovol2.CCUS))
       
       median_vec_PHYTOP.CCUS <-  median_vec_ALLBiovol2.CCUS 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       lines(  log10(median_vec_PHYTOP.CCUS) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = "forestgreen" , lwd = 5.5) # OK
       
       
       
       
       # Net Zoopl. , medians,  CCUS -------
       
       # all medians 
       
             zoo_NET.NBSSmatrix <- as.matrix(zoo_NET.CCUS.5[,39:93])
       names(as.data.frame (zoo_NET.NBSSmatrix))
       dim(zoo_NET.NBSSmatrix) # 83 samples
       rowSums(is.na(zoo_NET.NBSSmatrix))
       
       zoo_NET.NBSSmatrix_no_na <-  zoo_NET.NBSSmatrix[rowSums(is.na(zoo_NET.NBSSmatrix)) != ncol(zoo_NET.NBSSmatrix), ]
       rowSums(is.na(zoo_NET.NBSSmatrix_no_na))
       dim(zoo_NET.NBSSmatrix_no_na)
       # only 5 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_NET.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_NET.NBSSmatrix_no_na[,2])
       dim(zoo_NET.NBSSmatrix_no_na)   
       # only 5 samples with data
       
             # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[30:55]
       logmedians5_small <- logmedians5[1:30]
       logCvector_large <- logCvector[30:55]
       logCvector_small <- logCvector[1:30]
       
       logCvector_largezoo_NET <- logCvector_large
       logCvector_smallzoo_NET <- logCvector_small
       
       logmedians5_largezoo_NET <- logmedians5_large
       
       
       lines(   logmedians5_largezoo_NET ~ logCvector_largeNET,
                pch = 16, col = "dodgerblue" , lwd = 5.5)
       
       
       
       # UVP  Zoopl., medians,  CCUS -------
       
       zoo_UVP.NBSSmatrix <- as.matrix(zoo_UVP.CCUS.5[,39:93])
       names(as.data.frame (zoo_UVP.NBSSmatrix))
       dim(zoo_UVP.NBSSmatrix) # 34 samples
       rowSums(is.na(zoo_UVP.NBSSmatrix))
       
       zoo_UVP.NBSSmatrix_no_na <-  zoo_UVP.NBSSmatrix[rowSums(is.na(zoo_UVP.NBSSmatrix)) != ncol(zoo_UVP.NBSSmatrix), ]
       rowSums(is.na(zoo_UVP.NBSSmatrix_no_na))
       dim(zoo_UVP.NBSSmatrix_no_na)
       # 34 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_UVP.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_UVP.NBSSmatrix_no_na[,2])
       dim(zoo_UVP.NBSSmatrix_no_na)   
       # 40 samples with data
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       
       
       median_vec_UVPcarb_NA <- median_vec_ALLBiovol_NA
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, col = "darkmagenta" , lwd = 5.5)
       
       
       # Mesop. fish, medians,  CCUS -------
       
       meso_fi.NBSSmatrix <- as.matrix(Meso_fiALL.CCUS.5[,37:91])
       dim(meso_fi.NBSSmatrix)
       rowSums(is.na(meso_fi.NBSSmatrix))
       
       meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
       rowSums(is.na(meso_fi.NBSSmatrix_no_na))
       
       median_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, median_NA50) 
       mean_vec_ALLBiovol_NA50 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA50) 
       mean_vec_ALLBiovol_NA10 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA10) 
       
       
       length(meso_fi.NBSSmatrix_no_na[,2])
       dim(meso_fi.NBSSmatrix_no_na)   
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[41]) #  cutoff at -1.2 log10 C ind.-1
       logmedians5_large <- logmedians5[41:55]
       logmedians5_small <- logmedians5[1:41]
       logCvector_large <- logCvector[41:55]
       logCvector_small <- logCvector[1:41]
       
       logCvector_largeMESO_FI <- logCvector_large
       logCvector_smallMESO_FI <- logCvector_small
       
       logmedians5_largeMESO_FI <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col = "darkorange" , lwd = 5.5)
       
       
       # Micronekt, medians,  CCUS -------
       
       
       Micronekton.NBSSmatrix <- as.matrix(Micronekton.ALL.CCUS.5[,37:91])
       dim(Micronekton.NBSSmatrix) # 83 samples
       rowSums(is.na(Micronekton.NBSSmatrix))
       
       Micronekton.NBSSmatrix_no_na <-  Micronekton.NBSSmatrix[rowSums(is.na(Micronekton.NBSSmatrix)) != ncol(Micronekton.NBSSmatrix), ]
       rowSums(is.na(Micronekton.NBSSmatrix_no_na))
       dim(Micronekton.NBSSmatrix_no_na)
       # 40 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       max98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
       
       
       length(Micronekton.NBSSmatrix_no_na[,2])
       dim(Micronekton.NBSSmatrix_no_na)   
       # 40 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[42]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[42:55]
       logmedians5_small <- logmedians5[1:42]
       logCvector_large <- logCvector[42:55]
       logCvector_small <- logCvector[1:42]
       
       logCvector_largeMicronekton <- logCvector_large
       logCvector_smallMicronekton <- logCvector_small
       
       logmedians5_largeMicronekton <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekton ~ logCvector_largeMicronekt,
                pch = 16, col = "darkred" , lwd = 5.5)
       
      
       max98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
       lines(   log10(max98_vec_ALLBiovol_NA) ~ logCvector,
                pch = 16, col = "purple" , lwd = 5.5)
       
        
       # # medians with selection cutoff point  
       # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
       #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
       # 
       
       
       
       # Add legend -----------
       
       leg.names <- c( "CCUS", "ALL Areas",  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
       col.taxa  <- c(  "darkorange" , alpha("slateblue4", 0.4),  "forestgreen" , "dodgerblue" ,  "darkmagenta", "darkorange" , "darkred"  ) 
       linetypes <- c(  2,2, 1,1,1,1,1   ) 
       
       legend("topright", leg.names, lty = linetypes,
              lwd = 2.5, 
              col = col.taxa)
       
       ### end of grey plot with coloured medians (CCUS) ###
       
       ########### 
       ##########
       
       
       # For Paper ------------
       # EQU
       # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
       # FOR PAPER -  grey plot wit coloured medians 
       ### Plot EQU with medians and regression lines ------------------- 
       
       # plot EQU  NBSS Data (loop) --------------
     
       lnumb <-3
       
       plot (log10(Data.All.ALL.EQU.5.NBSSmatrix [lnumb,]) ~ log10(X_vector_gCind),
             ylim = c(-7, 13),col = "white",
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "EQU")
       
       for (i in 1: (nrow(Data.All.ALL.EQU.5.NBSSmatrix)) ){
         
         points (log10(Data.All.ALL.EQU.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                 pch = 16, col = alpha("darkgrey", 0.3))
         
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), EQU -----------------------  
       library (MASS)
       
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.EQU.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.EQU.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.955
       
       rlm1EQU.All.naive <- rlm1
       
       
       abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("slateblue4", 0.4) )
       
       abline(rlm1EQU.All.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("darkorange", 0.9) )
       
       
       # grey medians for ALL data (backround grey curved lines) ------
       
       phytopl.All.Atl.NBSSmatrix.1.EQU <- as.matrix(phytopl.unfiltr.ALL.EQU.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.EQU)
       # 110 samples
       
       median_vec_ALLBiovol2.EQU <- apply(phytopl.All.Atl.NBSSmatrix.1.EQU, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.EQU)
       
       median_vec_ALLBiovol2.EQU <- as.vector(unlist(median_vec_ALLBiovol2.EQU))
       
       median_vec_PHYTOP.EQU <-  median_vec_ALLBiovol2.EQU 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       
       
       lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = alpha("grey27", 0.5), lwd = 4.5)
       
       
       # Darkgrey median lines as background ---------- 
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col =alpha("grey27", 0.5), lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
                pch = 16, col = alpha("grey27", 0.5) , lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeNET ~ logCvector_largeNET,
                pch = 16, col = alpha("grey27", 0.6) , lwd = 4.7)
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, 
               col = alpha("grey27", 0.55), lwd = 4.6)
       
       
       # median lines for EQU area ----------
       
       phytopl.All.Atl.NBSSmatrix.1.EQU <- as.matrix(phytopl.unfiltr.ALL.EQU.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.EQU)
       # 110 samples
       
       median_vec_ALLBiovol2.EQU <- apply(phytopl.All.Atl.NBSSmatrix.1.EQU, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.EQU)
       
       median_vec_ALLBiovol2.EQU <- as.vector(unlist(median_vec_ALLBiovol2.EQU))
       
       median_vec_PHYTOP.EQU <-  median_vec_ALLBiovol2.EQU 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       lines(  log10(median_vec_PHYTOP.EQU) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = "forestgreen" , lwd = 5.5) # OK
       
       
       
       
       # Net Zoopl. , medians,  EQU -------
       
       # all medians 
       
     
       zoo_NET.NBSSmatrix <- as.matrix(zoo_NET.EQU.5[,39:93])
       names(as.data.frame (zoo_NET.NBSSmatrix))
       dim(zoo_NET.NBSSmatrix) # 83 samples
       rowSums(is.na(zoo_NET.NBSSmatrix))
       
       zoo_NET.NBSSmatrix_no_na <-  zoo_NET.NBSSmatrix[rowSums(is.na(zoo_NET.NBSSmatrix)) != ncol(zoo_NET.NBSSmatrix), ]
       rowSums(is.na(zoo_NET.NBSSmatrix_no_na))
       dim(zoo_NET.NBSSmatrix_no_na)
       # only 5 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_NET.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_NET.NBSSmatrix_no_na[,2])
       dim(zoo_NET.NBSSmatrix_no_na)   
       # only 5 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[30:55]
       logmedians5_small <- logmedians5[1:30]
       logCvector_large <- logCvector[30:55]
       logCvector_small <- logCvector[1:30]
       
       logCvector_largezoo_NET <- logCvector_large
       logCvector_smallzoo_NET <- logCvector_small
       
       logmedians5_largezoo_NET <- logmedians5_large
       
       
       lines(   logmedians5_largezoo_NET ~ logCvector_largeNET,
                pch = 16, col = "dodgerblue" , lwd = 5.5)
       
       
       
       # UVP  Zoopl., medians,  EQU -------
       
       zoo_UVP.NBSSmatrix <- as.matrix(zoo_UVP.EQU.5[,39:93])
       names(as.data.frame (zoo_UVP.NBSSmatrix))
       dim(zoo_UVP.NBSSmatrix) # 34 samples
       rowSums(is.na(zoo_UVP.NBSSmatrix))
       
       zoo_UVP.NBSSmatrix_no_na <-  zoo_UVP.NBSSmatrix[rowSums(is.na(zoo_UVP.NBSSmatrix)) != ncol(zoo_UVP.NBSSmatrix), ]
       rowSums(is.na(zoo_UVP.NBSSmatrix_no_na))
       dim(zoo_UVP.NBSSmatrix_no_na)
       # 34 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_UVP.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_UVP.NBSSmatrix_no_na[,2])
       dim(zoo_UVP.NBSSmatrix_no_na)   
       # 40 samples with data
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       
       
       median_vec_UVPcarb_NA <- median_vec_ALLBiovol_NA
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, col = "darkmagenta" , lwd = 5.5)
       
       
       # Mesop. fish, medians,  EQU -------
       
       
       
       meso_fi.NBSSmatrix <- as.matrix(Meso_fiALL.EQU.5[,37:91])
       dim(meso_fi.NBSSmatrix)
       rowSums(is.na(meso_fi.NBSSmatrix))
       
       meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
       rowSums(is.na(meso_fi.NBSSmatrix_no_na))
       
       median_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, median_NA50) 
       mean_vec_ALLBiovol_NA50 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA50) 
       mean_vec_ALLBiovol_NA10 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA10) 
       
       
       length(meso_fi.NBSSmatrix_no_na[,2])
       dim(meso_fi.NBSSmatrix_no_na)   
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[41]) #  cutoff at -1.2 log10 C ind.-1
       logmedians5_large <- logmedians5[41:55]
       logmedians5_small <- logmedians5[1:41]
       logCvector_large <- logCvector[41:55]
       logCvector_small <- logCvector[1:41]
       
       logCvector_largeMESO_FI <- logCvector_large
       logCvector_smallMESO_FI <- logCvector_small
       
       logmedians5_largeMESO_FI <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col = "darkorange" , lwd = 5.5)
       
       
       # Micronekt, medians,  EQU -------
       
       
       Micronekton.NBSSmatrix <- as.matrix(Micronekton.ALL.EQU.5[,37:91])
       dim(Micronekton.NBSSmatrix) # 83 samples
       rowSums(is.na(Micronekton.NBSSmatrix))
       
       Micronekton.NBSSmatrix_no_na <-  Micronekton.NBSSmatrix[rowSums(is.na(Micronekton.NBSSmatrix)) != ncol(Micronekton.NBSSmatrix), ]
       rowSums(is.na(Micronekton.NBSSmatrix_no_na))
       dim(Micronekton.NBSSmatrix_no_na)
       # 40 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(Micronekton.NBSSmatrix_no_na[,2])
       dim(Micronekton.NBSSmatrix_no_na)   
       # 40 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[42]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[42:55]
       logmedians5_small <- logmedians5[1:42]
       logCvector_large <- logCvector[42:55]
       logCvector_small <- logCvector[1:42]
       
       logCvector_largeMicronekton <- logCvector_large
       logCvector_smallMicronekton <- logCvector_small
       
       logmedians5_largeMicronekton <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekton ~ logCvector_largeMicronekt,
                pch = 16, col = "darkred" , lwd = 5.5)
       
       max98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
       lines(   log10(max98_vec_ALLBiovol_NA) ~ logCvector,
                pch = 16, col = "purple" , lwd = 5.5)
       
       
       
       # # medians with selection cutoff point  
       # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
       #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
       # 
       
       
       
       # Add legend -----------
       
       leg.names <- c( "EQU", "ALL Areas",  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
       col.taxa  <- c(  "darkorange" , alpha("slateblue4", 0.4),  "forestgreen" , "dodgerblue" ,  "darkmagenta", "darkorange" , "darkred"  ) 
       linetypes <- c(  2,2, 1,1,1,1,1   ) 
       
       legend("topright", leg.names, lty = linetypes,
              lwd = 2.5, 
              col = col.taxa)
       
       ### end of grey plot with coloured medians (EQU) ###
       
      
       ########### 
       ##########
       
       
       # For Paper ------------
       # BUS
       # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
       # FOR PAPER -  grey plot wit coloured medians 
       ### Plot BUS with medians and regression lines ------------------- 
       
       # plot BUS  NBSS Data (loop) --------------
       
      
       plot (log10(Data.All.ALL.BUS.5.NBSSmatrix [3,]) ~ log10(X_vector_gCind),
             ylim = c(-7, 13),col = "white",
             xlim = c(-14, 3),
             xlab = "log10(gC ind.-1)",
             ylab = "log10(gC m-3 / gC ind.-1)",
             main= "BUS")
       
              for (i in 1:nrow(Data.All.ALL.BUS.5.NBSSmatrix)){  
         
         points (log10(Data.All.ALL.BUS.5.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                 #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                 pch = 16, col = alpha("darkgrey", 0.3))
         
       }
       
       # Naive linear model (all data): ------------------
       # matrix to to vector (transpose), BUS -----------------------  
       library (MASS)
       
       
       
       x1 <- rep( logCvector, nrow(Data.All.ALL.BUS.5.NBSSmatrix) )
       y1 <- c ( t(Data.All.ALL.BUS.5.NBSSmatrix))
       
       df1 <- data.frame(x1, y1, ylog = log10(y1))
       names(df1)
       
       # -INF to NA!!!
       
       # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
       
       df1[df1=="-Inf"]<-NA
       
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       
       summary(df1)
       summary(df2)
       
       mean(df2$y1)
       
       summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.954
       summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
       # slope = -0.955
       
       rlm1BUS.All.naive <- rlm1
       
       
       abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("slateblue4", 0.4) )
       
       abline(rlm1BUS.All.naive, xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = alpha("darkorange", 0.9) )
       
       
       # grey medians for ALL data (backround grey curved lines) ------
       
       phytopl.All.Atl.NBSSmatrix.1.BUS <- as.matrix(phytopl.unfiltr.ALL.BUS.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.BUS)
       # 110 samples
       
       median_vec_ALLBiovol2.BUS <- apply(phytopl.All.Atl.NBSSmatrix.1.BUS, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.BUS)
       
       median_vec_ALLBiovol2.BUS <- as.vector(unlist(median_vec_ALLBiovol2.BUS))
       
       median_vec_PHYTOP.BUS <-  median_vec_ALLBiovol2.BUS 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       
       
       lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = alpha("grey27", 0.5), lwd = 4.5)
       
       
       # Darkgrey median lines as background ---------- 
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col =alpha("grey27", 0.5), lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
                pch = 16, col = alpha("grey27", 0.5) , lwd = 4.5)
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeNET ~ logCvector_largeNET,
                pch = 16, col = alpha("grey27", 0.6) , lwd = 4.7)
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, 
               col = alpha("grey27", 0.55), lwd = 4.6)
       
       
       # median lines for BUS area ----------
       
       phytopl.All.Atl.NBSSmatrix.1.BUS <- as.matrix(phytopl.unfiltr.ALL.BUS.5[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1.BUS)
       # 110 samples
       
       median_vec_ALLBiovol2.BUS <- apply(phytopl.All.Atl.NBSSmatrix.1.BUS, 2, median_NA50) 
       
       class(median_vec_ALLBiovol2.BUS)
       
       median_vec_ALLBiovol2.BUS <- as.vector(unlist(median_vec_ALLBiovol2.BUS))
       
       median_vec_PHYTOP.BUS <-  median_vec_ALLBiovol2.BUS 
       
       logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
       
       lines(  log10(median_vec_PHYTOP.BUS) ~ logX_vector_gCind_PHYTOPL,
               pch = 16, col = "forestgreen" , lwd = 5.5) # OK
       
       
       
       
       # Net Zoopl. , medians,  BUS -------
       
       # all medians 
       
       
       zoo_NET.NBSSmatrix <- as.matrix(zoo_NET.BUS.5[,39:93])
       names(as.data.frame (zoo_NET.NBSSmatrix))
       dim(zoo_NET.NBSSmatrix) # only 6 samples
       rowSums(is.na(zoo_NET.NBSSmatrix))
       
       zoo_NET.NBSSmatrix_no_na <-  zoo_NET.NBSSmatrix[rowSums(is.na(zoo_NET.NBSSmatrix)) != ncol(zoo_NET.NBSSmatrix), ]
       rowSums(is.na(zoo_NET.NBSSmatrix_no_na))
       dim(zoo_NET.NBSSmatrix_no_na)
       # only 5 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_NET.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_NET.NBSSmatrix_no_na[,2])
       dim(zoo_NET.NBSSmatrix_no_na)   
       # only 5 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[30:55]
       logmedians5_small <- logmedians5[1:30]
       logCvector_large <- logCvector[30:55]
       logCvector_small <- logCvector[1:30]
       
       logCvector_largezoo_NET <- logCvector_large
       logCvector_smallzoo_NET <- logCvector_small
       
       logmedians5_largezoo_NET <- logmedians5_large
       
       
       lines(   logmedians5_largezoo_NET ~ logCvector_largeNET,
                pch = 16, col = "dodgerblue" , lwd = 5.5)
       
       
       
       # UVP  Zoopl., medians,  BUS -------
       
       zoo_UVP.NBSSmatrix <- as.matrix(zoo_UVP.BUS.5[,39:93])
       names(as.data.frame (zoo_UVP.NBSSmatrix))
       dim(zoo_UVP.NBSSmatrix) # 34 samples
       rowSums(is.na(zoo_UVP.NBSSmatrix))
       
       zoo_UVP.NBSSmatrix_no_na <-  zoo_UVP.NBSSmatrix[rowSums(is.na(zoo_UVP.NBSSmatrix)) != ncol(zoo_UVP.NBSSmatrix), ]
       rowSums(is.na(zoo_UVP.NBSSmatrix_no_na))
       dim(zoo_UVP.NBSSmatrix_no_na)
       # 34 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( zoo_UVP.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(zoo_UVP.NBSSmatrix_no_na[,2])
       dim(zoo_UVP.NBSSmatrix_no_na)   
       # 40 samples with data
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       
       
       median_vec_UVPcarb_NA <- median_vec_ALLBiovol_NA
       
       lines(  log10(median_vec_UVPcarb_NA) ~ logCvector,pch = 16, col = "darkmagenta" , lwd = 5.5)
       
       
       # Mesop. fish, medians,  BUS -------
       
       
       
       meso_fi.NBSSmatrix <- as.matrix(Meso_fiALL.BUS.5[,37:91])
       dim(meso_fi.NBSSmatrix)
       rowSums(is.na(meso_fi.NBSSmatrix))
       
       meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
       rowSums(is.na(meso_fi.NBSSmatrix_no_na))
       
       median_vec_ALLBiovol_NA <- apply( meso_fi.NBSSmatrix_no_na, 2, median_NA50) 
       mean_vec_ALLBiovol_NA50 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA50) 
       mean_vec_ALLBiovol_NA10 <- apply(meso_fi.NBSSmatrix_no_na, 2, mean_NA10) 
       
       
       length(meso_fi.NBSSmatrix_no_na[,2])
       dim(meso_fi.NBSSmatrix_no_na)   
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[41]) #  cutoff at -1.2 log10 C ind.-1
       logmedians5_large <- logmedians5[41:55]
       logmedians5_small <- logmedians5[1:41]
       logCvector_large <- logCvector[41:55]
       logCvector_small <- logCvector[1:41]
       
       logCvector_largeMESO_FI <- logCvector_large
       logCvector_smallMESO_FI <- logCvector_small
       
       logmedians5_largeMESO_FI <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
                pch = 16, col = "darkorange" , lwd = 5.5)
       
       
       # Micronekt, medians,  BUS -------
       
       
       Micronekton.NBSSmatrix <- as.matrix(Micronekton.ALL.BUS.5[,37:91])
       dim(Micronekton.NBSSmatrix) # 83 samples
       rowSums(is.na(Micronekton.NBSSmatrix))
       
       Micronekton.NBSSmatrix_no_na <-  Micronekton.NBSSmatrix[rowSums(is.na(Micronekton.NBSSmatrix)) != ncol(Micronekton.NBSSmatrix), ]
       rowSums(is.na(Micronekton.NBSSmatrix_no_na))
       dim(Micronekton.NBSSmatrix_no_na)
       # 40 samples with real data (no NA's)
       
       median_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, median_NA50) 
       
       
       length(Micronekton.NBSSmatrix_no_na[,2])
       dim(Micronekton.NBSSmatrix_no_na)   
       # 40 samples with data
       
       
       # medians (curved line) ------
       logmedians5 <- log10(median_vec_ALLBiovol_NA)
       length(logmedians5) # 55 data points with medians
       (cut_point_Xval <- logCvector[42]) #  cutoff at -4.5 log10 C ind.-1
       logmedians5_large <- logmedians5[42:55]
       logmedians5_small <- logmedians5[1:42]
       logCvector_large <- logCvector[42:55]
       logCvector_small <- logCvector[1:42]
       
       logCvector_largeMicronekton <- logCvector_large
       logCvector_smallMicronekton <- logCvector_small
       
       logmedians5_largeMicronekton <- logmedians5_large
       
       
       # all medians (no cutoff selection), but olnly < 50% NA
       lines(   logmedians5_largeMicronekton ~ logCvector_largeMicronekt,
                pch = 16, col = "darkred" , lwd = 5.5)
       
       
       # # medians with selection cutoff point  
       # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
       #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
       # 
       
       
       max98_vec_ALLBiovol_NA <- apply( Micronekton.NBSSmatrix_no_na, 2, fun98_quant) 
       lines(   log10(max98_vec_ALLBiovol_NA) ~ logCvector,
                pch = 16, col = "purple" , lwd = 5.5)
       
       
       # Add legend -----------
       
       leg.names <- c( "BUS", "ALL Areas",  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
       col.taxa  <- c(  "darkorange" , alpha("slateblue4", 0.4),  "forestgreen" , "dodgerblue" ,  "darkmagenta", "darkorange" , "darkred"  ) 
       linetypes <- c(  2,2, 1,1,1,1,1   ) 
       
       legend("topright", leg.names, lty = linetypes,
              lwd = 2.5, 
              col = col.taxa)
       
       ### end of grey plot with coloured medians (BUS) ###
     
       
       # compare biomass UVP vs net-caught zoopl for a specific bin (BUS region) --------
       # abline(v = -4.6) 
       # abline(v = -4.4) 
        abline(v = -4.515700, col = "red") # bin nr. 30
#       abline(v = -4.214670 , col = "purple") 
       
    #    *bookmark
       
       X_vector_gCind
       log10(X_vector_gCind)
       log10(X_vector_gCind)[30]
       
       logmedians5_largezoo_NET
       logmedians5_largezoo_NET[1]# 2.59484 in log scale , 393 not log
       log10(median_vec_UVPcarb_NA)#  0.7997177  in lg scalle, 6.30 not log
       
       10^2.59484; log10(393.4051)# 2.5948
       10^0.7997177; log10(6.30)#  0.799
       # UVP / net zoopl ratio (BUS, smallest data-rich size bin):
       393 / 6.3 # 62.38
       # median net zoopl biomass was 62 times larger than  median UVP biomass (BUS, smallest data-rich size bin)
       
       
  
       
       ########### 
       ##########
       
       #### END OF GREY NBSS PLOTS BY REGIONS, unfiltered ######
       
      
       #####################################
       #####################
       # 7 PLOTS FOR PAPER ------------------
       # Phytoplanktn NBSS ploots by Temp, Chla -------------
       # Cold, Subtr. , Tropical ---------------------
       # Oligotrophic, Mesotrophic, Eutrophic ------------------
       # Composed Lattice plot with NBSS and RR Slope value
       ###########
       
       
       # Divide into 9 groups(7 actually exist with real samples) ----------
       
       
       plot(phytopl.All.Atl$Chlorophyll__mg_m_3_insitu ,
            phytopl.All.Atl$SST_C_insitu,
            #               xlim = c(-2,  1.5),
             ylim = c(3, 30),
            col = alpha ("darkgreen", 0.3), pch = 16)
       #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
       abline(h = c(20, 27.5))
       abline(v = c(0.2, 1))
       
       
       
       plot(log10(1+phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) ,
            phytopl.All.Atl$SST_C_insitu,
            #               xlim = c(-2,  1.5),
            ylim = c(3, 30),
            col = alpha ("darkgreen", 0.3), pch = 16)
       #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
       abline(h = c(20, 27.5))
       abline(v = c(  0.0792,  0.301))
       log10(1+1) # 0.477
       log10(1+0.2) # 0.0792
       
       library(colourvalues)
       phytopl.All.Atl.p.rsq.slope_ratio.ok$z_colours_fromTEMPinsitu   <-  colour_values(phytopl.All.Atl.p.rsq.slope_ratio.ok$Temperature__C_insitu, palette = "viridis")
       phytopl.All.Atl.p.rsq.slope_ratio.ok$z_colours_fromSSTinsitu   <-  colour_values(2+phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu^8, palette = "viridis")
       phytopl.All.Atl.p.rsq.slope_ratio.ok$z_colours_fromrob_reg__slopes   <-  colour_values(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes, palette = "viridis")
       
       
       # FOR PAPER colour points plot Chla SST col = Phtopl NBSS Slope --------------- 
       length(phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu)
       plot(log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) ,
            phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu,
            #               xlim = c(-2,  1.5),
            ylim = c(3, 30),
            main = "NBSS Slopes, n = 458",
            xlab =  "log10(1+Chla in situ (mg m-3))",
            ylab= "SST in situ (°C)",
            col =alpha(  phytopl.All.Atl.p.rsq.slope_ratio.ok$z_colours_fromrob_reg__slopes , 0.7), pch = 16 )
       #abline(rlm1.trop, lwd = 2.5, lty = 2, col = "darkorange")
       abline(h = c(20, 27.5))
       abline(v = c(  0.0792,  0.301))
       
       log10(1+0.2) # 0.0792   # 0.2 mg m-1, OK
       log10(1+1) # 0.301 # log10(1+ ( 1 mg m-3  )  ) # 1 mg m-3, OK
       log10(1+2) # 0.477 (old boundary (2 mg m-3))
       
       
       
       # Insert a Legend
       #?legend
       
       leg.slopes <- c(-2, -1.5, -1, -0.5, 0)
       leg.cols <- color_values(leg.slopes, palette = "viridis")
       
       legend("bottomright", legend = paste (leg.slopes), col =leg.cols,
              text.col = "darkgrey",  pch = 16,
              merge = TRUE, bg = "gray90", trace=TRUE)
       
       
       #### 7 NBSS PLOTS by temp & chla groups (cold, mesotrophic, etc) 
       ### NBSS plots by temperature and chla groups
       
       ### Calculate median with ALL data (not only the ones filtered for linearity)
       
       
       ### COLD waters (< 20 C° SST in situ ) Phyopl. NBSS PLOT ### --------------------
       
       Cold.phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(Cold.phytopl.All.Atl.unfiltr[,37:91])
       Cold.phytopl.NBSS.matrix.1 <- Cold.phytopl.All.Atl.NBSSmatrix.1
       dim(Cold.phytopl.All.Atl.NBSSmatrix.1)
       
       # PLot ALL NBSS Data, Carbon  -------
       
       #nrow(Cold.phytopl.NBSS.matrix.ok3) # 458 high-quality datasets, OK
       
       # plot all points and rlm models (blue lines), and red median line  ----------
       
       lnumb <- 34
       plot (log10(Cold.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind), 
             main = "Cold (SST < 20 ºC) unfiltered data",
             xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
             ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
             xlim = c(-15, -5), ylim = c(0, 14) ,
             pch = 16, col = "white")
       
       for (lnumb in 1 : nrow(Cold.phytopl.NBSS.matrix.1) ) 
       {    
         points (log10(Cold.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
                 pch = 16, col = alpha("forestgreen", 0.2))
       }
       
       (median(Cold.phytopl.All.Atl$rob_reg__slopes) ) # Cold Phytopl. = -1.04 
       length(Cold.phytopl.All.Atl$rob_reg__slopes) # 196 unfiltred COLD phytopl datasets 
       
       # median RR slope : -1.256 (ALL Atlantic)
       # (min(Cold.phytopl.All.Atl$rob_reg__slopes)) 
       # (max(Cold.phytopl.All.Atl$rob_reg__slopes) )
       # # min - max  RR slope : -0.44155
       
       
       
       # matrix to to vector (transpose), TSWA -----------------------  
       x1 <- rep( (log10(X_vector_gCind)), nrow(Cold.phytopl.NBSS.matrix.1) )
       y1 <- c ( t(Cold.phytopl.NBSS.matrix.1))
       
       # Regression based on ALL DATA ------------- 
       # regression line (Rob. regr.)
       #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
       #points(  log10(y1) ~ x1 )
       # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
       abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
                lwd = 2.5, lty = 2, col = "darkorange" )
       summary(rlm2)
       # slope = -1.146 , using all green points, rlm
       # means and medians without considering NAs  
       mean (c(4,NA, 18,34)) # NA
       median( c(4,NA, 18, 34)  ) #NA
       
       # apply median by columns
       
       median_vec_ALLBiovol <- apply(Cold.phytopl.NBSS.matrix.1, 2, median_) 
       mean_vec_ALLBiovol <- apply(Cold.phytopl.NBSS.matrix.1, 2, mean_) 
       
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)
       
       
       # apply mean by columns (only if less thn 50% NAs)
       
       # mean_vec_ALLBiovol2 <- apply(Cold.phytopl.NBSS.matrix.ok3, 2, mean_NA50) 
       
       median_vec_ALLBiovol2 <- apply(Cold.phytopl.NBSS.matrix.1, 2, median_NA50) 
       
       dim(median_vec_ALLBiovol2)
       
       class(median_vec_ALLBiovol2)
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
       
       #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
       
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
       
       dfok<- df1
       #dfok <- df1[complete.cases(df1$y1),]
       #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               col = "red" , lwd = 4.5)
       
       
       # rlm based on medians only
       
       abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
               lwd = 2.5, lty = 2, col = "darkgreen")
       
       summary(rlm_med)
       # slope: -1.1573 (based on medians < 50% NAs )
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               
               
               col = "red" , lwd = 4.5)
       
       # col = "darkmagenta" , lwd = 2.5)
       
       ### END OF COLD PHYTOPL NBSS PLOT ### ------------
       
       
       ### START OF SubTrop.phytOPICAL PHYTOPL NBSS PLOT ### ------------
       ### SubTrop waters ( 20 to 27.5 C°) Phyopl. NBSS PLOT ### --------------------
       
       SubTrop.phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(SubTrop.phytopl.All.Atl[,37:91])
       SubTrop.phytopl.NBSS.matrix.1 <- SubTrop.phytopl.All.Atl.NBSSmatrix.1
       dim(SubTrop.phytopl.All.Atl.NBSSmatrix.1)
       
       # PLot ALL NBSS Data, Carbon  -------
       
       #nrow(SubTrop.phytopl.NBSS.matrix.ok3) # 458 high-quality datasets, OK
       
       # plot all points and rlm models (blue lines), and red median line  ----------
       
       lnumb <- 34
       plot (log10(SubTrop.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind), 
             main = "Subtropical (SST 20 to 27.5 ºC), filtered data",
             xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
             ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
             xlim = c(-15, -5), ylim = c(0, 14) ,
             pch = 16, col = "white")
       
       for (lnumb in 1 : nrow(SubTrop.phytopl.NBSS.matrix.1) ) 
       {    
         points (log10(SubTrop.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
                 pch = 16, col = alpha("forestgreen", 0.2))
       }
       
       (median(SubTrop.phytopl.All.Atl$rob_reg__slopes) ) # SubTrop Phytopl. = -1.04 
       length(SubTrop.phytopl.All.Atl$rob_reg__slopes) # 196 unfiltred SubTrop phytopl datasets 
       
       # median RR slope : -1.256 (ALL Atlantic)
       # (min(SubTrop.phytopl.All.Atl$rob_reg__slopes)) 
       # (max(SubTrop.phytopl.All.Atl$rob_reg__slopes) )
       # # min - max  RR slope : -0.44155
       
       
       
       # matrix to to vector (transpose), TSWA -----------------------  
       x1 <- rep( (log10(X_vector_gCind)), nrow(SubTrop.phytopl.NBSS.matrix.1) )
       y1 <- c ( t(SubTrop.phytopl.NBSS.matrix.1))
       
       # regression line (Rob. regr.)
       #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
       #points(  log10(y1) ~ x1 )
       # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
       abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
                lwd = 2.5, lty = 2, col = "darkorange" )
       summary(rlm2)
       # slope = -1.146 , using all green points, rlm
       # means and medians without considering NAs  
       mean (c(4,NA, 18,34)) # NA
       median( c(4,NA, 18, 34)  ) #NA
       
       # apply median by columns
       
       median_vec_ALLBiovol <- apply(SubTrop.phytopl.NBSS.matrix.1, 2, median_) 
       mean_vec_ALLBiovol <- apply(SubTrop.phytopl.NBSS.matrix.1, 2, mean_) 
       
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)
       
       
       # apply mean by columns (only if less thn 50% NAs)
       
       # mean_vec_ALLBiovol2 <- apply(SubTrop.phytopl.NBSS.matrix.ok3, 2, mean_NA50) 
       
       median_vec_ALLBiovol2 <- apply(SubTrop.phytopl.NBSS.matrix.1, 2, median_NA50) 
       
       dim(median_vec_ALLBiovol2)
       
       class(median_vec_ALLBiovol2)
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
       
       #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
       
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
       
       dfok<- df1
       #dfok <- df1[complete.cases(df1$y1),]
       #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               col = "red" , lwd = 4.5)
       
       
       # rlm based on medians only
       
       abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
               lwd = 2.5, lty = 2, col = "darkgreen")
       
       summary(rlm_med)
       # slope: -1.1573 (based on medians < 50% NAs )
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               
               
               col = "red" , lwd = 4.5)
       
       # col = "darkmagenta" , lwd = 2.5)
       
       ### END OF SubTrop.phytOPICAL PHYTOPL NBSS PLOT ### ------------
       
       
       
       ### START OF Tropic.phytOPICAL PHYTOPL NBSS PLOT ### ------------
       ### Tropic waters ( 20 to 27.5 C°) Phyopl. NBSS PLOT ### --------------------
       
       Tropic.phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(Tropic.phytopl.All.Atl.unfiltr[,37:91])
       Tropic.phytopl.NBSS.matrix.1 <- Tropic.phytopl.All.Atl.NBSSmatrix.1
       dim(Tropic.phytopl.All.Atl.NBSSmatrix.1)
       
       # PLot ALL NBSS Data, Carbon  -------
       
       #nrow(Tropic.phytopl.NBSS.matrix.ok3) # 458 high-quality datasets, OK
       
       # plot all points and rlm models (blue lines), and red median line  ----------
       
       lnumb <- 34
       plot (log10(Tropic.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind), 
             main = "Tropical (SST > 27.5 ºC), unfiltered data",
             xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
             ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
             xlim = c(-15, -5), ylim = c(0, 14) ,
             pch = 16, col = "white")
       
       for (lnumb in 1 : nrow(Tropic.phytopl.NBSS.matrix.1) ) 
       {    
         points (log10(Tropic.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
                 pch = 16, col = alpha("forestgreen", 0.2))
       }
       
       (median(Tropic.phytopl.All.Atl$rob_reg__slopes) ) # Tropic Phytopl. = -1.04 
       length(Tropic.phytopl.All.Atl$rob_reg__slopes) # 196 unfiltred Tropic phytopl datasets 
       
       # median RR slope : -1.256 (ALL Atlantic)
       # (min(Tropic.phytopl.All.Atl$rob_reg__slopes)) 
       # (max(Tropic.phytopl.All.Atl$rob_reg__slopes) )
       # # min - max  RR slope : -0.44155
       
       
       
       # matrix to to vector (transpose), TSWA -----------------------  
       x1 <- rep( (log10(X_vector_gCind)), nrow(Tropic.phytopl.NBSS.matrix.1) )
       y1 <- c ( t(Tropic.phytopl.NBSS.matrix.1))
       
       # regression line (Rob. regr.)
       #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
       #points(  log10(y1) ~ x1 )
       # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
       abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
                lwd = 2.5, lty = 2, col = "darkorange" )
       summary(rlm2)
       # slope = -1.146 , using all green points, rlm
       # means and medians without considering NAs  
       mean (c(4,NA, 18,34)) # NA
       median( c(4,NA, 18, 34)  ) #NA
       
       # apply median by columns
       
       median_vec_ALLBiovol <- apply(Tropic.phytopl.NBSS.matrix.1, 2, median_) 
       mean_vec_ALLBiovol <- apply(Tropic.phytopl.NBSS.matrix.1, 2, mean_) 
       
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)
       
       
       # apply mean by columns (only if less thn 50% NAs)
       
       # mean_vec_ALLBiovol2 <- apply(Tropic.phytopl.NBSS.matrix.ok3, 2, mean_NA50) 
       
       median_vec_ALLBiovol2 <- apply(Tropic.phytopl.NBSS.matrix.1, 2, median_NA50) 
       
       dim(median_vec_ALLBiovol2)
       
       class(median_vec_ALLBiovol2)
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
       
       #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
       
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
       
       dfok<- df1
       #dfok <- df1[complete.cases(df1$y1),]
       #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               col = "red" , lwd = 4.5)
       
       
       # rlm based on medians only
       
       abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
               lwd = 2.5, lty = 2, col = "darkgreen")
       
       summary(rlm_med)
       # slope: -1.1573 (based on medians < 50% NAs )
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               
               
               col = "red" , lwd = 4.5)
       
       # col = "darkmagenta" , lwd = 2.5)
       
       ### END OF Tropic.phytOPICAL PHYTOPL NBSS PLOT ### ------------
       
       
       #### By Chla. classes -------------------
       
       #############
       # three Chla groups FOR PAPER -----------------------
       # classified by  Chla  in situ !
       
       # Separate OLigo-, meso-, and  eutrophic samples --------
       # boundaries: 0.2 and 2 mg Chla m-3 in situ
       
       # find boundaries -------------------
       length(phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 727 datasets
       # 458 total datasets
       
       
       oligotr.phytopl.All.Atl <- subset(phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 0.2 )
       summary(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
       length(oligotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) 
       #  n = 97 all Datasets,  (88 datasets onlt top linear)
       
       
       mesotr.phytopl.All.Atl <- subset(phytopl.All.Atl, Chlorophyll__mg_m_3_insitu > 0.2 )
       mesotr.phytopl.All.Atl <- subset(mesotr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 1 )#>2 mg m3 Chla 
       summary(mesotr.phytopl.All.Atl$SST_C_insitu) # 
       length(mesotr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       #  n = 299, all Datasets
       
       
       eutr.phytopl.All.Atl <-subset(phytopl.All.Atl, Chlorophyll__mg_m_3_insitu > 1 )
       eutr.phytopl.All.Atl <-subset(eutr.phytopl.All.Atl, Chlorophyll__mg_m_3_insitu < 30 )
       summary(eutr.phytopl.All.Atl$SST_C_insitu) # < 31 mg 
       length(eutr.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) 
       # 155, All datasets
       
      
       #########################
       # 7 groups, 9 groups
       # CREATE 7 GROUPS, (first create 9 groups, but 7 actually have enough data) ---------------   
       # UNFILTERED Data (linear & non-linear shaped NBSS)
       
       
       # Cold ----------
       
       Cold.oligotroph.phytopl.All.Atl <- subset(Cold.phytopl.All.Atl.unfiltr , Chlorophyll__mg_m_3_insitu < 0.2 )
       summary(Cold.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
       length(Cold.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) 
       # 12 datasets top-linear, 15 datasets ALL
       
       Cold.mesotroph.phytopl.All.Atl <- subset(Cold.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 0.2 & Chlorophyll__mg_m_3_insitu < 1)
       summary(Cold.mesotroph.phytopl.All.Atl$SST_C_insitu) # 
       length(Cold.mesotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # 83 datasets top-linear, 117 datasets ALL
       
       Cold.eutroph.phytopl.All.Atl <- subset(Cold.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 1 )
       summary(Cold.eutroph.phytopl.All.Atl$SST_C_insitu) # 
       length(Cold.eutroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # 99 datasets top-linear, 148 datasets ALL
       
       # Subtropical------------------
        
       SubTrop.oligotroph.phytopl.All.Atl <- subset(SubTrop.phytopl.All.Atl.unfiltr , Chlorophyll__mg_m_3_insitu < 0.2 )
       summary(SubTrop.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
       length(SubTrop.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) 
       # 48 datasets, 53 datasets ALL
       
       SubTrop.mesotroph.phytopl.All.Atl <- subset(SubTrop.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 0.2 & Chlorophyll__mg_m_3_insitu < 1)
       summary(SubTrop.mesotroph.phytopl.All.Atl$SST_C_insitu) # 
       summary(SubTrop.mesotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       length(SubTrop.mesotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # 107 datasets top-linear, 107 datasets ALL
       
       SubTrop.eutroph.phytopl.All.Atl <- subset(SubTrop.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 1 )
       summary(SubTrop.eutroph.phytopl.All.Atl$SST_C_insitu) # 
       length(SubTrop.eutroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # 5  datasets top-linear, 6 datasets ALL (not useful!)
       
       # Tropical (> 27.5 °C)----------- 
       
       Tropic.oligotroph.phytopl.All.Atl <- subset(Tropic.phytopl.All.Atl.unfiltr , Chlorophyll__mg_m_3_insitu < 0.2 )
       summary(Tropic.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 0 to 0.19 mg m3 Chla 
       length(Tropic.oligotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) 
       # 27 datasets, 28 datasets ALL
       
       Tropic.mesotroph.phytopl.All.Atl <- subset(Tropic.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 0.2 & Chlorophyll__mg_m_3_insitu < 1)
       summary(Tropic.mesotroph.phytopl.All.Atl$SST_C_insitu) # 
       summary(Tropic.mesotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       length(Tropic.mesotroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # 61 datasets top-linear, 63 datasets ALL
       
       # Tropic.eutroph.phytopl.All.Atl <- subset(Tropic.phytopl.All.Atl.unfiltr, Chlorophyll__mg_m_3_insitu > 1 )
       # summary(Tropic.eutroph.phytopl.All.Atl$SST_C_insitu) # 
       # length(Tropic.eutroph.phytopl.All.Atl$Chlorophyll__mg_m_3_insitu) # 
       # # 0  datasets top-linear, 0 datasets ALL (not useful!)
       
       ###########
       # 7 Plots by temperatore and Chl a groups  ------------
       
      
       # plot and matrix as function --------------
       
      
       
       ### START OF function  funplot.matrixlinesmedians ###
       
         funplot.matrixlinesmedians <- function (Y, Title) {
       
       Tropic.phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(Y[,37:91])
       Tropic.phytopl.NBSS.matrix.1 <- Tropic.phytopl.All.Atl.NBSSmatrix.1
       # plot all points and rlm models (blue lines), and red median line  ----------
       
       lnumb <- 3
       plot (log10(Tropic.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind), 
             main = paste(Title),
             xlab = "log10(gC ind. -1) ",
             ylab = "log10(gC m-3 gC-1 ind.) ",
             xlim = c(-15, -5), ylim = c(0, 14) ,
             pch = 16, col = "white")
       
       for (lnumb in 1 : nrow(Tropic.phytopl.NBSS.matrix.1) ) 
       {    
         points (log10(Tropic.phytopl.NBSS.matrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
                 pch = 16, col = alpha("forestgreen", 0.2))
       }
       
       x1 <- rep( (log10(X_vector_gCind)), nrow(Tropic.phytopl.NBSS.matrix.1) )
       y1 <- c ( t(Tropic.phytopl.NBSS.matrix.1))
       
       
       
       # Naive global model, All Atlantic, All Taxa ALL data
       
       #Slope 
       # -0.9548    +- 0.0026
       # Interc
       # -3.2915   +-  0.0191 
       
       
       abline(a = -3.2915,  b= -0.9548  , xlim = c(-11, -1), ylim = c(0, 14), 
              lwd = 2.5, lty = 2, col = "grey" )
       
       
       
       rlm987 <- rlm(log10(y1) ~ x1)
        abline(rlm987 ,  
                 lwd = 2.5, lty = 2, col = "darkorange" )
    # 
       this.slope <- as.numeric (coef(rlm987)[2])
    
       text( -14.2, 5.5, paste(round (this.slope, 3)),col = "darkorange"   )
       
          
       
       # Naive global model, All Atlantic, All Taxa ALL data
       
       #Slope 
       # -0.9548    +- 0.0026
       # Interc
       # -3.2915   +-  0.0191 
       
       text( -5.6, 5.5, paste(round (-0.9548, 3))  , col = "grey")
       
       
       # apply median by columns
       
       median_vec_ALLBiovol2 <- apply(Tropic.phytopl.NBSS.matrix.1, 2, median_NA50) 
       
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
       
       
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
       
       dfok<- df1
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               col = "red" , lwd = 4.5)
       
       
       # rlm based on medians only
       
       # abline( rlm_med <-rlm( log10(dfok$y1) ~ log10(dfok$x1)  ),
       #         lwd = 2.5, lty = 2, col = "darkgreen")
       # 
       # as.numeric (coef(rlm_med)[2])
       # 
      
       
       #summary(rlm_med)
       # slope: -1.9671 (based on medians < 50% NAs )
       
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
      
               
               col = "red" , lwd = 4.5)
       
      
      legend("topright", c("medians (>50% NA's)", "This subset", "All Taxa, All Regions"),
               lty = c(1, 2, 2),
               lwd = c(4.5, 2.5, 2.5),
               col = c( "red","darkorange",  "darkgrey" ) )
       
        
       ### END OF FUNCTION funplot.matrixlinesmedians ###
       }
       
       funplot.matrixlinesmedians(Tropic.phytopl.All.Atl.unfiltr, "Tropical, unfiltr")
       
       # plotdatalist <- list(       Tropic.oligotroph.phytopl.All.Atl,
       # Tropic.mesotroph.phytopl.All.Atl,
       # Tropic.eutroph.phytopl.All.Atl,
       #        SubTrop.oligotroph.phytopl.All.Atl,
       # SubTrop.mesotroph.phytopl.All.Atl,
       # SubTrop.eutroph.phytopl.All.Atl,
       # Cold.oligotroph.phytopl.All.Atl,
       # Cold.mesotroph.phytopl.All.Atl,
       # Cold.eutroph.phytopl.All.Atl)
       
      # funplot.matrixlinesmedians(plotdatalist[[1]], ", unfiltr")
       
      
       # funplot.matrixlinesmedians(Tropic.eutroph.phytopl.All.Atl, # 0 datasets!
       #                            "Tropic.eutroph.phytopl.All.Atl, unfiltr") # 0 datasets!
      
       funplot.matrixlinesmedians(Tropic.mesotroph.phytopl.All.Atl,
                                  "Tropic.mesotroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(Tropic.oligotroph.phytopl.All.Atl,
                                  "Tropic.oligotroph.phytopl.All.Atl, unfiltr")
        
       funplot.matrixlinesmedians(SubTrop.eutroph.phytopl.All.Atl,
                                  "SubTrop.eutroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(SubTrop.mesotroph.phytopl.All.Atl,
                                  "SubTrop.mesotroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(SubTrop.oligotroph.phytopl.All.Atl,
                                  "SubTrop.oligotroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(Cold.eutroph.phytopl.All.Atl,
                                  "Cold.eutroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(Cold.mesotroph.phytopl.All.Atl,
                                  "Cold.mesotroph.phytopl.All.Atl, unfiltr")
       
       funplot.matrixlinesmedians(Cold.oligotroph.phytopl.All.Atl,
                                  "Cold.oligotroph.phytopl.All.Atl, unfiltr")
       
       # Naive global model, All Atlantic, All Taxa ALL data
       
       #Slope 
       # -0.9548    +- 0.0026
       # Interc
       # -3.2915   +-  0.0191 
     
       
        
log10(X_vector_gCind)      

abline( v = -12.2)

abline( v = -10.6)

# -10.60206 log10 g C ind. for Emiliania huxleyi

       ################
       # NBSS Slopes by latitude --------------------------
       
       # Phytopl. -----------
      
       attach (phytopl.All.Atl.p.rsq.slope_ratio.ok)
       mod1 <- (rob_reg__slopes ~ abs(Latitude) )
       plot(mod1, data = phytopl.All.Atl.p.rsq.slope_ratio.ok,
             xlim = c(0 ,50),
            main = "Linear-shaped phytopl. NBSS, n = 458", 
            # ylim = c(-2.2, 0.1),
            col = alpha ("darkgreen", 0.3), pch = 16)

       abline(rlm1 <- rlm (mod1) , lwd = 2.5, lty = 2, col = "darkorange")
       dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
       median (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope)
       
       text ( 40,  -3.5 ,"Median slope: -1.26" )
       
       # lo <- loess(rob_reg__slopes ~ abs(Latitude) )
       # lines(predict(lo), col='red', lwd=2)
       smoothingSpline = smooth.spline( abs(Latitude) ,rob_reg__slopes,  spar=0.7)
       lines(smoothingSpline,  lwd = 2.5, lty = 1, col = alpha("navy", 0.7))
       
       symbols(  mod1, circles = 1.5^phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu )
       
       
       ### FOR PAPER
       # NBSS sloped vs SST, ALL phytopl data used for regression, rlm --------------
       # lm and RR ar highçly siognifocant, wether ALL data or only clean "p", R² OK data 
       
       # linear model is highly Significant (p < e-16, R²: 0.29,n = 552)
       # for ALL Phttopl. Data as well! 
       
       
       mod1 <- (phytopl.All.Atl$new_slopes ~ phytopl.All.Atl$SST_C_insitu)
       
       Ylab.NBBS <- "log10(gC ind.-1)"
       XlabNBSS  <-  "log10(gC m-3 / gC ind.-1)"
       
       
       getwd()
       
       #svg("phytopl_All_Atl_RRslopessSSTv2.svg" , width = 5.6, height =  4)
       
                plot(mod1, data = phytopl.All.Atl,
                     main = "Phytoplankton, ATLANTIC, All Data, n = 552",
                     xlab = XlabNBSS,
                     ylab = Ylab.NBBS, bty = "n",
                col = alpha("darkgreen", 0.3), pch = 16)
      text(28.5, 0.2,"p < e-16, R² = 0.29, n = 552" )
                  
  abline(rlm1 <- rlm (mod1) , lwd = 2.5, lty = 2, col = "darkorange")
    median (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope)
       
       summary ( lm1 <- lm(mod1)) # p-value: < 2.2e-16
       summary ( rlm1 <- rlm(mod1))
       #plot(mod1)
       
      # dev.off()
       
       ### FOR PAPER
       # NBSS slopes vs SST, Oly "TOP  phytopl model" ("p", R²) data used for regression, rlm --------------
       # lm and RR ar highly siognifocant, wether ALL data or only clean "p", R² OK data 
       
       # linear model is highly Significant (p < e-16, R²: 0.29,n = 552)
       # for ALL Phttopl. Data as well! 
       
       
       getwd()
       
      # svg("phytopl_Atl_FILTERED_RRslopessSSTv1.svg" , width = 5.6, height =  4)
        
       plot(rob_reg__slopes ~ SST_C_insitu, data = phytopl.All.Atl.p.rsqok,
            main = "Phytoplankton, ATLANTIC, All Data, n = 552",
            xlab = XlabNBSS,
            ylab = Ylab.NBBS, bty = "n",
            col = alpha("darkgreen", 0.3), pch = 16)
       
       text(28.5, 0.2,"p < e-16, R² = 0.29, n = 552" )
       
       abline(rlm1 <- rlm (mod1, data = phytopl.All.Atl.p.rsqok) ,
              lwd = 2.5, lty = 2, col = "darkorange")
      
        dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
       median (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope)
       
       summary ( lm1 <- lm(mod1 ,data = phytopl.All.Atl.p.rsqok) )# p-value: < 2.2e-16
       summary ( rlm1 <- rlm(mod1, data = phytopl.All.Atl.p.rsqok))
       #plot(mod1)
       
       # dev.off()
      
       
       
       
       
       
       
       
       
        
       ##### NBSS vs Chlorophyll (Top Models and ALL Data) -------------------------------
       
       
       
      ### FOR PAPER
       
       ##### NBSS vs Chlorophyll (Top Models and ALL Data) -------------------------------
       
       
       ### FOR PAPER
       # Unfiltered data
       # Phytoplankton NBSS slopes vs Chlorophyll, Only "TOP  phytopl model" ("p", R²) data used for regression, rlm --------------
       # lm and RR ar highly siognifocant, wether ALL data or only clean "p", R² OK data 
       
       # linear model is highly Significant (p < e-16, R²: 0.29,n = 552)
       # for ALL Phttopl. Data as well! 
       
       
       getwd()
       
        #svg("phytopl_Atl_ALL_UNFILTERED_RRslopesChla1.svg" , width = 5.6, height =  4)
       
        attach(phytopl.All.Atl)
       
       phytopl.All.Atl$log1Chlainsitu <- log(1+phytopl.All.Atl$Chlorophyll__mg_m_3_insitu)
       
       lm5 <- lm(new_slopes ~ log1Chlainsitu, data = phytopl.All.Atl)
       rlm5 <- rlm(new_slopes ~ log1Chlainsitu, data = phytopl.All.Atl)
       
       plot(new_slopes ~ log1Chlainsitu, data = phytopl.All.Atl,
            main = "Phytoplankton, ATLANTIC, All Data, n = 561 unfiltered data",
            ylab = "NBSS slope",
            xlab = "log10(1+Chla in situ)", bty = "n",
            col = alpha("darkgreen", 0.3), pch = 16)
       
       text(2, 0.49,"p < e-16, R² = 0.2, n = 561" )
       
       abline(rlm5, 
              lwd = 2.5, lty = 2, col = "darkorange")
       
       dim(phytopl.All.Atl)
       median (phytopl.All.Atl$new_slopes)
       
       summary ( lm5 )# p-value: < 2.2e-16
       summary ( rlm5)
       aovperm(rlm5) 
      
        
        #dev.off()
       
       
       
           # FILTERED Data
          # Phytoplankton NBSS slopes vs Chlorophyll, Only "TOP  phytopl model" ("p", R²) data used for regression, rlm --------------
       # lm and RR ar highly siognifocant, wether ALL data or only clean "p", R² OK data 
       
       # linear model is highly Significant (p < e-16, R²: 0.29,n = 552)
       # for ALL Phttopl. Data as well! 
       
       
       getwd()
       
       # svg("phytopl_Atl_FILTERED_RRslopessSSTv1.svg" , width = 5.6, height =  4)
       attach(phytopl.All.Atl.p.rsq.slope_ratio.ok)
       
       phytopl.All.Atl.p.rsq.slope_ratio.ok$log1Chlainsitu <- log(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu)
       
       lm3 <- lm(rob_reg__slopes ~ log1Chlainsitu, data = phytopl.All.Atl.p.rsq.slope_ratio.ok)
       rlm3 <- rlm(rob_reg__slopes ~ log1Chlainsitu, data = phytopl.All.Atl.p.rsq.slope_ratio.ok)
       
       plot(rob_reg__slopes ~ log1Chlainsitu, data = phytopl.All.Atl.p.rsq.slope_ratio.ok,
            main = "Phytoplankton, ATLANTIC, All Data, n = 454",
            xlab = "NBSS slope",
            ylab = "log10(+Chla in situ)", bty = "n",
            col = alpha("darkgreen", 0.3), pch = 16)
       
       text(28.5, 0.2,"p < e-16, R² = 0.29, n = 454" )
       
       abline(rlm3, 
              lwd = 2.5, lty = 2, col = "darkorange")
       
       dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
       median (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slope)
       
       summary ( lm3 )# p-value: < 2.2e-16, R2 = 0.35
       summary ( rlm3 )
       #plot(mod1)
       
       # dev.off()
       
       
       
       # Net- caught Zoopl. ----------- not so nice... few, badly formatted data ------
       # missing TSWA and BUS data (only 4 samples each!) 
       # Temperateure, chlorophyll, mesh size used!
       
       attach (zoopl_All_Atl4_Net)
       mod2 <- (rob_reg__slope ~ abs(Latitude))
       plot(mod2, data = zoopl_All_Atl4_Net,
            xlim = c(0 ,50), 
            main = "Linear-shaped net-caught zoopl. NBSS, n = 433", 
            # ylim = c(-2.2, 0.1),
            col = alpha ("navy", 0.3), pch = 16)
       
       abline(rlm2 <- rlm (mod2) , lwd = 2.5, lty = 2, col = "darkorange")
       summary(lm2 <- lm (mod2))
       abline(h = median (zoopl_All_Atl4_Net$rob_reg__slope), col = "darkgrey")
       median (zoopl_All_Atl4_Net$rob_reg__slope)
       # Median slope: -0.62, n = 433 linear-shaped net-caught zoopl. NBSS
       text ( 40,  -3.5 ,"Median slope: -0.62")
       
       dim(zoopl_All_Atl4_Net) # 134 samples, OK
     
       cor.test(rob_reg__slope , abs(Latitude), method = "spearman")
       #  NOT Significant!! Zoopl NBSS slope s INVARIABLE Avross Latitudes! p-value = 0.2052
       
       # Position of maximum ---------------
       summary(zoopl.All.Atl$x.value.of.maximum)
       
       attach (zoopl.All.Atl_p_rsq_OK)
       mod2 <- (log10(1+x.value.of.maximum) ~ abs(Latitude))
       plot(mod2, data = zoopl.All.Atl_p_rsq_OK,
            xlim = c(0 ,50), 
            main = "Linear-shaped net-caught zoopl. NBSS, n = 433", 
            # ylim = c(-2.2, 0.1),
            col = alpha ("navy", 0.3), pch = 16)
       
       #abline(rlm2 <- rlm (mod2) , lwd = 2.5, lty = 2, col = "darkorange")
       summary(lm2 <- lm (mod2))
       abline(h = median (zoopl.All.Atl_p_rsq_RR_OK4$rob_reg__slope), col = "darkgrey")
       median (zoopl.All.Atl_p_rsq_RR_OK4$rob_reg__slope)
       # Median slope: -0.62, n = 433 linear-shaped net-caught zoopl. NBSS
       text ( 40,  -3.5 ,"Median slope: -0.62")
       
       
       cor.test(x.value.of.maximum , abs(Latitude), method = "spearman")
       #  p = 0.024
       # Zoopl NBSS maximum position is related to Latitude! p-value = 0.024
       
       boxplot( zoo_NET.TSWA.5$x.value.of.maximum,
         zoo_NET.BUS.5$x.value.of.maximum,
         zoo_NET.EQU.5$x.value.of.maximum,
         zoo_NET.CCUS.5$x.value.of.maximum)
       
       boxplot( zoo_NET.TSWA.5$rob_reg__slope,
                zoo_NET.BUS.5$rob_reg__slope,
                zoo_NET.EQU.5$rob_reg__slope,
                zoo_NET.CCUS.5$rob_reg__slope)
       
       
       Data1 <- data.frame(
         Y<- c(zoo_NET.TSWA.5$rob_reg__slope,
             zoo_NET.BUS.5$rob_reg__slope,
             zoo_NET.EQU.5$rob_reg__slope,
             zoo_NET.CCUS.5$rob_reg__slope),
         Site =factor(rep(c("TSWA", "BUS", "EQU", "CCUS"),
                          times=c(length(zoo_NET.TSWA.5$rob_reg__slope), length(zoo_NET.BUS.5$rob_reg__slope),
                                  length(zoo_NET.EQU.5$rob_reg__slope), length(zoo_NET.CCUS.5$rob_reg__slope)
                                  )))
       )
       
       summary(aov( Y ~Site, data =  Data1)) # 6.19e-06 ***
       aovperm( Y ~ Site, data =  Data1) # p = 2e-04 
       #library(PMCMR)
       library(PMCMRplus)
       PMCMRplus::kwAllPairsNemenyiTest(Y ~Site, data =  Data1)
       attach(Data1)
     #   boxplot( Y ~ Site, main = "Zoopl. NBSS Slopes")
     #   text (2.5, -0.3,  "p<0.0001")
     #   text (3.5, -0.3,  "p= 0.036")
     #   dim(Data1)
     # #  dim()
       
       
       # UVP Zoopl. -----------
    
       attach (zoopl_All_Atl4_UVP)
       mod2 <- (rob_reg__slope ~ abs(Latitude))
       plot(mod2, data = zoopl_All_Atl4_UVP,
            xlim = c(0 ,50), 
            main = "Linear-shaped net-caught zoopl. NBSS, n = 433", 
            # ylim = c(-2.2, 0.1),
            col = alpha ("navy", 0.3), pch = 16)
       
       abline(rlm2 <- rlm (mod2) , lwd = 2.5, lty = 2, col = "darkorange")
       summary(lm2 <- lm (mod2))
       abline(h = median (zoopl_All_Atl4_UVP$rob_reg__slope), col = "darkgrey")
       median (zoopl_All_Atl4_UVP$rob_reg__slope)
       # Median slope: -0.62, n = 433 linear-shaped UVP zoopl. NBSS
       text ( 40,  -3.5 ,"Median slope: -0.531")
       
       dim(zoopl_All_Atl4_UVP) # 271 UVP samples with linear Shape
       
       cor.test(rob_reg__slope , abs(Latitude), method = "spearman")
       #  NOT Significant!! Zoopl NBSS slope s INVARIABLE Avross Latitudes! p-value = 0.2052
       
       # Position of maximum , zoopl_All_Atl.unfiltr_UVP , ALL UVP data ---------------
       summary(zoopl_All_Atl.unfiltr_UVP$x.value.of.maximum)
       
       attach (zoopl_All_Atl.unfiltr_UVP)
       mod2 <- (log10(x.value.of.maximum) ~ abs(Latitude))
       plot(mod2, data = zoopl_All_Atl.unfiltr_UVP,
            xlim = c(0 ,50), 
            main = "|All UVP zoopl. NBSS, n = 433", 
            # ylim = c(-2.2, 0.1),
            col = alpha ("navy", 0.3), pch = 16)
       
       abline(rlm2 <- rlm (mod2) , lwd = 2.5, lty = 2, col = "darkorange")
       summary(lm2 <- lm (mod2))
       #abline(h = median (zoopl.All.Atl_p_rsq_RR_OK4$rob_reg__slope), col = "darkgrey")
       median (zoopl_All_Atl.unfiltr_UVP$rob_reg__slope)
       # Median slope: -0.62, n = 433 linear-shaped net-caught zoopl. NBSS
       text ( 40,  -3.5 ,"Median slope: -0.62")
       
       
       cor.test(x.value.of.maximum , abs(Latitude), method = "spearman")
       #  p = 0.024
       # Zoopl NBSS maximum position is related to Latitude! p-value = 0.024
       
       boxplot( zoo_NET.TSWA.5$x.value.of.maximum,
                zoo_NET.BUS.5$x.value.of.maximum,
                zoo_NET.EQU.5$x.value.of.maximum,
                zoo_NET.CCUS.5$x.value.of.maximum)
       
       boxplot( zoo_NET.TSWA.5$rob_reg__slope,
                zoo_NET.BUS.5$rob_reg__slope,
                zoo_NET.EQU.5$rob_reg__slope,
                zoo_NET.CCUS.5$rob_reg__slope)
       
       
       Data1 <- data.frame(
         Y=c(zoo_NET.TSWA.5$rob_reg__slope,
             zoo_NET.BUS.5$rob_reg__slope,
             zoo_NET.EQU.5$rob_reg__slope,
             zoo_NET.CCUS.5$rob_reg__slope),
         Site =factor(rep(c("TSWA", "BUS", "EQU", "CCUS"),
                          times=c(length(zoo_NET.TSWA.5$rob_reg__slope), length(zoo_NET.BUS.5$rob_reg__slope),
                                  length(zoo_NET.EQU.5$rob_reg__slope), length(zoo_NET.CCUS.5$rob_reg__slope)
                          )))
       )
       
       summary(aov( Y ~Site, data =  Data1)) # 6.19e-06 ***
       aovperm( Y ~ Site, data =  Data1) # p = 2e-04 
       #library(PMCMR)
       library(PMCMRplus)
       PMCMRplus::kwAllPairsNemenyiTest(Y ~Site, data =  Data1)
       attach(Data1)
       # boxplot( Y ~ Site, main = "Zoopl. NBSS Slopes")
       # text (2.5, -0.3,  "p<0.0001")
       # text (3.5, -0.3,  "p= 0.036")
       # dim(Data1)
       # #dim()
       
       
       
       
       
       ### NON-Linear analysis ###########--------------  
       ### Kernel smoothing ###############
       ### Kernel smoothing and postion of the maximum estimation (minimum 5 data) ####
       
       
       
       ####################################################
       # # Functional analysis ####
       # General functional data analysis
       # The packages listed below provide “infrastructure” for representing and handling function-valued data and/or implement many widely applicable functional data methods:
       #   
       #   fda provides object-types for functional data with corresponding functions for smoothing, plotting and simple regression models, c.f. Ramsay et al. (2009, doi:10.1007/978-0-387-98185-7).
       # fdapace provides functional principal component based methods for sparsely or densely sampled random trajectories and time courses for functional regression and correlation, for longitudinal data analysis, the analysis of stochastic processes from samples of realized trajectories, and for the analysis of underlying dynamics.
       # fdasrvf performs alignment, PCA, and regression of multidimensional or unidimensional functions using the square-root velocity framework (Srivastava et al., 2011). This framework allows for elastic analysis of functional data through phase and amplitude separation.
       # fda.usc provides routines for exploratory and descriptive analysis of functional data such as depth measurements, outlier detection, as well as unsupervised and supervised classification, (univariate, nonparametric) regression models with a functional covariate and functional analysis of variance.
       # fds contains 19 data sets with functional data.
       # funData provides S4 classes for univariate and multivariate functional and image data and utility functions.
       # tf provides S3 vector data types for functional data with arithmetic and summary methods, derivation, integration and smoothing, plotting, data import and export, and data wrangling (subsetting, sub-assigning, zooming, extracting functional features) that allow including such vectors in data frames for joint analysis of functional and scalar variables.
       
       
     ### Figures actually used in the paper -------------
       
                # Fig. 1 (map)
      
        
         # Fig. 2 ---------------------------
       
       # Start of Figure 2 ------------------
       
       #############
       ## T-S diagram --- points with colours, FINAL PLOT FOR PAPER !!!! ----------------
       
       # colour plot with ggplot2 ----------------
       
       summary(PP.ATL.df4$Salinity_insitu)
       summary(Chlorophyll__mg_m_3_insitu)
       summary(SST_C_insitu)
       
       
       df4insitu <- data.frame ( SST = PP.ATL.df4$SST_C_insitu,
                                 Sal =  PP.ATL.df4$Salinity_insitu,
                                 Temp = PP.ATL.df4$Temperature__C_insitu,
                                 NBSS_slopes =  PP.ATL.df4$rob_reg__slopes, 
                                 logChl_a = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu))
       
       dim(df4insitu)# N = 458   datasets (with NAs)
       summary(df4insitu)
       df4insitu <- df4insitu[complete.cases(df4insitu),]
       dim(df4insitu) # 575 complete datsets
       summary(df4insitu)
       # N = 205   datasets (without NAs)
       df4insitu <- sort_by(df4insitu, df4insitu$SST )
       
       
       names(df4insitu)
       
       c <- ggplot(df4insitu, aes(x = Sal, y = SST))+ #+ xlim(31.5, 38) + ylim( 5,30 )
         geom_point(aes(color = logChl_a) , size = 3) +
         scale_color_gradientn(colors = colors4b) +
         theme(legend.position = "right")
       
       # End of Figure 2 ------------
       
       
       
       # Start  of Figure 3 ------------
    
          
       ### Calculate median with ALL data (not only the ones filtered for linearity)
       
       
       phytopl.All.Atl.NBSSmatrix.1 <- as.matrix(phytopl.All.Atl[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.1)
       
       phytopl.All.Atl.NBSSmatrix.OK3 <- as.matrix(phytopl.All.Atl.p.rsq.slope_ratio.ok[,37:91])
       dim(phytopl.All.Atl.NBSSmatrix.OK3)
       
       
    
       
       
       # PLot ALL NBSS Data, Carbon  -------
      
      nrow(phytopl.All.Atl.NBSSmatrix.1) # 721  ALL phytopl. datasets, OK
       
      nrow(phytopl.All.Atl.NBSSmatrix.OK3) # 458 high-quality phytopl. datasets, OK
       
       # plot all points and rlm moldes (blue lines), and red median line  ----------
       
       lnumb <- 34
       plot (log10( phytopl.All.Atl.NBSSmatrix.1[lnumb,]) ~ log10(X_vector_gCind), 
             main = "ATLANTIC, n = 727",
             xlab = "log10(gC ind. -1) ",
             ylab = "log10(gC m-3 gC-1 ind.) ",
             xlim = c(-15, -5), ylim = c(0, 14) ,
             pch = 16, col = "white")
       
       for (lnumb in 1 : nrow(phytopl.All.Atl.NBSSmatrix.OK3) )  {
         
       abline( rlm(log10( phytopl.All.Atl.NBSSmatrix.OK3[lnumb,]) ~ 
                     log10(X_vector_gCind)), 
               col =alpha ("grey", 0.3))

             } 
       
       for (lnumb in 1 : nrow(phytopl.All.Atl.NBSSmatrix.1) ) 
       {    
         points (log10( phytopl.All.Atl.NBSSmatrix.1[lnumb,]) ~ log10(X_vector_gCind) ,
                 pch = 16, col = alpha("forestgreen", 0.2))
       }
       
       # Carbon  (458 high-quality datasets, OK)
       (median(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
       # median RR slope : -1.256
       (min(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes)) 
       (max(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) )
       # min - max  RR slope : -0.44155
       
       
       
       # matrix to to vector (transpose),  -----------------------  
       x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.NBSSmatrix.1) )
       y1 <- c ( t(phytopl.All.Atl.NBSSmatrix.1))
       
       # regression line (Rob. regr.)
       #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
       #points(  log10(y1) ~ x1 )
       # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
       abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
                lwd = 2.5, lty = 2, col = "darkorange" )
       summary(rlm2)
    #   rlm2$residuals
       
       #rsq <- ( 1 - sum(rlm2$residuals^2)/ (sum(( log10(y1)  - mean( log10 (y1))^2))  ))
       
       # slope = -1.146 , using all green points, rlm
       # means and medians without considering NAs  
       mean (c(4,NA, 18,34)) # NA
       median( c(4,NA, 18, 34)  ) #NA
    
       
       
       # apply median by columns (only if less thn 50% NAs)
       
       
       median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.NBSSmatrix.1, 2,  median_NA50) 
       
       length(median_vec_ALLBiovol2)
       
       class(median_vec_ALLBiovol2)
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
    
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
 
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "yellow" , lwd = 2.5)
      
       
       x2 <- rep( (log10(X_vector_gCind)), nrow(phytopl.NBSS.matrix.1) )
       y2 <- c ( t(phytopl.NBSS.matrix.1))
       
       df2 <- data.frame (y2 = y2 ,x2 = x2, logy2 = log10(y2)  )
       
       dfok2 <- df2
       
       dfok2[dfok2=="-Inf"]<-NA
       dfok2[dfok2=="NaN"]<-NA
       
      
       dim(df2)
       dim(dfok2)
       
       summary(df2)
       summary(dfok2)
       
       
       abline(  (rlm2 <- rlm(logy2 ~ x2,data = dfok2)),  
                lwd = 2.5, lty = 2, col = "darkorange" )
       summary(rlm2)
       
      # permuco::aovperm(rlm2)
     
       lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),
               pch = 16, col = "pink3" , lwd = 2.7)
       
       lines(  log10(median_vec_ALLBiovol2) ~ log10(X_vector_gCind),
               pch = 16, col = "red" , lwd = 2.7)
       
       
       ### end of Figure 3 --------------
       
       # Fig. Sxxxx: boxplot
       
       
       # Start of Figure 4 --------------------
      
       # plot and matrix as function --------------
       
       
       ### START OF function  funplot.matrixlinesmedians ###
       
       fun98_calc_quant <- function(x) {quantile(x, probs = 0.98, na.rm =  TRUE) }
       
       
       
       funplot.matrixlinesmedians(phytopl.unfiltr.ALL.TSWA.5, "TSWA, unfiltr")
      
       fun95_calc_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
       
       
       fun98.add.98quant.line.max.vec <- function(y, xvec) { 
       
         mat.y <- as.matrix(y[,37:91])
         
       fun98_calc_quant <- function(x) {quantile(x, probs = 0.98, na.rm =  TRUE) }
       
       max_vec_ALL <- apply(mat.y, 2, fun98_calc_quant ) 
        # 
       median_vec_PHYTOP <- as.vector(unlist(max_vec_ALL))
       
      
       df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (xvec))
       
       dfok<- df1
       #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
       
       lines(  log10(dfok$y1) ~ log10(dfok$x1), 
               col = "red" , lwd = 4.5)
       
       dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (xvec))
       lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
               col = "purple" , lwd = 4.5)
       
       }
       
       funplot.matrixlinesmedians(phytopl.unfiltr.ALL.TSWA.5, "TSWA, unfiltr")
       fun98.add.98quant.line.max.vec(phytopl.unfiltr.ALL.TSWA.5, X_vector_gCind)
       text("98% quantiles", x = -9.3, y = 8.5, col = "purple")
       
        
       funplot.matrixlinesmedians(phytopl.unfiltr.ALL.EQU.5, "EQU, unfiltr")
       fun98.add.98quant.line.max.vec(phytopl.unfiltr.ALL.EQU.5, X_vector_gCind)
       text("98% quantiles", x = -9.3, y = 8.5, col = "purple")
       
        
       funplot.matrixlinesmedians(phytopl.unfiltr.ALL.CCUS.5, "CCUS, unfiltr")
       fun98.add.98quant.line.max.vec(phytopl.unfiltr.ALL.CCUS.5, X_vector_gCind)
       text("98% quantiles", x = -7.3, y = 8.5, col = "purple")
       
          
       funplot.matrixlinesmedians(phytopl.unfiltr.ALL.BUS.5, "BUS, unfiltr")
       fun98.add.98quant.line.max.vec(phytopl.unfiltr.ALL.BUS.5, X_vector_gCind)
       text("98% quantiles", x = -7.3, y = 8.5, col = "purple")
       
       
       #end of Figure 4 -----------------
    
       
       #######################   
       # Figure 5 Phytoplankton color plot (linear shaped only) -------------------------------------
       
       
       ### FINAL PLOT FOR PAPER ---------------------------------
       #  plot for paper OK, FOR PAPER ---------              
       ## NBSS_slopes
       
       
       df4insitu_v9_slop_chl_Temp <- data.frame (Temp = PP.ATL.df4$Temperature__C_insitu,
                                 NBSS_slopes =  PP.ATL.df4$rob_reg__slopes, 
                                 logChl_a = log10(1+PP.ATL.df4$Chlorophyll__mg_m_3_insitu))
       
              dim(df4insitu_v9_slop_chl_Temp)# N = 458   datasets (with NAs)
       summary(df4insitu_v9_slop_chl_Temp)
       df4insitu_v9_slop_chl_Temp <- df4insitu_v9_slop_chl_Temp[complete.cases(df4insitu_v9_slop_chl_Temp),]
       dim(df4insitu_v9_slop_chl_Temp) # 415 complete datsets
      
        summary(df4insitu_v9_slop_chl_Temp)
       # N = 415   datasets (without NAs)
       df4insitu_v9_slop_chl_Temp <- sort_by(df4insitu_v9_slop_chl_Temp, df4insitu_v9_slop_chl_Temp$Temp )
       
       
       c2 <- ggplot(df4insitu_v9_slop_chl_Temp, aes(x = logChl_a, y = Temp)) #+ xlim(31.5, 38) + ylim( 5,30 )
       
       
       c2 + geom_point(aes(color = NBSS_slopes) , size = 3) +
         scale_color_gradientn(colors = colors5) +
         theme(legend.position = "right")
       
       #end of Figure 5 -----------------
       
       
       #######################   
       # Figure 6 Phytoplankton NBSS Slopes scatter plots Slope vs Temperature and Chla -------------------------------------
       
       #######################   
       # Figure 7 Phytoplankton NBSS  plots Slope by Temperature and Chla groups -------------------------------------
       
       #######################   
       # Figure 8 Phytoplankton NBSS  plots Slope by Temperature and Chla groups -------------------------------------
       
       #######################   
       # Figure 9 NET Zooplankton NBSS  plot -------------------------------------
       
       lnumb <- 34
       plot (log10(zoopl.NBSSmatrix_4.NET[lnumb, ]) ~ logCvector, 
             main = "ATLANTIC, Net-caught Zooplankton",
             xlab = "log10( C Biomass (gC ind.-1)) ",
             ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
             xlim = c(-7, 1), 
             ylim = c(-4, 5) ,
             pch = 16, col = "white")
       for (lnumb in 1 : nrow(zoopl.NBSSmatrix_4.NET) ) 
       {    
         results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                               zoopl.NBSSmatrix_4.NET[lnumb,] )
         
         b =  results$results$slope_OLSR
         #[1]  -0.6459868
         
         a = results$results$intercept_OLSR
         #[1] -0.7477856
         
         abline( a= a , b= b, col =alpha ("grey65", 0.3))
         
         points (log10(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na[lnumb,]) ~ logCvector ,
                 pch = 16, col = alpha("navy", 0.2))
         
         
         # abline( rlm(log10( zoopl.NBSSmatrix_4.NET[lnumb,]) ~ logCvector),
         #         col =alpha ("grey65", 0.3))
       }
       
       
       names(zoopl_All_Atl4NET)
       median(zoopl_All_Atl4NET$intercept_OSLR_max.to.NA)
       # medioan lm intercept :-3.464819  carbon
       median(zoopl_All_Atl4NET$rob_reg_intercept)
       
       
       # Net zoopl., n = 111
       # Biovolume  (462 useful Data)
       (median(zoopl_All_Atl4NET$rob_reg__slope) )
       # median RR slope : -0.914 (bivol)
       # median RR slope : -0.7288432 (CARBON)
       
       (min(zoopl_All_Atl4NET$rob_reg__slope)) 
       (max(zoopl_All_Atl4NET$rob_reg__slope) )
       # min to max  RR slope 
       # -2.1248 to -0.6658 
       
       summary(zoopl_All_Atl4NET$rob_reg__slope)
       # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
       # -2.1248 -1.0195 -0.9147 -0.9821 -0.8328 -0.6658 
       
       # 
        # abline( a = mean(zoopl_All_Atl4NET$rob_reg_intercept),
        #         b =  mean(zoopl_All_Atl4NET$rob_reg__slope)   ,
        #   lwd = 2.5, lty = 2, col = "darkorange" )
       # Median RLM model fom all sample-by sample a  analysis
        abline( a = median(zoopl_All_Atl4NET$intercept_OSLR_max.to.NA),
                b =  median(zoopl_All_Atl4NET$slope_OSLR_max.to.NA)   ,
                lwd = 2.5, lty = 2, col = "darkorange" )
        
        # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
        
        abline(a = -3.2915, b = -0.9548 ,  
               lwd = 2.5, lty = 2, col = "grey44" )
        
        
        
               
       # apply median by columns (only if less thn 50% NAs)
       
       
       median_vec_ALLBiovol2 <- apply(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na, 2,  median_NA50) 
       
       length(median_vec_ALLBiovol2)
       
       class(median_vec_ALLBiovol2)
       
       median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
       
       df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
       
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "yellow" , lwd = 2.5)
       
       
       x2 <- rep( (log10(X_vector_gCind)), nrow(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na) )
       y2 <- c ( t(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na))
       
       df2 <- data.frame (y2 = y2 ,x2 = x2, logy2 = log10(y2)  )
       
       dfok2 <- df2
       
       dfok2[dfok2=="-Inf"]<-NA
       dfok2[dfok2=="NaN"]<-NA
       
       
       dim(df2)
       dim(dfok2)
       
       summary(df2)
       summary(dfok2)
            # matrix to to vector (transpose),  -----------------------  
       x1 <- rep( (log10(X_vector_gCind)), nrow(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na) )
       y1 <- c ( t(zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na))
       
       # regression line (Rob. regr.)
       #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
       #points(  log10(y1) ~ x1 )
       # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
      #  abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
      #           lwd = 2.5, lty = 2, col = "darkorange" )
      #  summary(rlm2)
      # # rlm2$residuals
       
       # 
       # 
       # abline(  (rlm2 <- rlm(logy2 ~ x2,data = dfok2)),  
       #          lwd = 2.5, lty = 2, col = "darkorange" )
       # summary(rlm2)
       # 
       # permuco::aovperm(rlm2)
       
       # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),
       #         pch = 16, col = "pink3" , lwd = 2.7)
       # 
       lines(  log10(median_vec_ALLBiovol2) ~ log10(X_vector_gCind),
               pch = 16, col = "salmon" , lwd = 2.7)
       
       lines(   logmedians5_largeNET ~ logCvector_largeNET,
                pch = 16, col = "red" , lwd = 5.5)
       
       
      summary(zoopl_All_Atl4NET$N.bins.used)
      # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      # 5.00    8.00    9.00    9.59   11.00   18.00 
             
# for Net-caught zooplankn five t 18 bins wer used for each linear model fit (median: 9 useful linear bins).
      
      summary(zoopl_All_Atl4NET$x.value.of.maximum)
      # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
      # 1.190e-07 1.910e-06 3.810e-06 1.166e-05 1.530e-05 6.100e-05 
            
    # The position of the maximum (start o the linear model fit) varied
      # considerably between samples, from 1.2 e-7 to 6.1 e-5 gC ind.-1 (median position of the maximum : 3.81 e-06 gC ind.1) 
      abline( v = log10(3.810e-06))
             
      ### end of Figure 9 (Net zoopl) --------------
      
    
      ### Figure 10 UVP Zoopl. --------------
      
     
      ### FOR PAPER!!
       #######################   
      # Figure 10 UVP  Zooplankton NBSS  plot -------------------------------------
      
     
      lnumb <- 34
      plot (log10(zoopl.NBSSmatrix_4.UVP[lnumb, ]) ~ logCvector, 
            main = "ATLANTIC, UVP Zooplankton",
            xlab = "log10( C Biomass (gC ind.-1)) ",
            ylab = "log10(Normalized C Biomass (gC m-3 / gC ind.-1) ",
            xlim = c(-7, 1), 
            ylim = c(-4, 5) ,
            pch = 16, col = "white")
      for (lnumb in 1 : nrow(zoopl.NBSSmatrix_4.UVP) ) 
      {    
        results <- NBSS.select.w.lm.MS.B.RR ( 10^logCvector , 
                                              zoopl.NBSSmatrix_4.UVP[lnumb,] )
        
        b =  results$results$slope_OLSR
        #[1]  -0.6459868
        
        a = results$results$intercept_OLSR
        #[1] -0.7477856
        
        abline( a= a , b= b, col =alpha ("grey65", 0.3))
        
        points (log10(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na[lnumb,]) ~ logCvector ,
                pch = 16, col = alpha("darkmagenta", 0.2))
        
        
        # abline( rlm(log10( zoopl.NBSSmatrix_4.UVP[lnumb,]) ~ logCvector),
        #         col =alpha ("grey65", 0.3))
      }
      
      names(zoopl_All_Atl4UVP)
      median(zoopl_All_Atl4UVP$intercept_OSLR_max.to.NA)
      # medioan lm intercept :-3.464819  carbon
      median(zoopl_All_Atl4UVP$rob_reg_intercept)
      
      
      # UVP zoopl., n = 111
      # Biovolume  (462 useful Data)
      (median(zoopl_All_Atl4UVP$rob_reg__slope) )
      # median RR slope : -0.914 (bivol)
      # median RR slope : -0.7288432 (CARBON)
      
      (min(zoopl_All_Atl4UVP$rob_reg__slope)) 
      (max(zoopl_All_Atl4UVP$rob_reg__slope) )
      # min to max  RR slope 
      # -2.1248 to -0.6658 
      
      summary(zoopl_All_Atl4UVP$rob_reg__slope)
      # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      # -2.1248 -1.0195 -0.9147 -0.9821 -0.8328 -0.6658 
      
      # 
      # abline( a = mean(zoopl_All_Atl4UVP$rob_reg_intercept),
      #         b =  mean(zoopl_All_Atl4UVP$rob_reg__slope)   ,
      #   lwd = 2.5, lty = 2, col = "darkorange" )
      # 
      abline( a = median(zoopl_All_Atl4UVP$intercept_OSLR_max.to.NA),
              b =  median(zoopl_All_Atl4UVP$slope_OSLR_max.to.NA)   ,
              lwd = 2.5, lty = 2, col = "darkorange" )
      
      
      
      # apply median by columns (only if less thn 50% NAs)
      
      
      median_vec_ALLBiovol2 <- apply(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na, 2,  median_NA50) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "yellow" , lwd = 2.5)
      
      
      x2 <- rep( (log10(X_vector_gCind)), nrow(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na) )
      y2 <- c ( t(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na))
      
      df2 <- data.frame (y2 = y2 ,x2 = x2, logy2 = log10(y2)  )
      
      dfok2 <- df2
      
      dfok2[dfok2=="-Inf"]<-NA
      dfok2[dfok2=="NaN"]<-NA
      
      
      dim(df2)
      dim(dfok2)
      
      summary(df2)
      summary(dfok2)
      # matrix to to vector (transpose),  -----------------------  
      x1 <- rep( (log10(X_vector_gCind)), nrow(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na) )
      y1 <- c ( t(zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na))
      
      # regression line (Rob. regr.)
      #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
      #points(  log10(y1) ~ x1 )
      # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
      #  abline(  (rlm2 <- rlm(log10(y1) ~ x1)),  
      #           lwd = 2.5, lty = 2, col = "darkorange" )
      #  summary(rlm2)
      # # rlm2$residuals
      
      # 
      # 
      # abline(  (rlm2 <- rlm(logy2 ~ x2,data = dfok2)),  
      #          lwd = 2.5, lty = 2, col = "darkorange" )
      # summary(rlm2)
      # 
      # permuco::aovperm(rlm2)
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),
      #         pch = 16, col = "pink3" , lwd = 2.7)
      # 
      lines(  log10(median_vec_ALLBiovol2) ~ log10(X_vector_gCind),
              pch = 16, col = "red" , lwd = 2.7)
      
      
      
      summary(zoopl_All_Atl4UVP$N.bins.used)
      # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      # 5.00    8.00    9.00    9.59   11.00   18.00 
      
      # for UVP-caught zooplankn five t 18 bins wer used for each linear model fit (median: 9 useful linear bins).
      
      summary(zoopl_All_Atl4UVP$x.value.of.maximum)
      # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
      # 1.190e-07 1.910e-06 3.810e-06 1.166e-05 1.530e-05 6.100e-05 
      
      # The position of the maximum (start o the linear model fit) varied
      # considerably between samples, from 1.2 e-7 to 6.1 e-5 gC ind.-1 (median position of the maximum : 3.81 e-06 gC ind.1) 
      #abline( v = log10(3.810e-06))
      
      
      # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
      
      abline(a = -3.2915, b = -0.9548 ,  
             lwd = 2.5, lty = 2, col = "grey44" )
      
      
      
      
      ### end of Figure 10 --------------
      
      
      ### Figure 11 - two graphs UVP Shapes Dome vs ...  and Slope vs temp.--------------
      
      
      
      ### Figure 12 Ecosystem NBSS ALL Atlantic--------------
    ### Two graphs in oe plot
      
      
      # ALL Data (colour clouds)
        
      
      
      # grey plot with  medians
    
      
      
      ### ECOSYSTEM  NBSS, ALL ATLANTIC --------------- 
      #### ALL TAXA ------------------
      # For Paper ALL TAxa ------------
      # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
      # FOR PAPER -  grey plot wit coloured medians 
      ### Plot All data (grey) with medians and regression lines ------------------- 
      
      # plot ALL NBSS Data (loop) --------------
      
      library(scales)
      
      lnumb <- length(Data.All.Atl$Longitude)
      lnumb# 2840 data (lines) in C format
      
      lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
      
      plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
            ylim = c(-7, 13),col = "white",
            xlim = c(-14, 3),
            xlab = "log10(gC ind.-1)",
            ylab = "log10(gC m-3 / gC ind.-1)",
            main= "NBSS, All Atlantic, All data, n = 2840, C units")
      
      for (i in 1: (nrow(Data.All.Atl.NBSSmatrix)) ){
        
        points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                pch = 16, col = alpha("darkgrey", 0.05))
        
      }
      
      # Naive linear model (all data): ------------------
      # matrix to to vector (transpose), TSWA -----------------------  
      library (MASS)
      
      x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
      y1 <- c ( t(Data.All.Atl.NBSSmatrix))
      
      df1 <- data.frame(x1, y1, ylog = log10(y1))
      names(df1)
      
      # -INF to NA!!!
      
      # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
      
      df1[df1=="-Inf"]<-NA
      
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      
      summary(df1)
      summary(df2)
      
      mean(df2$y1)
      
      summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.954
      summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.955
      
      rlm1.ALL.ATL.ALLTaxa.naive <- rlm1
      
      abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
             lwd = 2.5, lty = 2, col = "darkorange" )
      
      abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
             lwd = 2.5, lty = 2, col = "grey44" )
      
      
      # medians (curved line) ------
      
      # Phtopl - ALL Areas - unfiltered OK , only medians with < 50% NA
      median_vec_PHYTOP <-  median_vec_ALLBiovol2 
      
      logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
      
      
      
      
      median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.NBSSmatrix.1, 2,  median_NA50) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "yellow" , lwd = 2.5)
      
      
      x2 <- rep( (log10(X_vector_gCind)), nrow(phytopl.NBSS.matrix.1) )
      y2 <- c ( t(phytopl.NBSS.matrix.1))
      
      df2 <- data.frame (y2 = y2 ,x2 = x2, logy2 = log10(y2)  )
      
      dfok2 <- df2
      
      dfok2[dfok2=="-Inf"]<-NA
      dfok2[dfok2=="NaN"]<-NA
      
      
      dim(df2)
      dim(dfok2)
      
      summary(df2)
      summary(dfok2)
      
      
      abline(  (rlm2 <- rlm(logy2 ~ x2,data = dfok2)),  
               lwd = 2.5, lty = 2, col = "darkorange" )
      summary(rlm2)
      
      # permuco::aovperm(rlm2)
      
       median_vec_PHYTOP <- median_vec_ALLBiovol2
      
      lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
              pch = 16, col = "forestgreen" , lwd = 5.5)
      
      
      
      # Net Zoopl. , medians,  - ALL Areas -------
      
      
      #abline (v =  c (-4.9, -2)) 
      
      # all medians (no cutoff selection), but olnly < 50% NA
      
      # Filtered NET zoopl. samples only-------------
      # lines(   logmedians5_largeNET ~ logCvector_largeNET,
      #         pch = 16, col = "dodgerblue" , lwd = 5.5)
      
      
      # Unfiltered NET zoopl. samples -------
      
      median_vec_ALLNET.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na, 2, median_NA50) 
      
      
      logmedians5 <- log10(median_vec_ALLNET.NA_unfiltr)
      length(logmedians5) # 55 data points with medians
      (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
      logmedians5_large <- logmedians5[30:55]
      logmedians5_small <- logmedians5[1:30]
      logCvector_large <- logCvector[30:55]
      logCvector_small <- logCvector[1:30]
      
      logCvector_largeNET <- logCvector_large
      logCvector_smallNET <- logCvector_small
      
      logmedians5_largeNET <- logmedians5_large
      logmedians5_smallNET <- logmedians5_small
      
      
      
      
      # Unfiltered NET zoopl. samples, ALL --------
      
      lines(   logmedians5_largeNET ~ logCvector_largeNET,
               pch = 16, col = "dodgerblue" , lwd = 5.5)
      
      
      # UVP ----------
      
      # UVP  Zoopl. , medians,  - ALL Areas -------
      
      # Unfiltered UVP zoopl. samples -------
      
      median_vec_ALLUVP.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na, 2, median_NA50) 
      
      
      logmedians5 <- log10(median_vec_ALLUVP.NA_unfiltr)
      
      
      lines(  log10(median_vec_ALLUVP.NA_unfiltr) ~ logCvector,
              pch = 16, col = "magenta" , lwd = 5.5)
      
      
      
      # Mesop. fish, medians,  - ALL Areas -------
      
      
      # all medians (no cutoff selection), but olnly < 50% NA
      lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
               pch = 16, col = "darkorange" , lwd = 5.5)
      
      
      # Micronekt, medians, - ALL Areas -------
      
      
      # all medians (no cutoff selection), but olnly < 50% NA
      lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
               pch = 16, col = "darkred" , lwd = 5.5)
      
      
      # # medians with selection cutoff point  
      # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
      #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
      # 
      
      
      # Add legend -----------
      
      # leg.names <- c(  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
      # col.taxa  <- c(  "forestgreen" , "dodgerblue" ,  "magenta", "darkorange" , "darkred"  ) 
      # 
      # legend("topright", leg.names, lty = 1,
      #        lwd = 5.5, 
      #        col = col.taxa)
      # 
      
       
      ### end of grey plot with coloured medians (ALL Atlantic) ### ------------
     
      
      
      
      ### ECOSYSTEM  NBSS, ALL ATLANTIC --------------- 
      #### ALL TAXA ------------------
      # For Paper ALL TAxa ------------
      # GREY PLOT WITH COLORED MEDIANS and "naive" rlm line in dark grey ###
      # FOR PAPER -  grey plot wit coloured medians 
      ### Plot All data (grey) with medians and regression lines ------------------- 
      
      # plot ALL NBSS Data (loop) --------------
      
      library(scales)
      
      lnumb <- length(Data.All.Atl$Longitude)
      lnumb# 2840 data (lines) in C format
      
      lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
      
      plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
            ylim = c(-7, 13),col = "white",
            xlim = c(-14, 3),
            xlab = "log10(gC ind.-1)",
            ylab = "log10(gC m-3 / gC ind.-1)",
            main= "NBSS, All Atlantic, All data, n = 2840, C units")
      
      for (i in 1: (nrow(Data.All.Atl.NBSSmatrix)) ){
        
        points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                #            pch = 16, col = alpha(Data.All.Atl$target.org.colour[i], 0.1))
                pch = 16, col = alpha("darkgrey", 0.05))
        
      }
      
      # Naive linear model (all data): ------------------
      # matrix to to vector (transpose), TSWA -----------------------  
      library (MASS)
      
      x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
      y1 <- c ( t(Data.All.Atl.NBSSmatrix))
      
      df1 <- data.frame(x1, y1, ylog = log10(y1))
      names(df1)
      
      # -INF to NA!!!
      
      # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
      
      df1[df1=="-Inf"]<-NA
      
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      
      summary(df1)
      summary(df2)
      
      mean(df2$y1)
      
      summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.954
      summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.955
      
      rlm1.ALL.ATL.ALLTaxa.naive <- rlm1
      
      abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
             lwd = 2.5, lty = 2, col = "darkorange" )
      
      abline(rlm1.ALL.ATL.ALLTaxa.naive, xlim = c(-11, -1), ylim = c(0, 14), 
             lwd = 2.5, lty = 2, col = "grey44" )
      
      
      # medians (curved line) ------
      
      # Phtopl - ALL Areas - unfiltered OK , only medians with < 50% NA
      median_vec_PHYTOP <-  median_vec_ALLBiovol2 
      
      logX_vector_gCind_PHYTOPL <- log10(X_vector_gCind)
      
      
      
      
      median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.NBSSmatrix.1, 2,  median_NA50) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      df1 <- data.frame (y1 = ((median_vec_ALLBiovol2)) ,x1 =   (X_vector_gCind))
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "yellow" , lwd = 2.5)
      
      
      x2 <- rep( (log10(X_vector_gCind)), nrow(phytopl.NBSS.matrix.1) )
      y2 <- c ( t(phytopl.NBSS.matrix.1))
      
      df2 <- data.frame (y2 = y2 ,x2 = x2, logy2 = log10(y2)  )
      
      dfok2 <- df2
      
      dfok2[dfok2=="-Inf"]<-NA
      dfok2[dfok2=="NaN"]<-NA
      
      
      dim(df2)
      dim(dfok2)
      
      summary(df2)
      summary(dfok2)
      
      
      abline(  (rlm2 <- rlm(logy2 ~ x2,data = dfok2)),  
               lwd = 2.5, lty = 2, col = "darkorange" )
      summary(rlm2)
      
      # permuco::aovperm(rlm2)
      
      median_vec_PHYTOP <- median_vec_ALLBiovol2
      
      lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
              pch = 16, col = "forestgreen" , lwd = 5.5)
      
      
      
      # Net Zoopl. , medians,  - ALL Areas -------
      
      
      abline (v =  c (-4.9, -2)) 
      
      # all medians (no cutoff selection), but olnly < 50% NA
      
      # Filtered NET zoopl. samples only-------------
      # lines(   logmedians5_largeNET ~ logCvector_largeNET,
      #         pch = 16, col = "dodgerblue" , lwd = 5.5)
      
      
      # Unfiltered NET zoopl. samples -------
      
      median_vec_ALLNET.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_NET.MATRIX_no_na, 2, median_NA50) 
      
      
      logmedians5 <- log10(median_vec_ALLNET.NA_unfiltr)
      length(logmedians5) # 55 data points with medians
      (cut_point_Xval <- logCvector[30]) #  cutoff at -4.5 log10 C ind.-1
      logmedians5_large <- logmedians5[30:55]
      logmedians5_small <- logmedians5[1:30]
      logCvector_large <- logCvector[30:55]
      logCvector_small <- logCvector[1:30]
      
      logCvector_largeNET <- logCvector_large
      logCvector_smallNET <- logCvector_small
      
      logmedians5_largeNET <- logmedians5_large
      logmedians5_smallNET <- logmedians5_small
      
      
      
      
      # Unfiltered NET zoopl. samples, ALL --------
      
      lines(   logmedians5_largeNET ~ logCvector_largeNET,
               pch = 16, col = "dodgerblue" , lwd = 5.5)
      
      
      # UVP ----------
      
      # UVP  Zoopl. , medians,  - ALL Areas -------
      
      # Unfiltered UVP zoopl. samples -------
      
      median_vec_ALLUVP.NA_unfiltr <- apply( zoopl_All_Atl4.unfiltr_UVP.MATRIX_no_na, 2, median_NA50) 
      
      
      logmedians5 <- log10(median_vec_ALLUVP.NA_unfiltr)
      
      
      lines(  log10(median_vec_ALLUVP.NA_unfiltr) ~ logCvector,
              pch = 16, col = "magenta" , lwd = 5.5)
      
      
      
      # Mesop. fish, medians,  - ALL Areas -------
      
      
      # all medians (no cutoff selection), but olnly < 50% NA
      lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
               pch = 16, col = "darkorange" , lwd = 5.5)
      
      
      # Micronekt, medians, - ALL Areas -------
      
      
      # all medians (no cutoff selection), but olnly < 50% NA
      lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
               pch = 16, col = "darkred" , lwd = 5.5)
      
      
      # # medians with selection cutoff point  
      # lines(  logmedians5_smallMESO_FI ~ logCvector_smallMESO_FI,
      #         pch = 16, col = alpha ("salmon",1) , lwd = 3.5)
      # 
      
      
      # Add legend -----------
      
      leg.names <- c(  "Phytopl.",  "Net zoopl." , "UVP Zoopl.",  "Mesop. Fish" , "Micronekton"  ) 
      col.taxa  <- c(  "forestgreen" , "dodgerblue" ,  "magenta", "darkorange" , "darkred"  ) 
      
      legend("topright", leg.names, lty = 1,
             lwd = 5.5, 
             col = col.taxa)
      
      
      
      ### end of grey plot with coloured medians (ALL Atlantic) ### ------------
      
      # insert maximum (95 % quantile)
      
      fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
      
      max_vec_ALL <- apply(phytopl.NBSS.matrix.ok3, 2, fun95_quant ) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
      
      
      
      median_vec_PHYTOP <-  median_vec_ALLBiovol2 
      
      X_vector_gCind_PHYTOPL <- X_vector_gCind
      
      
      df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))
      
      dfok<- df1
      #dfok <- df1[complete.cases(df1$y1),]
      #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
      
      # lines(  log10(dfok$y1) ~ log10(dfok$x1), 
      #         col = "red" , lwd = 4.5)
      # 
      # dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
      # lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
      #         col = "purple" , lwd = 4.5)
      # 
      # 
      # 
      # 
      
      
      
      
       ###  start of colour NBSS  plot  all data, al taxa -------------
      # (nbss figs for Micronekton and mesop fish -> suppl mat )
     
      # end of figure 12 -------- 
       
      
      
      # TWO NBSS plots
      # Oligotrophic vs Mesotophic  linear  samples, ALL PhYTOPL , n = 311
      
      # Separate PHYTOPL into two goups by Chlorophyl ---------------------
      #  Chlorophyl < 1 mg m-3  vs  Chlorophyl > 1mg m-3
      # 1 mg m-3 = 0.3 1og10 (1+x)
      # 1 mg m-3 =  0 1og10 ()
      # oligotrophic from mestrophic
      # only for the 311 top quality PHYTOPL samples with top linear models
      
      phytopl.All.Atl.ok5 <- phytopl.All.Atl.p.rsq.slope_ratio.ok
      dim(phytopl.All.Atl.ok5)
      # 311 top quality samples with top linear models
      
      phytopl.All.Atl.ok5.oligotr  <- subset (phytopl.All.Atl.ok5, phytopl.All.Atl.ok5$Chlorophyll__mg_m_3_insitu < (1) )
      dim(phytopl.All.Atl.ok5.oligotr)
      # 214 oligotrophic top linear  samples
      
      phytopl.All.Atl.ok5.mesotr  <- subset (phytopl.All.Atl.ok5, phytopl.All.Atl.ok5$Chlorophyll__mg_m_3_insitu > (1) )
      dim(phytopl.All.Atl.ok5.mesotr)
      # 95 mesotrophic top linear  samples
      
      
      
      ######## 1. Oligotrophic Phytoplankton ########
      ############ NBSS PLOT ########
      # NBSS FIG. PICOPL. OK (p <0.05) COLD Waters ( SST in situ < 18.4 degrees calsius) only, n = 15 samples --------------
      
      P.DATA <- phytopl.All.Atl.ok5.oligotr # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED
      
      dim(phytopl.All.Atl.ok5.oligotr) # 214 samples
      # plot all points and rlm moldes (blue lines), and red median line  ----------
      
      
      phytopl.All.Atl.PICO.p_OK.MATRIX <- as.matrix(P.DATA[,37:91])
      
      lnumb <- 3
      
      mod1<-  lm( log10(as.numeric(P.DATA[lnumb, 38:91])) ~ log10(X_vector_gCind_PHYTOPL[38:91]))
      abline(mod1) 
      abline(v = -12) 
      
      plot (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind), 
            main = "ATLANTIC, Nano- and Microphytopl., n = 214, Oligotrophic waters",
            xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
            ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
            xlim = c(-15, -5), ylim = c(0, 14) ,
            pch = 16, col = "white")
      
      # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
      
      abline(a = -3.2915, b = -0.9548 ,  
             lwd = 2.5, lty = 2, col = "grey44" )
      
      
      # phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727b[,38:92])
      # dim(phytopl.NBSS.matrix.unfiltr727)
      # 
      # 
      # for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
      # {    
      #   points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
      #           pch = 16, col = alpha("lightblue", 0.5))
      # }
      
      names(P.DATA[lnumb, 37:91])
      
      for (lnumb in 1 : nrow(P.DATA) ) 
      { 
        # abline( rlm(log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind)), 
        #         col =alpha ("cornsilk4", 0.1))
        
        mod1<-  rlm( log10(as.numeric(P.DATA[lnumb, 37:91])) ~ log10(X_vector_gCind_PHYTOPL))
#                                                                     ,  na.action = "na.omit")
        abline( mod1, col =alpha ("cornsilk4", 0.2))
        
      }
      
      for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
      {    
        points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
                pch = 16, col = alpha("darkgrey", 0.2))
      }
      
      
      for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
      {    
        points (log10(phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind) ,
                pch = 16, col = alpha("pink", 0.7))
      }
      
      
      for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
      {    
        points (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,6:55]) ~ log10(X_vector_gCind[6:55]) ,
                pch = 16, col = alpha("forestgreen", 0.7))
      }
      
      abline (v = -12)
      
      
      # Carbon  (458 high-quality datasets, OK)
      (median( P.DATA$Pico.rr.slopes) )
      dim(P.DATA)
      summary(P.DATA$Pico.rr.slopes)
      
      dim(phytopl.All.Atl.PICO.p_OK.MATRIX)
      
      # matrix to to vector (transpose),  -----------------------  
      x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) )
      y1 <- c ( t(phytopl.All.Atl.PICO.p_OK.MATRIX))
      
      # regression line PICOPL, Rob. regr.
      
      med.slope.rr.PICO <- median(P.DATA$Pico.rr.slopes)
      med.interc.rr.PICO <- median(P.DATA$Pico.rr.interc)
      
      abline( a = med.interc.rr.PICO, b = med.slope.rr.PICO, lwd = 2.5, lty = 2, col = "forestgreen")
      
      # regression line (Rob. regr.)
      #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
      #points(  log10(y1) ~ x1 )
      # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
      # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
      #        lwd = 2.5, lty = 2, col = "darkorange" )
      # 

            abline(mod5<- rlm (log10(y1) ~ x1),  
             lwd = 2.5, lty = 2, col = "darkorange" )
      
               
            summary(mod5)
            # b =  -1.11
            
            # clean (remove NA)
            df6<- data.frame( x1 = x1 , y1 = y1)
            
            df6noNA <- na.omit(df6)
            
            
            
            # 95% CI for slope
            fun.bootstr.rr.slopeLOGY(  df6noNA$x1 , df6noNA$y1 , 200)
            #fun.bootstr.rr.slopeLOGY(  df6noNA$x1 , df6noNA$y1 , 20000)
            # 2.5%       50%     97.5% 
            # -1.167412 -1.113337 -1.059926 
            # 
            # "p" value
            permuco::aovperm(  rlm (log10(df6noNA$y1) ~ df6noNA$x1) )
            # p  2e-04
            
            text(x = -12, y = 2, "b = -1.11, 95%CI: -1.17 to -1.06", col = "darkorange"  )
            #text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )
            
            
            
            
                        
      
      # means and medians without considering NAs  
      mean (c(4,NA, 18,34)) # NA
      median( c(4,NA, 18, 34)  ) #NA
      
      mean_   <- function(...) mean(..., na.rm=T)
      median_ <- function(...) median(..., na.rm=T)
      
      mean_ (c(4,NA, 8, 34)) #15.33, OK
      median_ (c(4,NA, 8,34)) #8, OK
      
      
      # calculate median (or mean) only if less than 50% are NA ------------
      
      vec1 <- c(4, 5,NA,NA, NA,  34)
      
      sum(is.na(vec1))
      percNA <- sum(is.na(vec1))/ length(vec1)
      
      
      mean_NA50   <- function(...)  {
        
        percNA <- sum(is.na(...))/ length(...)
        
        if(percNA <= 0.5) {
          mean(..., na.rm=T) }
        
        else{ ;NA }
      }
      
      median_NA50   <- function(...)  {
        
        percNA <- sum(is.na(...))/ length(...)
        
        if(percNA <= 0.5) {
          median(..., na.rm=T) }
        
        else{ ;NA }
      }
      
      vec2 <- c(4, 5,NA,NA, NA,NA,  34)
      vec3 <- c(4, 5,NA,34)
      
      mean_NA50(vec2)
      mean_NA50(vec3)
      
      
      # apply median by columns
      
      median_vec_ALLBiovol <- apply( phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_) 
      mean_vec_ALLBiovol <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, mean_) 
      
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)
      
      
      # apply mean by columns (only if less thn 50% NAs)
      
      # mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 
      
      median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_NA50) 
      
      fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
      
      max_vec_ALL <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, fun95_quant ) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
      
      
      
      median_vec_PHYTOP <-  median_vec_ALLBiovol2 
      
      X_vector_gCind_PHYTOPL <- X_vector_gCind
      
      
      df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))
      
      dfok<- df1
      #dfok <- df1[complete.cases(df1$y1),]
      #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
      
      lines(  log10(dfok$y1) ~ log10(dfok$x1), 
              col = "red" , lwd = 4.5)
      
      dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
      lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
              col = "purple" , lwd = 4.5)
    
      
      # Add legend -----------
      
      leg.names <- c(   "Oligotr. Phytopl.", "All Taxa and Regions"  ) 
      col.taxa  <- c(  "darkorange" , "grey44"  ) 
      
      legend("topright", leg.names, lty = 2,lwd = 2.5,
             col = col.taxa)
      
      
        
      
      ############ END OF NBSS PLOT (oligotrohic Phytoplankton) ########
      
      
      
      ######## 2. Mesotrophic Phytoplankton ########
      ############ NBSS PLOT ########
       
      P.DATA <- phytopl.All.Atl.ok5.mesotr # DEFINE OBJECT (DATA.FRAME) TO BE ANALYZED AND PLOTTED
      
      dim(phytopl.All.Atl.ok5.mesotr) # 95 samples
      # plot all points and rlm moldes (blue lines), and red median line  ----------
      
      
      phytopl.All.Atl.PICO.p_OK.MATRIX <- as.matrix(P.DATA[,37:91])
      
      lnumb <- 3
      
      mod1<-  lm( log10(as.numeric(P.DATA[lnumb, 38:91])) ~ log10(X_vector_gCind_PHYTOPL[38:91]))
      abline(mod1) 
      abline(v = -12) 
      
      plot (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind), 
            main = "ATLANTIC, Nano- and Microphytopl., n = 95, Mesotrophic waters",
            xlab = "log10(Carbon Biomass (mgC ind. -1)) ",
            ylab = "log10(Normalized Carbon Biomass (mgC m-3 mgC-1) ",
            xlim = c(-15, -5), ylim = c(0, 14) ,
            pch = 16, col = "white")
      
      # Naive rlm model with ALL Data and All Taxa ,  a = -3.292,  b= -0.955 
      
      abline(a = -3.2915, b = -0.9548 ,  
             lwd = 2.5, lty = 2, col = "grey44" )
      
      
      # phytopl.NBSS.matrix.unfiltr727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727b[,38:92])
      # dim(phytopl.NBSS.matrix.unfiltr727)
      # 
      # 
      # for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
      # {    
      #   points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
      #           pch = 16, col = alpha("lightblue", 0.5))
      # }
      
      names(P.DATA[lnumb, 37:91])
      
      for (lnumb in 1 : nrow(P.DATA) ) 
      { 
        # abline( rlm(log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind)), 
        #         col =alpha ("cornsilk4", 0.1))
        
        mod1<-  rlm( log10(as.numeric(P.DATA[lnumb, 37:91])) ~ log10(X_vector_gCind_PHYTOPL))
        #                                                                     ,  na.action = "na.omit")
        abline( mod1, col =alpha ("cornsilk4", 0.2))
        
      }
      
      for (lnumb in 1 : nrow(phytopl.NBSS.matrix.unfiltr727) ) 
      {    
        points (log10(phytopl.NBSS.matrix.unfiltr727[lnumb,]) ~ log10(X_vector_gCind) ,
                pch = 16, col = alpha("darkgrey", 0.2))
      }
      
      
      for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
      {    
        points (log10(phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,]) ~ log10(X_vector_gCind) ,
                pch = 16, col = alpha("pink", 0.7))
      }
      
      
      for (lnumb in 1 : nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) ) 
      {    
        points (log10( phytopl.All.Atl.PICO.p_OK.MATRIX[lnumb,6:55]) ~ log10(X_vector_gCind[6:55]) ,
                pch = 16, col = alpha("forestgreen", 0.7))
      }
      
      abline (v = -12)
      
      
      # Carbon  (458 high-quality datasets, OK)
      (median( P.DATA$Pico.rr.slopes) )
      dim(P.DATA)
      summary(P.DATA$Pico.rr.slopes)
      
      dim(phytopl.All.Atl.PICO.p_OK.MATRIX)
      
      # matrix to to vector (transpose),  -----------------------  
      x1 <- rep( (log10(X_vector_gCind)), nrow(phytopl.All.Atl.PICO.p_OK.MATRIX) )
      y1 <- c ( t(phytopl.All.Atl.PICO.p_OK.MATRIX))
      
      # regression line PICOPL, Rob. regr.
      
      med.slope.rr.PICO <- median(P.DATA$Pico.rr.slopes)
      med.interc.rr.PICO <- median(P.DATA$Pico.rr.interc)
      
      abline( a = med.interc.rr.PICO, b = med.slope.rr.PICO, lwd = 2.5, lty = 2, col = "forestgreen")
      
      # regression line (Rob. regr.)
      #plot(log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14) )
      #points(  log10(y1) ~ x1 )
      # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)))
      # abline(rlm (log10(y1) ~ x1, xlim = c(-11, -1), ylim = c(0, 14)), 
      #        lwd = 2.5, lty = 2, col = "darkorange" )
      # 
      
      abline(mod105<- rlm (log10(y1) ~ x1),  
             lwd = 2.5, lty = 2, col = "darkorange" )
      
      summary(mod105)
      # b =  -0.9888
      
      # clean (remove NA)
      df6<- data.frame( x1 = x1 , y1 = y1)
      
      df6noNA <- na.omit(df6)
      
      
      
     # 95% CI for slope
             fun.bootstr.rr.slopeLOGY(  df6noNA$x1 , df6noNA$y1 , 200)
             #fun.bootstr.rr.slopeLOGY(  df6noNA$x1 , df6noNA$y1 , 20000)
             # 2.5%        50%      97.5% 
             # -1.0364601 -0.9884978 -0.9426193         
             
             
      # "p" value
      permuco::aovperm(  rlm (log10(df6noNA$y1) ~ df6noNA$x1) )
      # p  2e-04
      
      
      text(x = -12, y = 2, "b = -0.99, 95%CI: -1.03 to -0.94", col = "darkorange"  )
      #text(x = -10, y = 0.5, "based on all green dots, not sample-by-sample", col = "darkorange"  )
      
      
      
      # means and medians without considering NAs  
      mean (c(4,NA, 18,34)) # NA
      median( c(4,NA, 18, 34)  ) #NA
      
      mean_   <- function(...) mean(..., na.rm=T)
      median_ <- function(...) median(..., na.rm=T)
      
      mean_ (c(4,NA, 8, 34)) #15.33, OK
      median_ (c(4,NA, 8,34)) #8, OK
      
      
      # calculate median (or mean) only if less than 50% are NA ------------
      
      vec1 <- c(4, 5,NA,NA, NA,  34)
      
      sum(is.na(vec1))
      percNA <- sum(is.na(vec1))/ length(vec1)
      
      
      mean_NA50   <- function(...)  {
        
        percNA <- sum(is.na(...))/ length(...)
        
        if(percNA <= 0.5) {
          mean(..., na.rm=T) }
        
        else{ ;NA }
      }
      
      median_NA50   <- function(...)  {
        
        percNA <- sum(is.na(...))/ length(...)
        
        if(percNA <= 0.5) {
          median(..., na.rm=T) }
        
        else{ ;NA }
      }
      
      vec2 <- c(4, 5,NA,NA, NA,NA,  34)
      vec3 <- c(4, 5,NA,34)
      
      mean_NA50(vec2)
      mean_NA50(vec3)
      
      
      # apply median by columns
      
      median_vec_ALLBiovol <- apply( phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_) 
      mean_vec_ALLBiovol <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, mean_) 
      
      
      # lines(  log10(median_vec_ALLBiovol) ~ log10(X_vector_gCind),pch = 16, col = "red" , lwd = 2.5)
      
      
      # apply mean by columns (only if less thn 50% NAs)
      
      # mean_vec_ALLBiovol2 <- apply(phytopl.NBSS.matrix.ok3, 2, mean_NA50) 
      
      median_vec_ALLBiovol2 <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, median_NA50) 
      
      fun95_quant <- function(x) {quantile(x, probs = 0.95, na.rm =  TRUE) }
      
      max_vec_ALL <- apply(phytopl.All.Atl.PICO.p_OK.MATRIX, 2, fun95_quant ) 
      
      length(median_vec_ALLBiovol2)
      
      class(median_vec_ALLBiovol2)
      
      median_vec_ALLBiovol2 <- as.vector(unlist(median_vec_ALLBiovol2))
      
      #mean_vec_ALLBiovol2[mean_vec_ALLBiovol2 == "NULL"] <- NA
      
      
      
      median_vec_PHYTOP <-  median_vec_ALLBiovol2 
      
      X_vector_gCind_PHYTOPL <- X_vector_gCind
      
      
      df1 <- data.frame (y1 = ((median_vec_PHYTOP)) ,x1 =   (X_vector_gCind))
      
      dfok<- df1
      #dfok <- df1[complete.cases(df1$y1),]
      #plot(log10(dfok$y1) ~ log10(dfok$y1),  type = "l")
      
      lines(  log10(dfok$y1) ~ log10(dfok$x1), 
              col = "red" , lwd = 4.5)
      
      dfmax <- data.frame (y1 = ((max_vec_ALL)) ,x1 =   (X_vector_gCind))
      lines(  log10(dfmax$y1) ~ log10(dfmax$x1), 
              col = "purple" , lwd = 4.5)
      
      # Add legend -----------
      
      leg.names <- c(   "Mesotr. Phytopl.", "All Taxa and Regions"  ) 
      col.taxa  <- c(  "darkorange" , "grey44"  ) 
      
      legend("topright", leg.names, lty = 2,lwd = 2.5,
             col = col.taxa)
      
      
      
      ############ END OF NBSS PLOT (Mesotrophic Phytoplankton) ########
      
      
      
      
          
      #### Fig xx scheme (Discussion & Graphical Absract) ----------------
      
      ##########################
      #######################
      ####################
      # SUMMARY tables (p values 95%CI, etc) 
      
      
      
      # Function for complete RR linear model analysis
      # As complete function, with "p" value (SLOPE, rr), and 95CI%, and intercept 95%CI, and lm summaary -----------
      
      
      
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI <- function (x, y , n.runsBOOT = 1000, nruns.aovperm = 2000) 
      {
        
        
        df4 <- data.frame (x = x , y = y) 
        
        
        
        # confidence interval for INTERCEPT:
        
        # Bootstrap function ------------
        # creating vectors to store bootstrap regression coefficients
        # As fuuncion 
        fun.bootstr.rr.intercept.SIMPLE <- function (x,y, nruns) {
          
          # B = 2500 # n of runs
          B <- nruns # n of runs
          
          df <- data.frame (x = x , y = y) 
          df <- na.omit(df) # delete rows with NA
          
          boot.out <- rep(NA, B)
          vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
          ### LOOP ###
          for(i in 1:B){ #starting loop
            
            ##creating samples of observation IDs with replacement, of same size    as original sample
            boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
            
            #matching response and explanatory variable values to bootstrap sample   observation IDs
            ysam <- df$y[boot.id]
            xsam <- df$x[boot.id]
            
            #generating bootstrap SLR model for each bootstrap sample
            boot.mod <- rlm(  ysam ~  xsam )       
            
            #storing regression coefficient values for each bootstrap SLR
            boot.out[i] <- coef(boot.mod)[1]
          }
          #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
          (medianboot <-median(boot.out))
          # -0.883927
          ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
          # 2.5%      97.5% 
          #   -1.0279427 -0.7463216 
          
          
          ;  boot.CI 
        }
        
        
        
        ### 2. WITH  LOG10(Y)  --------------
        # confidence interval for slope:
        
        # Bootstrap function ------------
        # creating vectors to store bootstrap regression coefficients
        # As fuuncion 
        fun.bootstr.rr.slopeY <- function (x,y, nruns) {
          
          # B = 2500 # n of runs
          B <- nruns # n of runs
          
          df <- data.frame (x = x , y = y) 
          df <- na.omit(df) # delete rows with NA
          
          boot.out <- rep(NA, B)
          vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
          ### LOOP ###
          for(i in 1:B){ #starting loop
            
            ##creating samples of observation IDs with replacement, of same size    as original sample
            boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
            
            #matching response and explanatory variable values to bootstrap sample   observation IDs
            ysam <- df$y[boot.id]
            xsam <- df$x[boot.id]
            
            #generating bootstrap SLR model for each bootstrap sample
            boot.mod <-rlm(  ysam ~  xsam )       
            
            #storing regression coefficient values for each bootstrap SLR
            boot.out[i] <- coef(boot.mod)[2]
          }
          #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
          (medianboot <-median(boot.out))
          # -0.883927
          ( boot.CI <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
          # 2.5%      97.5% 
          #   -1.0279427 -0.7463216 
          
          
          ;  boot.CI 
        }
        
        
        
        
        # confidence interval for slope:
        
        # Bootstrap function ------------
        # creating vectors to store bootstrap regression coefficients
        # As fuuncion 
        fun.bootstr.rr.slope.SIMPLE <- function (x,y, nruns) {
          
          # B = 2500 # n of runs
          B <- nruns # n of runs
          
          df <- data.frame (x = x , y = y) 
          df <- na.omit(df) # delete rows with NA
          
          boot.out <- rep(NA, B)
          vector.id <- 1:length(df$x)   # vector of observation IDs (indices)
          ### LOOP ###
          for(i in 1:B){ #starting loop
            
            ##creating samples of observation IDs with replacement, of same size    as original sample
            boot.id <- sample(vector.id, length(df$y), replace=TRUE) 
            
            #matching response and explanatory variable values to bootstrap sample   observation IDs
            ysam <- df$y[boot.id]
            xsam <- df$x[boot.id]
            
            #generating bootstrap SLR model for each bootstrap sample
            boot.mod <-rlm( ( ysam) ~  xsam )       
            
            #storing regression coefficient values for each bootstrap SLR
            boot.out[i] <- coef(boot.mod)[2]
          }
          #determining 2.5% and 97.5% quantiles for both bootstrap regression coefficient samples and displaying
          (medianboot <-median(boot.out))
          # -0.883927
          ( boot.CI.slope <- quantile(boot.out, c(0.025, 0.5,  0.975), type = 2) )
          # 2.5%      97.5% 
          #   -1.0279427 -0.7463216 
          ;boot.CI.slope
          
        }
        
        
        
        boot.CI.intercept <-  fun.bootstr.rr.intercept.SIMPLE(df4$x, df4$y, n.runsBOOT)
        boot.CI.intercept <- round(boot.CI.intercept, 3)
        
        boot.CI.slope <-  fun.bootstr.rr.slope.SIMPLE(df4$x, df4$y, n.runsBOOT)
        boot.CI.slope <- round(boot.CI.slope, 3)
        
        
        
        AIC_res_rlm <- AIC (rlm( df4$y ~df4$x))
        
        summary_lm <- summary(lm( (df4$y) ~df4$x))
        
        summary_rlm <- summary(rlm( (df4$y) ~df4$x))
        
        
        N <- length( (df4$y))
        
        N_used <- length( (df4$y[ !is.na (df4$y)]))
        
        p_val_aovperm.rlm <- permuco::aovperm(  rlm( y ~ x )  , np = nruns.aovperm)
        
        # Calculate R-sq (rlm)
        y_pred.rlm1 <- predict(rlm1 <- rlm(y ~ x))
        rss <- sum((y - y_pred.rlm1)^2)
        tss <- sum((y - mean(y))^2)
        R_squared.rlm1 <- 1 - rss / tss
        
        
        interc.rlm.med.ci <-  paste("a =" , boot.CI.intercept[2] ,  "(95%CI:" , boot.CI.intercept[1],  "to" , boot.CI.intercept[3], ")"   )
        
        
        slope.rlm.med.ci <-   paste("b =" ,  boot.CI.slope[2],     "(95%CI:" , boot.CI.slope[1],      "to" , boot.CI.slope[3], ")"   )
        
        
        res_text<- ("AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,
,   boot.CI.intercept, boot.CI.slope")
        
        res <- list( AIC.rlm = AIC_res_rlm, summary.lm = summary_lm ,   summary.rlm = summary_rlm , 
                     p_val_aovperm.rlm = p_val_aovperm.rlm, N_used =  N_used, 
                     R_squared.rlm1 = R_squared.rlm1  ,  boot.CI.intercept= boot.CI.intercept,
                     boot.CI.slope = boot.CI.slope, paste(res_text), interc.rlm.med.ci, slope.rlm.med.ci )
        
        ; res
        
      }
      
      
      
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI (x= 1:10, y= c(1000,3,8,8,110,20,23,27,30,31), 
                                                        200, 2000)
      
      
      plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
            ylim = c(-10, 13),
            xlab = "log10(gC ind.-1)",
            ylab = "log10(gC m-3 / gC ind.-1)",
            main= "NBSS, All Atlantic, All data, n = 2840, C units")
      
      
      for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
        
        points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
                pch = 16, col = alpha("darkgrey", 0.1))
      }
      
      
      ######## TABLES FOR PAPER  ########
      
      # Table 1:  Summary of NBSS linear model analysis, based on pooled data. Results shown are from robust regression analysis with bootstrap (for 95%CI) and permutation test (for "p" values). 	--------------									
      #        ALL in CARBON NBSS units										
      #            N samples, 	n bins	"p" R2	"a" (95%CI)	"b" (95%CI)
      # 1. Phytoplankton (ALL sizes, all regions )										
      # 2. Phytoplankton (Nano- to Microphytopl.)										
      # 3. Phytoplankton (Picopl. only)										
      # 4. Net-caught Zoopl.										
      # 5. UVP Zoopl.										
      # 6. Micronekton										
      # 7.  Mesopelagic fishes										
      # 8. All Taxa	1,644	19.424	p < 0.0001	0.86		b = -0.955 (95%CI: -0.960 to -0.949).	pooled bins, RR	19.422		The overall RR NBSS linear model (all data, all taxa) slope was close to -1 ( b = -0.955, bootstrap 95% CI = -0.960 to -0.949, n = 19.424 bins, R²: 0.86, permutation p: < 0.0001)
      
      
      # 1. Phytoplankton (ALL sizes, all regions ) (pooled data)	---------	
      
       phytopl.All.Atl.MATRIX.unfiltr.nolims727 <- as.matrix(phytopl.All.Atl.unfiltr.nolims727[,38:92])
      # 
       dim(phytopl.All.Atl.MATRIX.unfiltr.nolims727) # 727
      # 
      # 
       x <- x1 <- rep( logCvector, nrow(phytopl.All.Atl.MATRIX.unfiltr.nolims727) )
       y <- y1 <- c ( t(phytopl.All.Atl.MATRIX.unfiltr.nolims727))
      # 
      df1 <- data.frame(x1, y1)
      names(df1)
      # abline(rlm( log10(y) ~ x))
      # 
             df2<- na.omit(df1)
       dim(df1)
       dim(df2)
      # 
      summary(df1)
       summary(df2)
       m1 <- rlm ( log10(df2$y1) ~ df2$x1)
       abline( m1, lty = 2, lwd = 3, col = "darkgreen")
       summary(m1)
      # # (Intercept)   -5.2800  
      # # df2$x1        -1.1460    
      
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 25 )
      # fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 
      #                                                    20000, 20000)
      # 
      #fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 
      #                                                                                                              2000, 2000)
      # $AIC.rlm
      # [1] 23697.21
      # 
      # $summary.lm
      # 
      # Call:
      #   lm(formula = (df4$y) ~ df4$x)
      # 
      # Residuals:
      #   Min      1Q  Median      3Q     Max 
      # -4.6341 -0.6480  0.1056  0.8778  3.9611 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) -5.08261    0.11010  -46.16   <2e-16 ***
      #   df4$x       -1.11894    0.01036 -108.02   <2e-16 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 1.252 on 7194 degrees of freedom
      # Multiple R-squared:  0.6186,	Adjusted R-squared:  0.6186 
      # F-statistic: 1.167e+04 on 1 and 7194 DF,  p-value: < 2.2e-16
      # 
      # 
      # $summary.rlm
      # 
      # Call: rlm(formula = (df4$y) ~ df4$x)
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -4.73812 -0.70835  0.01713  0.79651  3.92227 
      # 
      # Coefficients:
      #   Value     Std. Error t value  
      # (Intercept)   -5.2800    0.1022   -51.6789
      # df4$x         -1.1460    0.0096  -119.2206
      # 
      # Residual standard error: 1.122 on 7194 degrees of freedom
      # 
      # $p_val_aovperm.rlm
      # Anova Table
      # Resampling test using freedman_lane to handle nuisance variables and 2000 permutations.
      # SS   df     F parametric P(>F) resampled P(>F)
      # x         18279    1 11669                0           5e-04
      # Residuals 11269 7194                                       
      # 
      # $N_used
      # [1] 7196
      # 
      # $R_squared.rlm1
      # [1] 0.6163762
      # 
      # $boot.CI.intercept
      # 2.5%    50%  97.5% 
      # -5.470 -5.283 -5.084 
      # 
      # $boot.CI.slope
      # 2.5%    50%  97.5% 
      # -1.163 -1.146 -1.128 
      # 
      # [[9]]
      # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
      # 
      # [[10]]
      # [1] "a = -5.283 (95%CI: -5.47 to -5.084 )"
      # 
      # [[11]]
      # [1] "b = -1.146 (95%CI: -1.163 to -1.128 )"
      # 
      # 
      # 
      
      
      # 2. Phytoplankton (Nano- to Microphytopl.)	(pooled data)	 ----------		
      
      dim(phytopl.All.Atl.ok5) # 311 samples
      names(phytopl.All.Atl.ok5[,37:91])
      df.raw <- phytopl.All.Atl.ok5
      
      mat.raw <- as.matrix(df.raw[,37:91])
      # 
      dim(mat.raw) # 311 samples
      
      x <- x1 <- rep( logCvector, nrow(mat.raw) )
      y <- y1 <- c ( t(mat.raw))
      # 
      df1 <- data.frame(x1, y1)
      names(df1)
      #abline(rlm( log10(y) ~ x))
      # 
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      # 
      summary(df1)
      summary(df2)
      m2 <- rlm ( log10(df2$y1) ~ df2$x1)
      abline( m2, lty = 2, lwd = 3, col = "purple")
      summary(m2)
      # # (Intercept)   -4.2632 
      # # df2$x1      -1.0635    
      
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 200 )
      # fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 
      #                                                    20000, 20000 )
      # 
      # $AIC.rlm
      # [1] 7476.076
      # 
      # $summary.lm
      # 
      # Call:
      #   lm(formula = (df4$y) ~ df4$x)
      # 
      # Residuals:
      #   Min      1Q  Median      3Q     Max 
      # -3.5776 -0.6761 -0.0316  0.7004  3.0833 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) -4.24265    0.20903  -20.30   <2e-16 ***
      #   df4$x       -1.06127    0.02056  -51.62   <2e-16 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.9197 on 2796 degrees of freedom
      # Multiple R-squared:  0.488,	Adjusted R-squared:  0.4878 
      # F-statistic:  2664 on 1 and 2796 DF,  p-value: < 2.2e-16
      # 
      # 
      # $summary.rlm
      # 
      # Call: rlm(formula = (df4$y) ~ df4$x)
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -3.58198 -0.67656 -0.03275  0.70092  3.07889 
      # 
      # Coefficients:
      #   Value    Std. Error t value 
      # (Intercept)  -4.2632   0.2137   -19.9543
      # df4$x        -1.0635   0.0210   -50.6095
      # 
      # Residual standard error: 1.027 on 2796 degrees of freedom
      # 
      # $p_val_aovperm.rlm
      # Anova Table
      # Resampling test using freedman_lane to handle nuisance variables and 2000 permutations.
      # SS   df    F parametric P(>F) resampled P(>F)
      # x         2254    1 2664                0           5e-04
      # Residuals 2365 2796                                      
      # 
      # $N_used
      # [1] 2798
      # 
      # $R_squared.rlm1
      # [1] 0.4879453
      # 
      # $boot.CI.intercept
      # 2.5%    50%  97.5% 
      # -4.718 -4.268 -3.801 
      # 
      # $boot.CI.slope
      # 2.5%    50%  97.5% 
      # -1.113 -1.064 -1.020 
      # 
      # [[9]]
      # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
      # 
      # [[10]]
      # [1] "a = -4.268 (95%CI: -4.718 to -3.801 )"
      # 
      # [[11]]
      # [1] "b = -1.064 (95%CI: -1.113 to -1.02 )"
      # 
      # >  5e-04
      # [1] 5e-04
      # > 0.0005
      # [1] 5e-04
      
      
      
      # 3. Phytoplankton (Picopl. only)	(pooled data)	-----------
      
      dim(phytopl.All.Atl.atleast.4.PICO) # 279 samples
      names(phytopl.All.Atl.atleast.4.PICO[,38:92])
     
       df.raw <- phytopl.All.Atl.atleast.4.PICO
      
      mat.raw <- as.matrix(df.raw[,38:92])
      # 
      dim(mat.raw) # 311 samples
      
      x <- x1 <- rep( logCvector, nrow(mat.raw) )
      y <- y1 <- c ( t(mat.raw))
      # 
      df1 <- data.frame(x1, y1)
      names(df1)
      #abline(rlm( log10(y) ~ x))
      # 
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      # 
      summary(df1)
      summary(df2)
      m3 <- rlm ( log10(df2$y1) ~ df2$x1)
      abline( m3, lty = 2, lwd = 3, col = "red")
      summary(m3)
      # # (Intercept)   -6.0481
      # # df2$x1    -1.2119  
      
      
       fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 25 )
   
       # fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1),
       #                                                    20000, 20000 )
       
       # 
       # $AIC.rlm
       # [1] 13185.76
       # 
       # $summary.lm
       # 
       # Call:
       #   lm(formula = (df4$y) ~ df4$x)
       # 
       # Residuals:
       #   Min      1Q  Median      3Q     Max 
       # -4.3193 -0.5914  0.0356  0.8332  3.3505 
       # 
       # Coefficients:
       #   Estimate Std. Error t value Pr(>|t|)    
       # (Intercept) -5.92164    0.12647  -46.82   <2e-16 ***
       #   df4$x       -1.19289    0.01168 -102.11   <2e-16 ***
       #   ---
       #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
       # 
       # Residual standard error: 1.182 on 4147 degrees of freedom
       # Multiple R-squared:  0.7154,	Adjusted R-squared:  0.7154 
       # F-statistic: 1.043e+04 on 1 and 4147 DF,  p-value: < 2.2e-16
       # 
       # 
       # $summary.rlm
       # 
       # Call: rlm(formula = (df4$y) ~ df4$x)
       # Residuals:
       #   Min       1Q   Median       3Q      Max 
       # -4.40453 -0.67616 -0.03455  0.74800  3.31675 
       # 
       # Coefficients:
       #   Value     Std. Error t value  
       # (Intercept)   -6.0481    0.1180   -51.2618
       # df4$x         -1.2119    0.0109  -111.1893
       # 
       # Residual standard error: 1.051 on 4147 degrees of freedom
       # 
       # $p_val_aovperm.rlm
       # Anova Table
       # Resampling test using freedman_lane to handle nuisance variables and 2000 permutations.
       # SS   df     F parametric P(>F) resampled P(>F)
       # x         14565    1 10426                0           5e-04
       # Residuals  5793 4147                                       
       # 
       # $N_used
       # [1] 4149
       # 
       # $R_squared.rlm1
       # [1] 0.7140306
       # 
       # $boot.CI.intercept
       # 2.5%    50%  97.5% 
       # -6.315 -6.051 -5.814 
       # 
       # $boot.CI.slope
       # 2.5%    50%  97.5% 
       # -1.235 -1.212 -1.190 
       # 
       # [[9]]
       # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
       # 
       # [[10]]
       # [1] "a = -6.051 (95%CI: -6.315 to -5.814 )"
       # 
       # [[11]]
       # [1] "b = -1.212 (95%CI: -1.235 to -1.19 )"
       # 
       # 
       
       
       
      # 4. Net-caught Zoopl. (sample-wise data analysis only, no fixed size range)  ----------
      
       summary(zoopl_All_Atl4NET$rob_reg__slope)
       # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
       # -4.16424 -1.08440 -0.84045 -0.86508 -0.46466 -0.02725 
       # >      
       summary( zoopl_All_Atl4NET$rob_reg_intercept)
       length(zoopl_All_Atl4NET$rob_reg_intercept)
       # 117 top-linear samples (p<0.05)
         
       # sample-by-sample analysis 
      #  median of all sample-wise linear RR models 
      
       # Bootstrap function (95% CIs for median and mean)
       fun.boot.median.mean <- function(data, nboot = 10000) {
         medians <- numeric(nboot)
         means <- numeric(nboot)
         
         for (i in 1:nboot) {
           bootsample <- sample(data, replace = TRUE)
           medians[i] <- median(bootsample)
           means[i] <- mean(bootsample)
         }
         
         list(
           mean_ci = quantile(means, probs = c(0.025, 0.975)),
           median_ci = quantile(medians, probs = c(0.025, 0.975))
         )
       }
       
       # Run bootstrap (top-quality linear samples only)
         
       summary(zoopl_All_Atl4NET$rob_reg__slope)
       # Min.     1st Qu.   Median     Mean  3rd Qu.     Max. 
       # -4.16424 -1.08440 -0.84045 -0.86508 -0.46466 -0.02725 
      nboot = 200
     
        res <- fun.boot.median.mean(zoopl_All_Atl4NET$rob_reg__slope, nboot = nboot)
       
       res$median_ci # 95%CI for slope (Net-caught zoopl.) 
       # 2.5%      97.5% 
       #   -0.9225717 -0.7006903  
       # 
       
      # Slope  (b) median =  -0.8404 (95%CI: -0.9225 to -0.700 )
       
         summary( zoopl_All_Atl4NET$rob_reg_intercept)
         #   Min.    1st Qu.   Median     Mean  3rd Qu.     Max. 
         # -15.2226  -3.1585  -1.6380  -1.7897  -0.1874   1.8706 
       
         res <- fun.boot.median.mean(zoopl_All_Atl4NET$rob_reg_intercept, nboot = nboot)
         
         res$median_ci # 95%CI for intercept (Net-caught zoopl.) 
         # 2.5%      97.5% 
         #   -2.5427993 -0.7980364 
       
         # Interc.(a) median =  -1.638 (95%CI:   -2.5428 to -0.798 )
         
         
         
       ## Bootstrap for All Net zoopl. samples (unfiltered)
         
         # length(zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA[!is.na(zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA)]) 
         
        # res <- fun.boot.median.mean(zoopl_All_Atl4.unfiltr_NET$rob_reg__slope, nboot = nboot)
         
         #res$median_ci # 95%CI for slope (Net-caught zoopl.) 
          
         summary( zoopl_All_Atl4.unfiltr_NET$intercept_OSLR_max.to.NA)
         # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
         # -15.2226  -3.1585  -1.6380  -1.7897  -0.1874   1.8706 
         #       
       #  res <- fun.boot.median.mean(zoopl_All_Atl4.unfiltr_NET$rob_reg_intercept, nboot = nboot)
         
         #res$median_ci # 95%CI for intercept (Net-caught zoopl.) 
         # 
         
         summary( zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA)
#          Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -2.1337 -1.1315 -0.9591 -1.0036 -0.8407 -0.2018     773  
          
         
         
         
         
      # 5. UVP Zoopl.	(sample-wise data analysis only , no fixed size range)  -----------
      
           
         # run bootstrap (top-linear models , p < 0.05)
         summary(zoopl_All_Atl4UVP  $rob_reg__slope)
         #   Min.   1st Qu.  Median   Mean 3rd Qu.    Max. 
         # -3.3354 -1.1323 -0.7205 -0.8773 -0.4938 -0.0721 
                    
         length(zoopl_All_Atl4_UVP$rob_reg__slope) 
         # 202 top-quality UVP samples
         
         res <- fun.boot.median.mean(zoopl_All_Atl4UVP$rob_reg__slope, nboot = nboot)
         
         res$median_ci # 95%CI for slope (UVP zoopl.) 
         # 2.5%      97.5% 
         #   -0.7732587 -0.6492971 
         # # 
         
         # Slope  (b) median = -0.7205 (95%CI: -0.7732 to -0.6492)
         
         
         
         summary( zoopl_All_Atl4UVP$rob_reg_intercept)
         # Min.      1st Qu.   Median     Mean  3rd Qu.     Max. 
         # -10.6666  -3.1346  -1.5832  -2.2174  -0.8094   1.7318
         # 
         res <- fun.boot.median.mean(zoopl_All_Atl4UVP$rob_reg_intercept, nboot = nboot)
         
         res$median_ci # 95%CI for intercept (Net-caught zoopl.) 
         # 2.5%     97.5% 
         #   -1.857219 -1.338804 
         
         # Interc.(a) median =   -1.5832 (95%CI:   -1.857219 to -1.338804  )
         
         
         
         # COMPARE UVP  vs NET ZOOPL slopes (unfiltr. ALLL, OLSR)
         wilcox.test (zoopl_All_Atl4.unfiltr_NET$slope_OSLR_max.to.NA,
                      zoopl_All_Atl.unfiltr_UVP$slope_OSLR_max.to.NA)
         # p-value = 9.772e-05, OK
         
        
         # COMPARE UVP  vs NET ZOOPL slopes (filtr. top models, OLSR)
           wilcox.test (zoopl_All_Atl4NET$rob_reg__slope,
                       zoopl_All_Atl4UVP $rob_reg__slope)
         # #  p-value = 0.6836, n.s.
         
         
         # COMPARE UVP  vs NET ZOOPL slopes (filtr. top models, RR)
         #  wilcox.test (zoopl_All_Atl4NET$rob_reg__slope,
         #               zoopl_All_AtlUVP $rob_reg__slope)
         # #  p-value = 0.6836, n.s.
         
         
         
            
      
      # 6. Micronekton (pooled data) -----------
      
       mat.raw <- as.matrix(micronektCLEAN.All.Atl  [,37:91])
       names(micronektCLEAN.All.Atl[,37:91])
       
       dim(mat.raw)
       rowSums(is.na(mat.raw))
       
       mat._no_na <-  mat.raw[rowSums(is.na(mat.raw)) != ncol(mat.raw), ]
       rowSums(is.na(mat._no_na))
       
       dim(mat._no_na)# 86  samples with NBSS data, OK
       
       mat.2 <- mat._no_na
       # 
       dim(mat.2) # 86 samples
       
       x <- x1 <- rep( logCvector, nrow(mat.2) )
       y <- y1 <- c ( t(mat.2))
       # 
       df1 <- data.frame(x1, y1)
       names(df1)
       #abline(rlm( log10(y) ~ x))
       # 
       df2<- na.omit(df1)
       dim(df1)
       dim(df2)
       # 
       summary(df1)
       summary(df2)
       m21 <- rlm ( log10(df2$y1) ~ df2$x1)
       abline( m21, lty = 2, lwd = 3, col = "brown")
       summary(m21)
       # # (Intercept)   -4.006
       # # df2$x1      -0.921 
       
       
       fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 25 )
#       fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 20000, 20000 )
       
      
       # >        fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 20000, 20000 )
       # $AIC.rlm
       # [1] 455.9761
       # 
       # $summary.lm
       # 
       # Call:
       #   lm(formula = (df4$y) ~ df4$x)
       # 
       # Residuals:
       #   Min       1Q   Median       3Q      Max 
       # -3.07488 -0.41857  0.06642  0.44868  1.54437 
       # 
       # Coefficients:
       #   Estimate Std. Error t value Pr(>|t|)    
       # (Intercept) -4.01717    0.06629 -60.596  < 2e-16 ***
       #   df4$x       -0.94692    0.11295  -8.384 3.77e-15 ***
       #   ---
       #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
       # 
       # Residual standard error: 0.5927 on 250 degrees of freedom
       # Multiple R-squared:  0.2194,	Adjusted R-squared:  0.2163 
       # F-statistic: 70.28 on 1 and 250 DF,  p-value: 3.771e-15
       # 
       # 
       # $summary.rlm
       # 
       # Call: rlm(formula = (df4$y) ~ df4$x)
       # Residuals:
       #   Min       1Q   Median       3Q      Max 
       # -3.11549 -0.43385  0.03574  0.43335  1.50376 
       # 
       # Coefficients:
       #   Value    Std. Error t value 
       # (Intercept)  -3.9766   0.0647   -61.5058
       # df4$x        -0.9065   0.1102    -8.2292
       # 
       # Residual standard error: 0.6453 on 250 degrees of freedom
       # 
       # $p_val_aovperm.rlm
       # Anova Table
       # Resampling test using freedman_lane to handle nuisance variables and 20000 permutations.
       # SS  df     F parametric P(>F) resampled P(>F)
       # x         24.69   1 70.28        3.775e-15           5e-05
       # Residuals 87.83 250                                       
       # 
       # $N_used
       # [1] 252
       # 
       # $R_squared.rlm1
       # [1] 0.2180565
       # 
       # $boot.CI.intercept
       # 2.5%    50%  97.5% 
       # -4.113 -3.976 -3.846 
       # 
       # $boot.CI.slope
       # 2.5%    50%  97.5% 
       # -1.139 -0.910 -0.686 
       # 
       # [[9]]
       # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
       # 
       # [[10]]
       # [1] "a = -3.976 (95%CI: -4.113 to -3.846 )"
       # 
       # [[11]]
       # [1] "b = -0.91 (95%CI: -1.139 to -0.686 )"
       # 
       
       
       
      # 7. Mesopelagic fishes	(pooled data) ------------
      
      meso_fi.NBSSmatrix <- as.matrix(mesop_fish.CLEAN.All.Atl[,37:91])
      names(mesop_fish.All.Atl[,37:91])
     
      dim(meso_fi.NBSSmatrix)
      rowSums(is.na(meso_fi.NBSSmatrix))
      
      meso_fi.NBSSmatrix_no_na <-  meso_fi.NBSSmatrix[rowSums(is.na(meso_fi.NBSSmatrix)) != ncol(meso_fi.NBSSmatrix), ]
      rowSums(is.na(meso_fi.NBSSmatrix_no_na))
      
      dim(meso_fi.NBSSmatrix_no_na)# 86  samples with NBSS data, OK
     
      mat.raw <- meso_fi.NBSSmatrix_no_na
      # 
      dim(mat.raw) # 86 samples
      
      x <- x1 <- rep( logCvector, nrow(mat.raw) )
      y <- y1 <- c ( t(mat.raw))
      # 
      df1 <- data.frame(x1, y1)
      names(df1)
      #abline(rlm( log10(y) ~ x))
      # 
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      # 
      summary(df1)
      summary(df2)
      m22 <- rlm ( log10(df2$y1) ~ df2$x1)
      abline( m22, lty = 2, lwd = 3, col = "darkorange")
      summary(m22)
      # # (Intercept)   -4.7560 
      # # df2$x1       -0.8857   
      
      
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 25 )
      # 
      # fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1,
      #                                                    y = log10(df2$y1), 20000, 20000 )
      # 
      # 
      # >       fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1,
      #                                                            +                                                          y = log10(df2$y1), 20000, 20000 )
      # $AIC.rlm
      # [1] 890.4751
      # 
      # $summary.lm
      # 
      # Call:
      #   lm(formula = (df4$y) ~ df4$x)
      # 
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -2.14394 -0.40901  0.07778  0.47032  1.35095 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) -4.76259    0.04794  -99.34   <2e-16 ***
      #   df4$x       -0.84749    0.06547  -12.95   <2e-16 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.6664 on 435 degrees of freedom
      # Multiple R-squared:  0.2781,	Adjusted R-squared:  0.2764 
      # F-statistic: 167.6 on 1 and 435 DF,  p-value: < 2.2e-16
      # 
      # 
      # $summary.rlm
      # 
      # Call: rlm(formula = (df4$y) ~ df4$x)
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -2.19654 -0.42786  0.04817  0.44648  1.33218 
      # 
      # Coefficients:
      #   Value    Std. Error t value 
      # (Intercept)  -4.7560   0.0504   -94.4117
      # df4$x        -0.8857   0.0688   -12.8748
      # 
      # Residual standard error: 0.6499 on 435 degrees of freedom
      # 
      # $p_val_aovperm.rlm
      # Anova Table
      # Resampling test using freedman_lane to handle nuisance variables and 20000 permutations.
      # SS  df     F parametric P(>F) resampled P(>F)
      # x          74.42   1 167.6                0           5e-05
      # Residuals 193.16 435                                       
      # 
      # $N_used
      # [1] 437
      # 
      # $R_squared.rlm1
      # [1] 0.2763088
      # 
      # $boot.CI.intercept
      # 2.5%    50%  97.5% 
      # -4.847 -4.756 -4.661 
      # 
      # $boot.CI.slope
      # 2.5%    50%  97.5% 
      # -1.026 -0.885 -0.743 
      # 
      # [[9]]
      # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
      # 
      # [[10]]
      # [1] "a = -4.756 (95%CI: -4.847 to -4.661 )"
      # 
      # [[11]]
      # [1] "b = -0.885 (95%CI: -1.026 to -0.743 )"
      # 
      
      
        
      # 8. All Taxa -----------	
      #1,644	19.424	p < 0.0001	0.86		b = -0.955 (95%CI: -0.960 to -0.949).	pooled bins, RR	19.422		The overall RR NBSS linear model (all data, all taxa) slope was close to -1 ( b = -0.955, bootstrap 95% CI = -0.960 to -0.949, n = 19.424 bins, R²: 0.86, permutation p: < 0.0001)
      
      #  Naive linear model (All data): ------------------
      # matrix to vector (transpose), TSWA -----------------------  
      library (MASS)
      
      x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
      y1 <- c ( t(Data.All.Atl.NBSSmatrix))
      
      df1 <- data.frame(x1, y1, ylog = log10(y1))
      names(df1)
      
      # -INF to NA!!!
      
      # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
      
      df1[df1=="-Inf"]<-NA
      
      df2<- na.omit(df1)
      dim(df1)
      dim(df2)
      
      summary(df1)
      summary(df2)
      
      mean(df2$y1)
      
      summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.954
      summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
      # slope = -0.955
      
      # ALL samples and TAXA (naive analysis)
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 25, 25 )
#      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 50000, 50000 )

      
      
      # 
      # >       fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df2$x1, y = log10(df2$y1), 5000, 5000 )
      # $AIC.rlm
      # [1] 67498.87
      # 
      # $summary.lm
      # 
      # Call:
      #   lm(formula = (df4$y) ~ df4$x)
      # 
      # Residuals:
      #   Min      1Q  Median      3Q     Max 
      # -6.9150 -0.7542  0.1337  0.9793  4.3516 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) -3.377955   0.020157  -167.6   <2e-16 ***
      #   df4$x       -0.953726   0.002761  -345.5   <2e-16 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 1.372 on 19422 degrees of freedom
      # Multiple R-squared:  0.8601,	Adjusted R-squared:   0.86 
      # F-statistic: 1.194e+05 on 1 and 19422 DF,  p-value: < 2.2e-16
      # 
      # 
      # $summary.rlm
      # 
      # Call: rlm(formula = (df4$y) ~ df4$x)
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -7.00813 -0.84655  0.04045  0.88630  4.26449 
      # 
      # Coefficients:
      #   Value     Std. Error t value  
      # (Intercept)   -3.2915    0.0191  -172.5388
      # df4$x         -0.9548    0.0026  -365.4549
      # 
      # Residual standard error: 1.292 on 19422 degrees of freedom
      # 
      # $p_val_aovperm.rlm
      # Anova Table
      # Resampling test using freedman_lane to handle nuisance variables and 5000 permutations.
      # SS    df      F parametric P(>F) resampled P(>F)
      # x         224640     1 119363                0           2e-04
      # Residuals  36552 19422                                        
      # 
      # $N_used
      # [1] 19424
      # 
      # $R_squared.rlm1
      # [1] 0.8594104
      # 
      # $boot.CI.intercept
      # 2.5%       50%     97.5% 
      # -3.334752 -3.291729 -3.247749 
      # 
      # $boot.CI.slope
      # 2.5%        50%      97.5% 
      # -0.9605444 -0.9548234 -0.9491516 
      # 
      # [[9]]
      # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
      # 
      
      
      ########################################      
      ########### Table 2 ###############################
      ########################################      
      # 
      # Table 2: Relationships between abiotic drivers (SST and CHl a) and NBSS parameters (NBSS RR linear model slope and  NBSS median biomass). Analyses based on Sampl-by-sample data (one NBSS anlysis fpr eacch sample).  Linear moodel: Robust regression (RR) with subsequent bootstrap (for 95%CI) and permutation tests (for "p" values). ------------------				
      # ALL in CARBON NBSS	units
      
       # Table 2a : Chl a vs NBSS slope. Univariate robust regression linear model, Temperature (SST, °C, in situ) only.----------
     			
     
      
       # Table 2.a.1. Phytoplankton (Nano- to Microphytopl.), Chla 	-------------		
      
      #1a.   THE BEST Best univarate model,log10(1+ Chl a in situ ) ----------
      lmlogChla_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                               log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
      summary(lmlogChla_insitu) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
      # CARBON (new):  R-squared: 0.0285,  n = 311, p-value: 0.0029
      plot(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
             log10(1+ phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu), 
           col = alpha ("darkgreen", 0.3), pch = 16)
      abline(lmlogChla_insitu, lwd = 2.5, lty = 2, col = "darkorange")
      
      dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
      
      
      cor.test (phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ,
                log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu))
      # Pearson, p-value = 0.002865
   
      
      cor.test (method = "spearman", phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ,
                log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu))
      # spearman, p-value = 0.01624
      
      lmlogChla_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                               log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
      rlmlogChla_insitu <- rlm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                                 log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
      
      summary(rlmlogChla_insitu) #
      summary(lmlogChla_insitu) #
      #R-squared   lm:  0.0285
      
      library(permuco)
      aovperm(rlmlogChla_insitu ,np = 2000)
      #aovperm(rlmlogChla_insitu ,np = 20000)
      # p = 0.0026
      
      
      df <- data.frame (x = log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) , 
                        y = phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes) 
      df <- na.omit(df) # delete rows with NA
      
      # Calculate R-sq (rlm)
      y <- df$y; x <- df$x
      y_pred.rlm1 <- predict(rlm1 <- rlm(y ~ x))
      rss <- sum((y - y_pred.rlm1)^2)
      tss <- sum((y - mean(y))^2)
      (R_squared.rlm1 <- 1 - rss / tss)
      # #R-squared   rlm:  0.02515
      
      # bootrap 95%CI for slope and intercept
      
      
      # fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df$x, 
      #                                                    y = (df$y), 2500, 2500 )
      # 
      fun.bootstr.rr.summary.simple.complete.p.R2.95CI(  x  =  df$x, 
                                                         y = (df$y), 2500, 2500 )
      
       #                                                                                                                  y = (df$y), 2500, 2500 )
      # $AIC.rlm
      # [1] 95.48717
      # 
      # $summary.lm
      # 
      # Call:
      #   lm(formula = (df4$y) ~ df4$x)
      # 
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -1.19944 -0.15749  0.03222  0.19335  0.67920 
      # 
      # Coefficients:
      #   Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) -1.10773    0.02423 -45.719  < 2e-16 ***
      #   df4$x        0.20281    0.06747   3.006  0.00286 ** 
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.28 on 308 degrees of freedom
      # Multiple R-squared:  0.0285,	Adjusted R-squared:  0.02535 
      # F-statistic: 9.036 on 1 and 308 DF,  p-value: 0.002865
      # 
      # 
      # $summary.rlm
      # 
      # Call: rlm(formula = (df4$y) ~ df4$x)
      # Residuals:
      #   Min       1Q   Median       3Q      Max 
      # -1.21078 -0.17423  0.01416  0.17697  0.66248 
      # 
      # Coefficients:
      #   Value    Std. Error t value 
      # (Intercept)  -1.0880   0.0234   -46.5918
      # df4$x         0.1895   0.0650     2.9146
      # 
      # Residual standard error: 0.2617 on 308 degrees of freedom
      # 
      # $p_val_aovperm.rlm
      # Anova Table
      # Resampling test using freedman_lane to handle nuisance variables and 2500 permutations.
      # SS  df     F parametric P(>F) resampled P(>F)
      # x          0.7082   1 9.036         0.002865          0.0028
      # Residuals 24.1413 308                                       
      # 
      # $N_used
      # [1] 310
      # 
      # $R_squared.rlm1
      # [1] 0.02515238
      # 
      # $boot.CI.intercept
      # 2.5%    50%  97.5% 
      # -1.135 -1.088 -1.041 
      # 
      # $boot.CI.slope
      # 2.5%   50% 97.5% 
      # 0.050 0.189 0.334 
      # 
      # [[9]]
      # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
      # 
      # [[10]]
      # [1] "a = -1.088 (95%CI: -1.135 to -1.041 )"
      # 
      # [[11]]
      # [1] "b = 0.189 (95%CI: 0.05 to 0.334 )"
      # 
      # 
      
      
    
      	
      # Table 2.a.2.  Phytoplankton (Picopl. only) vs Chla  -----------
      
      
      # * L.4  only with  significant PICOPL. NBSS models (p <005) * -------------
      
      phytopl.All.Atl.PICO.p_OK <- subset (phytopl.All.Atl.atleast.4.PICO , Pico.olsr.p_val < 0.05)
      dim(phytopl.All.Atl.PICO.p_OK) 
      # 119 samples with at least 4 Picopl. bins and p < 0.05
      
      summary(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes)
      # Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
      # -3.851  -3.330  -2.708  -2.487  -2.171   1.068 
      
      lm27Chla <- lm ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
                           log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu))
 summary(lm27Chla) # p-value: < 2.2e-16   , R-squared:  0.55 
      
 rlm27Chla <- rlm ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
        log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu))
      
 plot ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
          log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu))
 abline(rlm27Chla)     
 
 aovperm(rlm27Chla) #  p = 2e-04

 x = log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu)
 y = phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes
 df27Chla <- data.frame(  x = x , y = y)
 df27Chla <- na.omit(df27Chla)
x = df27Chla$x  ;  y =  df27Chla$y
 
 # Calculate R-sq (rlm)
 y_pred.rlm27Chla <- predict(rlm27Chla <- rlm(y ~ x))
 rss27 <- sum((y - y_pred.rlm27Chla)^2)
 tss27 <- sum((y - mean(y))^2)
 R_squared.rlm27Chla <- 1 - rss27 / tss27
 R_squared.rlm27Chla # 0.54 R-squared RR
 
 
 fun.bootstr.rr.summary.simple.complete.p.R2.95CI (x, y, 200, 2000)
 #fun.bootstr.rr.summary.simple.complete.p.R2.95CI (x, y, 20000, 20000)
#  
#  [[10]]
#  [1] "a = -3.336 (95%CI: -3.519 to -3.16 )"
#  
#  [[11]]
#  [1] "b = 4.205 (95%CI: 3.266 to 5.317 )"
#  
#  $N_used
#  [1] 117
#  
#  $R_squared.rlm1
#  [1] 0.5431391
#       
# #  aovperm
# # p =  2e-04
# #  p = 0.0002
#  
 
      
      
      
      
      
      
      
      # 
      # Table 2b :SST vs NBSS slope. Univariate robust regression linear model, Chl a  (log10(mg m-3, in situ)) only.			-----------
      # CARBON NBSS, R² = xx, p = xx, n = xx 				
      # a	95%CI	b	95%CI
      
      # Phytoplankton (Nano- to Microphytopl.)				
      # Phytoplankton (Picopl. only)				
      # UVP Zoopl.				
      
      # Table 2.b.1.  Phytoplankton (Nano- to Microphytopl.) vs SST -----------
      
 lmSST_insitu <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                          (phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu) )
 summary(lmSST_insitu)
 #p-value: not sign.
 
 
 rlmSST_insitu <- rlm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~ 
                      (phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu) )
 summary(rlmSST_insitu)
 aovperm(rlmSST_insitu)
 #p-value: not sign. (aovperm)
 
 
 
      # Table 2.b.2.  Phytoplankton (Picopl. only) vs SST -----------
      
 
 phytopl.All.Atl.PICO.p_OK <- subset (phytopl.All.Atl.atleast.4.PICO , Pico.olsr.p_val < 0.05)
 dim(phytopl.All.Atl.PICO.p_OK) 
 # 119 samples with at least 4 Picopl. bins and p < 0.05
 
 summary(phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes)
 # Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
 # -3.851  -3.330  -2.708  -2.487  -2.171   1.068 
 
 lm27SST <- lm ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
                    (phytopl.All.Atl.PICO.p_OK$SST_C_insitu))
 summary(lm27SST) # p-value: < 2.2e-16   , R-squared:  0.66 
 
 rlm27SST <- rlm ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
                      (phytopl.All.Atl.PICO.p_OK$SST_C_insitu))
 
 plot ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
          (phytopl.All.Atl.PICO.p_OK$SST_C_insitu))
 abline(rlm27SST)     
 
 aovperm(rlm27SST) #  p = 2e-04
 
 x = (phytopl.All.Atl.PICO.p_OK$SST_C_insitu)
 y = phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes
 df27SST <- data.frame(  x = x , y = y)
 df27SST <- na.omit(df27SST)
 x = df27SST$x  ;  y =  df27SST$y
 
 # Calculate R-sq (rlm)
 y_pred.rlm27SST <- predict(rlm27SST)
 rss27 <- sum((y - y_pred.rlm27SST)^2)
 tss27 <- sum((y - mean(y))^2)
 R_squared.rlm27SST <- 1 - rss27 / tss27
 R_squared.rlm27SST # 0.66 R-squared RR
 
 
 fun.bootstr.rr.summary.simple.complete.p.R2.95CI (x, y, 200, 2000)
 # $AIC.rlm
 # [1] 225.8753
 # 
 # $summary.lm
 # 
 # Call:
 #   lm(formula = (df4$y) ~ df4$x)
 # 
 # Residuals:
 #   Min       1Q   Median       3Q      Max 
 # -1.14018 -0.39538 -0.09765  0.32781  1.81081 
 # 
 # Coefficients:
 #   Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)  1.89226    0.30307   6.244 8.41e-09 ***
 #   df4$x       -0.18023    0.01226 -14.704  < 2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Residual standard error: 0.6511 on 109 degrees of freedom
 # Multiple R-squared:  0.6648,	Adjusted R-squared:  0.6617 
 # F-statistic: 216.2 on 1 and 109 DF,  p-value: < 2.2e-16
 # 
 # 
 # $summary.rlm
 # 
 # Call: rlm(formula = (df4$y) ~ df4$x)
 # Residuals:
 #   Min       1Q   Median       3Q      Max 
 # -1.07781 -0.34890 -0.04297  0.35325  1.96327 
 # 
 # Coefficients:
 #   Value    Std. Error t value 
 # (Intercept)   1.5313   0.2750     5.5684
 # df4$x        -0.1680   0.0111   -15.1069
 # 
 # Residual standard error: 0.5299 on 109 degrees of freedom
 # 
 # $p_val_aovperm.rlm
 # Anova Table
 # Resampling test using freedman_lane to handle nuisance variables and 2000 permutations.
 # SS  df     F parametric P(>F) resampled P(>F)
 # x         91.66   1 216.2                0           5e-04
 # Residuals 46.21 109                                       
 # 
 # $N_used
 # [1] 111
 # 
 # $R_squared.rlm1
 # [1] 0.6583132
 # 
 # $boot.CI.intercept
 # 2.5%   50% 97.5% 
 # 0.943 1.493 2.429 
 # 
 # $boot.CI.slope
 # 2.5%    50%  97.5% 
 # -0.198 -0.167 -0.142 
 # 
 # [[9]]
 # [1] "AIC_res_rlm, summary_lm ,  summary_rlm,  p_val_aovperm, N_used,  R_squared.rlm1 ,\n,   boot.CI.intercept, boot.CI.slope"
 # 
 # [[10]]
 # [1] "a = 1.493 (95%CI: 0.943 to 2.429 )"
 # 
 # [[11]]
 # [1] "b = -0.167 (95%CI: -0.198 to -0.142 )"
 # 
 # 
 
 
 
 # CALCULATING the trophic equilibrium constant "c"
 # b =  c  - m  # Borgmann, 1987, Andersen et al., 2008, Barnes et al., 2010
# b = size spectrum slope (usually in units of abundance,  NNSS, b = -2), 
 # c =  trophic equilibrium constant c =  log(TE) / log (PPMR), 
 # m = 0.75  (metabolic scaling) 
 
 # if
 #b = -1.05 # #(Barnes et al., 2010) and
  # m = 0.75 # then
# c = b + m =  -1.05 + 0.75 #  c =  -0.3 #assumption: b = -2, m = 0.75
# c = b + m = -1 + 0.75 #  c = -0.25 #assumption: b = -1, m = 0.75
# c = b + m = -0.96 + 0.75 #  c = -0.21 #assumption: b = -0.96, m = 0.75
 
  
 # if m = 1 (Picoplankton)
#  b = -1.05 #(Barnes et al., 2010) # and
 # m = 1 # then
 # c = b + m = -1.05 + 1 #  c =  -1.25 #assumption: b = -2
 # c = b + m = -1 + 1 #  c = 0 #assumption: b = -1
 
 # c = -0.25; m = 1; b = c-m  #  c = 0 #assumption: b = -1.25
 # Trophic equilibrium theory predicts steeper NBSS for picoplankton
 
 
  
 # Is the trophy equilibrium constant? log(TE) / log (PPPMR)  ?
 
 
 
 
 
      # Table 2.b.3. UVP Zoopl. vs SST		 ---------------------		
 ### SST vs UVP NBSS Slope, # p-value = 0.02388, spearman rank cor
 
 # For paper
 
 library(ggplot2)
 df2 <- data.frame( NBSS_slope = zoopl_All_Atl4UVP$rob_reg__slope , SST = zoopl_All_Atl4UVP$SST_C_insitu)
 
 ggplot(df2, aes(SST, NBSS_slope ) ) +
   scale_x_continuous(name ="SST (°C)", label = comma,
                      limits=c(12, 30))+
   theme(axis.text.y = element_text(face="bold", color="#993333", 
                                    size=10) ,
         axis.text.x = element_text(face="bold", color="#993333", 
                                    size=10) )+ 
   geom_point( alpha = 0.3, col = "navy", size = 2.5) +
   geom_smooth(method = "lm")
 
 
 
 cor.test( zoopl_All_Atl4UVP$rob_reg__slope , 
           zoopl_All_Atl4UVP$SST_C_insitu, method = "spearman")
 # p-value = 0.02388, weak but significant relationship
 # slighly flaater (more large-sized organisms) in colder waters
 
 
 summary( lm(zoopl_All_Atl4UVP$rob_reg__slope ~ zoopl_All_Atl4UVP$SST_C_insitu))
 # n.s.

  ### SST vs UVP NBSS Slope, # p-value = 0.02388, spearman rank cor
 
 
 
 # FOR PAPER!!---------------
 # UVP PLOT!!! --------------
 # NBSS maximum vs STT ! -----------------
 
 
 plot(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu,
      col = alpha ("navy", 0.3), pch = 16,
      xlab = "SST (°C, in situ)"  ,
      ylab = "ind.biomass at NBSS maximum (log10gC ind-1)"   )
 abline (v = c(20, 27.8), lwd = 2.5, lty = 2, col = "darkgrey")
 cor.test(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) , zoopl_All_Atl4UVP$SST_C_insitu,  
          method = "spearman") #  p-value = 0.009701
 
 # "spearman") #  p-value = 0.00970, pos of maximum vs SST
 ### SST vs UVP NBSS pos of maximum , # p-value = 0.0097, spearman rank cor
 # 196 Data pairs
 
 length(zoopl_All_Atl4UVP$x.value.of.maximum) # 202 UVP profiles (filtered for log-linear sections), 6 without SST data
 zoopl_All_Atl4UVP$SST_C_insitu
 hist(zoopl_All_Atl4UVP$x.value.of.maximum)
 
 summary(( lm(log10(zoopl_All_Atl4UVP$x.value.of.maximum ) ~ zoopl_All_Atl4UVP$SST_C_insitu)))    
 
 
 library(ggplot2)
 df1 <- data.frame( UVP_posMAX = log10(zoopl_All_Atl4UVP$x.value.of.maximum) , SST = zoopl_All_Atl4UVP$SST_C_insitu)
 
 ggplot(df1, aes(SST, UVP_posMAX ) ) +
   geom_point( alpha = 0.3, col = "navy", size = 2.5) +
   geom_smooth(span = 0.7)
 
 ggplot(df1, aes(SST, UVP_posMAX ) ) +
   theme(axis.text.y = element_text(face="bold", color="#993333", 
                                    size=10), 
         axis.text.x = element_text(face="bold", color="#993333", 
                                    size=10))     +
   scale_x_continuous(name ="SST (°C)", label = comma,
                      limits=c(12, 30))+
   scale_y_continuous( name ="ind. biomass at NBSS maximum (log10 gC ind-1) ", 
                       label = comma)+
   geom_point( alpha = 0.3, col = "navy", size = 2.5) +
   geom_smooth(span = 0.7)
 
 ### SST vs UVP NBSS Slope, # p-value = 0.02388, spearman rank cor
 
     
 
 
 
 
      
      
      
      # Table 2c : Multivariate model: SST and Chlorophyll a vs NBSS slope. Multivariate robust regression linear model, SST (°C) and Chl a  (log10(mg m-3, in situ)).----------			
      # CARBON NBSS				
      # Model summary			
      # R2	p	n	full equation
      
      
      # Phytoplankton (Nano- to Microphytopl.)				
   
 
 lmlogChla_SST <- lm(phytopl.All.Atl.p.rsq.slope_ratio.ok$rob_reg__slopes ~
                       phytopl.All.Atl.p.rsq.slope_ratio.ok$SST_C_insitu +
                          log10(1+phytopl.All.Atl.p.rsq.slope_ratio.ok$Chlorophyll__mg_m_3_insitu) )
 summary(lmlogChla_SST) # biovol: R-squared:  0.36, 456 degrees of freedom, p-value: < 2.2e-16
 #p-value: 0.01169
 
 dim(phytopl.All.Atl.p.rsq.slope_ratio.ok)
 
 
 
 # hytopl.All.p.rsq.slo_ratio.ok$SST_C_insi                          0.611    
 # log10(1 + phytopl.All.Atl.p.rsq.slope.ok$Chlorophyll__mg_m_3_in)  0.00497 ** 
 
 
   # Phytoplankton (Picopl. only)				
     
 lm27SST_Chla <- lm ( phytopl.All.Atl.PICO.p_OK$Pico.rr.slopes ~
                   log10(1+phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) +
                   phytopl.All.Atl.PICO.p_OK$SST_C_insitu)
 summary(lm27SST_Chla) # p-value: < 2.2e-16   , R-squared:  0.66 
 
 # 
 # Coefficients:
 #   Estimate Std. Error t value
 # (Intercept)                                                      0.27104    0.43755   0.619
 # log10(1 + phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu)  1.86299    0.38791   4.803
 # phytopl.All.Atl.PICO.p_OK$SST_C_insitu                          -0.12852    0.01561  -8.235
 # 
 # Pr(>|t|)    
 # (Intercept)                                                        0.537    
 # log10(1 + phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_ins) 5.17e-06 ***
 #   phytopl.All.Atl.PICO.p_OK$SST_C_insi                      5.06e-13 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Residual standard error: 0.5975 on 106 degrees of freedom
 # (10 observations deleted due to missingness)
 # Multiple R-squared:  0.7236,	Adjusted R-squared:  0.7184 
 # F-statistic: 138.7 on 2 and 106 DF,  p-value: < 2.2e-16
 # 
      
 
 relaimpo::calc.relimp(lm27SST_Chla)
 
 # Relative importance metrics: 
 #   
 #   lmg
 # log10(1 + phytopl.All.Atl.PICO.p_OK$Chlorophyll__mg_m_3_insitu) 0.3034522
 # phytopl.All.Atl.PICO.p_OK$SST_C_insitu                          0.4201314    
 #      
 # Chlorophyll__mg_m_3_insitu: 30% variab. explained
 # SST_C_insitu : 42%variab. explained
 
 
##################    
###########   
### Naive vs non-naive global linear model ###   ---------   
 
 ###########  
 # Alpha 1 (old) Naive linear model (all data): ------------------
 # matrix to to vector (transpose), TSWA -----------------------  
 library (MASS)
# plot ALL NBSS Data (loop) --------------
  
 # dev.off ()
   library(scales)
 
 Ylab.NBBS <- "log10(gC ind.-1)"
 XlabNBSS  <-  "log10(gC m-3 / gC ind.-1)"
 
 
 lnumb.all <- length(Data.All.Atl$Longitude)
 lnumb.all# 2840 data (lines) in C format
 
 lnumb <- round( runif( n = 1 ,min = 180, max = 2800), 0)
 
 plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
       ylim = c(-10, 13),
       xlab = "log10(gC ind.-1)",
       ylab = "log10(gC m-3 / gC ind.-1)",
       main= "NBSS, All Atlantic, All data, n = 2840, C units")
 
 
 for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
   
   points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
           pch = 16, col = alpha("darkgrey", 0.1))
 }
 
 # First Naive linear model (all data): ------------------
 
 x1 <- rep( logCvector, nrow(Data.All.Atl.NBSSmatrix) )
 y1 <- c ( t(Data.All.Atl.NBSSmatrix))
 
 df1 <- data.frame(x1, y1, ylog = log10(y1))
 names(df1)
 
 # -INF to NA!!!
 
 # df1$ylog <- replace(df1$ylog, "-Inf", NA) 
 
 df1[df1=="-Inf"]<-NA
 
 df2<- na.omit(df1)
 dim(df1)
 dim(df2)
 
 summary(df1)
 summary(df2)
 
 mean(df2$y1)
 
 summary(lm1 <- lm (log10(df2$y1) ~ df2$x1, na.action = NULL))
 # slope = -0.954
 summary(rlm1 <- rlm (log10(df2$y1) ~ df2$x1, na.action = NULL))
 # slope = -0.955
 
 
 abline(rlm1, xlim = c(-11, -1), ylim = c(0, 14), 
        lwd = 2.5, lty = 2, col = "darkorange" )
 
 
 # summary(log10(X_vector_gCind))
 # summary(as.vector((log10(Data.All.Atl.NBSSmatrix))))
 
 
 
 ###########  
 # Alpha - 2 -  non-naive global Model  using only data-rich size classes
 
 # 1. define data-rich size classes -------
 # Rationale 
 # For phytoplankton NBSS, linear models were adjusted using all nano- and microphytoplankton 
 # size bins in the range from -12 to -9 log10(gC cell-1), thus encompassing 10 size bins
 # from -11.74 to -9.03  log10(gC cell-1) ). 
 #   This most data-rich size range was selected to exclude picoplankton (which was not analysed consistently in all 
 #   cruises and regions) and large-sized microplankton (wich were rare in subsampling 
 #   procedures and thus often had zero counts, i.e., empty bins) p
 #   rior to linear model fitting. 
 
# dev.off()
 
 logCvector
 
 
 # Indices to keep
 DRS.no.overlap <- c(7:16, 29:37, 41:45)
 
 # Replace everything else with NA
 logCvector.with.NAs <- rep(NA, length(logCvector))
 logCvector.with.NAs[DRS.no.overlap] <- logCvector[DRS.no.overlap]
 
 logCvector.with.NAs
 
 # dev.off ()
 library(scales)
 
 Ylab.NBBS <- "log10(gC ind.-1)"
 XlabNBSS  <-  "log10(gC m-3 / gC ind.-1)"
 
 
 lnumb.all <- length(Data.All.Atl$Longitude)
 lnumb.all# 2840 data (lines) in C format
 
 lnumb <- 2
 
 plot (log10(Data.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
       ylim = c(-10, 13), col= "white",
       xlab = "log10(gC ind.-1)",
       ylab = "log10(gC m-3 / gC ind.-1)",
       main= "NBSS, All Atlantic, All data, n = 2840, C units")
 
 
 for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
   
   points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
           pch = 16, col = alpha("darkgrey", 0.1))
 }
 
 # 
 # 
 # for (i in 1:nrow(Data.All.Atl.NBSSmatrix)){
 #   
 #   points (log10(Data.All.Atl.NBSSmatrix[i,]) ~ (logCvector.with.NAs),
 #           pch = 16, col = alpha("darkblue", 0.1))
 # }
 # 
 
 
 
 
 logCvector[6:15] #  10 selected bins, phytoplankton
 
 median_vec_PHYTOP <- apply( phytopl.NBSS.matrix.ok3, 2, median_) 
 
 lines(  log10(median_vec_PHYTOP) ~ logX_vector_gCind_PHYTOPL,
         pch = 16, col = "forestgreen" , lwd = 5.5)
 
 
 rect( xleft = logCvector[7],  xright  = logCvector[16], 
        ybottom = 2, ytop = 3, col = "forestgreen", border =  "black")
 
 logCvector[7:16] #  20 selected bins, phytoplankton
 
     
 
 lines(   logmedians5_largeNET ~ logCvector_largeNET,
          pch = 16, col = "dodgerblue" , lwd = 5.5)
 
 rect( xleft = logCvector[30],  xright  = logCvector[37], 
       ybottom = -5, ytop = -4, col = "dodgerblue", border =  "black")
 
 logCvector[30:37] #  8 selected bins, net-caught zooplankton
 
 
 
 lines(  log10(median_vec_ALLUVP.NA_unfiltr) ~ logCvector,
         pch = 16, col = "magenta" , lwd = 5.5)
 
 median_vec_ALLUVP.NA_unfiltr
 
 v1 <- median_vec_ALLUVP.NA_unfiltr
 
 # Position of the first non-NA
 first_non_na <- which(!is.na(v1))[1]
 
 # Position of the last non-NA
 last_non_na  <- tail(which(!is.na(v1)), 1)
 
 first_non_na # 29
 last_non_na # 35
 
 logCvector[29:35] #  6 selected bins, UVP zooplankton
 
 rect( xleft = logCvector[29],  xright  = logCvector[35], 
       ybottom = -3, ytop = -2, col = "magenta", border =  "black")
 
 
 # Mesop. fish, medians,  - ALL Areas -------
 
 
 # all medians (no cutoff selection), but olnly < 50% NA
 lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
          pch = 16, col = "darkorange" , lwd = 5.5)
 
 v1 <- median_vec_ALLUVP.NA_unfiltr
 
 # Position of the first non-NA
 first_non_na <- which(!is.na(v1))[1]
 
 # Position of the last non-NA
 last_non_na  <- tail(which(!is.na(v1)), 1)
 
 first_non_na # 29
 last_non_na # 35
 
 logCvector[41:44] #  8 selected bins, Mesop. fish, 
 
 rect( xleft = logCvector[41],  xright  = logCvector[44], 
       ybottom = -6, ytop = -5, col = "darkorange", border =  "black")
 
 
 
 
 # Micronekt, medians, - ALL Areas -------
 
 
 # all medians (no cutoff selection), but olnly < 50% NA
 lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
          pch = 16, col = "darkred" , lwd = 5.5)
 
 
 
 v1 <- logmedians5_largeMicronekt
 
 # Position of the first non-NA
 first_non_na <- which(!is.na(v1))[1]
 
 # Position of the last non-NA
 last_non_na  <- tail(which(!is.na(v1)), 1)
 
 first_non_na # 29
 last_non_na # 35
 
 logCvector[42:45] #  8 selected bins, Micronekt
 
 rect( xleft = logCvector[42],  xright  = logCvector[45], 
       ybottom = -7, ytop = -6, col = "darkred", border =  "black")
 
 
 
 
 
 
 # complete vector of data-rich size classes
 
 logCvector[7:16] #  20 selected bins, phytoplankton
 
 logCvector[30:37] #  8 selected bins, net-caught zooplankton
 
 logCvector[29:35] #  6 selected bins, UVP zooplankton
 
 logCvector[41:44] #  4 selected bins, Mesop. fish, 
 
 logCvector[42:45] #  4 selected bins, Mesop. Micronekt
 
 DRS <- c( 7:16, 30:37, 29:35,  41:44,  42:45) # data-rich size classes
 DRS.no.overlap <- c( 7:16, 29:37,  41:45) # data-rich size classes, no overlap
 
 
 logCvector[DRS]
 
 logCvector[DRS.no.overlap]
 
 
 # Indices to keep
 DRS.no.overlap <- c(7:16, 29:37, 41:45)
 
 # Replace everything else with NA
 logCvector.with.NAs <- rep(NA, length(logCvector))
 logCvector.with.NAs[DRS.no.overlap] <- logCvector[DRS.no.overlap]
 
 logCvector.with.NAs
 
 
 #[bookm] 
 
 # LM based on medians (data-rich, selected ranges  only) ----------

 
  logmedians5_largeMicronekt ~ logCvector_largeMicronekt
  logCvector[41:44] #  4 selected bins, Mesop. fish, 
  logmedians6lengthOK_largeMicronekt <- rep(NA, length(logCvector))
  logmedians6lengthOK_largeMicronekt[42:45] <- logmedians5_largeMicronekt[1:4]

    # lines(   logmedians5_largeMicronekt ~ logCvector_largeMicronekt,
    #        pch = 16, col = "darkred" , lwd = 5.5)
     

   lines(   logmedians6lengthOK_largeMicronekt ~ logCvector,
           pch = 16, col = "darkred" , lwd = 5.5)
  
 
   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI
  
    logCvector[40:43] #  4 selected bins, Mesop. fish, 
   logmedians6lengthOK_largeMESO_FI <- rep(NA, length(logCvector))
   logmedians6lengthOK_largeMESO_FI[41:44] <- logmedians5_largeMESO_FI[1:4]
  
   
  
   # lines(   logmedians5_largeMESO_FI ~ logCvector_largeMESO_FI,
   #          pch = 16, col = "darkorange" , lwd = 5.5)
   # 
    lines(   logmedians6lengthOK_largeMESO_FI ~ logCvector,
            pch = 16, col = "darkorange" , lwd = 5.5)
   
    
    logmedians5_largeNET
    logCvector[30:37] #  4 selected bins, net-caught zoopl 
    logmedians6lengthOK__largeNET <- rep(NA, length(logCvector)) 
    logmedians6lengthOK__largeNET[30:37] <- logmedians5_largeNET[1:8]
    
    lines(   logmedians6lengthOK__largeNET ~ logCvector,
             pch = 16, col = "dodgerblue" , lwd = 5.5)
    
 
 Xvec = rep(logCvector, 5) 
   Yvec =  c(  log10(median_vec_PHYTOP),
      logmedians6lengthOK__largeNET ,
       log10(median_vec_ALLUVP.NA_unfiltr) , 
                logmedians6lengthOK_largeMicronekt,
 logmedians6lengthOK_largeMESO_FI)
  
    lines(   Yvec ~ Xvec,
            pch = 16, col = "black" , lwd = 5.5)
   
    # LM and RLM based on medians, slope = -1 (both lm and rlm), intercept =  -3.7 for rlm   and 3.5 for LM
   lm.all.medians <-  lm(Yvec ~ Xvec)  # LM based on medians (data-rich, selected ranges  only)
   rlm.all.medians <-  rlm(Yvec ~ Xvec)  # RLM based on medians (data-rich, selected ranges  only)
   
   summary(lm.all.medians) # slope = -1.00419 +- 0.0284, Intercept: -3.533 +- 0.174
   summary(rlm.all.medians)#  -1.0079  +- 0.0213, Intercept:  -3.686 +-  0.1303
   
   
   
   abline(( lm.all.medians ), col = ("red"))
   abline(rlm.all.medians, col = "orange")
   
#   [bookmark]
       
 # LM based on ALL data  (data-rich, selected ranges only) ---------------
 
      # compose a new, complete data matrix with cleaned data -------------
 
 
 #phytopl
   
   phytopl.All.Atl.NBSSmatrix 
   phytopl.All.Atl
   
   
       
   lnumb <- 2
   # 
   # plot (log10(phytopl.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
   #       ylim = c(-10, 13),
   #       xlab = "log10(gC ind.-1)",
   #       ylab = "log10(gC m-3 / gC ind.-1)",
   #       main= "NBSS, All Atlantic, All data, n = 2840, C units")
   # 
   # 
   for (i in 1:nrow(phytopl.All.Atl.NBSSmatrix)){
     
     points (log10(phytopl.All.Atl.NBSSmatrix[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("green", 0.1))
   }
   
   
   # Net-caught zoopl
   
   zoo_NET.9 <- zoo_NET.5
   dim(zoo_NET.9)
   names(zoo_NET.9) # NBSS data from 39:93
   summary(zoo_NET.9)
   
   zoo_NET.9.only.selcted.size.range <- zoo_NET.9
   zoo_NET.9.only.selcted.size.range[,39:93] <- NA 
   zoo_NET.9.only.selcted.size.range[(38+30):(38+37)] <- zoo_NET.9[(38+30):(38+37)] 
   zoo_NET.9.MATRIX <-   as.matrix(zoo_NET.9.only.selcted.size.range[,39:93])
   
   
   
   lnumb <- 2
   # 
   # plot (log10(phytopl.All.Atl.NBSSmatrix[lnumb,]) ~ log10(X_vector_gCind),
   #       ylim = c(-10, 13),
   #       xlab = "log10(gC ind.-1)",
   #       ylab = "log10(gC m-3 / gC ind.-1)",
   #       main= "NBSS, All Atlantic, All data, n = 2840, C units")
   # 
   # 
   for (i in 1:nrow(zoo_NET.9.MATRIX)){
     
     points (log10(zoo_NET.9.MATRIX[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("blue", 0.1))
   }
   
   
   # UVP zoopl
   
   zoo_UVP.5
  
   zoo_UVP.9 <- zoo_UVP.5
   dim(zoo_UVP.9)
   names(zoo_UVP.9) # NBSS data from 39:93
   summary(zoo_UVP.9)
   
   zoo_UVP.9.only.selcted.size.range <- zoo_UVP.9
   zoo_UVP.9.only.selcted.size.range[,39:93] <- NA 
   zoo_UVP.9.only.selcted.size.range[(38+29):(38+35)] <- zoo_UVP.9[(38+29):(38+35)] 
   zoo_UVP.9.MATRIX <-   as.matrix(zoo_UVP.9.only.selcted.size.range[,39:93])
   
   
     
   for (i in 1:nrow(zoo_UVP.9.MATRIX)){
     
     points (log10(zoo_UVP.9.MATRIX[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("magenta", 0.1))
   }
   
   
   
   # 
   
   zoo_UVP.5
   
   zoo_UVP.9 <- zoo_UVP.5
   dim(zoo_UVP.9)
   names(zoo_UVP.9) # NBSS data from 39:93
   summary(zoo_UVP.9)
   
   zoo_UVP.9.only.selcted.size.range <- zoo_UVP.9
   zoo_UVP.9.only.selcted.size.range[,39:93] <- NA 
   zoo_UVP.9.only.selcted.size.range[(38+29):(38+35)] <- zoo_UVP.9[(38+29):(38+35)] 
   zoo_UVP.9.MATRIX <-   as.matrix(zoo_UVP.9.only.selcted.size.range[,39:93])
   
   
   
   for (i in 1:nrow(zoo_UVP.9.MATRIX)){
     
     points (log10(zoo_UVP.9.MATRIX[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("magenta", 0.1))
   }
   
   # Mesop. fish
   
   meso_fi.NBSSmatrix <- as.matrix(mesop_fish.CLEAN.All.Atl[,37:91])
   names(mesop_fish.All.Atl[,37:91])
   
   
   mesop_fi.9 <- mesop_fish.CLEAN.All.Atl
   dim(mesop_fi.9)
   names(mesop_fi.9) # NBSS data from 39:93
   summary(mesop_fi.9)
   
   mesop_fi.9.only.selcted.size.range <- mesop_fi.9
   mesop_fi.9.only.selcted.size.range[,37:91] <- NA 
   mesop_fi.9.only.selcted.size.range[(36+41):(36+44)] <- mesop_fi.9[(36+41):(36+44)] 
   mesop_fi.9.MATRIX <-   as.matrix(mesop_fi.9.only.selcted.size.range[,37:91])
   
   
   
   for (i in 1:nrow(mesop_fi.9.MATRIX)){
     
     points (log10(mesop_fi.9.MATRIX[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("darkorange", 0.1))
   }
   
   
   
   
   # Mesop. Micronekton
   
   
   micronektCLEAN.All.Atl
   
    names(micronektCLEAN.All.Atl[,37:91])
   
   
    micronekt.9 <- micronektCLEAN.All.Atl
   dim(micronekt.9)
   names(micronekt.9) # NBSS data from 39:93
   summary(micronekt.9)
   
   micronekt.9.only.selcted.size.range <- micronekt.9
   micronekt.9.only.selcted.size.range[,37:91] <- NA 
   micronekt.9.only.selcted.size.range[(36+42):(36+45)] <- micronekt.9[(36+42):(36+45)] 
   micronekt.9.MATRIX <-   as.matrix(micronekt.9.only.selcted.size.range[,37:91])
   
   
   
   for (i in 1:nrow(micronekt.9.MATRIX)){
     
     points (log10(micronekt.9.MATRIX[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("darkred", 0.1))
   }
   
   
   
   
   
   
   Data.MATRIX.All.data.rich.all.recombined <- rbind(phytopl.All.Atl.NBSSmatrix,
                                                     zoo_NET.9.MATRIX,
                                                     zoo_UVP.9.MATRIX,
                                                    micronekt.9.MATRIX, 
                                                     mesop_fi.9.MATRIX ) 
   
      
 
   lnumb <- 2
 
#   dev.off()  
   plot (log10(Data.MATRIX.All.data.rich.all.recombined[lnumb,]) ~ log10(X_vector_gCind),
         ylim = c(-10, 13), col= "white",
         xlab = "log10(gC ind.-1)",
         ylab = "log10(gC m-3 / gC ind.-1)",
         main= "NBSS, All Atlantic, All data, n = 7635 bins, C units")
   
   
   for (i in 1:nrow(Data.MATRIX.All.data.rich.all.recombined)){
     
     points (log10(Data.MATRIX.All.data.rich.all.recombined[i,]) ~ log10(X_vector_gCind),
             pch = 16, col = alpha("darkgreen", 0.1))
   }
   
   
   
   log.x1.data.rich.sizes <- rep( logCvector, nrow(Data.MATRIX.All.data.rich.all.recombined) )
   y1.data.rich.sizes <- c ( t(Data.MATRIX.All.data.rich.all.recombined))
   logy1.data.rich.sizes <- log10(y1.data.rich.sizes)
   
   df1.data.rich.sizes <- data.frame(log.x1.data.rich.sizes, logy1.data.rich.sizes)
   names(df1.data.rich.sizes)
   
   lm99.all.data.rich.sizes <- lm( logy1.data.rich.sizes ~ log.x1.data.rich.sizes) 
   abline(lm99.all.data.rich.sizes, col = "red")
   rlm99.all.data.rich.sizes <- rlm( logy1.data.rich.sizes ~ log.x1.data.rich.sizes) 
   abline(rlm99.all.data.rich.sizes, col = "darkorange")
 
   lines(   Yvec ~ Xvec,
            pch = 16, col = "darkgrey" , lwd = 5.5)
   
  
 # med.fun <-  function(x) {quantile(x, probs = 0.5, na.rm =  TRUE) }
  
     fun98_quant <-  function(x) {quantile(x, probs = 0.98, na.rm =  TRUE) }
   
     
   max98_vec_ALL.micronekt.all <- apply(micronekt.9.MATRIX, 2, fun98_quant ) 
   
   lines(   log10(max98_vec_ALL.micronekt.all) ~ log10(X_vector_gCind),
            pch = 16, col = "purple" , lwd = 5.5)
 
   # max_vec_ALL.mesop.fish.all <- apply(mesop_fi.9.MATRIX, 2, fun95_quant ) 
   # 
   # lines(   log10(max_vec_ALL.mesop.fish.all) ~ log10(X_vector_gCind),
   #          pch = 16, col = "lightblue" , lwd = 5.5)
   # 
   
      
   summary(lm99.all.data.rich.sizes) # slope:  -0.924 +- 0.004, Intercept: -3.16 +- 0.033
   summary(rlm99.all.data.rich.sizes)# slope:  -0.9320 +- 0.0038, Intercept:  -3.1438 +-  0.032
   
   
   # Delete columns with zeros
   
   
   # Suppose your data frame is named df
   df <- df[df[[ncol(df)]] != 0, ]
   
   df <- read.csv("~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/Data/Stations overview Atlantic size spectraV7_y.csv", sep=";")
   
   dim(df)
   # View(df)
   
   
   df2 <- df[df[[ncol(df)]] != 0, ]
   
   
   dim(df2)
   
   # View(df2)
   
   summary(df$SOMA.NBSS)
   summary(df2$SOMA.NBSS)
   
   hist(df$SOMA.NBSS)
   hist(df2$SOMA.NBSS)
   
   
   write.table(df2, file = "~/Papers/0000 - Paper_Fish_Zoopl_Size_Spectra_Brazil_Africa/Data/SizespectraV7_no_zeros.csv", sep=";")

   
   
   ### Calculate median biomass of microplankton () as an indicator of trophic state
   
   
   
   
   
      #################################      #################################      #################################
      
      ### END OF SCRIPT ###
      
      #################################      #################################      #################################
      
      
      # Going Back and Forward
      # When navigating through code (especially when navigating through a sequence of function calls) you often want to quickly return to the previous editing location. RStudio maintains a list of active navigations and allows you to traverse them using the Back and Forward commands (available on the Edit menu and on the far left of the source editor toolbar).
      # 
      # Back and Forward apply to the following navigation gestures:
      #   
      #   Opening a document (or switching tabs)
      # Going to a function definition
      # Jumping to a line
      # Jumping to a function using the function menu
      # You can also invoke Back and Forward using the Ctrl+F9/Ctrl+F10 (Cmd+F9/Cmd+F10 on the Mac) keyboard shortcuts.
      
     
      
      
   