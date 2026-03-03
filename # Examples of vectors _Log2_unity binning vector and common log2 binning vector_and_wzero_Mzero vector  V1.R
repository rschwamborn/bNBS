
# Examples of binning vectors

# Content: -----------
# 1. Intuitive standard log2 vector (not centered at the geometric mean M = 1) ----------
# 2. "Log2_unity vector" (center of the vector:  M_zero = 1, log (M_zero) = 0) ----------
# 3. Log(golden_ratio+1)_double_unity vector (center of the vector: log (M_zero) = (0), log (w_zero) = 0) ----------


# 1. Intuitive standard log2 vector (not centered at the geometric mean M = 1) ----------

# Standard log2-binning vector (breaks for binning) ---------------
# length of the vector: 80 elements 
# n <- 80
half <- floor(80 / 2)
exponents <- seq(-half, half - 1)   # symmetric around 0
log2vec <-   2^exponents  # the standard log2-binning vector, these are the breaks!
breaks <- log2vec

breaks [38:44]
# 0.125 0.250 0.500 1.000 2.000 4.000 8.000 16.000 # intuitive, simple breaks

##########################



# 2 "Log2_unity vector" (at the center of the vector:  M_zero = 1,  log(M_zero) = 0 ) ----------
# geometric progression with ratio  r=2.


generate_log2_unity_bins <- function(n_bins = 9) {
  # Central bin
  i0 <- ceiling(n_bins / 2)
  
  # Midpoints: powers of 2 around central bin
  midpoints <- 2^((1:n_bins) - i0)
  
  # Boundaries: geometric mean of consecutive midpoints
  boundaries <- numeric(n_bins + 1)
  for(i in 1:n_bins) {
    boundaries[i + 1] <- sqrt(midpoints[i] * midpoints[min(i + 1, n_bins)])
  }
  boundaries[1] <- midpoints[1] / sqrt(2) # lower edge of first bin
  
  list(
    midpoints = midpoints,
    boundaries = boundaries,
    central_bin_index = i0
  )
}

# Example usage
bins <- generate_log2_unity_bins(n_bins = 80)
bins$midpoints
bins$boundaries

length(bins$midpoints) #80

length(bins$boundaries) #81

# Inspect central bin
i0 <- bins$central_bin_index
cat("Central bin midpoint:", bins$midpoints[i0], "\n")
cat("Central bin edges:", bins$boundaries[i0], "-", bins$boundaries[i0 + 1], "\n")
cat("Central bin linear width:", bins$boundaries[i0 + 1] - bins$boundaries[i0], "\n")


breaks <- bins$boundaries


breaks [38:44]
#   0.176  0.353  0.707  1.414  2.828  5.656 11.313 # intuitive, simple breaks, centered at M = 1

# bins are log2-scaled,  OK

# geometric mean of i0
x <- as.vector(breaks [40:41])
# Geometric mean
geo_mean_zero <- exp(mean(log(x)))
geo_mean_zero # 1, OK

w_zero <- diff((breaks [40:41]))
w_zero # is 0.7071068, OK< but should ideally be w = 1 (double unity)


breaks[46] / breaks[45] # log2-scaled,OK ("octave-scaled")



# 3.Log(golden_ratio+1)-double-unity vector (at the center of the vector: logM = (0), log (w) = 0) ----------
# geometric progression with ratio  r=  1 + 1.618... (golden ratio + 1).
# logM = (0), log (w) = 0

generate_logGR1_unity_bins <- function(n_bins = 80, w0 = NULL) {
  # n_bins: total number of bins
  # w0: optional linear width of central bin (if NULL, use geometric ratio)
  
  # Central bin index
  i0 <- ceiling(n_bins / 2)
  
  # Number of bins above and below central bin
  n_up <- n_bins - i0
  n_down <- i0 - 1
  
  # Step 1: Choose ratio r between midpoints
  # Default: golden ratio (approx 1.618) to make midpoint geometric with reasonable spread
  r <- (1 + sqrt(5)) / 2  # default ratio
  if (!is.null(w0)) {
    # Solve for r so that central bin width = w0
    # width = sqrt(M0 * r) - sqrt(M0 / r) = sqrt(r) - 1/sqrt(r)
    # Solve: sqrt(r) - 1/sqrt(r) = w0
    f <- function(x) sqrt(x) - 1/sqrt(x) - w0
    r <- uniroot(f, lower = 1, upper = 10)$root
  }
  
  # Step 2: Generate midpoints
  midpoints <- numeric(n_bins)
  midpoints[i0] <- 1  # central bin midpoint
  
  # bins above central
  for (i in (i0 + 1):n_bins) {
    midpoints[i] <- midpoints[i - 1] * r
  }
  
  # bins below central
  for (i in (i0 - 1):1) {
    midpoints[i] <- midpoints[i + 1] / r
  }
  
  # Step 3: Compute boundaries from midpoints
  boundaries <- numeric(n_bins + 1)
  boundaries[1] <- midpoints[1] / sqrt(r)  # lower edge of first bin
  for (i in 1:n_bins) {
    boundaries[i + 1] <- sqrt(midpoints[i] * midpoints[min(i + 1, n_bins)])
  }
  
  # Step 4: Return as a list
  list(
    midpoints = midpoints,
    boundaries = boundaries,
    ratio = r,
    central_bin_index = i0
  )
}

# Example usage:
bins <- generate_logGR1_unity_bins(n_bins = 80, w0 = 1)

length(bins$midpoints) #80

length(bins$boundaries) #81

# Inspect central bin
i0 <- bins$central_bin_index
cat("Central bin midpoint:", bins$midpoints[i0], "\n")
cat("Central bin edges:", bins$boundaries[i0], "-", bins$boundaries[i0 + 1], "\n")
cat("Central bin linear width:", bins$boundaries[i0 + 1] - bins$boundaries[i0], "\n")

breaks <- bins$boundaries

breaks[38:44]

breaks[45] / breaks[44]#  2.618 scaled,OK ( 1 + golden ratio)



# geometric mean of i0
x <- as.vector(breaks [40:41])
# Geometric mean
geo_mean_zero <- exp(mean(log(x)))
geo_mean_zero # 1, OK

w_zero <- diff((breaks [40:41]))
w_zero #  1,  OK




# 4. #  another way to obtain the golden-ratio+1  - scaled vector with w_zero = 1, and M_zero = 1
# To make central bin width w_zero  = 1 while keeping midpoint = 1, you must increase the ratio r → midpoints are no longer powers of 2.
# not log-2 scaled 

generate_.B.logGR1_unity_bins <- function(n_bins = 9, w0_target = 1) {
  i0 <- ceiling(n_bins / 2)
  
  # Step 1: solve ratio r_new for central bin width = w0_target
  # central bin midpoint M0 = 1
  M0 <- 1
  # solve x^2 - (w0/M0) x - 1 = 0
  a <- 1
  b <- -w0_target / M0
  c <- -1
  x <- (-b + sqrt(b^2 - 4*a*c)) / (2*a) # positive root
  r_new <- x^2
  
  # Step 2: midpoints: powers of r_new around center
  midpoints <- numeric(n_bins)
  midpoints[i0] <- M0
  for(i in (i0+1):n_bins) midpoints[i] <- midpoints[i-1] * r_new
  for(i in (i0-1):1) midpoints[i] <- midpoints[i+1] / r_new
  
  # Step 3: boundaries from geometric mean
  boundaries <- numeric(n_bins+1)
  for(i in 1:n_bins) boundaries[i+1] <- sqrt(midpoints[i] * midpoints[min(i+1,n_bins)])
  boundaries[1] <- midpoints[1] / sqrt(r_new)
  
  list(
    midpoints = midpoints,
    boundaries = boundaries,
    central_bin_index = i0,
    r = r_new
  )
}

# Example
bins <- generate_.B.logGR1_unity_bins(n_bins = 80, w0_target = 1)

i0 <- bins$central_bin_index
cat("Central bin midpoint:", bins$midpoints[i0], "\n")
cat("Central bin edges:", bins$boundaries[i0], "-", bins$boundaries[i0 + 1], "\n")
cat("Central bin width:", bins$boundaries[i0 + 1] - bins$boundaries[i0], "\n")
cat("Ratio r used:", bins$r, "\n")



length(bins$midpoints) #80

length(bins$boundaries) #81

# Inspect central bin
i0 <- bins$central_bin_index
cat("Central bin midpoint:", bins$midpoints[i0], "\n")
cat("Central bin edges:", bins$boundaries[i0], "-", bins$boundaries[i0 + 1], "\n")
cat("Central bin linear width:", bins$boundaries[i0 + 1] - bins$boundaries[i0], "\n")


breaks <- bins$boundaries


breaks [38:44]
#    0.09016994  0.23606798  0.61803399  1.61803399  4.23606798 11.09016994 29.03444185
# bins are NOT log2-scaled!

# geometric mean of i0
x <- as.vector(breaks [40:41])
# Geometric mean
geo_mean_zero <- exp(mean(log(x)))
geo_mean_zero # 1, OK

w_zero <- diff((breaks [40:41]))
w_zero # is 0.7071068, should be 1

breaks [38:44]


breaks[45] / breaks[44]

# 2.618 scaled,OK ( 1 + golden ratio)


