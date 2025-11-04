start <- Sys.time()
# PROJECT 3 - EXTENDED STATISTICAL PROGRAMMING =================================

# Group 12 
# Aseel Alghamdi : S2901228
# Fenanda Dwitha Kurniasari : S2744048
# Nurmawiya : S2822251

# Aseel     :
#
# Fenanda   : 
#
# Nurmawiya : 
#

# Repository Link
# https://github.com/fenandadwithak/Team12/tree/main/Proj3

# ==============================================================================
##                                   Outline
## =============================================================================
## What you get:
## (1) 
## (2) 
## (3) 
## (4) 
## (5)
## (6)

## =============================================================================

# 
# 
#         

##==============================================================================jb
library(splines)
library(stats)
# Open data
#df = read.table('engcov.txt')

##============================Function Matrix===================================
spline_matrix = function (d, K=80) {
  # Function to make spline
  # 1. Define internal and full knot sequences
  # Internal knots (K - 2)
  internal_knots <- seq(min(d), max(d), length.out = K-2)
  
  # Full knot sequence with 2 extra boundary knots at each end (K + 4 total)
  full_knots <- c(rep(min(d), 4), internal_knots, rep(max(d), 4))
  
  # 2. Construct B-spline basis matrix (X_tilde)
  X_tilde <- splineDesign(knots = full_knots, x = d, outer.ok = TRUE)
  
  # 3. Create X = X_tilde (if no penalty transformation)
  X <- X_tilde
  
  # 4. Construct S, the penalty matrix
  #    S penalizes the roughness of the spline (e.g., 2nd derivative)
  #    Here we use finite differences to approximate the roughness penalty
  S <- crossprod(diff(diag(K),diff=2))
  
  # 5. Return everything as a list
  list(X_tilde = X_tilde, X = X, S = S, knots = full_knots)
}

# Testing
spline_matrix(t)

#########################STEP 1 & 2#################################
#Load data and basic setup
##########################################################

eng <- read.table("engcov.txt", header = TRUE)  # expects a column named "deaths"
y <- eng$deaths                                # response variable (counts per time unit)
Tn <- length(y)                                # number of time points (e.g., days)

# y: numeric vector of length Tn
# Tn: integer, number of observations

##########################################################
#Set up distributed-lag kernel (log-normal shape)
##########################################################

# The kernel defines how much influence past values have on today's predictor
d <- 1:80                   # allow up to 80 days of lag
edur <- 3.151               # log-normal mean (on log scale)
sdur <- 0.469               # log-normal standard deviation (on log scale)
pd <- dlnorm(d, meanlog = edur, sdlog = sdur)  # weight for each lag day
pd <- pd / sum(pd)          # normalize so weights sum to 1

# d: integer vector from 1 to 80
# pd: numeric vector of length 80, positive and sums to 1
















































































































































































































































































































































































































































































































































































# Fenanda ----

# ==============================================================================
##                          NON-PARAMETRIC BOOTSTRAPING
## =============================================================================

# Initialize number of replicate sample
n_bootstrap = 200
# Initialize matrix to store 



end <- Sys.time()
end-start