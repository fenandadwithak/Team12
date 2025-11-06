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

##########################################################
#Build cubic B-spline basis over time (X_tilde)
##########################################################

K <- 20     # number of spline basis functions (more = more flexible)
ord <- 4    # spline order (4 = cubic)

# Choose equally spaced knots across time
breaks <- seq(1, Tn, length.out = K - ord + 2)

# Build full knot vector:
#  •	repeat the first and last breaks "ord" times to stabilize behavior at edges
#  •	include interior breaks once
all_knots <- c(
  rep(breaks[1], ord),                       # left boundary replicated 'ord' times
  breaks[2:(length(breaks) - 1)],            # interior knots
  rep(breaks[length(breaks)], ord)           # right boundary replicated 'ord' times
)
stopifnot(length(all_knots) == K + ord)      # standard B-spline requirement
t_grid  <- 1:Tn                               # evaluation grid: time = 1..Tn
X_tilde <- splineDesign(                      # evaluate spline basis at each time
  all_knots, x = t_grid, ord = ord, outer.ok = TRUE
)  # dimension: Tn x K  (rows=time, cols=spline basis functions)
cat("dim(X_tilde):", paste(dim(X_tilde), collapse = " x "), "\n")

##########################################################
#Convolve spline basis with lag kernel to get final design matrix X
##########################################################
X <- matrix(0, nrow = Tn, ncol = K)  # initialize final design (Tn x K)

# For each lag dlag:
# •	weight w = pd[dlag]
# •	shift the "source" rows of X_tilde up by dlag days
# •	add w * shifted X_tilde into the appropriate rows of X

for (dlag in seq_along(pd)) {
  w <- pd[dlag]                       # scalar weight for this lag
  t_idx <- dlag:Tn                    # destination rows that have past data available
  src   <- 1:(Tn - dlag + 1)          # source rows from X_tilde (shifted up by dlag)
  X[t_idx, ] <- X[t_idx, ] + w * X_tilde[src, ]
}

# After the loop, X is a lag-augmented spline design matrix (Tn x K)
# If gamma are coefficients, the linear predictor is eta = X %*% gamma
cat("dim(X):", paste(dim(X), collapse = " x "), "\n")

##########################################################
#Smoothness penalty matrix S0 (second-difference penalty)
##########################################################

# Penalizes roughness in gamma; encourages smooth spline coefficients
D <- matrix(0, nrow = K - 2, ncol = K)        # (K-2) x K operator
for (i in 1:(K - 2)) D[i, i:(i + 2)] <- c(1, -2, 1)
S0 <- t(D) %*% D                               # penalty matrix (K x K), PSD

##########################################################
#Penalized Poisson objective and its gradient
##########################################################

# Model: y_t ~ Poisson(mu_t), log(mu_t) = eta_t = (X %*% gamma)_t
# Objective to MINIMIZE = Negative Log-Likelihood + smoothness penalty
#   NLL = sum(exp(eta) - y*eta)   (dropping constants)
#   Penalty = 0.5 * lambda * gamma' S0 gamma

pen_nll <- function(gamma, X, y, S0, lambda = 1e-1) {
  # gamma: K-vector of spline coefficients we are optimizing
  eta <- as.vector(X %*% gamma)                 # linear predictor (Tn-vector)
  mu  <- exp(eta)                               # Poisson mean (positive)
  nll <- sum(mu - y * eta)                      # Poisson negative log-likelihood
  # IMPORTANT: matrix multiply (%*%), NOT modulus (%%)
  pen <- 0.5 * lambda * as.numeric(t(gamma) %*% (S0 %*% gamma))  # smoothness cost
  nll + pen
}
pen_grad <- function(gamma, X, y, S0, lambda = 1e-1) {
  # Analytic gradient: faster and more accurate than numerical
  eta <- as.vector(X %*% gamma)
  mu  <- exp(eta)
  score <- t(X) %*% (mu - y)                    # gradient of NLL part
  as.vector(score + lambda * (S0 %*% gamma))    # add gradient of penalty
}

# Numerical gradient checker (slow but simple):
# verifies pen_grad matches finite differences -> catches bugs.

fd_grad <- function(fun, gamma, eps = 1e-6, ...) {
  g <- numeric(length(gamma))
  for (j in seq_along(gamma)) {
    e <- rep(0, length(gamma)); e[j] <- eps
    g[j] <- (fun(gamma + e, ...) - fun(gamma - e, ...)) / (2 * eps)
  }
  g
}

##########################################################
#Fit the model using BFGS optimization
##########################################################
set.seed(1)                  # reproducibility (not essential for optim, but fine)
gamma0 <- rep(0, K)          # start from all zeros
lambda <- 1e-1               # smoothing strength: larger -> smoother fitted curve

# Check that our gradient matches finite differences at the start point

g_anal <- pen_grad(gamma0, X, y, S0, lambda)
g_fd   <- fd_grad(pen_nll, gamma0, eps = 1e-6, X = X, y = y, S0 = S0, lambda = lambda)
max_rel_err <- max(abs(g_anal - g_fd) / pmax(1, abs(g_fd)))
cat("Max relative FD gradient error:", signif(max_rel_err, 4), "\n")

# Expect a very small error (e.g., ~1e-6 to 1e-8). If it's big, something is off.
# Use BFGS (a standard quasi-Newton optimizer) to minimize the objective

fit <- optim(
  par     = gamma0,                 # initial gamma
  fn      = pen_nll,                # objective function
  gr      = pen_grad,               # analytic gradient
  method  = "BFGS",                 # efficient for smooth problems
  control = list(reltol = 1e-10,    # tight tolerance (can relax if needed)
                 maxit = 1000,      # cap iterations
                 trace = 0),        # no console trace
  X = X, y = y, S0 = S0, lambda = lambda
)

##########################################################
#Extract fitted quantities
##########################################################

gamma_hat <- fit$par                      # estimated spline coefficients (K-vector)
eta_hat   <- as.vector(X %*% gamma_hat)   # fitted linear predictor (log-scale)
mu_hat    <- exp(eta_hat)                 # fitted mean counts (always > 0)

cat("Converged:", fit$convergence == 0, "\nFinal penalized NLL:", fit$value, "\n")



































































































































































































































































































































































































































































































































































# Fenanda ----

# ==============================================================================
##                          NON-PARAMETRIC BOOTSTRAPING
## =============================================================================

# Initialize number of replicate sample
n_bootstrap = 200
# Initialize matrix to store 



end <- Sys.time()
end-start