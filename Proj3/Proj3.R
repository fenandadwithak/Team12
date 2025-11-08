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
library(splines)

##=========================Data Preparation===================================
engcov <- read.table("engcov.txt", header = TRUE)
y  <- engcov$nhs
t  <- engcov$julian
n  <- length(y)

## Interval distribution from infection to death
d  <- 1:80
edur <- 3.151; sdur <- 0.469
pi_j <- dlnorm(d, edur, sdur)
pi_j <- pi_j / sum(pi_j) 

##============ EVALUATE MATRIX X, X_tilde, and S =============================
evaluate_mat <- function(t_seq, K = 80, ord = 4, poly_deg = 2) {
  # Function to construct matrix X, X_tilde, and S as main component to 
  #    predict f(t)/number of new infections occurring on day t
  
  # Input argument : 1. time/days of the year
  #                  2. number of splines basis
  #                     f(t) = β1(b1)(t) + β1(b1)(t) + ... + βK(bK)(t)
  #                     larger K, smoother line
  #                  3. ord, order spline (cubic spline = 4)
  #                  4. poly_deg, polynomial degree is polynomial for matrix X
  #                     poly_deg = 2
  #                     f(t) = β0 + β1(t) + β2(t^2)
  
  # Output Argument: 1. X = model matrix from polynomial linear model
  #                  2. X_tilde = model matrix from splines
  #                  2. S = penalty matrix
  #         
  
  # Define t min and t max 
  t_min <- min(t_seq)
  t_max <- max(t_seq)
  
  # middle K-2 knots = internal knots, 
  # K+4 evenly spaced knots, at the end, we add each 2 knots before and after internal knots
  middle_knots <- seq(t_min, t_max, length.out = K-2)
  knots <- c(rep(t_min, 2), middle_knots, rep(t_max, 2))
  
  # 2. Xtilde: Basis splines (model matrix)
  #    model f(t) = Xtilde x matrix(Beta)
  Xtilde <- splineDesign(knots, t_seq, ord = ord, outer.ok = TRUE)
  colnames(Xtilde) <- paste0("b", seq_len(ncol(Xtilde)))
  
  # 3. X: Model matrix polynomial 
  #    if we set poly deg 3, the model would be f(t) = Bo + B1(t) + B2(t^2)
  X <- as.matrix(sapply(0:poly_deg, function(p) t_seq^p))
  colnames(X) <- paste0("poly", 0:poly_deg)
  
  # 4. S: Penalty matrix (second difference)
  dK <- diff(diag(K), differences = 2)
  S <- crossprod(dK)
  
  return(list(Xtilde = Xtilde, X = X, S = S, knots = knots))
}

##============ FUNCTION PENALIZED NLL =============================
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