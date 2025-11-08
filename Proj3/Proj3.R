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
##                                   OUTLINE
# ==============================================================================
## What you get:
## (1) 
## (2) 
## (3) 
## (4) 
## (5)
## (6)

##============ EVALUATE MATRIX X, Xtilde, and S ===============================
library(splines)
engcov <- read.table("engcov.txt", header = TRUE)
y <- engcov$nhs; t <- engcov$julian; n <- length(y)

d <- 1:80; edur <- 3.151; sdur <- 0.469
pd <- dlnorm(d, edur, sdur)
pd <- pd / sum(pd)

make_matrices <- function(t, K=80) {
  knots <- seq(min(t)-30, max(t), length.out=K+4)
  #define a sequence starting from the first day of death recorded minus 30
  #(min(t)-30) because for the earliest records, deaths might be caused by the
  #infection happened around 30 days ago, until maximum day observed (max(t))
  # for as many as K+4 evenly spaced, so there will be a sequence of 84 numbers
  
  Xtilde <- splineDesign(knots=knots, x=(min(t)-30):max(t), ord=4,
                         outer.ok= TRUE)
  #construct matrix X-tilde using splineDesign with 84 knots as defined before
  #starting from min(t)-30 until max(t) with order of the spline = 4
  
  X <- matrix(0, nrow=length(t), ncol=ncol(Xtilde)) #the model matrix of deaths
  #for every value in length(t) = 150
  for (i in 1:length(t)) {
    j_max <- min(80, 29 + i)
    for (j in 1:j_max) {
      if ((30+i-j) >= 1)
        X[i,] <- X[i,] + Xtilde[30 + i - j,] * pd[j]
    }
  }
  S <- crossprod(diff(diag(ncol(X)), diff=2))
  list(Xtilde=Xtilde, X=X, S=S)
}

mats <- make_matrices(t)
X <- mats$X
S <- mats$S
Xtilde <- mats$Xtilde

##============ FUNCTION PENALIZED NLL ==========================================
# Model: y_t ~ Poisson(mu_t), log(mu_t) = eta_t = (X %*% gamma)_t
# Objective to MINIMIZE = Negative Log-Likelihood + smoothness penalty
#   NLL = sum(exp(eta) - y*eta)   (dropping constants)
#   Penalty = 0.5 * lambda * gamma' S0 gamma

pen_nll <- function(gamma, X, y, S0, lambda = 1e-1) {
  # gamma: K-vector of spline coefficients we are optimizing
  eta <- as.vector(X %*% gamma)               # linear predictor (Tn-vector)
  mu  <- exp(eta)                             # Poisson mean (positive)
  nll <- sum(mu - y * eta)                    # Poisson negative log-likelihood
  # IMPORTANT: matrix multiply (%*%), NOT modulus (%%)
  pen <- 0.5 * lambda * as.numeric(t(gamma) %*% (S0 %*% gamma))# smoothness cost
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
set.seed(1)               # reproducibility (not essential for optim, but fine)
gamma0 <- rep(0, K)       # start from all zeros
lambda <- 1e-1            # smoothing strength: larger -> smoother fitted curve

# Check that our gradient matches finite differences at the start point

g_anal <- pen_grad(gamma0, X, y, S0, lambda)
g_fd  <- fd_grad(pen_nll, gamma0, eps=1e-6, X=X, y=y, S0=S0, lambda=lambda)
max_rel_err <- max(abs(g_anal - g_fd) / pmax(1, abs(g_fd)))
cat("Max relative FD gradient error:", signif(max_rel_err, 4), "\n")

# Expect a very small error (e.g., ~1e-6 to 1e-8). If it's big, something is off
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

gamma_hat <- fit$par               # estimated spline coefficients (K-vector)
eta_hat   <- as.vector(X %*% gamma_hat) # fitted linear predictor (log-scale)
mu_hat    <- exp(eta_hat)  # fitted mean counts (always > 0)

cat("Converged:", fit$convergence==0, "\nFinal penalized NLL:", fit$value, "\n")



































































































































































































































































































































































































































































































































































# Fenanda ----

# ==============================================================================
##                          NON-PARAMETRIC BOOTSTRAPING
## =============================================================================

# Initialize number of replicate sample
n_bootstrap = 200
# Initialize matrix to store 



end <- Sys.time()
end-start