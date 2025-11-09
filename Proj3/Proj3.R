start <- Sys.time()
##=============== PROJECT 3 - EXTENDED STATISTICAL PROGRAMMING =================

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
##                                  OUTLINE
# ==============================================================================
## What you get:
## (1) 
## (2) 
## (3) 
## (4) 
## (5)
## (6)


##===================== Data Preparation & Load library ========================
library(splines) #load library splines

df <- read.table("engcov.txt", header = TRUE) # import datasets
y <- df$nhs; t <- df$julian
n <- length(y) # define total observation 
d <- 1:80; edur <- 3.151; sdur <- .469
pd <- dlnorm(d, edur, sdur) # density (PDF) function of death
pd <- pd / sum(pd) # probability of death

##==================== (1) Construct X, Xtilde, and S ==========================
make_matrices <- function(t, K=80) {
  # Function to construct matrices (Xtilde, X, S) required for predicting 
  #  number of new infections occuring on day t which expressed by an infection 
  #  rate function f(t)
  
  # Input/Arguments :(1) t, a vector of day since 2020. In this dataset, we use 
  #                      julian
  #                  (2) K, number of basis function, set to 80
  #                      larger K would result smoother model
  #
  # Output/Return   :(1) Xtilde, a model matrix (spline basis) to predict f(t) 
  #                      f(t) = b1(t)β1 + b2(t)β2 + ... + bk(t)βk
  #                      row : days; colomn : B-spline basis function
  #                  (2) X, model matrix for the expected death count (the 
  #                      observed outcome), after convolution with Xtilde and 
  #                      π(j)/probability func. for days from infection to death 
  #                  (3) S, penalised matrix. Matrix defining which aspect of β
  #                      are penalised
  
  # Construct X_tilde
  # Define a sequence starting from the first day of death recorded minus 30
  # (min(t)-30) because for the earliest records, deaths might be caused by the
  # infection happened around 30 days ago, until maximum day observed (max(t))
  # for as many as K+4 evenly spaced, so there will be a sequence of 84 numbers
  knots <- seq(min(t)-30, max(t), length.out=K+4)
  
  Xtilde <- splineDesign(knots=knots, 
                         x=(min(t)-30):max(t), 
                         ord=4,
                         outer.ok= TRUE)#construct matrix X-tilde using 
                                        #splineDesign with 84 knots as defined 
                                        #before starting from min(t)-30 until 
                                        #max(t) with order of the spline = 4
  
  # Construct X
  X <- matrix(0, nrow=length(t), ncol=ncol(Xtilde)) # initialize matrix to store
  #for every value in length(t) = 150
  # Execute for each observation of time (1:150)
  for (i in 1:length(t)) {
    j_max <- min(80, 29 + i) ## define upper bound for summation
    for (j in 1:j_max) { ## calculate for each days until death 
      # fill corresponding cells with corresponding element X_tilde and prob 
      # death time function
      X[i,] <- X[i,] + Xtilde[30 + i - j,] * pd[j] 
    }
  }
  
  # Construct S (Penalised Matrix)
  S <- crossprod(diff(diag(ncol(X)), diff=2))
  
  list(Xtilde=Xtilde, X=X, S=S) ## return the output into a list 
}##make_matrices

mats <- make_matrices(t) #
X <- mats$X ## X
S <- mats$S ## S
Xtilde <- mats$Xtilde ## X tilde

##============= (2) Function Penalised NLL and its Objective Func. =============
# Model has the structure of a Poisson GLM, then we need to compute penalised
# function. In this case, we want to define objective function PNLL(β)
# 
# Formula:
# NLL = −ℓ(γ)= ∑[e^(ηi)−yi.ηi] where μi=e^(ηi)=e^(Xi.γ)
# Penalty (P) = (λ/2)β'Sβ
#
# PNLL(β) = NLL(β) + P
#   where NLL = sum(exp(eta) - y*eta) (dropping constants)
#         P =  0.5 * lambda * gamma' S gamma 

pen_nll <- function(gamma, X, y, S, lambda) {
  # Function to compute penalised negative log likehood
  # Input/Argument : (1) Gamma : K-vector of spline coefficients
  #                  (2) y :  the deaths on day of the year ti
  #                  (3) S : penalised matrix
  #                  (4) lambda : smoothing parameter
  
  # Output/Return  : single numeric value of pen_nll
  eta <- as.vector(X %*% gamma)               # linear predictor 
  mu  <- exp(eta)                             # Poisson mean (positive)
  nll <- sum(mu - y * eta)                    # Poisson negative log-likelihood
  pen <- 0.5 * lambda * as.numeric(t(gamma) %*% (S %*% gamma))# penalty
  nll + pen # pnll
}##pen_nll


# Define Gradient vector of Objective Function/ its derivative vector w.r.t γ
pen_grad <- function(gamma, X, y, S, lambda) {
  # Analytic gradient: faster and more accurate than numerical
  eta <- as.vector(X %*% gamma)
  mu  <- exp(eta)
  score <- t(X) %*% (mu - y)                    # gradient of NLL part
  as.vector(score + lambda * (S %*% gamma))     # add gradient of penalty
}

# Checking the derivative (sp notes 74)
K <- 80
fd <- gamma0 <- rep(0, K)       # start from all zeros
lambda <- 5e-5
pen_nll0 <- pen_nll(gamma0, X, y, S, lambda) 
eps <- 1e-7
for (i in 1:length(gamma0)) {
  gamma1 <- gamma0; gamma1[i] <- eps
  pen_nll1 <- pen_nll(gamma1,X,y,S,lambda)
  fd[i] <- (pen_nll1 - pen_nll0)/eps
}
fd; pen_grad(gamma0, X, y, S, lambda) ## already same
range(fd - pen_grad(gamma0, X, y, S, lambda)) ## apx zero

##================ (3) Fit the model using BFGS optimization ===================
gamma0 <- rep(0, K)       # start from all zeros
lambda <- 5e-5            # smoothing strength: larger -> smoother fitted curve

# Use BFGS (a standard quasi-Newton optimizer) to minimize the objective
fit <- optim(
  par     = gamma0,                 # initial gamma
  fn      = pen_nll,                # objective function
  gr      = pen_grad,               # analytic gradient
  method  = "BFGS",                 # efficient for smooth problems
  control = list(maxit = 1000),     # cap iterations
  X=X, y=y, S=S, lambda=lambda
)

gamma_hat <- fit$par # estimated spline coefficients (K-vector)
eta_hat  <- as.vector(X %*% gamma_hat) # fitted linear predictor (log-scale)
mu_hat <- exp(eta_hat)  # fitted mean counts (always > 0)
f_hat <- Xtilde %*% gamma_hat #

fit$convergence # 0 means successful
fit$value #final objective value for reference

plot(t, y, pch=16, col="black",
     xlab="Day of year / Julian Day (Observed)",
     ylab="Daily Deaths (Fitted)",
     main=expression(paste("Observed vs Fitted,  λ=5e-5")))
lines(t, mu_hat, col="red", lwd=2)
legend("topright", legend=c("Observed", "Fitted"),
       col=c("black", "red"), pch=c(16, NA), lty=c(NA,1), bty="n")

##================ (4) Fit the model using BFGS optimization ===================
lambdas <- exp(seq(-13, -7, length=50))
BIC_vals <- numeric(length(lambdas))
best_fit <- NULL

for (i in seq_along(lambdas)) {
  fit <- optim(par=gamma0, fn=pen_nll, gr=pen_grad, method="BFGS",
               X=X, y=y, S=S, lambda=lambdas[i], control=list(maxit=1000))
  beta_hat <- exp(fit$par)
  mu_hat <- as.vector(X %*% beta_hat)
  W <- diag(as.vector(y / mu_hat^2))
  H0 <- t(X) %*% W %*% X
  H_lambda <- H0 + lambdas[i] * S
  EDF <- sum(diag(solve(H_lambda, H0)))
  ll <- sum(y * log(mu_hat) - mu_hat)
  BIC_vals[i] <- -2 * ll + log(length(y)) * EDF
  if (is.null(best_fit) || BIC_vals[i] < min(BIC_vals[1:i])) best_fit <- fit
}

best_lambda <- lambdas[which.min(BIC_vals)]
cat("Best lambda:", best_lambda, "\n")

##======================== (5) Bootstrap Uncertainty ===========================





































































































































































































# Fenanda ----

# ==============================================================================
##                          NON-PARAMETRIC BOOTSTRAPING
## =============================================================================

# Initialize number of replicate sample
n_bootstrap = 200
# Initialize matrix to store 



end <- Sys.time()
end-start