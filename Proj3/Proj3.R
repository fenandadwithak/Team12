start <- Sys.time()
##=============== PROJECT 3 - EXTENDED STATISTICAL PROGRAMMING =================

# Group 12 
# Aseel Alghamdi : S2901228
# Fenanda Dwitha Kurniasari : S2744048
# Nurmawiya : S2822251

# Aseel     : Make function to contruct X, X tilde, and S (make_matrice)
#             Make objective function and its derivative
#             Revise and put comment on code
# Fenanda   : Bootstrap Uncertainty
#             Make final plot
#             Revise and put comment on code
# Nurmawiya : Preliminary Sanity Check with initial lambda, making plot
#             Fit the model using BFGS
#             Revise and put comment on code

# Repository Link
# https://github.com/fenandadwithak/Team12/tree/main/Proj3

# ==============================================================================
##                                  OUTLINE
# ==============================================================================
# This script implements a full deconvolution model to infer the daily number
# of new Covid-19 infections in England during 2020 from observed hospital
# death counts. The work follows these steps:
#
# 1. Construct the matrices required for the model: 
#       - X̃ : the B-spline basis matrix used to represent f(t)
#       - X  : the corresponding model matrix for expected deaths, obtained by
#              convolving X̃ with the infection-to-death delay distribution
#       - S  : the second-difference penalty matrix controlling smoothness.
#
# 2. Write functions for the penalised negative log-likelihood and its gradient
#    with respect to γ, where β = exp(γ) ensures positivity. Verify the gradient
#    numerically using finite differences.
#
# 3. Fit the model with a fixed smoothing parameter λ = 5e-5 as an initial
#    diagnostic check, and plot the observed deaths, fitted deaths, and the
#    corresponding infection curve f(t) to confirm sensible behaviour.
#
# 4. Select the smoothing parameter by minimising the Bayesian Information
#    Criterion (BIC). A grid of λ values is explored over log-scale
#    seq(-13, -7, length = 50), and the effective degrees of freedom are
#    computed from the Hessian at each fit.
#
# 5. Quantify uncertainty in the estimated infection curve using
#    non-parametric bootstrap resampling. Bootstrap weights are generated and
#    the model is refitted for 200 iterations to obtain a distribution of f(t).
#
# 6. Produce the final visual outputs: 
#       - observed vs fitted hospital deaths, and 
#       - the estimated infection trajectory f(t) with 95% bootstrap intervals

##===================== Data Preparation & Load library ========================
library(splines) #provides functions for constructing B-spline basis matrices
library(ggplot2) #used later for plotting fitted deaths and infection curves

df <- read.table("engcov.txt", header = TRUE) #load covid-19 deaths datasets
y <- df$nhs #observed deaths on each day
t <- df$julian #corresponding day-of-year in 2020
n <- length(y) #total observation

#probability distribution for the delay from infection to death (1–80 days)
d <- 1:80; edur <- 3.151 #mean (on log-scale) of log-normal delay distribution
sdur <- .469 #sd (on log-scale)
pd <- dlnorm(d, edur, sdur) #unnormalised probabilities
pd <- pd / sum(pd) #normalise to sum to 1

##==================== (1) Construct X, Xtilde, and S ==========================
make_matrices <- function(t, K=80) {
  # Function to construct the matrices required for the deconvolution model:
  #   - Xtilde : B-spline basis for representing f(t)
  #   - X      : model matrix mapping beta to expected deaths (via convolution)
  #   - S      : second-difference penalty matrix for smoothing beta
  #
  # Arguments:
  #   t : vector of observation days (Julian day of 2020)
  #   K : number of B-spline basis functions (default 80)
  #
  # Returns:
  #   A list containing Xtilde, X, and S.
  
  ## ----- Construct Xtilde (spline basis for f(t)) -----
  # Internal knots run from (min death day - 30) to max death day.
  # The subtraction of 30 allows f(t) to include infections up to ~30 days
  # before the first recorded death.
  
  # internal knots: K-2 points spanning the modelling range of f(t)
  internal_knots <- seq(min(t) - 30, max(t), length.out = K - 2)
  step <- diff(internal_knots[2:3]) #calculating space between adjacent knots
  #constructing full knots sequence up to total of K + 4 = 84
  full_knots <- c(internal_knots[1] - step * 3:1, #3 exterior knots before
                  internal_knots, #internal knots
                  internal_knots[length(internal_knots)] + step * 1:3)
                  #3 exterior knots after
  
  #construct matrix X-tilde using splineDesign with 84 knots as defined 
  #before starting from min(t)-30 until max(t) with order of the spline = 4
  Xtilde <- splineDesign(knots=full_knots, 
                         x=(min(t)-30):max(t), 
                         ord=4, # cubic splines
                         outer.ok= TRUE)
  
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

mats <- make_matrices(t) #construct matrices for the observed days
X <- mats$X ## X
S <- mats$S ## S (second-difference penalty)
Xtilde <- mats$Xtilde ## X tilde

##============= (2) Function Penalised NLL and its Objective Func. =============

# In this part, we define objective function to be optimised. For this case, 
# it's already given that f(t) has possion GLM distribution

# First, we know that (log likehood funct of poisson) formula for poisson 
#        l = ∑[yi.log(μi)−μi−log(yi!)]
#        where  μ = Xβ
#        Also, we drop log(yi!) as it does not depends on parameter
#   
# NLL = - ∑(Xβ−y⋅log(Xβ))[Negative loglikelihood function], where μ = X β and 
#       β=exp(γ). The exponential reparameterisation ensures β>0.
#
# Then, the penalty is applied to β (control smoothness of the function), 
#       P = 0.5λβ'Sβ

# Penalised NLL = NLL + P

# Objective Function
pen_nll <- function(gamma, X, y, S, lambda, weight=1) {
  # Function to compute penalised negative log likehood
  # Input/Argument : (1) Gamma : K-vector of spline coefficients
  #                  (2) y : the deaths on day of the year ti
  #                  (3) S : penalised matrix
  #                  (4) lambda : smoothing parameter
  #                  (5) weight : weight for bootstraping sampled data
  #                      adjusted for more general use in this case
  
  # Output/Return  : single numeric value of pen_nll
  beta <- exp(gamma) # β = exp(γ)
  mu <- as.vector(X %*% beta) # μ = Xβ
  #ll formulae dropping the constant of lgamma(y+1)
  ll <- sum(weight*((y * log(mu)) - mu))
  # likelihood funct of possion dist
  penalty <- 0.5 * lambda * t(beta) %*% (S %*% beta) # penalty
  return(-ll + penalty)
}##pen_nll

# Define Gradient vector of Objective Function/ its derivative vector 
pen_grad <- function(gamma, X, y, S, lambda, weight=1) {
  # Function to compute derivative vector
  # Input/Argument : (1) Gamma : K-vector of spline coefficients
  #                  (2) y : the deaths on day of the year ti
  #                  (3) S : penalised matrix
  #                  (4) lambda : smoothing parameter
  #                  (5) weight : weight for bootstraping sampled data
  #                      adjusted for more general use in this case
  
  beta <- exp(gamma) # β = exp(γ), ensures positivity
  mu <- X %*% beta # μ = Xβ
  F <- t(X) %*% (weight*(y / mu - 1)) # ∂li/∂γj; or F =  diag(yi/µi−1).X.diag(β) 
  grad_ll <- -F * beta # total derivative for all data
  grad_pen <- lambda * (beta * (S %*% beta)) # total gradient for penalty 
  return(as.vector(grad_ll + grad_pen))
}

# Checking the derivative (From SP Notes pg 74)
K <- 80
fd <- gamma0 <- rep(0, K) #start from all zeros
lambda <- 5e-5 
pen_nll0 <- pen_nll(gamma0, X, y, S, lambda) 
eps <- 1e-7
for (i in 1:length(gamma0)) {
  gamma1 <- gamma0; gamma1[i] <- eps
  pen_nll1 <- pen_nll(gamma1,X,y,S,lambda)
  fd[i] <- (pen_nll1 - pen_nll0)/eps
}
range(fd - pen_grad(gamma0, X, y, S, lambda)) ## apx zero, correct gradient

##================ (3) Finding the sane starting values for gamma ===================
gamma0 <- rep(0, K)       # initial values for gamma
lambda <- 5e-5            # fixed smoothing parameter for sanity check

#minimise the penalised negative log-likelihood using BFGS.
#the analytic gradient (pen_grad) is supplied to improve optimisation stability
fit <- optim(
  par     = gamma0,                 # initial gamma
  fn      = pen_nll,                # objective function
  gr      = pen_grad,               # analytic gradient
  method  = "BFGS",                 # efficient for smooth problems
  X=X, y=y, S=S, lambda=lambda
)

beta_hat <- exp(fit$par) #estimated spline coefficients - beta hat
mu_hat <- as.vector(X %*% beta_hat) #fitted deaths - mu hat

# Construct data frames for plotting:
#   deaths_df : observed deaths and fitted deaths
#   infect_df : estimated infection curve f̂(t) = Xtilde*beta hat
deaths_df <- data.frame(day = t, deaths = y, fitted = mu_hat)
infect_df <- data.frame(day = (min(t)-30):max(t),
                        f_hat = as.vector(Xtilde %*% beta_hat))

# Plot the observed deaths, fitted deaths, and estimated infection curve f(t)
windows()
ggplot() +
  geom_point(data = deaths_df,
             aes(x = day, y = deaths, 
                 color = "Observed Deaths"), size=1.5) +
  geom_line(data = deaths_df,
            aes(x = day, y = fitted, 
                color = "Fitted Deaths"), size=1) +
  geom_line(data = infect_df,
            aes(x = day, y = f_hat, 
                color = "New Infection f(t)"), size = 1) +
  labs(
    x = "Day of year / Julian Day Observed",
    y = "Daily Deaths (Counts)",
    title = expression(paste("Daily death and estimated number of infection,",
                             lambda, "=5e-5")),
    color = "Legend"
  ) +
  scale_color_manual(
    values = c(
      "Observed Deaths" = "black",
      "Fitted Deaths" = "red",
      "New Infection f(t)" = "blue"
    )) +
  theme_bw()
gamma2 <- fit$par
##================ (4) Fit the model using BFGS optimization ===================
lambdas <- exp(seq(-13, -7, length=50))
BIC_vals <- numeric(length(lambdas))

for (i in seq_along(lambdas)) {
  fit <- optim(par=gamma2, fn=pen_nll, gr=pen_grad, method="BFGS",
               X=X, y=y, S=S, lambda=lambdas[i])
  beta_hat <- exp(fit$par)
  mu_hat <- as.vector(X %*% beta_hat)
  W <- diag(as.vector(y / mu_hat^2))
  H0 <- t(X) %*% W %*% X
  H_lambda <- H0 + lambdas[i] * S
  EDF <- sum(diag(solve(H_lambda, H0)))
  ll <- sum(y * log(mu_hat) - mu_hat)
  BIC_vals[i] <- -2 * ll + log(length(y)) * EDF
}

# Plot BIC against lambda
min_BIC_index = which.min(BIC_vals)
best_lambda <- lambdas[min_BIC_index] ## Optimum Lambda
plot(seq_along(lambdas), BIC_vals, type = "o",
     xlab = expression(lambda), ylab = "BIC",
     main = "BIC vs. Lambda")
points(min_BIC_index, BIC_vals[min_BIC_index],
       col = "red", pch = 19, cex = 1.5)

# Parameter (µ) when Lambda optimum
fit <- optim(par=gamma2, fn=pen_nll, gr=pen_grad, method="BFGS",
             X=X, y=y, S=S, lambda=lambdas[min_BIC_index])
mu_hat <- fit$par


##=============== (5) Non Parametric Bootstrap Uncertainty =====================
# Initialize number of replicate sample
n_bootstrap <- 200 

# Initialize matrix to store bootstrap replicates
mat_boots <- matrix(NA, nrow=nrow(Xtilde), ncol=n_bootstrap)

for (b in 1:n_bootstrap) {
  wb <- tabulate(sample(n, replace=TRUE), n)
  fit_b <- optim(gamma2, 
                 pen_nll, gr=pen_grad, 
                 method="BFGS",
                 y=y, X=X, S=S, 
                 lambda=best_lambda, 
                 weight=wb)
  beta_b <- exp(fit_b$par)
  mat_boots[,b] <- Xtilde %*% beta_b ## estimate number of new infection
}

##=====================(6) Final Plot===========================================

# Estimate the daily new infection rate f(t) with its 95% confidence limits
infect <- data.frame(time_ft=(min(t)-30):max(t), 
                    mean_ft=rowMeans(mat_boots), #estimated mean of f(t)
                    sd_ft=apply(mat_boots,1,sd), #estimated sd of f(t)
                    
                    # apply is used for calculate quantile per row/observed days
                    lb_ft = apply(mat_boots, 1, quantile, probs=0.025),
                    ub_ft = apply(mat_boots, 1, quantile, probs=0.975)
                    )

# Estimate death and observed death
beta_hat <- exp(gamma2) #
mu_hat <- as.vector(X %*% beta_hat) # mean deaths per day (from based fit model)
deaths <- data.frame(day = t, deaths = y, model_fit = mu_hat)

# Final Plot (Bootstrap)
ggplot() +
  ## Add 95% CI to plot
  geom_ribbon(
    data = infect,
    aes(x = time_ft, ymin = lb_ft, ymax = ub_ft, fill = "CI")
  ) +
  ## Add point represents death recorded from nhs data
  geom_point(
    data = deaths,
    aes(x = day, y = deaths, color = "Observed Deaths"),
    size = 1.5
  ) +
  ## Add line represents estimated death
  geom_line(
    data = deaths,
    aes(x = day, y = mu_hat, color = "Fitted Deaths"),
    size = 1
  ) +
  ## Add line represents estimated f(t) based on 200 replicates data
  geom_line(
    data = infect,
    aes(x = time_ft, y = mean_ft, color = "New Infection f(t)"),
    size = 1
  ) +
  ## Add label
  labs(
    x = "Day of year / Julian Day Observed",
    y = "Daily Deaths (Counts)",
    title = "The daily death and estimated number of daily infection",
    color = "Legend",
    fill = "Legend"
  ) +
  ## Add legend for line and points
  scale_color_manual(
    values = c(
      "Observed Deaths" = "black",
      "Fitted Deaths" = "red",
      "New Infection f(t)" = "blue"
    ),
    guide = "legend"
  ) +
  ## Add legend for CI
  scale_fill_manual(
    values = c("CI" = rgb(0.53, 0.81, 0.98, 0.25)),
    guide = "legend"
  ) +
  guides(
    color = guide_legend(
      title = "Legend",
      order = 1,
      override.aes = list(fill = NA)
    ),
    fill = guide_legend(
      title = NULL, ## no titles for CI's legend
      order = 1,
      override.aes = list(
        linetype = 0,
        shape = 22,
        size = 5,
        color = NA,
        fill = rgb(0.53, 0.81, 0.98, 0.25)
      )
    )
  ) +
  theme_bw()

end <- Sys.time()
end-start