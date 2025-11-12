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
## What you get:
## (1) 
## (2) 
## (3) 
## (4) 
## (5)
## (6)


##===================== Data Preparation & Load library ========================
library(splines) #load library splines
library(ggplot2)

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
  #                      larger K gives a more flexible/wigglier basis
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
  
  # Define Knots
  internal_knots <- seq(min(t)-30, max(t), length.out = K-2)
  range_knot = internal_knots[3]-internal_knots[2]
  before_knots = internal_knots[1] - 3:1*range_knot
  after_knots = internal_knots[length(internal_knots)] + 1:3*range_knot
  full_knots = c(before_knots, internal_knots, after_knots)
  
  Xtilde <- splineDesign(knots=full_knots, 
                         x=(min(t)-30):max(t), 
                         ord=4, # cubic splines
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

# In this part, we define objective function to be optimised. For this case, 
# it's already given that f(t) has possion GLM distribution

# First, we know that (log likehood funct of poisson) formula for poisson 
#        l = ∑[yi.log(μi)−μi−log(yi!)]
#        where  μ = Xβ
#        Also, we drop log(yi!) as it does not depends on parameter
#   
# NLL = ∑(Xβ−y⋅log(Xβ))[Negative loglikelihood function]
#       where β=exp(γ) (imposed for positivity)
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
  
  
  beta <- exp(gamma) # β = exp(γ)
  mu <- X %*% beta # μ = Xβ
  F <- t(X) %*% (weight*(y / mu - 1)) # ∂li/∂γj; or F =  diag(yi/µi−1).X.diag(β) 
  grad_ll <- -F * beta # total derivative for all data
  grad_pen <- lambda * (beta * (S %*% beta)) # total gradient for penalty 
  return(as.vector(grad_ll + grad_pen))
}

# Checking the derivative (From SP Notes pg 74)
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
  X=X, y=y, S=S, lambda=lambda
)

beta_hat <- exp(fit$par) # estimated spline coefficients (K-vector)
mu_hat <- as.vector(X %*% beta_hat)

deaths_df <- data.frame(day = t, deaths = y, fitted = mu_hat)
infect_df <- data.frame(day = (min(t)-30):max(t),
                        f_hat = as.vector(Xtilde %*% beta_hat))

# Combined Plot
windows()
ggplot() +
  geom_point(data = deaths_df,
             aes(x = day, y = deaths, color = "Observed Deaths"), size=1.5) +
  geom_line(data = deaths_df,
            aes(x = day, y = fitted, color = "Fitted Deaths"), size=1) +
  geom_line(data = infect_df,
            aes(x = day, y = f_hat, color = "New Infection f(t)"), size = 1) +
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
  theme_bw() +
  theme(
    axis.title.y.left = element_text(color="black"),
    axis.title.y.right = element_text(color="blue"),
    legend.position = c(0.85, 0.75)
  ) +
  scale_y_continuous(sec.axis = sec_axis(~., name = expression(hat(f)(t))))


##================ (4) Fit the model using BFGS optimization ===================
gamma2 <- fit$par
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

#plot BIC against lambda
best_lambda <- lambdas[which.min(BIC_vals)]
plot(seq_along(lambdas), BIC_vals)

cat("Best lambda:", best_lambda, "\n")

##=============== (5) Non Parametric Bootstrap Uncertainty =====================
# Initialize number of replicate sample
n_bootstrap <- 200 

# Initialize matrix to store bootstrap replicates
mat_boots <- matrix(NA, nrow=nrow(mats$Xtilde), ncol=n_bootstrap)

for (b in 1:n_bootstrap) {
  wb <- tabulate(sample(n, replace=TRUE), n)
  fit_b <- optim(gamma2, 
                 pen_nll, gr=pen_grad, 
                 method="BFGS",
                 y=y, X=X, S=S, 
                 lambda=best_lambda, 
                 weight=wb)
  beta_b <- exp(fit_b$par)
  mat_boots[,b] <- mats$Xtilde %*% beta_b ## estimate number of new infection
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
  theme_bw() +
  # Adjust position
  theme(
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.box.background = element_rect(color = "gray70", linetype = "dotted"),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.4, "lines"),
    legend.key.width  = unit(1, "lines"),
    legend.spacing.y  = unit(0.05, "lines"),
    legend.margin = margin(2, 2, 2, 2)
  )

end <- Sys.time()
end-start