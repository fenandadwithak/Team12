## ===============================================================
##        Practical 2 — Social Structure in SEIR Models
## ===============================================================
## What you get:
## (1) make_households()  : build household vector h
## (2) get_net()          : build regular-contact network alink
## (3) nseir()            : SEIR simulator with household/network/random mixing
## (4) plot_seir()        : tidy plotting for S/E/I/R counts
## (5) compare_four()     : runs the four required scenarios + 2x2 plot

## ===============================================================

## -------------------------
## (1) Households generator
## -------------------------
make_households <- function(n, hmax = 5, seed = NULL) {
  stopifnot(n >= 1, hmax >= 1)
  if (!is.null(seed)) set.seed(seed)
  
  ## Build random household sizes ~ Uniform{1,...,hmax} until we reach n
  sizes <- integer(0)
  while (sum(sizes) < n) sizes <- c(sizes, sample(1:hmax, 1))
  if (sum(sizes) > n) sizes[length(sizes)] <- sizes[length(sizes)] - (sum(sizes) - n)
  
  ## Assign people to households (household IDs 1..H)
  H   <- length(sizes)
  ids <- rep(seq_len(H), times = sizes)
  h   <- sample(ids, length(ids), replace = FALSE)  # shuffle people across IDs
  h
}

## ----------------------------------
## (2) Regular-contact network (alink)
## ----------------------------------
## Returns a list alink where alink[[i]] is the vector of i's regular (non-household) contacts.
## Two methods:
##   - method="exact": create an undirected edge (i,j) with probability
##         p_ij = min(1, nc * beta_i * beta_j / (mean(beta)^2 * (n-1)))
##     for non-household pairs only. This matches the brief exactly but is O(n^2).
##   - method="approx" (default): for each i, sample ~nc partners from non-household
##     candidates with probability ∝ beta_j, adding symmetric edges. This is fast (≈O(n*nc))
##     and yields the intended expected degrees; good for n up to ~10k.
get_net <- function(beta, h, nc = 15, method = c("approx", "exact"), seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  n <- length(beta); stopifnot(length(h) == n, nc >= 0)
  alink <- vector("list", n)
  beta_bar <- mean(beta)
  
  if (method == "approx") {
    for (i in seq_len(n)) {
      cand <- which(h != h[i])
      cand <- cand[cand != i]
      if (!length(cand)) next
      target <- min(nc, length(cand))
      if (target == 0) next
      pick <- sample(cand, size = target, replace = FALSE, prob = beta[cand])
      for (j in pick) {
        if (!(j %in% alink[[i]])) alink[[i]] <- c(alink[[i]], j)
        if (!(i %in% alink[[j]])) alink[[j]] <- c(alink[[j]], i)
      }
    }
  } else {  # exact Chung–Lu style probability, O(n^2)
    c0 <- nc / (beta_bar^2 * (n - 1))
    for (i in seq_len(n - 1)) {
      cand <- (i + 1):n
      cand <- cand[h[cand] != h[i]]             # exclude same-household
      if (!length(cand)) next
      pij <- pmin(1, c0 * beta[i] * beta[cand]) # edge probs
      hits <- runif(length(cand)) < pij
      js <- cand[hits]
      if (!length(js)) next
      alink[[i]] <- c(alink[[i]], js)
      for (j in js) alink[[j]] <- c(alink[[j]], i)
    }
  }
  alink
}

## -------------------------------
## (3) SEIR with social structure
## -------------------------------
## states: 1=S, 2=E, 3=I, 4=R
## alpha = c(alpha_h, alpha_c, alpha_r)  (household, network, random mixing)
## delta: P(I -> R) per day
## gamma: P(E -> I) per day
## pinf : initial infected proportion (start directly in I)
nseir <- function(beta, h, alink,
                  alpha = c(0.1, 0.01, 0.01),
                  delta = 0.2, gamma = 0.4,
                  nc = 15, nt = 100, pinf = 0.005,
                  seed = NULL,
                  exact_random_mix = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(beta)
  stopifnot(length(h) == n, length(alink) == n, length(alpha) == 3)
  
  a_h <- alpha[1]; a_c <- alpha[2]; a_r <- alpha[3]
  beta_bar <- mean(beta)
  
  ## Initial states
  S <- rep(1L, n)
  I0 <- sample.int(n, size = max(1L, floor(n * pinf)))
  S[I0] <- 3L  # Infectious at t=1
  
  out <- matrix(NA_integer_, nrow = nt, ncol = 4)
  colnames(out) <- c("S","E","I","R")
  
  for (t in seq_len(nt)) {
    isS <- (S == 1L); isE <- (S == 2L); isI <- (S == 3L)
    
    idxS <- which(isS)
    p_hh  <- rep(0, length(idxS))
    p_net <- rep(0, length(idxS))
    p_mix <- rep(0, length(idxS))
    
    ## (A) Household infections: 1 - (1 - a_h)^(# infectious in household)
    if (a_h > 0) {
      I_counts <- tapply(isI, h, sum)                   # infectious per household
      hh_counts <- I_counts[as.character(h[idxS])]
      hh_counts[is.na(hh_counts)] <- 0
      p_hh <- 1 - (1 - a_h)^(hh_counts)
    }
    
    ## (B) Regular-network infections: 1 - (1 - a_c)^(# infectious neighbors)
    if (a_c > 0) {
      isI_num <- as.integer(isI)
      kinf <- vapply(alink[idxS], function(v) if (length(v)) sum(isI_num[v]) else 0L, integer(1L))
      p_net <- 1 - (1 - a_c)^(kinf)
    }
    
    ## (C) Random mixing:
    ## Exact definition: per pair (i in I, j in S), P(i infects j) =
    ##   a_r * nc * beta_i * beta_j / (beta_bar^2 * (n-1)).
    ## Combine across all infectives i to get j's overall probability.
    if (a_r > 0) {
      c0 <- a_r * nc / (beta_bar^2 * (n - 1))
      if (!exact_random_mix) {
        ## Fast hazard approximation: p = 1 - exp(- c0 * beta_j * sum(beta_i over I))
        sum_beta_I <- sum(beta[isI])
        lam <- c0 * beta[idxS] * sum_beta_I
        p_mix <- 1 - exp(-lam)
      } else {
        ## Exact product over infectious set: p = 1 - Π_i (1 - c0 * beta_j * beta_i)
        bi <- beta[isI]
        for (k in seq_along(idxS)) {
          bj <- beta[idxS[k]]
          p_mix[k] <- 1 - exp(sum(log1p(-c0 * bj * bi)))
        }
      }
    }
    
    ## Combine independent channels
    p_inf <- 1 - (1 - p_hh) * (1 - p_net) * (1 - p_mix)
    p_inf[p_inf < 0] <- 0; p_inf[p_inf > 1] <- 1
    
    ## Transitions (synchronous updates)
    newE_idx <- idxS[ runif(length(idxS)) < p_inf ]
    idxE <- which(isE)
    newI_idx <- idxE[ runif(length(idxE)) < gamma ]
    idxI <- which(isI)
    newR_idx <- idxI[ runif(length(idxI)) < delta ]
    
    if (length(newE_idx)) S[newE_idx] <- 2L
    if (length(newI_idx)) S[newI_idx] <- 3L
    if (length(newR_idx)) S[newR_idx] <- 4L
    
    out[t,] <- c(sum(S==1L), sum(S==2L), sum(S==3L), sum(S==4L))
  }
  
  data.frame(day = 1:nt, out, row.names = NULL)
}

## ----------------
## (4) Nice plotting
## ----------------
plot_seir <- function(res, main = "", legend_pos = "right") {
  matplot(res$day, res[,c("S","E","I","R")],
          type = "l", lty = 1, lwd = 2, xlab = "Day", ylab = "Population",
          main = main, col = 1:4)
  legend(legend_pos, bty = "n", legend = c("S","E","I","R"),
         lty = 1, lwd = 2, col = 1:4)
}

## -------------------------------
## (5) The four required scenarios
## -------------------------------
## Scenarios:
##   1) Full model with default parameters.
##   2) Remove household + regular network; keep equal average contacts via random mixing:
##        alpha_h = alpha_c = 0, alpha_r = 0.04
##   3) Full model but set beta to its mean (constant sociability).
##   4) Combine (2) and (3): constant beta + random mixing only.
compare_four <- function(n = 3000, hmax = 5, nc = 15, nt = 120,
                         alpha = c(0.1, 0.01, 0.01),
                         delta = 0.2, gamma = 0.4, pinf = 0.005,
                         seed = 1,
                         net_method = "approx",
                         mix_exact = FALSE) {
  set.seed(seed)
  beta <- runif(n, min = 0, max = 1)           # U(0,1) as specified
  h    <- make_households(n, hmax)             # households
  net  <- get_net(beta, h, nc, method = net_method)  # regular contacts
  
  ## 1) Full model
  res1 <- nseir(beta, h, net, alpha=alpha, delta=delta, gamma=gamma,
                nc=nc, nt=nt, pinf=pinf, exact_random_mix = mix_exact)
  
  ## 2) Random mixing only (keep average contacts via a_r = 0.04)
  alpha2 <- c(0, 0, 0.04)
  res2 <- nseir(beta, h, vector("list", n), alpha=alpha2, delta=delta, gamma=gamma,
                nc=nc, nt=nt, pinf=pinf, exact_random_mix = mix_exact)
  
  ## 3) Full model, constant beta
  beta_bar <- mean(beta)
  beta3 <- rep(beta_bar, n)
  res3 <- nseir(beta3, h, net, alpha=alpha, delta=delta, gamma=gamma,
                nc=nc, nt=nt, pinf=pinf, exact_random_mix = mix_exact)
  
  ## 4) Constant beta + random mixing only
  res4 <- nseir(beta3, h, vector("list", n), alpha=alpha2, delta=delta, gamma=gamma,
                nc=nc, nt=nt, pinf=pinf, exact_random_mix = mix_exact)
  
  ## 2x2 plots
  op <- par(mfrow=c(2,2), mar=c(4,4,2,1))
  on.exit(par(op))
  plot_seir(res1, main = "Full model")
  plot_seir(res2, main = "Random mixing only (a_r = 0.04)")
  plot_seir(res3, main = "Full model, constant beta")
  plot_seir(res4, main = "Random mixing + constant beta")
  
  invisible(list(full=res1, rand_only=res2, full_constbeta=res3, rand_constbeta=res4))
}

## --------------------------
## (Optional) quick smoke test
## --------------------------
# source("practical2_seir.R")
# out <- compare_four(n=2000, hmax=5, nc=15, nt=100, seed=42,
#                     net_method="approx", mix_exact=FALSE)
# str(out$full)

## ----------
## Team note:
## ----------
## - Aseel Alghamdi(S2901228): households + network +nseir core
## - Fenanda Dwitha Kurniasari(S2744048):  vectorization checks  +plotting + scenarios + report
## - Nurmawiya(S2822251): plotting + scenarios + report
## ===============================================================