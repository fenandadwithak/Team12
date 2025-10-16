# PROJECT 2 - EXTENDED STATISTICAL PROGRAMMING =================================
# Group 12 
# Aseel Alghamdi : S2901228
# Fenanda Dwitha Kurniasari : S2744048
# Nurmawiya : S2822251

# Aseel : nseir function, cross-check code, cross-check the entire 
#         code and revise code & comments
# Fenanda : plot function, model comparison, cross-check the entire code
#           and revise code & comments
# Nurmawiya : n people and household, contact network, cross-check the entire
#             code and revise code & comments

# 1) Create n people assigned to corresponding household where maximum size is 5
n <- 10000 #population size
hmax <- 5 #maximum household size
#creating n vector of which the number refers to household
#-and the occurrences refers to size of corresponding household
set.seed(1)
sizes <- sample(1:hmax, n, replace = TRUE)
h <- sample(rep(1:length(sizes), times = sizes)[1:n]) #shuffle people accross hh

# 2) Contact network
get.net <- function (beta, nc=15, h) {
  # Function to create links (daily contacts) for each person
  
  # nc: average number of contacts per person, set to 15
  # beta: sociability parameter, higher means more likely to form links
  # OUTPUT : list of connected networks from each person (output.net)
  
  n <- length(beta) ## calculate total obs of beta, equal to total population
  beta_bar <- mean(beta) ## calculate mean of sociability parameter
  
  # Allocate possible connection formed in each individual
  
  output.net <- vector("list", n) ## list to store the result
  output.net0 <- integer(n) ## vector to store number of contact made per person
  
  # Make list of n elements which in each element initially has nc subelement
  # Later, we would allocate (drop/add sub-element and store possible contact
  # Initialise contact list for each person with nc=15 placeholders
  for (k in 1:n) output.net[[k]] <- integer(nc) ## set initial value output.net
  
  # Loop over each individual
  for (i in 1:(n-1)) {
    j <- (i+1):n # generating candidate contacts
    flag <- h[i] != h[j]  # flag to filter same household 
    j <- j[flag] # filter candidates contacts from different household
    
    if (length(j) == 0) next # skip iteration when there are no candidate
    
    # Probability of  a link between persons i and j 
    p_link <- nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
    
    # Select the candidates contacts based on p-link
    linked <- j[runif(length(j)) < p_link] #vector possible candidate contacts
    
    # If there is chosen candidate(s), store it into output.net and output.net0
    if (length(linked) > 0) {
      # Position to fill
      pos = (output.net0[i]+1):(output.net0[i]+length(linked))
      
      # Assign link[[i]]to candidates from "linked" at given position
      output.net[[i]][pos] <- linked
      
      # Assign number of chosen candidates at i index in output.net0
      output.net0[i] <- output.net0[i] + length(linked)
      
      # Vice versa, assign "i" to their linked candidates
      for (j1 in linked) {
        # Assign i to output.net's of its chosen candidates
        output.net[[j1]][output.net0[j1]+1] <- i 
        # Assign the number of candidates (add by 1) to its chosen candidates
        output.net0[j1] <- output.net0[j1] + 1 
      }
    }
    ## to the next value of i
  }
  
  # Drop sub-element/rejected candidates and replace person with no contact 
  for (k in 1:n) {
    if (output.net0[k] > 0) { 
      # Filter assigned links when chosen number of chosen candidates < nc
      output.net[[k]] <- output.net[[k]][1:output.net0[k]]
    } else {
      # if there is no chosen candidate, replace sub-element into NULL
      output.net[[k]] <- integer(0)
    }
  }
  return(output.net)
} ## get.net

# 3) nseir function
nseir <- function(beta, h, alink,
                  alpha = c(0.1, 0.01, 0.01),
                  delta = 0.2, gamma = 0.4,
                  nc = 15, nt = 100, pinf = 0.005,
                  seed = NULL,
                  exact_random_mix = FALSE) {
  
  # SEIR stochastic simulation model
  # beta: sociability parameter of each person
  # h: household list
  # alink: list of regular contacts of each person returned by get.net
  # alpha[1]: daily prob i infecting j if in same household (a_h)
  # alpha[2]: daily prob i infecting j if in regular contact (a_c)
  # alpha[3]: daily prob i infecting j through random mixing (a_r)
  # delta: daily prob I -> R
  # gamma: daily prob E -> I
  # nc: average number of daily contacts for each person
  # nt: number of days to simulate
  # pinf: proportion of the initial population to randomly start in the I state
  
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
} ##nseir

# 4) Function Plot
dyn.plot <- function(seir, title = "SEIR Dynamics") {
  # Function to create dynamics plot of the simulated population in each states
  
  # x-axis : days to simulate (from 1 to 100)
  # y-axis : number of individuals per states
  #          explanation of y-axis
  #          ---------------------|
  #          at the initial day, number of susceptible individual is the largest
  #          total individuals in all states S+E+I+R = 1000
  #          these numbers would dynamically change day by day in each states
  #          depend on the defined models
  
  # input : nseir model from function nseir 
  # output : Dynamic Plot days vs number of individuals per states
  
  # Plot the states in S(susceptible)
  # Create the base plot
  plot(seir$day, seir$S, type = "p", ##plot day vs number of susceptible 
       ylim = c(0, max(seir$S)), xlim = c(0, max(seir$day)), ## set maximum range
       xlab = "Days", ylab = "Number of individuals", ##labelling the plot
       col = "blue", pch = 1, main = title) ##color and type of points in plot
  
  # Add points for each states
  points(seir$day, seir$E, col = "orange", pch = 1) ##add states 
  points(seir$day, seir$I, col = "red", pch = 1)
  points(seir$day, seir$R, col = "darkgreen", pch = 1)
  
  # Add legend
  legend(
    "right",        
    legend = c("Susceptible", "Exposed", "Infectious", "Recovered"),
    col = c("blue", "orange", "red", "darkgreen"),
    pch = 1, bty = "n", ##type of points and no box legend 
    cex = 0.8,       ## smaller text
    pt.cex = 0.8,    ## smaller symbols
    y.intersp = 0.6, ## spacing horizontal
    x.intersp = 0.4  ## spacing between symbol and text 
  )
} ## dyn.plot


# 5) Comparing Model
# Setting up beta vector and alink pairs
# beta is vector of n values (n = 1000) ~ U(0,1)
# alink contains list of network made per each person 
set.seed(345) ## to ensure we get same vector values of beta for all model
beta <- runif(n); alink <- get.net(beta, nc = 15, h = h)

# Scenario 1: Full model with default parameters
seir1 <- nseir(beta, h, alink)

# Scenario 2: Random mixing - Setting alpha[1] = alpha[2] = 0 and 
#             alpha[3] = 0.04
seir2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))

# Scenario 3: Full model with constant beta (beta vector = average beta)
seir3 <- nseir(rep(mean(beta), n), h, alink)

# Scenario 4: Constant beta + random mixing
seir4 <- nseir(rep(mean(beta), n), h, alink, alpha = c(0, 0, 0.04))

# Plot all scenarios 
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)) 
dyn.plot(seir1, "Full Model Default Parameter")
dyn.plot(seir2, "Random Mixing Only")
dyn.plot(seir3, "Constant Beta")
dyn.plot(seir4, "Constant Beta and Random Mixing")