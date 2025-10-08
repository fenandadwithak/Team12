# 1) Create n people assigned to corresponding household where maximum size is 5
n <- 1000 #population size
hmax <- 5 #maximum household size
#creating n vector of which the number refers to household
#and the occurences refers to size of corresponding household
prob <- rep (1/hmax, hmax) #probability of discrete uniform for each size
set.seed(3)
h <- sample(rep(1:n, times=sample(1:hmax, n, replace=TRUE, prob = prob))[1:n])

# 2) Contact network
#nc = average number of contacts per person
#beta = sociability parameter
get.net <- function (beta, nc=15, h) {
  n <- length(beta)
  beta_bar <- mean(beta)
  links <- vector("list", n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (h[i] != h[j]) { # exclude household members
        p_link <- nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
        if (runif(1) < p_link) {
          links[[i]] <- c(links[[i]], j)
          links[[j]] <- c(links[[j]], i)
        }
      }
    }
  }
  return(links)
}

# 3) nseir function
nseir<-function(beta, h, alink, alpha=c(.1,.01,.01),
                delta=.2, gamma=.4, nc=15, nt=100, pinf=.005){
  n <- length(beta)
  state <- rep("S", n)
  state[sample(1:n, pinf * n)] <- "I"
  
  res <- matrix(0, nrow = nt, ncol = 4)
  colnames(res) <- c("S", "E", "I", "R")
  
  beta_bar <- mean(beta)
  
  for (t in 1:nt) {
    newE <- rep(FALSE, n)
    newI <- (state == "E") & (runif(n) < gamma)
    newR <- (state == "I") & (runif(n) < delta)
    
    for (i in which(state == "I")) {
      # household infections
      hh_members <- which(h == h[i] & state == "S")
      newE[hh_members]<-newE[hh_members]|(runif(length(hh_members))< alpha[1])
      
      # network infections
      contacts <- alink[[i]]
      contacts <- contacts[state[contacts] == "S"]
      newE[contacts] <- newE[contacts] | (runif(length(contacts)) < alpha[2])
      
      # random mixing
      targets <- which(state == "S")
      p_random <- alpha[3]*nc*beta[i]*beta[targets]/(beta_bar^2 * (n - 1))
      newE[targets] <- newE[targets] | (runif(length(targets)) < p_random)
    }
    
    state[newE & state == "S"] <- "E"
    state[newI] <- "I"
    state[newR] <- "R"
    
    res[t, ] <- table(factor(state, levels = c("S", "E", "I", "R")))
  }
  
  return(data.frame(t = 1:nt, res))
}

# 4) Function Plot
plot.seir <- function(res, title = "SEIR Dynamics") {
  plot(res$t, res$S, type = "l", ylim = c(0, max(res$S)),
       xlab = "Days", ylab = "Number of individuals",
       col = "blue", lwd = 2, main = title)
  lines(res$t, res$E, col = "orange", lwd = 2)
  lines(res$t, res$I, col = "red", lwd = 2)
  lines(res$t, res$R, col = "darkgreen", lwd = 2)
  legend("right", legend = c("S", "E", "I", "R"),
         col = c("blue", "orange", "red", "darkgreen"), lwd = 2)
}

# 5) Compare 
set.seed(123)

beta <- runif(n)
alink <- get.net(beta, nc = 15, h = h)

# Scenario 1: Full model
res1 <- nseir(beta, h, alink)

# Scenario 2: Random mixing only
res2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))

# Scenario 3: Constant beta
res3 <- nseir(rep(mean(beta), n), h, alink)

# Scenario 4: Constant beta + random mixing
res4 <- nseir(rep(mean(beta), n), h, alink, alpha = c(0, 0, 0.04))

# Plot all scenarios side by side
par(mfrow = c(2, 2))
plot.seir(res1, "Full Model")
plot.seir(res2, "Random Mixing Only")
plot.seir(res3, "Constant Beta")
plot.seir(res4, "Constant Beta + Random Mixing")