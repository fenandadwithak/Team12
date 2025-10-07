# 1) Create household assignment
n <- 1000 #population size
hmax <- 5 #maximum household size
#number of households with size uniformly distributed of max = 5
hh_n <- ceiling(n/((hmax+1)/2))
#creating n vector of which the number refers to household
#and the occurence of the number refers to size of corresponding household
h <- rep(1:hh_n, times = sample(1:hmax, hh_n, replace = TRUE))[1:n]

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

#3 nseir function
nseir<-function(beta,h,alink,alpha=c(.1,.01,.01),
                delta=.2,gamma=.4,nc=15,nt=100,pinf=.005){
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