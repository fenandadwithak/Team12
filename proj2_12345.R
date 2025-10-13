# PROJECT 2 - EXTENDED STATISTICAL PROGRAMMING =================================

# Group 12 
# Aseel Alghamdi : S2901228
# Fenanda Dwitha Kurniasari : S2744048
# Nurmawiya : S2822251

# Aseel : create n people and its household, contact network, cross-check 
#         the entire code & revise comments
# Fenanda : create plot function, compare the model, cross-check the entire code
#           and revise comments
# Nurmawiya : create nseir function, cross-check code, cross-check the entire 
#             code and revise comments


# SEIR (Suspectible, Exposed, Infections, Recovered) Model Simulation
# 1) Create n people assigned to corresponding household where maximum size is 5
  n <- 1000 #population size
  hmax <- 5 #maximum household size
  #creating n vector of which the number refers to household
  #-and the occurrences refers to size of corresponding household
  set.seed(3)
  h <- rep(1:n, times = sample(1:hmax, n, replace=TRUE))[1:n]

# 2) Contact network
get.net <- function (beta, nc=15, h) {
  # Function to create links (daily contacts) for each person
  # nc: average number of contacts per person
  # beta: sociability parameter, higher means more likely to form links
  n <- length(beta)
  beta_bar <- mean(beta) #mean sociability parameter
  links <- vector("list", n) #set null vector for links of each person
  
  for (i in 1:(n-1)) { #loop over person i from 1 to n-1
    for (j in (i+1):n) { #loop over person j from i+1 to n to avoid duplicates
      if (h[i] != h[j]) { #only connect if they are not in the same household
        #prob that i and j have contact (daily link)
        p_link <- nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
        
        #If daily link of i and j person > random value U(0,1),
        #then create a link between i and j
        if (runif(1) < p_link) {
          links[[i]] <- c(links[[i]], j) #Store j into list who connects with i
          links[[j]] <- c(links[[j]], i) #Store i in j's list of contacts 
        }
      }
    }
  }
  return(links) #return the list of daily contacts for each person
} ## get.net

# 3) nseir function
nseir <- function(beta, h, alink, alpha=c(.1, .01, .01),
                  delta=.2, gamma=.4, nc=15, nt=100, pinf=.005){
  # SEIR stochastic simulation model
  # beta: sociability parameter of each person
  # h: household list
  # alink: list of regular contacts of each person returned by get.net
  # alpha[1]: daily prob i infecting j if in same household (alpha h)
  # alpha[2]: daily prob i infecting j if in regular contact (alpha c)
  # alpha[3]: daily prob i infecting j through random mixing (alpha r)
  # delta: daily prob I -> R
  # gamma: daily prob E -> I
  # nc: average number of daily contacts for each person
  # nt: number of days to simulate
  # pinf: proportion of the initial population to randomly start in the I state
  
  # output : nt x 5 matrix contains the number of people across states (S,E,I,R) 
  #          in each days (from 1st day to nt(100th) day). It depicts the 
  #          dynamic change in the number of people in each states based on 
  #          their daily network (from same household and daily contact) and 
  #          random mixing (Irrespective of household or regular network 
  #          relations)
  
  n <- length(beta) #population size
  state <- rep("S", n) #initialize susceptible
  state[sample(1:n, pinf * n)] <- "I" #randomly choose initial in the I state
  
  res <- matrix(0, nrow = nt, ncol = 4) #set up storage for pop in each state
  colnames(res) <- c("S", "E", "I", "R") #naming the column
  
  beta_bar <- mean(beta) #mean sociability parameter
  
  for (t in 1:nt) { #simulate over nt days
    newE <- rep(FALSE, n) #to record who becomes exposed
    #using random deviates: runif(n)
    newI <- (state == "E") & (runif(n) < gamma) #E -> I with prob gamma
    newR <- (state == "I") & (runif(n) < delta) #I -> R with prob delta
    
    #Loop over each infectious person to spread infection
    for (i in which(state == "I")) {
      # household infections: infect susceptible in same household
      hh_members <- which(h == h[i] & state == "S")
      newE[hh_members] <- newE[hh_members]|(runif(length(hh_members))< alpha[1])
      
      # network infections: infect regular contacts from alink
      contacts <- alink[[i]]
      contacts <- contacts[state[contacts] == "S"]
      newE[contacts] <- newE[contacts] | (runif(length(contacts)) < alpha[2])
      
      # random mixing: infect random susceptible
      targets <- which(state == "S")
      #S -> E
      p_random <- alpha[3]*nc*beta[i]*beta[targets]/(beta_bar^2 * (n - 1))
      newE[targets] <- newE[targets] | (runif(length(targets)) < p_random)
    }
    
    #update states based on new infections
    state[newE & state == "S"] <- "E"
    state[newI] <- "I"
    state[newR] <- "R"
    
    #record daily counts of each state
    res[t, ] <- table(factor(state, levels = c("S", "E", "I", "R")))
  }
  
  #returns a list of elements S, E, I, R
  return(data.frame(t = 1:nt, res))
} ## nseir

# 4) Function Plot
dyn.plot <- function(res, title = "SEIR Dynamics") {
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
  
  # ylim = c(0, max(res$S, res$E, res$I, res$R))
  # Plot the states in S(susceptible)
  # Create the base plot
  plot(res$t, res$S, type = "p", ##plot day vs number of suspectible 
       ylim = c(0, max(res$S)), xlim = c(0, max(res$t)), ## set maximum range
       xlab = "Days", ylab = "Number of individuals", ##labelling the plot
       col = "blue", pch = 1, main = title) ##color and type of points in plot
  
  # Add points for each states
  points(res$t, res$E, col = "orange", pch = 1) ##add states 
  points(res$t, res$I, col = "red", pch = 1)
  points(res$t, res$R, col = "darkgreen", pch = 1)
  
  # Add legend
  legend(
    "right",        
    legend = c("Suspectible", "Exposed", "Infectious", "Recovered"),
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
  res1 <- nseir(beta, h, alink)

  # Scenario 2: Random mixing - Setting alpha[1] = alpha[2] = 0 and alpha[3] = 0.04
  res2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))

  # Scenario 3: Full model with constant beta (beta vector = average beta)
  res3 <- nseir(rep(mean(beta), n), h, alink)

  # Scenario 4: Constant beta + random mixing
  res4 <- nseir(rep(mean(beta), n), h, alink, alpha = c(0, 0, 0.04))

  # Plot all scenarios 
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)) 
  dyn.plot(res1, "Full Model Default Parameter")
  dyn.plot(res2, "Random Mixing Only")
  dyn.plot(res3, "Constant Beta")
  dyn.plot(res4, "Constant Beta + Random Mixing")