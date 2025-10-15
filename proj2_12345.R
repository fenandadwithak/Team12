# PROJECT 2 - EXTENDED STATISTICAL PROGRAMMING =================================

# Group 12 
# Aseel Alghamdi : S2901228
# Fenanda Dwitha Kurniasari : S2744048
# Nurmawiya : S2822251

# Aseel : create n people and its household, contact network, cross-check 
#         the entire code and revise code & comments
# Fenanda : create plot function, compare the model, cross-check the entire code
#           and revise code & comments
# Nurmawiya : create nseir function, cross-check code, cross-check the entire 
#             code and revise code & comments

# SEIR (Susceptible, Exposed, Infections, Recovered) Model Simulation
# Goals : create a model and how its use to investigate the role of household
#         and network structure on epidemic dynamics.


# 1) Create n people assigned to corresponding household where maximum size is 5
n <- 1000 #population size
hmax <- 5 #maximum household size
#creating n vector of which the number refers to household
#-and the occurrences refers to size of corresponding household
set.seed(3)
sizes <- sample(1:hmax, n, replace = TRUE)
h <- rep(1:length(sizes), times = sizes)[1:n]

# 2) Contact network
get.net <- function (beta, nc=15, h) {
  # Function to create links (daily contacts) for each person
    
  # nc: average number of contacts per person, set to 15
  # beta: sociability parameter, higher means more likely to form links
  # OUTPUT : list of connected networks from each person (output.net)
    
  n <- length(beta) ## calculate total obs of beta, or equal to total population
  beta_bar <- mean(beta) ## calculate mean(beta)/"sociability parameter"
    
  # Allocate the possible connection formed in each individual
    
  output.net <- vector("list", n) ## list to store the result
  output.net0 <- integer(n) ## vector to store number of contact made per person
    
  # Create list of n elements which in each element has "ub_nc" sub-element
  # Later, we would to allocate (drop/add sub-element and store possible 
  # contacts until overall number of contacts achieved (nc=15)
  for (k in 1:n) output.net[[k]] <- integer(nc) ## set initial value output.net
    
  # Loop over each individual
  for (i in 1:(n-1)) {
    j <- (i+1):n # candidates contacts were generated
    flag <- h[i] != h[j]  # "FALSE" if i & j are belong to same household
      # means that pairing wouldn't be generated
    j <- j[flag] # filter candidates contacts from different household
      
    if (length(j) == 0) next # skip iteration when there are no candidate
      
    # Probability of  a link between persons i and j 
    p_link <- nc * beta[i] * beta[j] / (beta_bar^2 * (n - 1))
      
    # Select the candidates contacts based on p-link
    # linked is a vector containing possible candidate contacts
    linked <- j[runif(length(j)) < p_link]
      
    # If there is chosen candidate(s), store it into output.net and output.net0
    if (length(linked) > 0) {
      # Position to fill
      pos = (output.net0[i]+1):(output.net0[i]+length(linked))
        
      # Assign link[[i]] with number of candidates stored in "linked" at the
      # defined position
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
      # if person has "chosen candidates" to link and 
      # number of chosen candidates < nc, filter "assigned" sub-element 
      # in "output.net"
      output.net[[k]] <- output.net[[k]][1:output.net0[k]]
    } else {
      # if there is no chosen candidate, replace sub-element into NULL
      output.net[[k]] <- integer(0)
    }
  }
  return(output.net)
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
  state[sample(1:n, round(pinf * n))] <- "I" #randomly choose initial in the I
  
  seir <- matrix(0, nrow = nt, ncol = 4) #set up storage for pop in each state
  colnames(seir) <- c("S", "E", "I", "R") #naming the column
  
  beta_bar <- mean(beta) #mean sociability parameter
  
  for (t in 1:nt) {
    #simulate over nt days
    newE <- logical(n) #to record who becomes exposed
    #using random deviates: runif(n)
    newI <- (state == "E") & (runif(n) < gamma) #E -> I with prob gamma
    newR <- (state == "I") & (runif(n) < delta) #I -> R with prob delta
    inf <- which(state == "I") #locating infectious
    
    #Loop over each infectious person to spread infection
    for (i in inf) {
      # household infections: infect susceptible in same household
      hh <- which(h == h[i]) #find everyone in the same household
      sus_hh <- hh[state[hh] == "S"] #keep only those who are susceptible
      if (length(sus_hh)) { #will skip if there's no susceptible
        #store newE for sus_hh if random value < alpha[1]
        newE[sus_hh] <- newE[sus_hh] | (runif(length(sus_hh)) < alpha[1])
        #using OR logical to accumulate infections, not to reset them
      }
      
      # network infections: infect regular contacts from alink
      if (length(alink[[i]])) { #will skip if there's no regular contacts
        #store only the member(s) with regular contacts in S
        sus_net <- alink[[i]][state[alink[[i]]] == "S"]
        if (length(sus_net)) { #will skip if there's no sus_net
          #store newE for sus_net if random value < alpha[2]
          newE[sus_net] <- newE[sus_net] | (runif(length(sus_net)) < alpha[2])
        }
      }
      
      # random mixing: infect random susceptible
      sus <- which(state == "S") #locating susceptible
      if (length(sus)) { #will skip if there is no sus
        #probability irrespective of household or regular contacts
        p_random <- alpha[3]*nc*beta[i]*beta[sus] / (beta_bar^2*(n-1))
        #store newE for sus if random value < p_random
        newE[sus] <- newE[sus] | (runif(length(sus)) < p_random)
      }
    }
    
    #update states based on new infections
    state[newE & state == "S"] <- "E" #store susceptible who get exposed
    state[newI] <- "I" #store exposed to infectious if random value < gamma
    state[newR] <- "R" #store infectious to recovered if random value < delta
    
    #record daily counts of each state
    seir[t, ] <- tabulate(factor(state, levels=c("S", "E", "I", "R")), nbins=4)
  }
  
  #returns a list of elements S, E, I, R
  return(data.frame(t = 1:nt, seir))
} ## nseir

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
  plot(seir$t, seir$S, type = "p", ##plot day vs number of susceptible 
       ylim = c(0, max(seir$S)), xlim = c(0, max(seir$t)), ## set maximum range
       xlab = "Days", ylab = "Number of individuals", ##labelling the plot
       col = "blue", pch = 1, main = title) ##color and type of points in plot
  
  # Add points for each states
  points(seir$t, seir$E, col = "orange", pch = 1) ##add states 
  points(seir$t, seir$I, col = "red", pch = 1)
  points(seir$t, seir$R, col = "darkgreen", pch = 1)
  
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
