# Function Plot
dyn.plot <- function(res, title = "SEIR Dynamics") {
  # Function to create dynamics plot of the simulated population in each states
  
  # x-axis : days to simulate (from 1 to 100)
  # y-axis : number of individuals per states
  #          explanation of y-axis
  #          ---------------------|
  #          at the initial day,number of suspectible individual is the largest
  #          total individuals in all states S+E+I+R = 1000
  #          these numbers would dynamically change day by day in each states
  #            depend on the defined models
  
  # input : nseir model from function nseir 
  # output : Dynamic Plot days vs number of individuals per states
  
  #ylim = c(0, max(res$S, res$E, res$I, res$R))
  # Plot the states in S(suspectible)
  # plot the initial day, where 
  plot(res$t, res$S, type = "p", ylim = c(0, max(res$S, res$E, res$I, res$R)),
       xlab = "Days", ylab = "Number of individuals",
       col = "blue", pch = 1, main = title)
  
  # Add points for each compartment
  points(res$t, res$E, col = "orange", pch = 1)
  points(res$t, res$I, col = "red", pch = 1)
  points(res$t, res$R, col = "darkgreen", pch = 1)
  
  # Add legend
  legend("right", legend = c("S", "E", "I", "R"),
         col = c("blue", "orange", "red", "darkgreen"), pch =1)
}


# Compare 
set.seed(123)
beta <- runif(n)
alink <- get.net(beta, nc = 15, h = h)

# Scenario 1: Full model
res1 <- nseir(beta, h, alink)

# Scenario 2: Random mixing only
res2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))

# Scenario 3: Constant beta only
res3 <- nseir(rep(mean(beta), n), h, alink)

# Scenario 4: Constant beta + random mixing
res4 <- nseir(rep(mean(beta), n), h, alink, alpha = c(0, 0, 0.04))

# --- Plot all scenarios side by side
par(mfrow = c(2, 2))
dyn.plot(res1, "Full Model")
dyn.plot(res2, "Random Mixing Only")
dyn.plot(res3, "Constant Beta")
dyn.plot(res4, "Constant Beta + Random Mixing")