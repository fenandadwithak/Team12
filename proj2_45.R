# Function Plot
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

# Compare 
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

# --- Plot all scenarios side by side
par(mfrow = c(2, 2))
plot.seir(res1, "Full Model")
plot.seir(res2, "Random Mixing Only")
plot.seir(res3, "Constant Beta")
plot.seir(res4, "Constant Beta + Random Mixing")