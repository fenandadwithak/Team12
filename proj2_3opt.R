nseir_fast <- function(beta, h, alink, alpha = c(.1, .01, .01),
                       delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
  n <- length(beta)
  state <- rep("S", n)
  state[sample.int(n, round(pinf * n))] <- "I"
  
  seir <- matrix(0, nrow = nt, ncol = 4)
  colnames(seir) <- c("S", "E", "I", "R")
  
  beta_bar <- mean(beta)
  
  # Precompute household membership and contact list lookups
  h_groups <- split(seq_len(n), h) #convert h vector into list, per element hh 
  
  # loop over days (1 to 100)
  for (t in seq_len(nt)) {
    # Transition updates (vectorized)
    newI <- (state == "E") & (runif(n) < gamma) 
    newR <- (state == "I") & (runif(n) < delta)
    inf <- which(state == "I") #location infectious person
    sus <- which(state == "S") #location susceptible
    
    newE <- logical(n)
    
    ## 1. Household infections
    # Only check households that contain at least one infected member
    inf_hh <- unique(h[inf]) #ID Household that contain inf. person
    # loop over each elements in defined id household 
    for (hh in inf_hh) {
      members <- h_groups[[hh]] #get all household member from id household
      sus_hh <- members[state[members] == "S"] #filter only "S"
      if (length(sus_hh) > 0) {
        newE[sus_hh] <- (runif(length(sus_hh)) < alpha[1])
      }
    }
    
    ## 2. Network infections from daily contact 
    if (length(inf) > 0) { 
      #search candidate contacts of "I" people with state "S"
      sus_net_all <- unlist(lapply(inf, function(i) {
        contacts <- alink[[i]] #get candidate person to be exposed (from alink)
        if (length(contacts)>0) contacts[state[contacts] == "S"] else integer(0)
      })) ## save into vector sus_net_all, all eligible person to exposed from daily contact
      
      if (length(sus_net_all) > 0) {
        sus_net_all <- unique(sus_net_all) #possibility 2 "I" person have contacted with same person
        newE[sus_net_all] <- newE[sus_net_all] | (runif(length(sus_net_all)) < alpha[2])
      }
    }
    
    ## 3. Random mixing (fully vectorized) ----
    if (length(inf) > 0 && length(sus) > 0) {
      p_infect_sus <- alpha[3] * nc * sum(beta[inf]) * beta[sus] /
        (length(inf) * beta_bar^2 * (n - 1))
      newE[sus] <- newE[sus] | (runif(length(sus)) < p_infect_sus)
    }
    
    ## 4. State updates ----
    state[newE & state == "S"] <- "E"
    state[newI] <- "I"
    state[newR] <- "R"
    
    ## ---- 5. Record counts ----
    seir[t, ] <- tabulate(factor(state, levels = c("S", "E", "I", "R")), nbins = 4)
  }
  
  data.frame(t = 1:nt, seir)
}
