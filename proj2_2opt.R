get.net <- function (beta, nc=15, h) {
  # Function to create links (daily contacts) for each person
  
  # nc: average number of contacts per person, set to 15
  # beta: sociability parameter, higher means more likely to form links
  # OUTPUT : list of connected networks from each person (output.net)
  
  n <- length(beta) ## calculate total obs of beta, or equal to total population
  beta_bar <- mean(beta) ## calculate mean(beta)/"sociability parameter"
  
  # Allocate the possible connection formed in each individual
  
  # For each individual, we assume an upper bound on the number of contact is 
  # equal to 3*nc (nc=15).This value is selected because the actual number of
  # contacts generated probabilistically (via the link probability, p_link)
  # and may fall short of the target mean number of contact (nc). So, to ensure
  # that the expected number of contact achieved, we use "oversample" 3x target
  # (as reccomended 2-4x target)
 
  ub_nc <- ceiling(nc * 3)  ## the upper bound number of contact per individual
  output.net <- vector("list", n) ## list to store the result
  output.net0 <- integer(n) ## vector to store number of contact made per idv
  
  # Create list of n elements which in each element has "ub_nc" sub-element
  # Later, we would to allocate (drop sub-element and store possible contacts)
  # until overall number of contacts achieved (nc=15)
  for (k in 1:n) output.net[[k]] <- integer(ub_nc) 
  
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
      # if person has "chosen candidates" to link, filter "assigned" sub-element 
      # in "output.net"
      output.net[[k]] <- output.net[[k]][1:output.net0[k]]
    } else {
      # if there is no chosen candidate, replace sub-element into NULL
      output.net[[k]] <- integer(0)
    }
  }
  return(output.net)
} ## get.net
  

alink = get.net(beta,nc=15,h)
# check the actual nc 
mean(sapply(alink, length))
