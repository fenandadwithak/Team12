start.time <- Sys.time()

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

# STEP 4.a =====================================================================
loc_1 <- grep("[",a,fixed=TRUE) #location/index of [
length_a <- length(a)

all_loc <- NULL
for (i in loc_1){
  loc_2 <- grep("]",a[i:min((i+100),length_a)],fixed=TRUE)
  loc_3 <- grep(".",a[i:min((i+100),length_a)],fixed=TRUE)
  
  if(length(loc_2)>0){
    loc_exc <- c(i:(loc_2[1]+(i)-(1)))
  } else{loc_exc <- c(i:(loc_3[1]+(i)-(1)))}
  
  ifelse (i==loc_1[1], 
          all_loc <- loc_exc, 
          all_loc <- append(all_loc,loc_exc))
}

a <- a[-(all_loc)]

# STEP 4.b & 4.c================================================================
#remove fully uppercase letter exclude I and A, and remove numbers
upnum_loc <- which(toupper(a)==a & !(a %in% c("I", "A")) | grepl("[0-9]", a))
a <- a[-(upnum_loc)]

#remove all underscore, dash, parentheses, asterisk
a <- gsub("[*()_-]", "",a)

# STEP 4.d - 4.f================================================================
#punctuation split and lower case function
split_punct <- function (x){
  punct <- c(",", ".", ";", "!", ":", "?")
  for (i in punct) {
    x <- gsub(paste0('[', i, ']'), paste0("#", i), x)
  }
  x <- tolower(unlist(strsplit(x, "[#]", perl = TRUE)))
}

a <- split_punct(a)

#write.table(a,"cleaned_a.txt",sep="\t",row.names=FALSE) #result check

# STEP 5========================================================================
b <- unique(a)
freq <- tabulate(match(a,b))
b <- which(rank(-freq) <= 1000) #average ties method, rank 1 = words
#final dataset of b contains indices of top ~1000 from the unique words

# STEP 6========================================================================
b_word <- a[b]
n <- length(a)
mlag <- 4
mrow <- n - mlag
mcol <- mlag + 1
M1 <- match(a,b_word) #M1 = token

M <- matrix(NA, mrow, mcol)
for (i in 0:4) {
  M[,i+1] <- M1[(i+1):(mrow+i)]
}

# STEP 7-9======================================================================
next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  
  k.match <- match(key,b_word)
  loc.key <- which(is.finite(k.match))
  key.n <- k.match[loc.key]
  ##set.seed(1)
  ## If key is too long, use only the last mlag elements
  if (length(key.n) > mlag) key.n <- tail(key.n, mlag)
  
  u_all <- c()          # will store candidate next-word tokens
  p_all <- c()          # will store corresponding probabilities
  
  ## Loop over i = 1...length(key)
  ## Each iteration uses shorter and shorter suffix of the key
  for (i in seq_along(key.n)) {
    ## Columns of M to match: from (mlag - length(key) + i) to mlag
    mc <- mlag - (length(key.n) - i)  # start column
    ## Find rows where M[,cols] exactly matches the current suffix of key
    ii <- colSums(!(t(M[, mc:mlag, drop = FALSE]) == key.n[i:length(key.n)]))
    row.match <- which(ii == 0 & is.finite(ii))
    
    if (length(row.match) > 0) {
      ## Get the (mlag+1)-th column for matching rows: the "next token"
      u <- M[row.match, mlag + 1]
      u <- u[!is.na(u)]
      ## Compute probability weight for this suffix
      prob <- rep(w[length(key.n) - i + 1] / length(u), length(u))
      
      ## Store results
      u_all <- c(u_all, u)
      p_all <- c(p_all, prob)
    }
  }
  
  ## If no matches found, sample a random common word according to freq in M1
  if (length(u_all) == 0) {
    tab <- table(M1[!is.na(M1)])
    return(sample(as.integer(names(tab)), size = 1, prob = tab))
  }
  
  ## Otherwise sample next token according to combined probabilities
  return (sample(u_all, size = 1, prob = p_all))
}


femael.predict <- function(M, M1) {
  repeat {
    key <- readline(prompt = "Please input the key: ")
    
    if (length(key) > 0 && 
        is.na(suppressWarnings(as.numeric(key))) &&
        key!="") {
      
      # Generate words until we reach 5 tokens total
      while (length(key) < 5) {
        key <- unlist(strsplit(key, " "))
        key <- split_punct(key)
        nxt <- next.word(key, M, M1)
        key <- c(key, a[nxt])  # append predicted word
      }
      
      cat("The result is:\n")
      print(paste(key, collapse = " "))
      break #Exit loop if condition is satisfied
    } else {
      cat("Invalid input. Please input another key.\n")
    }
  }
}

femael.predict(M,M1)
romeo
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken