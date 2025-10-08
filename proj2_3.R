####1.a########################################################################
#
n <- 1000
hmax <- 5 

set.seed(3) # Fix the randomness so we get the same result each time
house_sizes <- sample(1:hmax, size = ceiling(n / ((hmax + 1) / 2)), replace = TRUE)
#create a sample contain a numbers from 1 to hmax (1:hmax)
# Size since our sample is normal we have the average ((hmax + 1) / 2) and we did that to know how many time we need to repeat the sample that has mean = 3
#OR
# Generate household sizes: randomly pick sizes between 1 and hmax
# The total number of households is estimated by dividing n by the average household size
# Average of uniform(1:hmax) is (1 + hmax)/2, so we divide n by that

sum(house_sizes) #to make sure my sample >= n

# Each household number is repeated according to its size
h <- rep(seq_along(house_sizes), times = house_sizes)[1:n]

# Shuffle the household assignments randomly
h <- sample(h, n, replace = FALSE)
