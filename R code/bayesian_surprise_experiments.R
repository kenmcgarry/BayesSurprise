# bayesian_surprise_experiments.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 24/06/2021; 
# https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.4/topics/KLD
# https://memosisland.blogspot.com/2015/08/practical-kullback-leibler-kl.html

# Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact would
# be very boring!

# Estimate P(M|D)

## Kullback-Leibler (KL) Divergence: Discrete Case 
# Two discrete random variables X and Y, having realizations X=xk and Y=yl, over m and n singletons 
# k=1,...,n and l=1,...,m respectively, are given. Two-way contingency table of realizations of X 
# and Y along the same measured data points, N observations, can be constructed by counting co-occurrences, 
# or cross-classifying factors. Normalizing the resulting frequency table will produce joint and 
# marginal probabilities.


# data for testing KL
freqs1 = c(1/5, 1/5, 3/5)
freqs2 = c(1/10, 4/10, 1/2)

## entropy, joint entropy, conditional entropy
# entropy
cat("entropy using entropy.plugin ",entropy.plugin(freqs1, unit="log2"),"\n")
cat("entropy using entropy.empirical ",entropy.empirical(freqs1, unit="log2"),"\n")
cat("entropy using entropy packagel ",entropy(freqs1, unit="log2"),"\n")


# Joint Entropy
y2d <- discretize2d(freqs1, freqs2, numBins1=3, numBins2=3)
sum(y2d)
entropy(y2d )

# Conditional Entropy
# H(X|Y)
freqsdata <- cbind(freqs1,freqs2)
H <- discretize2d(freqs1,freqs2,numBins1=nrow(freqsdata)^(1/3), numBins2=nrow(freqsdata)^(1/3))
condentropy(H[,1], H[,2], method = "emp")

##------------- Kullback-Liebler experiments -------------
px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
cat("KL using kullback from laplacedemon ", kullback(px,py),"\n")
cat("KL using entropy package ",KL.plugin(px,py),"\n")

##-------------- Jim Albert version ---------------- 
prior <- c(0.25,0.25,0.25,0.25)
names(prior) <- c(0.2,0.25,0.3,0.35)

prior <- rep(0.25, 11)
names(prior) <- 20:30
y <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)

results <- discrete.bayes(dpois,prior,y)
print(results)

# Jim Albert PDF Version
# In this example, suppose one wishes to learn about a baseball player's probability of getting
# a hit p. I believe that a reasonable set of probabilities are 0.20, 0.21, ..., 0.36 and I assign
# these probabilities the corresponding weights 1, 1, 1, 2, 2, 2, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1. 
# Create a vector prior with values 1, 1, 1, 2, 2, 2, ..., then name the probability entries
# with the proportion values. The probability vector is normalized by dividing by its sum.

options(digits = 4) # avoids too many digits
options(width = 60)  # breaks up long outputs into shorter easy to see sections

# observe the player's hitting performance for four periods of 80 at-bats (opportunities) - for
# these four periods, he is 20 for 80, 22 for 80, 19 for 80, and 29 for 80. I place the hit counts
# in the vector y and the sample sizes in the vector n
y <- c(20, 22, 19, 29)
n <- c(80, 80, 80, 80)

prior <- c(1, 1, 1, 2, 2, 2, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1)
names(prior) <- seq(0.2, 0.36, by = 0.01)
prior <- prior/sum(prior)
prior

# obtain the posterior probabilities and the predictive probability by using discrete.bayes
# with arguments dbinom, prior, and y. We add the additional argument size=n which is the
# vector of fixed sample sizes used in the function dbinom. The sampling density is the binomial 
# density dbinom that corresponds to the sampling density
out <- discrete.bayes(dbinom, prior, y, size = n)
posterior <- out$prob
print(posterior)

par(mfrow = c(2, 1))
barplot(prior, main = "Prior Probabilities", xlab = "p", ylim = c(0, 0.2))
barplot(posterior, main = "Posterior Probabilities", xlab = "p")

#-----------------------------------------------------



