# francesco.R
# https://www.r-bloggers.com/2019/05/bayesian-models-in-r-2/
# Example of Bayesian updating of a belief given new data.

# Likelihood refers to the probability of observing the data that has been observed assuming that 
# the data came from a specific scenario.

# The posterior can be computed from three key values:
# 1. A likelihood distribution, P(D|M) ;
# 2. A prior distribution, P(M) ;
# 3. The 'average likelihood',
# Example of playing roulette in the casino. 
# Among other things, you can bet on hitting either black (B) or red (r) with supposedly equal 
# probability. For simplification, we assume P(B)=P(r)=0.5 and that we recall the ten past draws 
# before we placed a bet: B, B, r, B, B, B, r, B, B, B 
# So based on past history we have 8/10 coming up f(B)

# Step 1. All possible ways (likelihood distribution)
rangeP <- seq(0, 1, length.out = 100)
plot(rangeP, dbinom(x = 8, prob = rangeP, size = 10),type = "l", xlab = "P(Black)", ylab = "Density")

# Step 2. Update your belief (prior distribution)
# We expressed our belief that P(M) must be 0.5 or close, e.g. P(B) \sim Normal(0.5, 0.1) . 
# For comparison, overlay this prior distribution with the likelihood from the previous step.
lines(rangeP, dnorm(x = rangeP, mean = .5, sd = .1) / 15, col = "red")

# The prior is now shown in red. In the code above, I divided the prior by a constant solely for 
# scaling purposes. Keep in mind that distribution density only matters for the posterior.
# Computing the product between the likelihood and the prior is straightforward, and gives us the 
# numerator from the theorem. The next bit will compute and overlay the unstandardised posterior 
# of P(B) , P(D|M) X P(M) . The usage of a sequence of estimates for P(B) to 
# reconstruct probability distributions is called grid approximation.

likelihood <- dbinom(x = 8, prob = rangeP, size = 10)
prior <- dnorm(x = rangeP, mean = .5, sd = .1)
lines(rangeP, likelihood * prior, col = "green")

# we have successfully used the ten roulette draws (black) to updated our prior (red) into the 
# unstandardised posterior (green). Why is it called 'unstandardised'? The answer comes with 
# the denominator from the theorem.

# Step 3. Make it sum up to one (standardising the posterior)
# An important property of any probability density or mass function is that it integrates to one. 
# P(M|D) = P(D|M) x P(M)  

unstdPost <- likelihood * prior
stdPost <- unstdPost / sum(unstdPost)
lines(rangeP, stdPost, col = "blue")
legend("topleft",legend = c("Likelihood", "Prior", "Unstd Posterior", "Posterior"),text.col = 1:4, bty = "n")



