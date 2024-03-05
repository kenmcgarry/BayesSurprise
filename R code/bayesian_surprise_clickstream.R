# bayesian_surprise_clickstream.R

library(clickstream)

# fitting a simple Markov chain and predicting the next click
clickstreams <- c("User1,h,i,c,p,p,c,p,c,h,c,p,p,c,p,p,o",
                  "User2,i,c,i,c,c,c,d",
                  "User3,h,i,c,p,p,c,i,c,p,c,c,p,c,c,i,d",
                  "User4,c,c,p,c,d",
                  "User5,h,i,c,p,p,c,p,p,c,p,p,p,i,p,o",
                  "User6,i,h,i,c,p,p,c,c,p,p,c,p,c,d")
cls <- as.clickstreams(clickstreams, header = TRUE)
mc <- fitMarkovChain(cls,order=2)
startPattern <- new("Pattern", sequence = c("h", "i", "c"))
predict(mc, startPattern ,5)

#
# predict with predefined absorbing probabilities
#
startPattern <- new("Pattern", sequence = c("h", "c"),
                    absorbingProbabilities = data.frame(d = 0.2, o = 0.8))

predict(mc, startPattern, 3)
        
