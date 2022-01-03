# bayesian_surprise_experiments-1.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 23/09/2021; 


## EXPERIMENT 2 - simulated data, motifs already in place. Main purpose is to build PST with a larger
# dataset of several similar sequences and test with several new sequences of motifs.

seq.motifs.train2 <- TraMineR::seqdef(c("M6-M2-M3-M4-M6-M1-M2-M3-M4",
                                        "M6-M5-M2-M3-M4-M5-M1-M2-M3-M4-M5-M2-M3-M4-M5-M1-M2-M6",
                                        "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M3-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M3-M5-M6-M2-M3-M4",
                                        "M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                                        "M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M6",
                                        "M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M6",
                                        "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6"))

seq.motifs.test2 <- TraMineR::seqdef(c("M3-M6-M4-M5-M1-M2-M3-M4-M5-M6-M1-M5-M6-M2-M3-M4-M1",
                                       "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                                       "M6-M5-M6",
                                       "M4-M5-M1-M2-M3-M4-M5-M1",
                                       "M5-M6-M6-M6-M1-M2-M1-M2-M1-M2-M1-M2-M1-M2"))

# Build the PST on training data
model2.pst <- pstree(seq.motifs.train2,nmin=1,with.missing = TRUE)
summary(model2.pst)

# Test the PST
whichone <- 1
testseq <- seqdef(seq.motifs.test2[whichone])
probtest <- PST::predict(model2.pst, seq.motifs.test2, decomp=TRUE,output = "prob")
probtest  

probtrain <- cprob(seq.motifs.train2 , L = 0, prob = TRUE)  # overall probs for training dataset
probtrain

####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION ----------------------------------
# A lot of processing just to attach sequence items (names) to the probababiliies. 
# The 1st record is the probabilities for these test data symbols as deemed by the PST.
# The 2nd record is the probabilities for the train data. 
probtestnames <- strsplit(seq.motifs.test2[1,], "-")
colnames(probtest) <- unlist(probtestnames)

probcompare <- as.data.frame(probtest)
probcompare <- rbind(probcompare,rep(1:length(probtest)))
rownames(probcompare) <- c("testseq","alphabet")

for(i in 1:length(probtest)){
  sym <- colnames(probtest)[i]
  
  indx <- grepl(sym, colnames(probtest))
  valindex <- grepl(sym, colnames(probtrain))
  probcompare[2,indx] <- probtrain[valindex]
  
}

### ------------------------------------------------------------------------------------------------------------

# Compare probabilities from output from learned PST on new test data and train data probs.
# These reside in the probcompare dataframe.

print(probcompare)

tempprob <- probcompare  # make a temp structure just in case of errors
tempprob[tempprob > 0.5] <- 1
tempprob[tempprob < 0.5] <- 0
print(tempprob)

smat <- as.matrix(tempprob)  # convert to matrix to access contents
success <- sum(smat[1,])  # add up the 1's in testseq probabilities
scount <- length(tempprob)

# Now produce freq counts and then cast problem as per advice from Francesco...
# https://www.r-bloggers.com/2019/05/bayesian-models-in-r-2/
# Step 1. All possible ways (likelihood distribution)
rangeP <- seq(0, 1, length.out = 100)
#plot(rangeP, dbinom(x = success, prob = rangeP, size = scount),type = "l", xlab = "P(Black)", ylab = "Density")

# Step 2. State your belief (prior distribution)
#lines(rangeP, dnorm(x = rangeP, mean = .3, sd = .1) / 15, col = "red")
likelihood <- dbinom(x = success, prob = rangeP, size = scount)
prior <- dnorm(x = rangeP, mean = .3, sd = .1)
#lines(rangeP, likelihood * prior, col = "green")

# Step 3. Make it sum up to one (standardising the posterior)
# An important property of any probability density or mass function is that it integrates to one. 
# P(M|D) = P(D|M) x P(M)  

unstdPost <- likelihood * prior
stdPost <- unstdPost / sum(unstdPost)
stdPrior <- prior / sum(prior)
#lines(rangeP, stdPost, col = "blue")
#legend("topleft",legend = c("Likelihood", "Prior", "unstdPost","Posterior"),text.col = 1:4, bty = "n")

# calculate the difference between prior and posterior using kullback-leibler divergence
cat("\nKullback-Leibler distance of ",
    round(kullback_leibler_distance(stdPost,stdPrior,testNA=F, unit="log2"),digits=2), 
    " bits (log2) between prior and posterior distributions")




