# bayesian_surprise_experiments-1.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 23/11/2021; 

# Note: Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact could
# be very boring!

##### EXPERIMENT 5 - SEPSIS SEQUENCE ANOMALIES ------------------

data(sepsis)  # from seqDetect package but use my RData file

load("sepsis.seq.RData")  # The data in sepsis took a lot of preprocessing to get in a fit state for PST !!
                         # See bayesian_surprise_make-data-experiment5.R

set.seed(101) # Set Seed so that same sample can be reproduced in future
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(sepsis.seq), size = floor(.75*nrow(sepsis.seq)), replace = FALSE)
sepsis.seq.train <- sepsis.seq[seqsample, ]
sepsis.seq.test  <- sepsis.seq[-seqsample, ]

# Build the PST on training data
model5.pst <- pstree(sepsis.seq.train, L = 10, nmin = 2, ymin = 0.001)


# Test the PST
probtest <- PST::predict(model5.pst, sepsis.seq.test , decomp=TRUE,output = "prob")
head(probtest)  

probtrain <- cprob(sepsis.seq.train , L = 0, prob = TRUE)  # overall probs for training dataset
head(probtrain)

kl <- rep(0,nrow(sepsis.seq.test))
ent <- rep(0,nrow(sepsis.seq.test))
surp <- rep(0,nrow(sepsis.seq.test))
decay <- rep(0,nrow(sepsis.seq.test))
interest <- rep(0,nrow(sepsis.seq.test))
x <- seq(1, length.out = nrow(sepsis.seq.test))  # number of test samples for x-axis


for(i in 1:nrow(sepsis.seq.test)){
  #probtest <- PST::predict(model5.pst, sepsis.seq.test[i,], decomp=TRUE,output = "prob")  
  
  ####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION -----------------------------------------------
  # A lot of processing just to attach sequence items (names) to the probababiliies. MAke a dataframe with the
  # columns as the test data sequence. 
  # The 1st record is the probabilities for these test data symbols as deemed by the PST.
  # The 2nd record is the probabilities for the train data. 
  probtestnames <- seqconc(sepsis.seq.test[i,])
  probtestnames <- seqdecomp(probtestnames)
  
  data1 <- t(na.omit(as.data.frame(probtest[i,]) ))
  
  colnames(data1) <- probtestnames
  probcompare <- as.data.frame(data1)
  probcompare <- rbind(probcompare,rep(1:length(data1)))
  rownames(probcompare) <- c("testseq","trainseq")
  
  for(j in 1:length(data1)){
    sym <- colnames(data1)[j]
    
    indx <- which(colnames(data1) == sym) 
    valindex <- which(colnames(data1) == sym)
    
    probcompare[2,indx] <- probtrain[valindex]  # get the overall prob for those syms
    
  }
  
  probcompare[is.na(probcompare)] <- 0
  
  ### ------------------------------------------------------------------------------------------------------------
  # Compare probabilities from output from learned PST on new test data and train data probs.
  # These reside in the probcompare dataframe.
  
  #print(probcompare)
  
  tempprob <- probcompare  # make a temp structure just in case of errors
  tempprob[tempprob > 0.5] <- 1
  tempprob[tempprob < 0.5] <- 0
  #print(tempprob)
  
  smat <- as.matrix(tempprob)  # convert to matrix to access contents
  success <- sum(smat[1,])  # add up the 1's in testseq probabilities
  scount <- length(tempprob)
  
  #print(success)
  
  # Now produce freq counts and then cast problem as per advice from Francesco...
  # https://www.r-bloggers.com/2019/05/bayesian-models-in-r-2/
  # Step 1. All possible ways (likelihood distribution)
  rangeP <- seq(0, 1, length.out = 100)
  #plot(rangeP, dbinom(x = success, prob = rangeP, size = scount),type = "l", xlab = "P(Black)", ylab = "Density")
  
  # Step 2. State your belief (prior distribution)
  #lines(rangeP, dnorm(x = rangeP, mean = .3, sd = .1) / 15, col = "red")
  likelihood <- dbinom(x = success, prob = rangeP, size = scount)
  prior <- dnorm(x = rangeP, mean = .5, sd = .1)
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
  
  #cat("\nKullback-Leibler distance of ",round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, unit="log2"),digits=2), 
  #    " nats between prior and posterior distributions")
  
  #kl[i] <- KL.plugin(stdPost,stdPrior,unit="log2")
  kl[i] <- round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, unit="log10"),digits=2)
  ent[i] <- shannon.conditional.entropy(stdPost,stdPrior)  
  #shannon.conditional.entropy(stdPost,stdPrior)  
  #ent[i] <- entropy.Dirichlet(stdPost,stdPrior,unit="log10")
  ###surp[i] <- kl[i]/2   ###
  surp[i] <- mean(unstdPost) - mean(stdPrior)  ###
  decay[i] <- 1/(mean(stdPrior)*i) * log(1 - 1/((mean(stdPrior)*i) + (mean(unstdPost)*i)))
  
}   # end of my massive outer for loop


# Initially, how many patterns are interesting in the test set without decay?
interest[surp > mean(surp)] <- 1

# Using decay, how many patterns are genuinely interesting prior to reaching the zero value cut-off?
# This is a process based on similar patterns reappearing over time, as they are presented.

# Gather these decay identified patterns and examine why they are deemed interesting



## PLOTS AND STATS FOR PAPER

# original data in a 'wide' format
df <- data.frame(x, ent, kl, surp)

# melt the data to a long format
df2 <- melt(data = df, id.vars = "x")

# do the state sequence plot for the paper
seqiplot(sepsis.seq,with.legend="right")

# do the entropy histogram for the paper
hist(seqient(sepsis.seq.test),xlab="entropy",ylim=c(0, 70),main="")

# do the surprise histogram for the paper
hist(surp,xlab="surprise",ylim=c(0, 70),main="")

# do the decay plot of surprise for the paper
df3 <- data.frame(x,decay)
df3 <- na.omit(df3)


# wacky plot for decay
p <- ggplot(data=df3, aes(x = x, y = decay))  +
  geom_line() +
  theme_ipsum() +
  xlab("N") +
  ylab("Surprise decay") 
show(p) # plot it!

# A wacky way of plotting, but RStudio has trouble with any ggplot2 plots using "source" commands
# they will simply not appear unless you do it this way!
p <- ggplot(data=df2, aes(x = x, y = value, colour = variable))  +
       geom_line() +
       theme_ipsum() +
       xlab("Time") +
       ylab("Information value") +
       geom_hline(yintercept=mean(surp), linetype='dotted', col = 'black',size=2) +
       annotate("text", x = 100, y = mean(surp), label = "Bayesian 'Wow' Level", vjust = -0.5, size=6)

show(p)  # plot it!

# xtable for latex tables
xtable(head(sepsis.seq.train[c(20,19,10,36,37,39),]))
seqlength(sepsis.seq.train[70:80,])
seqlength(sepsis.seq.train[c(20,19,10,36,37,39),])

### interestingness - via differences, distances and support
# So what are the differences if any between those new sequences flagged as interesting and the others?


# 1. SEQUENCE ENTROPY ----------------------------------------------------
sepsis.ient <- seqient(sepsis.seq)
summary(sepsis.ient)

# 2. SEQUENCE COMPLEXITY (TURBULENCE) -------------------------------------
range(seqST(sepsis.seq))
summary(seqST(sepsis.seq))

# 3. LONGEST COMMON PREFIX (LCP) -------------------------------------------
lcp <- seqdist(sepsis.seq.test, method = "LCP",norm=FALSE)
range(lcp)
median(lcp)

# 4. LONGEST COMMON SUBSEQUENCE (LCS) --------------------------------------
## Normalized LCS distances to the most frequent sequence
sepsis.dref1 <- seqdist(sepsis.seq.test, method = "LCS",refseq = 0, norm = "auto")
# LCS distances between two subsets of sequences
set1 <- 1:125
set2 <- 126:251
sepsis.dref2 <- seqdist(sepsis.seq.test, method = "LCS",refseq = list(set1,set2))
median(sepsis.dref2)


# 5. OPTIMAL MATCHING DISTANCE (OMD) and clustering -------------------------
ccost <- seqsubm(sepsis.seq, method = "CONSTANT", cval = 2)
# compute the distances using the matrix and the default indel cost of 1
sepsis.OM <- seqdist(sepsis.seq, method = "OM", indel=3, sm = ccost)
sepsis.clust <- agnes(sepsis.OM, diss = TRUE, method = "ward")
plot(sepsis.clust, which.plots = 2)

cluster2 <- cutree(sepsis.clust, k = 2)
cluster2 <- factor(sepsis.clust, labels = c("Type 1", "Type 2"))

table(cluster2)
seqfplot(sepsis.seq, group = cluster2, pbarw = T)
round(sepsis.OM[1:5, 1:5], 1)

# sequence events
seqmtplot(sepsis.seq, group = cluster2, ylim=c(0,4))

# creates a 16 x 16 cost matrix (we have 16 symbols)
couts <- seqsubm(sepsis.seq, method = "TRATE")
round(couts, 2)
range(couts)


# 6. EVENT SEQUENCE DISCRIMINATION --------------------------------------------------------
#    create EVENT sequences and discriminate between surprising and the not-surprising sequences
sep.seqeperiod <- seqecreate(sepsis.seq.test, tevent = "period")
sep.seqe <- seqecreate(sepsis.seq.test)
sep.seqestate <- seqecreate(sepsis.seq.test, tevent = "state")

fsubseq <- seqefsub(sep.seqe, min.support = 20)
fsubseq[1:20]

#fsubseq <- seqefsub(sep.seqe, pmin.support = 0.01)

sepsisNEW <- data.frame(surp,decay,kl)  #
cohort <- factor(sepsisNEW$surp > mean(surp), labels = c("surprising","not-surprising")) # divide into two groups
#cohort <- factor(sepsis$birthyr > 1945, labels = c("<=1945",">1945"))

discrimcohort <- seqecmpgroup(fsubseq, group = cohort, method = "bonferroni")
discrimcohort[1:10]

seqdf <- cbind(as.character(discrimcohort$subseq), discrimcohort$data)
seqdf <- seqdf[-c(5,8,9)]  # drop the index, and two residuals
xtable(seqdf[1:10])

plot(discrimcohort[1:10],ylim=c(0,1.0))






