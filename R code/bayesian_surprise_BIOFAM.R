# bayesian_surprise_BIOFAM.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 01/03/2024; 

# Note: Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact could
# be very boring!

## EXPERIMENT 4 - Swiss household data
data("biofam", package = "TraMineR")

biofam.labels <- c("Parent","Left","Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced" )
biofam.seq <- seqdef(biofam, 10:25, labels = biofam.labels)
#seqdplot(biofam.seq)

biofam.short.labs <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
biofam.seq <- seqdef(biofam, 10:25, alphabet = 0:7,states = biofam.short.labs, labels = biofam.labels)

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(biofam.seq), size = floor(.75*nrow(biofam.seq)), replace = FALSE)
biofam.seq.train <- biofam.seq[seqsample, ]
biofam.seq.test  <- biofam.seq[-seqsample, ]

biofam.pst <- pstree(biofam.seq.train, L = 10, nmin = 2, ymin = 0.001)

biofam.seq.test <- impute(biofam.pst, biofam.seq.test, method="prob")  # get overall probs for each symbol

# Test the PST
probtest <- PST::predict(biofam.pst, biofam.seq.test, decomp=TRUE,output = "prob") # 
probtrain <- cprob(biofam.seq.train, L = 0, prob = TRUE)  # overall probs for training dataset

kl <- rep(0,nrow(biofam.seq.test))
ent <- rep(0,nrow(biofam.seq.test))
surp <- rep(0,nrow(biofam.seq.test))
decay <- rep(0,nrow(biofam.seq.test))
interest <- rep(0,nrow(biofam.seq.test))
x <- seq(1, length.out = nrow(biofam.seq.test))  # number of test samples for x-axis
  
for(i in 1:nrow(biofam.seq.test)){
  probtest <- PST::predict(biofam.pst, biofam.seq.test[i,], decomp=TRUE,output = "prob")  
  
####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION -----------------------------------------------
# A lot of processing just to attach sequence items (names) to the probababiliies. MAke a dataframe with the
# columns as the test data sequence. 
# The 1st record is the probabilities for these test data symbols as deemed by the PST.
# The 2nd record is the probabilities for the train data. 
probtestnames <- seqconc(biofam.seq.test[i,])
probtestnames <- seqdecomp(probtestnames)
colnames(probtest) <- probtestnames
probcompare <- as.data.frame(probtest)
probcompare <- rbind(probcompare,rep(1:length(probtest)))
rownames(probcompare) <- c("testseq","alphabet")

for(j in 1:length(probtest)){
  sym <- colnames(probtest)[j]
  
  indx <- which(colnames(probtest) == sym) 
  valindex <- which(colnames(probtest) == sym)
  
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
kl[i] <- round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, unit="log10",epsilon=0.001),digits=2)
ent[i] <- shannon.conditional.entropy(stdPost,stdPrior)  
#shannon.conditional.entropy(stdPost,stdPrior)  
#ent[i] <- entropy.Dirichlet(stdPost,stdPrior,unit="log10")
#surp2[i] <- kl[i]/2
surp[i] <- mean(unstdPost) - mean(stdPrior)  ###
decay[i] <- 1/(mean(stdPrior)*i) * log(1 - 1/((mean(stdPrior)*i) + (mean(unstdPost)*i)))

}   # end of my massive outer for loop

## PLOTS AND STATS FOR PAPER

# original data in a 'wide' format
df <- data.frame(x, ent, kl, surp)

# melt the data to a long format
df2 <- melt(data = df, id.vars = "x")

# do the state sequence plot for the paper
seqiplot(biofam.seq,with.legend="right")

# do the entropy histogram for the paper
hist(seqient(biofam.seq.train),xlab="entropy",ylim=c(0, 500),main="")

# do the surprise histogram for the paper
hist(surp,xlab="surprise",ylim=c(0, 70),main="")

# do the decay plot of surprise for the paper
df3 <- data.frame(x,decay)
df3 <- na.omit(df3)

df_ent <- df2[1:50,]
df_kl <- df2[501:550,]
df_sur <- df2[1451:1500,]
df_kl$x <- 1:nrow(df_sur)
df_sur$x <- 1:nrow(df_sur)

df_all <- rbind(df_ent,df_kl)
df_all <- rbind(df_all,df_sur)

# A wacky way of plotting, but RStudio has trouble with any ggplot2 plots using "source" commands
# they will simply not appear unless you do it this way!
p <- ggplot(data=df_all, aes(x = x, y = value, colour = variable))  +
  geom_line() +
  theme_ipsum() +
  xlab("Records") +
  ylab("Information value") +
  geom_hline(yintercept=mean(surp), linetype='dotted', col = 'black',size=2) +
  annotate("text", x = 20, y = mean(surp), label = "Bayesian 'Wow' Level", vjust = -0.5, size=6)

show(p)

# wacky plot for decay
p <- ggplot(data=df3, aes(x = x, y = decay))  +
  geom_line() +
  theme_ipsum() +
  xlab("N") +
  ylab("Surprise decay") 
  #geom_hline(yintercept=mean(surp), linetype='dotted', col = 'black',size=2) 
  #annotate("text", x = 100, y = mean(surp), label = "Bayesian 'Wow' Level", vjust = -0.5, size=6)

show(p)

# xtable function for generating LaTex tables
xtable(head(biofam.seq.train))

### interestingness - via differences, distances and support
# So what are the differences if any between those new sequences flagged as interesting and the others?


# 1. SEQUENCE ENTROPY ----------------------------------------------------
biofam.ient <- seqient(biofam.seq)
summary(biofam.ient)

# 2. SEQUENCE COMPLEXITY (TURBULENCE) -------------------------------------
range(seqST(biofam.seq))
summary(seqST(biofam.seq))

# 3. LONGEST COMMON PREFIX (LCP) -------------------------------------------
lcp <- seqdist(biofam.seq.test, method = "LCP",norm=FALSE)
range(lcp)

# 4. LONGEST COMMON SUBSEQUENCE (LCS) --------------------------------------
## Normalized LCS distances to the most frequent sequence
biofam.dref1 <- seqdist(biofam.seq.test, method = "LCS",refseq = 0, norm = "auto")
# LCS distances between two subsets of sequences
set1 <- 1:250
set2 <- 251:500
biofam.dref2 <- seqdist(biofam.seq.test, method = "LCS",refseq = list(set1,set2))

# 5. OPTIMAL MATCHING DISTANCE (OMD) and CLUSTERING -------------------------

# 6. EVENT SEQUENCE DISCRIMINATION --------------------------------------------------------
#    create EVENT sequences and discriminate between surprising and the not-surprising sequences
bf.seqeperiod <- seqecreate(biofam.seq.test, tevent = "period")
bf.seqe <- seqecreate(biofam.seq.test)
bf.seqestate <- seqecreate(biofam.seq.test, tevent = "state")

fsubseq <- seqefsub(bf.seqe, min.support = 20)
fsubseq[1:10]

fsubseq <- seqefsub(bf.seqe, pmin.support = 0.01)
biofam$surp <- surp  #
biofam$decay <- decay #
biofam$kl <- kl  #

biofamNEW <- data.frame(surp,decay,kl)  #

cohort <- factor(biofamNEW$surp > mean(surp), labels = c("surprising","not-surprising")) # divide into two groups

discrimcohort <- seqecmpgroup(fsubseq, group = cohort, method = "bonferroni")
discrimcohort[1:10]

# create the table for the paper
seqdf <- cbind(as.character(discrimcohort$subseq), discrimcohort$data)
seqdf <- seqdf[-c(5,8,9)]  # drop the index, and two residuals
xtable(seqdf)

# plot the discrimination plot between surprising and not-surprising sequences
#plot(discrimcohort[1:10],ylim=c(0,1.0),legend.cex=0.0001,legend.title="")
plot(discrimcohort[1:5,],with.legend=FALSE,horiz=FALSE)

