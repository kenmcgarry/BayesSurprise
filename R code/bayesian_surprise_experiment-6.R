# bayesian_surprise_experiments-6.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 29/12/2021; 

# Note: Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact could
# be very boring!

##### EXPERIMENT 6 - MVAD SEQUENCE ANOMALIES ------------------

data(mvad)
mvad.labels <- c("employment", "further education", "higher education", "joblessness", "school", "training")
mvad.seq <- seqdef(mvad, 15:86, labels = mvad.labels)


set.seed(101) # Set Seed so that same sample can be reproduced in future
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(mvad.seq), size = floor(.75*nrow(mvad.seq)), replace = FALSE)
mvad.seq.train <- mvad.seq[seqsample, ]
mvad.seq.test  <- mvad.seq[-seqsample, ]

# Build the PST on training data
model6.pst <- pstree(mvad.seq.train, L = 10, nmin = 2, ymin = 0.001)

# Test the PST
probtest <- PST::predict(model6.pst, mvad.seq.test , decomp=TRUE,output = "prob")
head(probtest)  

probtrain <- cprob(mvad.seq.train , L = 0, prob = TRUE)  # overall probs for training dataset
head(probtrain)

kl <- rep(0,nrow(mvad.seq.test))
ent <- rep(0,nrow(mvad.seq.test))
surp <- rep(0,nrow(mvad.seq.test))
decay <- rep(0,nrow(mvad.seq.test))
interest <- rep(0,nrow(mvad.seq.test))
x <- seq(1, length.out = nrow(mvad.seq.test))  # number of test samples for x-axis


for(i in 1:nrow(mvad.seq.test)){
  #probtest <- PST::predict(model6.pst, mvad.seq.test[i,], decomp=TRUE,output = "prob")  
  
  ####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION -----------------------------------------------
  # A lot of processing just to attach sequence items (names) to the probababiliies. MAke a dataframe with the
  # columns as the test data sequence. 
  # The 1st record is the probabilities for these test data symbols as deemed by the PST.
  # The 2nd record is the probabilities for the train data. 
  probtestnames <- seqconc(mvad.seq.test[i,])
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
  kl[i] <- round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, epsilon = 0.00001, unit="log10"),digits=2)
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
seqiplot(mvad.seq,with.legend="right")

# do the entropy histogram for the paper
hist(seqient(mvad.seq.test),xlab="entropy",ylim=c(0, 70),main="")

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
xtable(head(mvad.seq.train[c(20,19,10,36,37,39),]))
seqlength(mvad.seq.train[70:80,])
seqlength(mvad.seq.train[c(20,19,10,36,37,39),])


