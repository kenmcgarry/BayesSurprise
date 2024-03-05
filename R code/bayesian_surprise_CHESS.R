# bayesian_surprise_CHESS.R
## Wisconsin Card Test

setwd("C:/common_laptop/R-files/BayesianSurprise")

library(PST)
library(TraMineR)
library(stringr)
library(dplyr)

################# load in chess data and transform it for traminer ########## 
chess <- read.csv("chess2.csv")
# chess <- read.csv("chess1.csv")

chess_moves <- as.character(chess$moves[1:280])
chess_wins <- as.character(chess$winner[1:280])
#chess_moves <- as.character(chess$moves)
#chess_wins <- as.character(chess$winner)
chess <- paste0(chess_moves,sep=" ",chess_wins)  #paste0(a, b)
chess <- str_replace_all(chess," ","-")
chess.seq <- TraMineR::seqdef(chess)
head(chess.seq)

seqsample <- sample.int(n = nrow(chess.seq),size=floor(.95*nrow(chess.seq)), replace = FALSE)
chess.seq.train <- chess.seq[seqsample, ]
chess.seq.test  <- chess.seq[-seqsample, ]

# load in rather than recalculate!
load("chessworking1.RData") 
############################################
chess.pst <- pstree(chess.seq.train, L = 5, nmin = 1, ymin = 0.001, with.missing=TRUE)

probtest <- PST::predict(chess.pst, chess.seq.test, decomp=TRUE,output = "prob")
#probtest <- PST::predict(chess.pst, chess.seq.test, decomp=TRUE,output = "prob")

head(probtest)  
# Build the PST on training data
#model6.pst <- pstree(chess.seq.train, L = 10, nmin = 2, ymin = 0.001)

probtrain <- cprob(chess.seq.train , L = 0, prob = TRUE)  # overall probs for training dataset
head(probtrain)

kl <- rep(0,nrow(chess.seq.test))
ent <- rep(0,nrow(chess.seq.test))
surp <- rep(0,nrow(chess.seq.test))
decay <- rep(0,nrow(chess.seq.test))
interest <- rep(0,nrow(chess.seq.test))
x <- seq(1, length.out = nrow(chess.seq.test))  # number of test samples for x-axis


for(i in 1:nrow(chess.seq.test)){
  #probtest <- PST::predict(model6.pst, chess.seq.test[i,], decomp=TRUE,output = "prob")  
  
  ####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION -----------------------------------------------
  # A lot of processing just to attach sequence items (names) to the probababiliies. MAke a dataframe with the
  # columns as the test data sequence. 
  # The 1st record is the probabilities for these test data symbols as deemed by the PST.
  # The 2nd record is the probabilities for the train data. 
  probtestnames <- seqconc(chess.seq.test[i,])
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


#################### now plot sequences in traminer style ################
seqiplot(chess.seq,with.legend="right")
#########################################################################
# Initially, how many patterns are interesting in the test set without decay?
interest[surp > mean(surp)] <- 1

# Using decay, how many patterns are genuinely interesting prior to reaching the
# zero value cut-off? This is a process based on similar patterns reappearing 
# over time, as they are presented. Gather these decay identified patterns and 
# examine why they are deemed interesting

## PLOTS AND STATS FOR PAPER

# original data in a 'wide' format
df <- data.frame(x, ent, kl, surp)

# melt the data to a long format
df2 <- melt(data = df, id.vars = "x")

# do the state sequence plot for the paper
seqiplot(chess.seq,with.legend="right")

# do the entropy histogram for the paper
hist(seqient(chess.seq.test),xlab="entropy",ylim=c(0, 80),main="")

# do the surprise histogram for the paper
hist(surp,xlab="surprise",ylim=c(0, 80),main="")

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
# df2 is a long format contains ent, kl, sur. The magic numbers will be different for each data set
df_ent <- df2[1:50,]
df_kl <- df2[141:190,]
df_sur <- df2[281:330,]
df_sur$x <- 1:nrow(df_sur)

df_all <- rbind(df_ent,df_kl)
df_all <- rbind(df_all,df_sur)


p <- ggplot(data=df_all, aes(x = x, y = value, colour = variable))  +
  geom_line() +
  theme_ipsum() +
  xlab("Records") +
  ylab("Information value") +
  geom_hline(yintercept=mean(surp), linetype='dotted', col = 'black',size=2) +
  annotate("text", x = 15, y = mean(surp), label = "Bayesian 'Wow' Level", vjust = -0.5, size=6)

show(p)  # plot it!

# xtable for latex tables
xtable(head(chess.seq.train[c(20,19,10,36,37,39),]))
seqlength(chess.seq.train[70:80,])
seqlength(chess.seq.train[c(20,19,10,36,37,39),])

### interestingness - via differences, distances and support
# So what are the differences if any between those new sequences flagged as 
# interesting and the others?


# 1. SEQUENCE ENTROPY ----------------------------------------------------
chess.ient <- seqient(chess.seq)
summary(chess.ient)

# 2. SEQUENCE COMPLEXITY (TURBULENCE) -------------------------------------
range(seqST(chess.seq))
summary(seqST(chess.seq))

# 3. LONGEST COMMON PREFIX (LCP) -------------------------------------------
lcp <- seqdist(chess.seq.test, method = "LCP",norm=FALSE)
range(lcp)
median(lcp)

# 4. LONGEST COMMON SUBSEQUENCE (LCS) --------------------------------------
## Normalized LCS distances to the most frequent sequence
chess.dref1 <- seqdist(chess.seq.test, method = "LCS",refseq = 0, norm = "auto")
# LCS distances between two subsets of sequences
set1 <- 1:70
set2 <- 71:140
chess.dref2 <- seqdist(chess.seq.test, method = "LCS",refseq = list(set1,set2))
median(chess.dref2)


# 5. OPTIMAL MATCHING DISTANCE (OMD) and clustering -------------------------
ccost <- seqsubm(chess.seq, method = "CONSTANT", cval = 2)
# compute the distances using the matrix and the default indel cost of 1
chess.OM <- seqdist(chess.seq, method = "OM", indel=3, sm = ccost)
chess.clust <- agnes(chess.OM, diss = TRUE, method = "ward")
plot(chess.clust, which.plots = 2)

cluster2 <- cutree(chess.clust, k = 2)
#cluster2 <- factor(chess.clust, labels = c("Type 1", "Type 2"))

table(cluster2)
seqfplot(chess.seq, group = cluster2, pbarw = T)
round(chess.OM[1:5, 1:5], 1)

# sequence events
seqmtplot(chess.seq, group = cluster2, ylim=c(0,4))

# creates a 16 x 16 cost matrix (we have 16 symbols)
couts <- seqsubm(chess.seq, method = "TRATE")
round(couts, 2)
range(couts)


# 6. EVENT SEQUENCE DISCRIMINATION --------------------------------------------------------
#    create EVENT sequences and discriminate between surprising and the not-surprising sequences
chess.seqeperiod <- seqecreate(chess.seq.test, tevent = "period")
chess.seqe <- seqecreate(chess.seq.test)
chess.seqestate <- seqecreate(chess.seq.test, tevent = "state")

xss <- 10
fsubseq <- seqefsub(chess.seqe, min.support = xss)
fsubseq[1:14]

#fsubseq <- seqefsub(chess.seqe, pmin.support = 0.01)

chessNEW <- data.frame(surp,decay,kl)  #
cohort <- factor(chessNEW$surp > mean(surp), labels = c("surprising","not-surprising")) # divide into two groups
#cohort <- factor(sepsis$birthyr > 1945, labels = c("<=1945",">1945")) # example from traminer

discrimcohort <- seqecmpgroup(fsubseq, group = cohort, method = "bonferroni")
discrimcohort[1:10]

seqdf <- cbind(as.character(discrimcohort$subseq), discrimcohort$data)
seqdf <- seqdf[-c(5,8,9)]  # drop the index, and two residuals
xtable(seqdf[1:10,])

#plot(discrimcohort[1:10,],ylim=c(0,2.0),with.legend=FALSE)

# modify plots
plot(discrimcohort[1:5,],with.legend=FALSE,horiz=FALSE)



