# bayesian_surprise_WCST.R
## Wisconsin Card Sorting Test
## data generated by https://www.psytoolkit.org/experiment-library/wcst.html
## each participant receives 60 cards varying by rules (colour, shape or number). 
## The participants have to figure out the rules which change over time.

setwd("C:/common_laptop/R-files/BayesianSurprise")

library(PST)
library(TraMineR)
library(readxl)
library(dplyr)
library(TraMineR)
library(stringr)

############# read in all excel files ############################ 
file.list <- list.files(pattern='*.xlsx')
df.list <- lapply(file.list, read_excel)

############ join them by rows ##################################
df.wcst <- bind_rows(df.list, .id = "id")
wcst <- as.data.frame(df.wcst)  # conversion from tibble to dataframe makes it a bit easier

############ select the key variables and transform strings ###########
##
## columns 10 (status), 11 (actual card clicked), 12 Error
## Status 1=correct, 2=wrong card, 3=too slow
## card clicked [between 1-4]
## Error (1=error 0=correct)
##
## columns 6 (card description: shape)
## columns 7 (number of symbols)
## columns 8 (card colour)

### convert numerics to strings
wcst$Status <- ifelse(wcst$Status==1,"correct",ifelse(wcst$Status==2,"wrong","tooslow"))

### e.g. star;yellow;2TwoSymbols;outcome
test <- wcst$`Number of symbols`
test <-ifelse(test==1,"OneSymbol",ifelse(test==2,"TwoSymbols",ifelse(test==3,"ThreeSymbols",ifelse(test==4,"FourSymbols","FourSymbols"))))
wcst$`Number of symbols` <- test # update the main dataframe

#save(wcst,file = "wcst.RData")  ##save again after adding more data
##################################################################################
## Now form sequences
load("wcst.RData")

tempstr<- vector(mode = "list", length = nrow(wcst))

for(i in 1:nrow(wcst)){
  seqraw <- wcst[i,c(6,7,8,9,11)]  
  seqraw <- toString(seqraw)
  seqraw <- paste0(seqraw,sep="") 
  seqraw <- str_replace_all(seqraw," ","-")
  seqraw <- str_replace_all(seqraw,",","")
  tempstr[i] <- seqraw
}

wcst.seq <- TraMineR::seqdef(unlist(tempstr))

head(wcst.seq)
#wcst.seq.sample <- TraMineR::seqdef(wcst.seq.sample) # convert it into a seq object

seqsample <- sample.int(n = nrow(wcst.seq), size = floor(.9*nrow(wcst.seq)), replace = FALSE)
wcst.seq.train <- wcst.seq[seqsample, ]
wcst.seq.test  <- wcst.seq[-seqsample, ]

wcst.pst <- pstree(wcst.seq.train,L=4)#,L = 5, nmin = 1, ymin = 0.01)
# Test the PST
probtest <- PST::predict(wcst.pst, wcst.seq.test, decomp=TRUE,output = "prob")
head(probtest)  

################################### THE LOOP FROM HELL ########################
probtrain <- cprob(wcst.seq.train , L = 0, prob = TRUE)  # overall probs for training dataset
head(probtrain)

kl <- rep(0,nrow(wcst.seq.test))
ent <- rep(0,nrow(wcst.seq.test))
surp <- rep(0,nrow(wcst.seq.test))
decay <- rep(0,nrow(wcst.seq.test))
interest <- rep(0,nrow(wcst.seq.test))
x <- seq(1, length.out = nrow(wcst.seq.test))  # number of test samples for x-axis


for(i in 1:nrow(wcst.seq.test)){
  #probtest <- PST::predict(model5.pst, sepsis.seq.test[i,], decomp=TRUE,output = "prob")  
  
  ####  THIS BLOCK OF CODE BLOCK NEEDS TO RECODED AS A FUNCTION -----------------------------------------------
  # A lot of processing just to attach sequence items (names) to the probababiliies. MAke a dataframe with the
  # columns as the test data sequence. 
  # The 1st record is the probabilities for these test data symbols as deemed by the PST.
  # The 2nd record is the probabilities for the train data. 
  probtestnames <- seqconc(wcst.seq.test[i,])
  probtestnames <- seqdecomp(probtestnames)
  
  data1 <- t(na.omit(as.data.frame(probtest[i,]) ))
  #data1 <- t((as.data.frame(probtest[i,]) ))
  
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
  
  unstdPost <- likelihood * prior + (runif(n=length(likelihood), min=0.001, max=0.05))
  stdPost <- unstdPost / sum(unstdPost)
  stdPrior <- prior / sum(prior)
  #lines(rangeP, stdPost, col = "blue")
  #legend("topleft",legend = c("Likelihood", "Prior", "unstdPost","Posterior"),text.col = 1:4, bty = "n")
  # calculate the difference between prior and posterior using kullback-leibler divergence
  
  #cat("\nKullback-Leibler distance of ",round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, unit="log2"),digits=2), 
  #    " nats between prior and posterior distributions")
  
  #kl[i] <- KL.plugin(stdPost,stdPrior,unit="log2")
  kl[i] <- round(kullback_leibler_distance(stdPrior,stdPost,testNA=F, unit="log10",epsilon = 0.00001),digits=2)
  ent[i] <- shannon.conditional.entropy(stdPost,stdPrior)  
  #shannon.conditional.entropy(stdPost,stdPrior)  
  #ent[i] <- entropy.Dirichlet(stdPost,stdPrior,unit="log10")
  ###surp[i] <- kl[i]/2   ###
  surp[i] <- mean(unstdPost) - mean(stdPrior)  ###
  decay[i] <- 1/(mean(stdPrior)*i) * log(1 - 1/((mean(stdPrior)*i) + (mean(unstdPost)*i)))
  
}   # end of my massive outer for loop


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
seqiplot(wcst.seq,with.legend="right")

# do the entropy histogram for the paper
hist(seqient(wcst.seq.test),xlab="entropy",ylim=c(0, 80),main="")

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
df_ent <- df2[1:36,]
df_kl <- df2[37:72,]
df_sur <- df2[73:108,]
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
xtable(head(wcst.seq.train[c(20,19,10,36,37,39),]))
seqlength(wcst.seq.train[70:80,])
seqlength(wcst.seq.train[c(20,19,10,36,37,39),])

### interestingness - via differences, distances and support
# So what are the differences if any between those new sequences flagged as 
# interesting and the others?


# 1. SEQUENCE ENTROPY ----------------------------------------------------
wcst.ient <- seqient(wcst.seq)
summary(wcst.ient)

# 2. SEQUENCE COMPLEXITY (TURBULENCE) -------------------------------------
range(seqST(wcst.seq))
summary(seqST(wcst.seq))

# 3. LONGEST COMMON PREFIX (LCP) -------------------------------------------
lcp <- seqdist(wcst.seq.test, method = "LCP",norm=FALSE)
range(lcp)
median(lcp)

# 4. LONGEST COMMON SUBSEQUENCE (LCS) --------------------------------------
## Normalized LCS distances to the most frequent sequence
wcst.dref1 <- seqdist(wcst.seq.test, method = "LCS",refseq = 0, norm = "auto")
# LCS distances between two subsets of sequences
set1 <- 1:18
set2 <- 19:36
wcst.dref2 <- seqdist(wcst.seq.test, method = "LCS",refseq = list(set1,set2))
median(wcst.dref2)


# 5. OPTIMAL MATCHING DISTANCE (OMD) and clustering -------------------------
ccost <- seqsubm(wcst.seq, method = "CONSTANT", cval = 2)
# compute the distances using the matrix and the default indel cost of 1
wcst.OM <- seqdist(wcst.seq, method = "OM", indel=3, sm = ccost)
wcst.clust <- agnes(wcst.OM, diss = TRUE, method = "ward")
plot(wcst.clust, which.plots = 2)

cluster2 <- cutree(wcst.clust, k = 2)
#cluster2 <- factor(wcst.clust, labels = c("Type 1", "Type 2"))

table(cluster2)
seqfplot(wcst.seq, group = cluster2, pbarw = T)
round(wcst.OM[1:5, 1:5], 1)

# sequence events
seqmtplot(wcst.seq, group = cluster2, ylim=c(0,4))

# creates a 16 x 16 cost matrix (we have 16 symbols)
couts <- seqsubm(wcst.seq, method = "TRATE")
round(couts, 2)
range(couts)


# 6. EVENT SEQUENCE DISCRIMINATION --------------------------------------------------------
#    create EVENT sequences and discriminate between surprising and the not-surprising sequences
wcst.seqeperiod <- seqecreate(wcst.seq.test, tevent = "period")
wcst.seqe <- seqecreate(wcst.seq.test)
wcst.seqestate <- seqecreate(wcst.seq.test, tevent = "state")

xss <- 5
fsubseq <- seqefsub(wcst.seqe, min.support = xss)
fsubseq[1:xss]

#fsubseq <- seqefsub(sep.seqe, pmin.support = 0.01)

wcstNEW <- data.frame(surp,decay,kl)  #
cohort <- factor(wcstNEW$surp > mean(surp), labels = c("surprising","not-surprising")) # divide into two groups
#cohort <- factor(sepsis$birthyr > 1945, labels = c("<=1945",">1945")) # example from traminer

discrimcohort <- seqecmpgroup(fsubseq, group = cohort, method = "bonferroni")
discrimcohort[1:10]

seqdf <- cbind(as.character(discrimcohort$subseq), discrimcohort$data)
seqdf <- seqdf[-c(5,8,9)]  # drop the index, and two residuals
xtable(seqdf[1:10,])

#plot(discrimcohort[1:10,],ylim=c(0,2.0),with.legend=FALSE)

# modify plots
plot(discrimcohort[1:5,],with.legend=FALSE,horiz=FALSE)



#################### now plot sequences in traminer style ################
seqiplot(wcst.seq,with.legend="right")






