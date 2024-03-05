# bayesian_surprise_new.R
# changes based on reviewers comments for COGSYS-D-22-00205

setwd("C:/common_laptop/R-files/BayesianSurprise")

library(PST)
library(TraMineR)

# two new datasets
data(mvad)
data(actcal)

actcal.seq <- seqdef(actcal, 13:24, labels = c("FullTime", "PartTime", "LowPartTime", "NoWork"))
transition <- seqetm(actcal.seq, method = "transition")
transition

#####################################################
data("SRH", package = "PST")
set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(SRH.seq), size = floor(.75*nrow(SRH.seq)), replace = FALSE)
SRH.seq.train <- SRH.seq[seqsample, ]
SRH.seq.test  <- SRH.seq[-seqsample, ]

data(s1)
s2.seq <- seqdef(s1)
srh.pst <- pstree(s2.seq, L=3)
plot(S2)
plot(S2, horiz=TRUE,withlegend=FALSE)
plot(S2, nodePar=list(node.type="path", lab.type="prob", lab.pos=1, lab.offset=2, lab.cex=0.7), 
     edgePar=list(type="triangle"), withlegend=FALSE)


# New PST diagram for paper
# 
s1 <- "h-e-l-l-o-h-e-l-l-o"
s1.seq <- seqdef(s1)
PST1.pst <- pstree(s1.seq, L=3)#nmin=1, ymin=0.001)
#cmine(PST1.pst, pmin = 0.05, state = "H")
plot(PST1.pst, horiz=TRUE,withlegend=FALSE,axis=FALSE)

# new PST plot for paper
s3 <- "c-c-a-a-a-a-b-b-b-b-a-a-a-a-c-c"
s3.seq <- seqdef(s3)
S3 <- pstree(s3.seq, L=3)
plot(S3,horiz=FALSE,withlegend=FALSE,axis=FALSE)

################ experiment 1, do a facet plot###################
txt <- "    Seq     H   E   T   C
1   Seq_1   2   1   5   4
2   Seq_2   2   1   5   4
3   Seq_3   2   1   5   4
4   Seq_4   0   0   6   6
5   Seq_5   0   4   2   6 "

dat <- read.table(textConnection(txt), header = TRUE)
dat.m <- melt(dat)

ggplot(dat.m, aes(variable, value, group = Seq)) + 
  geom_line() + 
  facet_wrap(~Seq)

################ experiment 2, possibly ###################
# https://cran.r-project.org/web/packages/ggseqplot/vignettes/ggseqplot.html
data(biofam)
biofam <- biofam[1:200,]
biofam.seq <- seqdef(biofam[,10:25], xtlab=as.character(15:30), xtstep=3)
## Plotting transversal entropies by sex
seqplot.tentrop(biofam.seq, group=biofam$sex, legend.pos="bottomright")
slist <- list(woman = biofam.seq[biofam$sex=="woman",],
              man = biofam.seq[biofam$sex=="man",])
seqplot.tentrop.m(slist, legend.pos="bottomright")
## Plotting transversal entropies for women
## by father's social status
group <- biofam$cspfaj[biofam$sex=="woman"]
seqplot.tentrop(biofam.seq[biofam$sex=="woman",], group=group,
                main="Women, by father's social status", legend.pos="bottomright")

### also a possibility #####
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
                "Child", "Left+Child", "Left+Marr+Child", "Divorced")
## Here, we use only 100 cases selected such that all elements
## of the alphabet be present.
## (More cases and a larger k would be necessary to get a meaningful example.)
biofam.seq <- seqdef(biofam[501:600, ], 10:25, labels=biofam.lab)
diss <- seqdist(biofam.seq, method="LCS")
## Using 12 groups and default MDS sorting
seqplot.rf(biofam.seq, diss=diss, k=12,
           main="Non meaningful example (n=100)")
## With a user specified sorting variable
## Here time spent in parental home: there are ties
## We set a seed because of random order in ties
set.seed(123)
parentTime <- seqistatd(biofam.seq)[, 1]
seqplot.rf(biofam.seq, diss=diss, k=12, sortv=parentTime,main="Sorted by parent time")

###############################################################################
#### another possibility ####
data(ex1)
#s.ex1 <- seqdef(ex1[,1:13],weights=ex1[,"weights"])
s.ex1 <- seqdef(chess.seq[1:13,],weights=ex1[,"weights"])

seqlength(s.ex1)
seqlength(s.ex1, with.missing=FALSE)
group <- c(1,1,1,2,2,2,2)
ind.d <- seqindic.dyn(s.ex1, indic='cplx', with.missing=FALSE)
plot(ind.d, group=group, fstat=weighted.mean, na.rm=TRUE, conf=TRUE, ret=TRUE)
## Treating 'missing' as a regular state
ind.dm <- seqindic.dyn(s.ex1, indic='cplx', with.missing=TRUE)
plot(ind.dm, group=group, fstat=weighted.mean, na.rm=TRUE, conf=TRUE, ret=TRUE)


######################################    #########################################
data(mvad)
## Building a state sequence object
mvad.seq <- seqdef(mvad, 17:86)
## Sequence of typical states
mvad.si.gcse5eq <- seqimplic(mvad.seq, group=mvad$gcse5eq)
##Plotting the typical states
plot(mvad.si.gcse5eq, lwd=3, conf.level=c(0.95, 0.99))
## Printing the results
print(mvad.si.gcse5eq, xtstep=12)
##### motif version #####
msseq <- seqimplic(motifs.seq)#, group=mvad$gcse5eq)
##Plotting the typical states
plot(msseq, lwd=3, conf.level=c(0.95, 0.99))
## Printing the results
print(msseq, xtstep=12)

#################################### combined plots ############################
## compute dissimilarity matrix required for plot
diss <- seqdist(chess.seq, method = "LCS")
ggseqrfplot(chess.seq, diss = diss, k = 10) 

## adjusted version
ggseqrfplot(chess.seq, diss = diss, k = 10) &
  theme_ipsum(base_family = "") &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 10)) &
  plot_annotation(title = "Bayesian decay Plot")


################################## REMEDIAL PLOTS ##############################
motifs.seq <- TraMineR::seqdef(c("M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
  "M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M2-M3-M4-M5-M1-M2-M6",
  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M4-M2-M2-M3-M4-M5-M1-M6",
  "M4-M5-M1-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M1-M6",
  "M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
  "M1-M2-M3-M4-M5-M1-M1-M2-M3-M4-M5-M1-M2-M3",
  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6"))

seqsample <- sample.int(n = nrow(motifs.seq),size=floor(.8*nrow(motifs.seq)), replace = FALSE)
motifs.seq.train <- motifs.seq[seqsample, ]
motifs.seq.test  <- motifs.seq[-seqsample, ]

motifs.pst <- pstree(motifs.seq,L=5,ymin = 0.001)

z<-seqdef("M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6")
pqplot(motifs.pst,z, measure="logloss", plotseq=TRUE, seqscale=TRUE)

z<-seqdef("d4-e6-Nc3-white")
PST::pqplot(chess.pst, z, measure="logloss", plotseq=TRUE, seqscale=TRUE)

PST::pqplot(wcst.pst,  wcst.seq.test[10,], measure="logloss", plotseq=TRUE, seqscale=TRUE,ptype="b")

PST::pqplot(biofam.pst, biofam.seq.test[1,], measure="logloss", plotseq=TRUE, seqscale=TRUE)

z <- seqdef("ERRegistration-ERTriage-ERSepsisTriage-CRP-LacticAcid-Leucocytes-AdmissionIC-AdmissionNC-CRP-Leucocytes-CRP-AdmissionNC-ReleaseB")
PST::pqplot(sepsis.pst, z, measure="logloss", plotseq=TRUE, seqscale=TRUE)

z <- seqdef("B2-M-G2-B2-M-M-G2-B2-M-M-M")
PST::pqplot(srh.pst, z, measure="logloss", plotseq=TRUE, seqscale=TRUE)

##############################
data(mvad)
mvad.labels <- c("employment", "further education", "higher education","joblessness", "school", "training")
mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode,labels = mvad.labels)
seqiplot(mvad.seq, withlegend = F, title = "Index plot (10 first sequences)")
seqfplot(mvad.seq, pbarw = T, withlegend = F, title = "Sequence frequency plot")
seqdplot(mvad.seq, withlegend = F, title = "State distribution plot")
seqlegend(mvad.seq, fontsize = 1.3)  # just the legend

Entropy <- seqstatd(mvad.seq)$Entropy
plot(Entropy, main = "Decay of Surprise Value",col = "blue", 
     xlab = "Records", ylab = "Decay",type = "l")

### WCST
Entropy <- seqstatd(wcst.seq)$Entropy
plot(Entropy, main = "Decay of Surprise Value",col = "blue", 
     xlab = "Records", ylab = "Decay",type = "Ht")

##################### LOG-LOSS for entire data sets #############
## BIOFAM
logloss <- predict(biofam.pst, biofam.seq.test, output = "logloss")
avelogloss <- mean(logloss)
hist(logloss,col = "cyan", main="",freq=FALSE,ylim=c(0,5))  #ylim=c(0,200)
#abline(h=avelogloss,col="red",lty = 2)

## SEPSIS
logloss <- predict(sepsis.pst, sepsis.seq.test, output = "logloss")
avelogloss <- mean(logloss)
hist(logloss,col = "cyan", main="",freq=FALSE,ylim=c(0,1))  #ylim=c(0,200)
#abline(h=avelogloss,col="red",lty = 2)

## CHESS
logloss <- predict(chess.pst, chess.seq.test, output = "logloss")
avelogloss <- mean(logloss)
hist(logloss,col = "cyan", main="",freq=FALSE,ylim=c(0,.5))  #ylim=c(0,200)
#abline(h=avelogloss,col="red",lty = 2)

## WCST
logloss <- predict(wcst.pst, wcst.seq.test, output = "logloss")
logloss[is.infinite(logloss)] <- 0
avelogloss <- mean(logloss)
hist(logloss,col = "cyan", main="",freq=FALSE,ylim=c(0,1.5))  #ylim=c(0,200)
#abline(h=avelogloss,col="red",lty = 2)



#################################################################################

### LOGLOSS WITH CUTOFF POINT OF 80 ########
Turbulence <- seqST(mvad.seq)
summary(Turbulence)
Records <- Turbulence
hist(Records, col = "cyan", main = "Log-loss")
abline(h=80,col="red")

### LOGLOSS WITH CUTOFF POINT OF 15 ########
Turbulence <- seqST(chess.seq)
summary(Turbulence)
Records <- Turbulence
hist(Records, col = "cyan", main = "Log-loss")
abline(h=15,col="red")

### LOGLOSS WITH CUTOFF POINT OF 80 ########
Turbulence <- seqST(wcst.seq.test)
summary(Turbulence)
Records <- Turbulence
hist(Records, col = "cyan", main = "Log-loss")
abline(h=20,col="red")

### LOGLOSS WITH CUTOFF POINT OF 180 ########
Turbulence <- seqST(biofam.seq)
summary(Turbulence)
Records <- Turbulence
hist(Records, col = "cyan", main = "Log-loss")
abline(h=188,col="red")

### LOGLOSS WITH CUTOFF POINT OF 150 ########
Turbulence <- seqST(sepsis.seq)
summary(Turbulence)
Records <- Turbulence
hist(Records, col = "cyan", main = "Log-loss")
abline(h=150,col="red")



mvad.seqe <- seqecreate(mvad.seq)
fsubseq <- seqefsub(mvad.seqe, pmin.support = 0.05)
plot(fsubseq[1:15], col = "cyan")
discr <- seqecmpgroup(fsubseq, group = cl1.3fac)
plot(discr[1:6])


submat <- seqsubm(mvad.seq, method = "TRATE")
dist.om1 <- seqdist(mvad.seq, method = "OM", indel = 1,sm = submat)



cm2 <- cmine(motifs.pst, pmin = 0.5, state = c("M1", "M2"))
plot(cmine(motifs.pst, pmin = 0.5, state = c("M1", "M2")))
cprob(motifs.seq, L = 1, prob = TRUE, with.missing = TRUE)

