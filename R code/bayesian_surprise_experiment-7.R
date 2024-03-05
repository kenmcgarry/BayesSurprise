# bayesian_surprise_experiment-7.R
## Minecraft data by Samaneh Saadat
## https://github.com/SamanehSaadat/ContrastMotifDiscovery/blob/master/data/labeled_seqs.csv
## Samaneh Saadat and Gita Sukthankar, “Contrast Motif Discovery in Minecraft”, 
## Proceedings of the AAAI Conference on Artificial Intelligence 
## and Interactive Digital Entertainment (AIIDE), Oct 2020

setwd("C:/common_laptop/R-files/BayesianSurprise")

library(PST)
library(TraMineR)

minecraft <- read.csv("minecraft_seqs.csv")
#minecraft <- as.character(minecraft$seq)  # data needs to be a series strings
#minecraft <- as.factor(minecraft$seq)  # data needs to be a series strings

mc.seq.sample <- seqsep(as.character(minecraft$seq))   # stick the "-" symbol between all letters
mc.seq.sample <- TraMineR::seqdef(mc.seq.sample) # convert it into a seq object

seqsample <- sample.int(n = nrow(mc.seq.sample), size = floor(.5*nrow(mc.seq.sample)), replace = FALSE)
mc.seq.train <- mc.seq.sample[seqsample, ]
mc.seq.test  <- mc.seq.sample[-seqsample, ]


model7.pst <- pstree(mc.seq.train,  L=5, nmin = 1, ymin=0.001)
probtest <- PST::predict(model7.pst, mc.seq.test[1:5,], decomp=FALSE,output = "prob")
head(probtest)  







