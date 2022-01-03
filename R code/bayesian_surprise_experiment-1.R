# bayesian_surprise_experiments-1.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 22/09/2021; 

# Note: Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact could
# be very boring!

# This file reads in the sequence data for a given experiment, and will perform motif generation if the
# sequence is particularly long and/or has many letters in its alphabet. Some data sets will not require this
# process if they have a limited alphabet or are short sequences - we then go straight to building probabilistic 
# suffix trees and searching for interesting patterns via Bayesian surprise. 

## EXPERIMENT 1 - simulated data, find largest repeating sequence that occurs at least twice (motifs)
seq1 <- "kvuasvclhihijhijhhhihhjhikjhhjjjihjhijklhisvclfvbsvclfvbsvclfvbsvclfvbjjjkhihhijhijhijhijjjsvclfvbsvclfvbkhihkhihkvuasvclkvuasvclkvuasvcl"
motifs <- find_motifs(seq1)

# converts the sequence motif into a series of identifiers prefixed by "M" e.g. "M7" etc
# motifs[[2]] contains the list of motif names and order in which they appear
motiforder <- motif_sequence(seq1,motifs[[2]])  

# join the motifs in one string but separate them by "-" as required by seqdef() function
# e.g. "M2-M3-M5-M5-M6-M6-M6-M7-M5-M1"
motiforder <- paste(motiforder, collapse="-")  # convert into seqdef friendly format
seq1.motif.train <- TraMineR::seqdef(motiforder)
seq1.motif.test <-  TraMineR::seqdef("M6-M7-M5-M1-M1-M7-M4-M6")  # simulated test sequence

# build the PST
model1.pst <- pstree(seq1.motif.train,nmin=1,with.missing = TRUE)
summary(model1.pst)

# pass new data through PST and determine similarity or not
probs <- PST::predict(model1.pst, seq1.motif.test, decomp=TRUE,output = "prob")
probs

# Example plot of a PST for paper
s1.seq <- seqdef("A-B-C-B-B-B-C-C-A-A-C-A")
S1 <- pstree(s1.seq, L=3, nmin=1)  # maximal depth of 3, minimal samples is 1.
plot(S1, nodePar = list(node.type = "path", 
     lab.type = "prob",lab.pos = 1, lab.offset = 3, lab.cex = 0.9),
     edgePar = list(type = "triangle"), withlegend = "FALSE")

#------------------------------------------------------------------------------------------
## just testing random sequence generator
randdata <- randomClickstreams(states = c("P1","P2","P3","P4"),
    startProbabilities = c(0.3, 0.2, 0.2, 0.3),
    transitionMatrix = matrix(c(0.2, 0.2, 0.4, 0.2, 
                                0.1, 0.4, 0.4, 0.1,  
                                0.2, 0.1, 0.2, 0.5, 
                                0.4, 0.2, 0.2, 0.2), nrow = 4),
                                meanLength = 15, n = 20)
seq.list <- list()
for(i in 1:length(randdata)){
  seq.list[i] <- stri_join_list(randdata[i], sep = "-", collapse = NULL)
}

seq.list <- unlist(seq.list)
seq.rand.train <- TraMineR::seqdef(seq.list)
                 
modelrand.pst <- pstree(seq.rand.train, L = 10, nmin = 1, ymin = 0.001)

#----------------------------------------------------------------------------------------------------
## dengue virus DNA sequence, this virus is composed of 10,735 DNA symbols, coding for 10 proteins
load("dengvirus.RData")

seq <- toupper(seq)
denseq <- paste(seq,collapse ="")
#head(denseq)

# for DNA important motifs are : start codon:"ATG"; stop codons: are "TGA", "TAA", and "TAG"
motifs <- find_motifs(substr(denseq,1,1000))

# converts the sequence motif into a series of identifiers prefixed by "M" e.g. "M7" etc
# motifs[[2]] contains the list of motif names and order in which they appear
motiforder <- motif_sequence(denseq,motifs[[2]])  

# join the motifs in one string but separate them by "-" as required by seqdef() function
# e.g. "M2-M3-M5-M5-M6-M6-M6-M7-M5-M1"
motiforder <- paste(motiforder, collapse="-")  # convert into seqdef friendly format

seq.deng.train <- TraMineR::seqdef(motiforder)

modeldeng.pst <- pstree(seq.deng.train, L = 10, nmin = 1, ymin = 0.001)

probtest <- PST::predict(modeldeng.pst, seq.deng.train[1], decomp=TRUE,output = "prob") # 




