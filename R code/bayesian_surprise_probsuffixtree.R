# bayesian_surprise_probsuffixtree.R
# 16/07/2021
# Provides a framework for analysing state sequences with probabilistic suffix trees (PST), 
# the construction that stores variable length Markov chains (VLMC).

setwd("C:/common_laptop/R-files/BayesianSurprise")
# simulate data sequence, now we are getting somewhere
#s1 <- "ac-abc-abc-aa-aa-ac-ab-aa-aa-aa-bbb-bbb-abc-bbb"
#s2 <- "abc-abc-ab-c-cc-abb-bba-b-a-b-bcca-bbb-abc-bbb-bbb"
#s3 <- "ab-ab-ba-ba"

testdata1.seq <- TraMineR::seqdef(c("M2-M3-M5-M5"))
testdata2.seq <- TraMineR::seqdef(c("M2-M5-M3-M5"))
testdata3.seq <- TraMineR::seqdef(c("M1-M6-M3-M4"))
testdata4.seq <- TraMineR::seqdef(c("M2-M3-M5-M5-M2-M5-M3-M5"))
testdata5.seq <- TraMineR::seqdef(c("M6-M3"))

# experiment for detecting differences in repeating pattern M1...M6
motiftrain <- c("M1-M6-M3-M4-M5-M1-M2-M3-M4-M5-M5-M2-M6-M4",
                "M5-M1-M2-M3-M4-M5-M5-M2-M3-M4-M5-M1-M2",
                "M2-M3-M5-M5",
                "M2-M3-M5-M5",
                "M2-M3-M5-M5",
                "M2-M5-M2-M5",
                "M1-M6-M3-M4-M5-M1-M2",
                "M2-M3-M5-M5-M1-M2-M3-M4-M5-M1-M2")
seq.motif.train <- TraMineR::seqdef(motiftrain)
model.pst <- pstree(seq.motif.train,nmin=1,with.missing = TRUE)
summary(model.pst)
probs <- PST::predict(model.pst, testdata5.seq, decomp=TRUE,output = "prob")
probs

stop("Ken stops execution at this point giving an error message!")

# -------------------OK UNTIL THIS POINT --------------------------------------------------------
# Generate data - As a probabilistic suffix tree (PST) represents a generating model, it can 
# be used to generate artificial sequence data sets. Sequences are built by generating the 
# states at each successive position. The process is similar to sequence prediction (see predict), 
# except that the retrieved conditional probability distributions provided by the PST are used 
# to generate a symbol instead of computing the probability of an existing state. 
gendat.seq <- generate(model.pst, n=10, l=10, method="prob")


sim2.pst <- pstree(gendat.seq,nmin=1,with.missing = TRUE)
summary(sim2.pst)

gendat2.seq <- generate(sim2.pst, n=20, l=26, method="prob")
probs2 <- PST::predict(sim2.pst, gendat2.seq, decomp=FALSE,output = "prob")
probs2


sbig.seq <- seqdef(sbig)
sbig.seq.event <- seqecreate(sbig.seq)
s1.seq.event <- seqecreate(s1.seq)
s2.seq.event <- seqecreate(s2.seq)
s3.seq.event <- seqecreate(s3.seq)

seqpm(sbig.seq,"abc",sep="-")   # Find substring patterns in sequence
seqistatd(sbig.seq)    # State frequencies in each individual sequence
seqefsub(sbig.seq.event, min.support= 1)   # Search for frequent subsequences

# calculate relative frequencies and hence probabilities
seqlength <- seqelength(s2.seq.event)
statseq <- seqstatf(s2.seq)
statseq[,3] <- statseq[,1]/seqlength; 
colnames(statseq)[3] <- "RelativeFreqs"
print(statseq,digits=3)

# get the probabilities of a sequence, same results as relative frequencies
cprob(s1.seq, L = 0, prob = TRUE)
cprob(s2.seq, L = 0, prob = TRUE)
cprob(s3.seq, L = 0, prob = TRUE)

# query will retrieve counts or next symbol probability distribution
query(sim1.pst, "M7")

# cmine looks for all contexts yielding a probability of at least pmin = 0.5 for the state "ab"
cmine(sim1.pst, pmin = 0.5, state = "M1")

# generate plausible data from the PST model
gen1.seq <- generate(sim1.pst, l = 13, L = 10, n = 10)

## model selection
alpha <- c(1.0, 0.5, 0.1, 0.02, 0.01, 0.001)
C.list <- qchisq(1 - alpha, df = 6 - 1) / 2
C.tuned <- tune(sim1.pst, gain = "G2", C = C.list, output = "stats", criterion = "AIC")

