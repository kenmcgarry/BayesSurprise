# bayesian_surprise.R
# 20/01/2021
# https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.4/topics/KLD

# The difference between past, current and future patterns or events:
# How surprised should we be? If we are surprised by the new pattern is it valid, is it a novel pattern?
# Critera:distance, complexity, validity, statistical strength, unusual variables or variable associations 

setwd("C:/common_laptop/R-files/BayesianSurprise")

library(LaplacesDemon)
library(LearnBayes)
library(bnlearn)
library(gRain)
library(lattice)
library(gridExtra)

isachs <- read.table("C:/common_laptop/R-files/hypothesis/sachs.interventional.txt",header = TRUE, colClasses = "factor")

val.str1 <- paste("[PKC][PKA|PKC][praf|PKC:PKA]",
                  "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
                  "[pakts473|p44.42:PKA][P38|PKC:PKA]",
                  "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
                  "[PIP2|plcg:PIP3]",sep="")

val.str2 <- paste("[PKC][PKA|PKC][praf|PKC:PKA]",
                  "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
                  "[pakts473|p44.42:PKA][P38|PKC:PKA]",
                  "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
                  "[PIP2|plcg:PIP3]",sep="")

val.str3 <- paste("[PKC][PKA|PKC][praf|PKC:PKA]",
                  "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
                  "[pakts473|p44.42:PKA][P38|PKC:PKA]",
                  "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
                  "[PIP2|plcg:PIP3]",sep="")

val <- model2network(val.str1)

isachs = isachs[, 1:11]   # preprocess Sachs data from numeric to strings 
for (i in names(isachs))
  levels(isachs[, i]) = c("LOW", "AVG", "HIGH")
fitted = bn.fit(val, isachs, method = "bayes")
jtree = compile(as.grain(fitted))
jprop = setFinding(jtree, nodes = "p44.42",states = "LOW")

querygrain(jtree, nodes = "pakts473")$pakts473
querygrain(jprop, nodes = "pakts473")$pakts473
querygrain(jtree, nodes = "PKA")$PKA
querygrain(jprop, nodes = "PKA")$PKA
names(which.max(querygrain(jprop,nodes = c("PKA"))$PKA)) 

pats473 <- cpdist(fitted, nodes = "pakts473",evidence = (p44.42 == "LOW"))
withevid <- prop.table(table(pats473))
PKA <- cpdist(fitted, nodes = "PKA",evidence = (p44.42 == "LOW"))
noevid <- prop.table(table(PKA))

cpquery(fitted,
        event = (pakts473 == "LOW") & (PKA != "HIGH"),
        evidence = (p44.42 == "LOW") | (praf == "LOW"))

# barchart: with and without evidence
counts <- rbind(as.vector(noevid),as.vector(withevid))
colnames(counts) <- as.list(names(noevid))

barplot(counts, main="Gene activity", horiz = TRUE, ylab="Categories",
        xlab="probablity", col=c("white","grey"),xlim=c(0,1),
        legend = c("without surprise","with surprise"), beside=TRUE)

break

#--------------    The bn.fit.barchart() function is a bit limited  ----------------------------
bn.fit.barchart(fitted$PKA)
bn.fit.barchart(fitted$Raf)
bn.fit.barchart(fitted$pakts473)

#----------------   Marco Scutari version (personal communication) --------------------------------------#
isachs = isachs[, 1:11]
for (i in names(isachs))
  levels(isachs[, i]) = c("LOW", "AVG", "HIGH")

fitted = bn.fit(val, isachs, method = "bayes")
jtree = compile(as.grain(fitted))
jprop = setFinding(jtree, nodes = "p44.42", states = c("LOW"))


a = noevid#querygrain(jtree, nodes = "pakts473")$pakts473
b = withevid#querygrain(jprop, nodes = "pakts473")$pakts473
d = data.frame(A = a, B = b, C = factor(c("LOW", "AVG", "HIGH"),levels = c("LOW", "AVG", "HIGH")))


lattice.options(default.theme = canonical.theme(color = FALSE)) 

p1 = barchart(C ~ A + B, data = d, horizontal = TRUE, xlim=c(0,1),
              ylab = "pakts473", xlab = "probability", main = expression(P(pakts473)),
              auto.key = list(corner = c(0.9, 0.9), reverse.rows = TRUE,
              text = c("without evidence", "with evidence")))

a = querygrain(jtree, nodes = "PKA")$PKA
b = querygrain(jprop, nodes = "PKA")$PKA

d = data.frame(A = a, B = b, C = factor(c("LOW", "AVG", "HIGH"),levels = c("LOW", "AVG", "HIGH")))

p2 = barchart(C ~ A + B, data = d, horizontal = TRUE, xlim=c(0,1),
              ylab = "PKA", xlab = "probability", main = expression(P(PKA)),
              auto.key = list(corner=c(0.9,0.9), reverse.rows = TRUE,
                         text = c("without evidence", "with evidence")))

#postscript("cpquery.eps", width = 10, height = 4,horizontal = FALSE, paper = "special")
pdf(file="cpquery.pdf",width = 10, height = 4)

grid.arrange(p1, p2, ncol = 2)
dev.off()



#-----------------------------------------------------------------------------#

