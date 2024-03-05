# bayesian_surprise_main.R
# Surprising patterns could be novel and interesting, worthy of investigation. Surprise is the difference 
# between our expectations and actual experiences. How surprised should we be based on the size of differences? 
# ...and if we are surprised by this new pattern is it valid? 

# commenced: 20/01/2021;
# submitted: 07/01/2022;
# updated:   29/02/2024

library(dplyr)
library(ggplot2)
library(reshape2)
library(infotheo)
library(philentropy)
library(stringr)
library(stringi)
library(stringdist)
library(TraMineR)
library(TraMineRextras)
library(ggseqplot)
library(PST)
library(SeqDetect)
library(tidyverse)
library(tidyr)
library(hrbrthemes)
library(xtable)
library(cluster)
library(patchwork)

# my methods might be of use to:
# https://www.frontiersin.org/articles/10.3389/fnhum.2013.00492/full

setwd("C:/common_laptop/R-files/BayesianSurprise")
rm(list = ls()) # remove any legacy variables and data from previous analysis

#save.image(file='March3rd2024.RData')
#load('2ndMarch2024.RData')

source("bayesian_surprise_functions.R")
source("bayesian_surprise_strfunctions.R")

#source("bayesian_surprise_WCST.R") # Wisconsin Card sorting test
#source("bayesian_surprise_CHESS.R")   # CHESS moves dataset
#source("bayesian_surprise_BIOFAM.R") # Biofam
#source("bayesian_surprise_SEPSIS.R")  # Sepsis


