# bayesian_surprise_main.R
# Surprising patterns could be novel and interesting, worthy of investigation. Surprise is the difference 
# between our expectations and actual experiences. How surprised should we be based on the size of differences? 
# ...and if we are surprised by this new pattern is it valid? 

# 20/01/2021; 28/10/2021

library(dplyr)
library(ggplot2)
library(reshape2)
library(infotheo)
library(philentropy)
library(stringr)
library(stringi)
library(stringdist)
library(TraMineR)
library(PST)
library(SeqDetect)
library(tidyverse)
library(tidyr)
library(hrbrthemes)
library(xtable)

setwd("C:/common_laptop/R-files/BayesianSurprise")
rm(list = ls()) # remove any legacy variables and data from previous analysis

source("bayesian_surprise_functions.R")
source("bayesian_surprise_strfunctions.R")
#source("bayesian_surprise_sequences.R")
#source("bayesian_surprise_simple_example.R")
#source("bayesian_surprise_experiments.R")
#source("bayesian_surprise_experiment-1.R")
#source("bayesian_surprise_experiment-2.R")
#source("bayesian_surprise_experiment-3.R") # SRH
source("bayesian_surprise_experiment-4.R") # Biofam
#source("bayesian_surprise_experiment-5.R")  # Sepsis




