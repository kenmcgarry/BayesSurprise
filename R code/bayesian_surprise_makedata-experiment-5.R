# bayesian_surprise_experiments-1.R
# Surprising patterns could be novel and interesting, worthy of investigation.
# 23/11/2021; 

# Note: Bayesian surprise is the KL divergence from prior to posterior i.e. between the distribution of model 
# hypothesis before and after observing the data. We should not use entropy or other widely used information 
# bearing calculations because data such as outliers or random data, though containing a lot of information 
# because they are highly improbable, are not necessarily surprising in the human context and in fact could
# be very boring!

##### EXPERIMENT 5 - SEPSIS SEQUENCE ANOMALIES ------------------

data(sepsis)  # from seqDetect package

# Build the array of lists
cases <- distinct(sepsis, case_id)
casearray <- array(list(), nrow(cases))
patientarray <- array(list(), nrow(cases))

activity <- sepsis %>%  # just keep case_id and activity
  select(case_id,activity)

activity <- data.frame(lapply(activity, as.character), stringsAsFactors=FALSE) # convert to strings

for(i in 1:nrow(cases)){
  casearray[[i]] <- activity %>%
    filter(case_id == cases[i,1])
 
}

data1 <- casearray[[1]]  # convert and keep first record
data1 <- na.omit(data1) #some NA's have crept in, removing them reduces from 15214 to 15190 records

cases <- distinct(data1, case_id)

for(i in 1:nrow(cases)){
  casearray[[i]] <- activity %>%
    filter(case_id == cases[i,1])
}

# Now format into a strings of activities for each patient
for(i in 1:nrow(cases)){
  tempstr <- unlist(casearray[[i]]$activity)
  tempstr <- str_replace_all(tempstr, "[[:punct:]]", " ")
  tempstr <- str_replace_all(tempstr, fixed(" "), "") # remove whitespace
  tempstr <- paste(tempstr,collapse="-")
  tempstr <- str_replace_all(tempstr, "[\r\n\t\v\f]" , "")
  
  patientarray[[i]]$activity <- tempstr
  #cat("\n",i," tempstr",tempstr)
  if(length(casearray[[i]]$activity) > 30){cat(i," ")}  # 118, 347, 870 are huge!
  #if(length(casearray[[i]]$activity) > 30){patientarray[[i]] <- patientarray[[-(i)]]}
}

# These records are causing problesm for the PST::predict() function because of their size, delete them
patientarray <- patientarray[-c(98,108,111,118,124,162,193,199,257,271,291,294,303,345,347,348,355,369,375,408,497,500,502,524,560,562,647,651,677,678,688,715,743,776,801,818,870,879,910,920,924,930,933,973,985,1011,1037)]

# Now convert into seqdef() format  # M6-M2-M3-M4-M6-M1-M2-M3-M4"
rownames(patientarray) <- NULL
crappy <- cbind(unlist(patientarray))
sepsis.seq <- TraMineR::seqdef(crappy)
rownames(sepsis.seq) <- NULL

set.seed(101) # Set Seed so that same sample can be reproduced in future
# Now randomly select 75% of data as sample from total 'n' rows of the data  
seqsample <- sample.int(n = nrow(sepsis.seq), size = floor(.75*nrow(sepsis.seq)), replace = FALSE)
sepsis.seq.train <- sepsis.seq[seqsample, ]
sepsis.seq.test  <- sepsis.seq[-seqsample, ]

# Build the PST on training data
model5.pst <- pstree(sepsis.seq.train, L = 10, nmin = 2, ymin = 0.001)


# Test the PST
probtest <- PST::predict(model5.pst, sepsis.seq.test , decomp=TRUE,output = "prob")
probtest  

probtrain <- cprob(sepsis.seq.train , L = 0, prob = TRUE)  # overall probs for training dataset
probtrain




