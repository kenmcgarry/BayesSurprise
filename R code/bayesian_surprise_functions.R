# bayesian_surprise_functions.R
# 08/07/2021
# Some functions written by myself, others taken from packages and stackoverflow contributors.

## discrete.bayes() from Jim Albert, Computes the posterior distribution for an arbitrary one parameter 
# distribution for a discrete prior distribution.
discrete.bayes<- function(df, prior, y, ...) 
{
  param <- as.numeric(names(prior))
  lk <- function(j) prod(df(y, param[j], ...))
  likelihood <- sapply(1:length(param), lk)
  pred <- sum(prior * likelihood)
  prob <- prior * likelihood/pred
  obj <- list(prob = prob, pred = pred)
  class(obj) <- "bayes"
  obj
}


## Kullback() from LaplacesDemon ------------------------------------------------------------------
kullback <- function (px, py, base = exp(1)) 
{
  if (!is.vector(px)) 
    px <- as.vector(px)
  if (!is.vector(py)) 
    py <- as.vector(py)
  n1 <- length(px)
  n2 <- length(py)
  if (!identical(n1, n2)) 
    stop("px and py must have the same length.")
  if (any(!is.finite(px)) || any(!is.finite(py))) 
    stop("px and py must have finite values.")
  if (any(px <= 0)) 
    px <- exp(px)
  if (any(py <= 0)) 
    py <- exp(py)
  px[which(px < .Machine$double.xmin)] <- .Machine$double.xmin
  py[which(py < .Machine$double.xmin)] <- .Machine$double.xmin
  px <- px/sum(px)
  py <- py/sum(py)
  KLD.px.py <- px * (log(px, base = base) - log(py, base = base))
  KLD.py.px <- py * (log(py, base = base) - log(px, base = base))
  sum.KLD.px.py <- sum(KLD.px.py)
  sum.KLD.py.px <- sum(KLD.py.px)
  mean.KLD <- (KLD.px.py + KLD.py.px)/2
  mean.sum.KLD <- (sum.KLD.px.py + sum.KLD.py.px)/2
  out <- list(KLD.px.py = KLD.px.py, KLD.py.px = KLD.py.px, 
              mean.KLD = mean.KLD, sum.KLD.px.py = sum.KLD.px.py, sum.KLD.py.px = sum.KLD.py.px, 
              mean.sum.KLD = mean.sum.KLD, intrinsic.discrepancy = min(sum.KLD.px.py,sum.KLD.py.px))
  #return(out)
  return(min(sum.KLD.px.py,sum.KLD.py.px))
}


# ### discretize2d ::entropy package -----------------
# usage : discretize2d(rbind(freqs1,freqs2))
discretize2d <- function (x1, x2, numBins1, numBins2, r1 = range(x1), r2 = range(x2)) 
{
  b1 = seq(from = r1[1], to = r1[2], length.out = numBins1 + 1)
  b2 = seq(from = r2[1], to = r2[2], length.out = numBins2 + 1)
  y2d = table(cut(x1, breaks = b1, include.lowest = TRUE),  cut(x2, breaks = b2, include.lowest = TRUE))
  return(y2d)
}

# ### entropy.empirical ::entropy package -----------------
# usage : entropy.empirical(rbind(freqs1,freqs2))
entropy.empirical <- function (y, unit = c("log", "log2", "log10")) 
{
  return(entropy.plugin(freqs.empirical(y), unit = unit))
}

# ### freqs.empirical ::entropy package -----------------
# usage : freqs.empirical(rbind(freqs1,freqs2))
freqs.empirical <- function(y) 
{
  return(y/sum(y))
}


# ### entropy.plugin ::entropy package -----------------
# usage : entropy((freqs1),unit="log2")
entropy <- function (y, lambda.freqs, method = c("ML", "MM","Jeffreys", "Laplace", "SG", "minimax", 
                          "CS", "NSB", "shrink"), unit = c("log","log2", "log10"), verbose = TRUE, ...) 
{
  method = match.arg(method)
  if (method == "ML") 
    H = entropy.empirical(y, unit = unit)
  if (method == "MM") 
    H = entropy.MillerMadow(y, unit = unit)
  if (method == "NSB") 
    H = entropy.NSB(y, unit = unit, ...)
  if (method == "CS") 
    H = entropy.ChaoShen(y, unit = unit)
  if (method == "Jeffreys") 
    H = entropy.Dirichlet(y, a = 1/2, unit = unit)
  if (method == "Laplace") 
    H = entropy.Dirichlet(y, a = 1, unit = unit)
  if (method == "SG") 
    H = entropy.Dirichlet(y, a = 1/length(y), unit = unit)
  if (method == "minimax") 
    H = entropy.Dirichlet(y, a = sqrt(sum(y))/length(y), 
                          unit = unit)
  if (method == "shrink") 
    H = entropy.shrink(y, lambda.freqs = lambda.freqs, unit = unit, 
                       verbose = verbose)
  return(H)
}

### entropy.shrink ::entropy package -----------------
# usage : entropy.shrink(rbind(freqs1,freqs2))
entropy.shrink <- function (y, lambda.freqs, unit = c("log", "log2","log10"), verbose = TRUE) 
{
  f = freqs.shrink(y, lambda.freqs = lambda.freqs, verbose = verbose)
  h = entropy.plugin(f, unit = unit)
  attr(h, "lambda.freqs") = attr(f, "lambda.freqs")
  return(h)
}


### entropy.plugin ::entropy package -----------------
# usage : entropy.plugin(rbind(freqs1,freqs2))
entropy.plugin <- function (freqs, unit = c("log", "log2", "log10")) 
{
  unit = match.arg(unit)
  freqs = freqs/sum(freqs)
  H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") 
    H = H/log(2)
  if (unit == "log10") 
    H = H/log(10)
  return(H)
}

####  freqs.dirichlet ::entropy package ----------------
# usage : freqs.Dirichlet(freqs1,freqs2)
freqs.Dirichlet <- function (y, a) 
{
  ya = y + a
  na = sum(ya)
  pa = ya/na
  return(pa)
}

#### entropy dirichlet ::entropy package ----------------
# usage : entropy.Dirichlet(freqs1,freqs2)
entropy.Dirichlet <- function (y, a, unit = c("log", "log2", "log10")) 
{
  return(entropy.plugin(freqs.Dirichlet(y, a), unit = unit))
}

### KL.plugin ::entropy package ----------------
# usage : KL.plugin(freqs1,freqs2)
KL.plugin <- function (freqs1, freqs2, unit = c("log", "log2","log10")) 
{
  unit = match.arg(unit)
  freqs1 = freqs1/sum(freqs1)
  freqs2 = freqs2/sum(freqs2)
  if (any(!(freqs2 > 0))) 
    warning("Vanishing value(s) in argument freqs2!")
  LR = ifelse(freqs1 > 0, log(freqs1/freqs2), 0)
  KL = sum(freqs1 * LR)
  if (unit == "log2") 
    KL = KL/log(2)
  if (unit == "log10") 
    KL = KL/log(10)
  return(KL)
}

### MI mutual information ::entropy package 
### MI = H(X) + H(Y) - H(X,Y)
# usage : mi.plugin(rbind(freqs1,freqs2))
mi.plugin <- function (freqs2d, unit = c("log", "log2", "log10")) 
{
  unit = match.arg(unit)
  freqs2d = as.matrix(freqs2d/sum(freqs2d))
  freqs.x = rowSums(freqs2d)
  freqs.y = colSums(freqs2d)
  freqs.null = freqs.x %o% freqs.y
  MI = KL.plugin(freqs2d, freqs.null, unit = unit)
  return(MI)
}

### from nowhere --------------------
# returns entropy of a vector in bits
shannon.entropy <- function(p)   # defined as H
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  h <- (-sum(log2(p.norm)*p.norm))
  return(h)
}


### from infotheo --------------------
# defined as H(X|Y)
shannon.conditional.entropy <- function (X, Y = NULL) 
{
  if (is.null(Y)) 
    Hres <- shannon.entropy(X)
  else {
    Hyx <- shannon.entropy(data.frame(Y, X))
    Hy <- shannon.entropy(Y)
    Hres <- Hyx - Hy
  }
  return(Hres)
}

## learnBayes
pdisc<- function (p, prior, data) 
{
  s = data[1]
  f = data[2]
  p1 = p + 0.5 * (p == 0) - 0.5 * (p == 1)
  like = s * log(p1) + f * log(1 - p1)
  like = like * (p > 0) * (p < 1) - 999 * ((p == 0) * (s > 
                                                         0) + (p == 1) * (f > 0))
  like = exp(like - max(like))
  product = like * prior
  post = product/sum(product)
  return(post)
}

## learnBayes
discint <- function (dist, prob) 
{
  x = dist[, 1]
  p = dist[, 2]
  n = length(x)
  sp = sort(p, index.return = TRUE)
  ps = sp$x
  i = sp$ix[seq(n, 1, -1)]
  ps = p[i]
  xs = x[i]
  cp = cumsum(ps)
  ii = 1:n
  j = ii[cp >= prob]
  j = j[1]
  eprob = cp[j]
  set = sort(xs[1:j])
  v = list(prob = eprob, set = set)
  return(v)
}

## learnBayes
pdiscp <- function(p, probs, n, s) 
{
  pred = 0 * s
  for (i in 1:length(p)) {
    pred = pred + probs[i] * dbinom(s, n, p[i])
  }
  return(pred)
}

####------------------------------------------------------------------------------------------------


## Adds new valid sequence(s) to an existing database
addseq.database <- function(olddatabase,newsym){
  if(!is.null(newsym)){
    tmp <- cbind(unlist(newsym))
    tmpstr <- str_c(tmp, collapse = "-")
    newdatabase <- paste(olddatabase,tmpstr,sep="-")
    return(newdatabase)
  }
  else{return(olddatabase)}
  
}

## Identifies unique sequences to add to the alphabet
# Any novel sequence/pattern will have NA for probability value when used by predict(). 
# Find these and return them.
novel.sequences <- function(newseq,newprobs){
  
  anyNA <- which(is.na(newprobs))  # get index of NA's
  
  if(!is.null(anyNA)){
    allstrings <- unlist(strsplit(newseq, "-", fixed = TRUE))
    contexts <- allstrings[anyNA]
    cat(bold$green("\nFound",length(anyNA),"novel patterns in new sequence(",paste(contexts,collapse = ','),")\n"))
  }else{contexts <- NULL}
  
  return(contexts)
}





##-------------------------- BAYESIAN SURPRISE FUNCTIONS -----------------------

# Compare probablities between the original PST model and the test sequence PST probs
# What are the differences and where? This is the starting point of Bayesian Surprise.
# Use the entropy, BS and mi variables for plotting.
compare_surprise <- function(pstmodel, traindata, testdata){
  
  differences <- data.frame(entropy=, BS=, mi=)
  probtrain <- cprob(traindata, L = 0, prob = TRUE)  # overall probs for training dataset
  
  probtest <- PST::predict(model3.pst, testdata, decomp=TRUE,output = "prob")
  
  
  return(differences)
}


## Bayesian interest (decay parameter)  
## As new test data arrives and is not seen in the training data set, it may be identified as 
## interesting. But if we see more of the same then it will stop being interesting. So overall, for 
## the test data the Bayesian Interest value will eventually decay to zero. Without the decay of 
## Bayesian surprise each test record will always be interesting no matter how many times it appears.  

bayes_decay <- function(surp){
  
  
  
  return(decay_para)
}

logLoss = function(pred, actual){
  -1*mean(log(pred[model.matrix(~ actual + 0) - pred > 0]))
}

### Correct formula in native R 
logLoss2 <- function(pred, actual){
  -mean(actual * log(pred) + (1 - actual) * log(1 - pred))
}



