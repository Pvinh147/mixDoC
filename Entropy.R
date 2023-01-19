# Entropy Source File
# Requirements: mxAlgebra for class weights is standardized and named "stdweights"

### Entropy Function -----------------------------------------------------------
# Function for individual class probabilities
indClassProbs <- function(model){
  cp <- mxEval(stdweights, model)
  cp2 <- as.vector(cp)
  cps <- diag(length(cp2))
  diag(cps) <- cp2
  subs <- model@submodels
  of <- function(num) {
    return(mxEval(objective, subs[[num]]))
  }
  rl <- sapply(1:length(names(subs)), of)
  raw <- (rl %*% cps)
  tot <- 1 / apply(raw, 1, sum)
  div <- matrix(rep(tot, length(cp2)), ncol=length(cp2))
  icp <- raw * div
  return(icp)
  
}



# Function for entropy
entropy <- function(classProbs){
  # this function takes a matrix of class probabilities from a
  # mixture model with people in rows and classes in columns
  # and returns the entropy statistic, as given by
  # Ramaswamy, Desarbo, Reibstein, and Robinson 1993
  n <- dim(classProbs)[1]
  k <- dim(classProbs)[2]
  e <- 1-(sum(na.omit(-classProbs*log(classProbs)))/(n*log(k)))
  return(e)
}