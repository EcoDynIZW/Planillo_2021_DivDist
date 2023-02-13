## 2018.08.24

### FUNCTIONS FOR MULTIVARIATE MODEL SELECTION ###

## Function to compute AICc values of multivariate models 
# Function based on Nakamura et al. 2015 to get AICc values of the model ALTHOUGH CORRECTED: average by species number (univariate models)
# Included in package mglmn, but here applied to specific model

get.AICc <- function(mymodel, myresponse){
  n.param <- length(mymodel$coefficients)/ncol(myresponse) + 1 # +1 because scaling parameter negative binomial
  log.L.temp <- logLik(mymodel)
  AICc.value <- sum(-2*log.L.temp + 2*nrow(myresponse)*(n.param)/(nrow(myresponse) - (n.param)-1))
  return(AICc.value)
}

get.avgAICc <- function(mymodel, myresponse){
  n.param <- length(mymodel$coefficients)/ncol(myresponse) + 1 # +1 because scaling parameter negative binomial
  log.L.temp <- logLik(mymodel)
  AICc.value <- sum(-2*log.L.temp + 2*nrow(myresponse)*(n.param)/(nrow(myresponse) - (n.param)-1))/ncol(myresponse)
  return(AICc.value)
}

# A different way with more parameters but it should be the dame
# get.avgAICc <- function(mymodel, mydata, myresponse){
#   n.param <- length(mymodel$coefficients)/ncol(myresponse) + 1 # +1 because scaling parameter negative binomial
#   log.L.temp <- logLik(mymodel)
#   AICc.value <- sum(-2*log.L.temp + 2*nrow(mydata)*(n.param)/(nrow(mydata) - (n.param)-1))/ncol(myresponse)
#   print(AICc.value)
# }


### For gaussian models

get.AICc.gaus <- function(mymodel, myresponse){
  n.param <- length(mymodel$coefficients)/ncol(myresponse) 
  log.L.temp <- logLik(mymodel)
  AICc.value <- sum(-2*log.L.temp + 2*nrow(myresponse)*(n.param)/(nrow(myresponse) - (n.param)-1))
  return(AICc.value)
}

get.avgAICc.gaus <- function(mymodel, myresponse){
  n.param <- length(mymodel$coefficients)/ncol(myresponse) 
  log.L.temp <- logLik(mymodel)
  AICc.value <- sum(-2*log.L.temp + 2*nrow(myresponse)*(n.param)/(nrow(myresponse) - (n.param)-1))/ncol(myresponse)
  return(AICc.value)
}