# caPathsFacsFunctions.R
# Functions to help with the normalization and handling of data


## FUNCTION 1: TAKE THE RMA NORMALIZED ESETS FOR CANCER PATHWAY PERTURBATION 
## DATA FROM SYNAPSE, ADJUST THE DATA TO REMOVE TREATMENT EFFECTS, PERFORM A
## FAST SVD AND RETURN THE SVD OBJECT

removeTreatment <- function(synEnt){
  require(mGenomics) # Source from depot.sagebase.org
  Z <- model.matrix(~ factor(synEnt$objects$treatment))
  Y <- exprs(synEnt$objects$rmaEset)
  tmpMat <- Y - t(Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% t(Y))
  svdObj <- fs(tmpMat)
  return(svdObj)
}

findOutlier <- function(svdObj, pc){
  require(car)
  outlierObj <- outlierTest(lm(svdObj$v[ , pc] ~ 1))
  return(outlierObj)
}
