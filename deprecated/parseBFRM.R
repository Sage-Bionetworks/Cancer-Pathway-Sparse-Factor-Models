# parseBFRM.R
# Functions to parse BFRM output into a list object of Factors by Gene Symbols
# that meet a 0.99 posterior probability threshold

parseBFRM <- function(bfrmResult){
  print('Loading Posterior Probabilities')
  mPostPib <- bfrmResult@results$mPostPib
  numFacs <- dim(mPostPib)[2]
  facParse <- function(x){
    incFeatureLogical <- x > 0.99
    incFeatureInd <- grep('TRUE', incFeatureLogical)
    return(incFeatureInd)
  }
  print('Creating list of Factors by Features with Posterior Probability > 0.99')
  facList <- apply(mPostPib[ , 2:numFacs], 2, facParse)
  return(facList)
} 