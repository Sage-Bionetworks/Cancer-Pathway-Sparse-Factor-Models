## bfrmHelperFunctions.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

# parseBFRM.R
# Functions to parse BFRM output into a list object of Factors by feature indices
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
  
  ## Annotate these indices with GENE SYMBOLS
  platform <- annotValue(bfrmResult, 'assayPlatform')
  annotQueryResult <- 
      synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn308579"')
  platformEntID <- annotQueryResult[platform, 2]
  platformEnt <- loadEntity(platformEntID)
  platformAnnotations <- platformEnt$objects[[1]]
  
  annotateList <- function(Indices, platformAnnotations){
    geneSymbolsMappings <- platformAnnotations[Indices, 'Symbol']
  }
  annotFacList <- lapply(facList, )
  
  return(facList)
} 
