## bfrmHelperFunctions.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

# parseBFRM.R
# Functions to parse BFRM output into a list object of Factors by feature indices
# that meet a 0.99 posterior probability threshold

parseBFRM <- function(bfrmResultEnt){
  require(bfrm)
  require(synapseClient)
  print('Loading Posterior Probabilities')
  bfrmResult <- bfrmResultEnt$objects$evolveRes
  mVarIn <- bfrmResultEnt$objects$evolveRes@results$mVariablesIn
  allFeatNames <- bfrmResultEnt$objects$featureNames
  incFeatNames <- allFeatNames[mVarIn]
  mPostPib <- bfrmResult@results$mPostPib
  numCols <- dim(mPostPib)[2]
  facParse <- function(x){
    incFeatureLogical <- x > 0.99
    incFeatureInd <- grep('TRUE', incFeatureLogical)
    facFeatureNames <- incFeatNames[incFeatureInd]
    return(facFeatureNames)
  }
  print('Creating list of Factors by Features with Posterior Probability > 0.99')
  facList <- apply(mPostPib[ , 2:numCols], 2, facParse)
  
  ## Annotate these indices with GENE SYMBOLS
  print('Annotating the Factors by Gene Symbol')
  platform <- annotValue(bfrmResultEnt, 'assayPlatform')
  annotQueryResult <- 
      synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn308579"')
  platformEntID <- annotQueryResult[grep(platform, annotQueryResult[ , 1]), 2]
  platformEnt <- loadEntity(platformEntID)
  platformDir <- platformEnt$cacheDir
  platformAnnotations <- read.delim(paste(platformDir, 'probeSetAnnotations.txt', sep = '/'),
                                    header = TRUE, sep = '\t')
  platAnnot <- platformAnnotations[ , 2:3]
  rownames(platAnnot) <- platformAnnotations[ , 1]
  annotateList <- function(Indices, platAnnot){
    geneSymbolsMappings <- platAnnot[Indices, 'Symbol']
  }
  annotFacList <- lapply(facList, annotateList, platAnnot)
  names(annotFacList) <- paste('Factor', 1:(numCols - 1), sep = '')
  return(annotFacList)
} 
