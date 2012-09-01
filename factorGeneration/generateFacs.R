## generateFacs.R

## Brian M. Bot
## Sage Bionetworks
## Seattle, Washington
## brian.bot@sagebase.org

generateFacs <- function(esetId){
  
  require(synapseClient)
  require(bfrm)
  require(Biobase)
  
  ## READ IN EXPRESSION DATA
  normEnt <- loadEntity(esetId)
  normMat <- exprs(normEnt$objects[[1]])
  
  trtInd <- as.numeric(grepl(tolower(annotValue(normEnt, "treatmentString")), 
                             tolower(colnames(normMat))))
  
  ## BFRM PRIMER USING THE PERTUBATION STATUS AS DESIGN VARIABLE
  bfrmRes <- bfrm(normMat, design=trtInd)
  topProbeInd <- which(bfrmRes@results$mPostPib[, 2] >= 0.99)
  
  ## SWITCH INTO EVOLUTIONARY MODE SEEDING THE SEARCH WITH RESULTS FROM BFRM RESULTS ABOVE
  evolveRes <- evolve(normMat, 
                      init = topProbeInd,
                      priorpsia=2,
                      priorpsib=0.005,
                      varThreshold=0.85, 
                      facThreshold=0.95,
                      maxVarIter=30,
                      minFacVars=10,
                      maxFacVars=length(topProbeInd),
                      maxFacs=50,
                      maxVars=length(topProbeInd))
  
  ## UPLOAD BACK INTO FACTOR LIBRARY PROJECT -- bfrmResult objects STUDY
  newEnt <- 
    createEntity(Data(list(name=paste(annotValue(normEnt, "treatmentString"), 
                                              " perturbation - bfrmResult object", sep=""), parentId="syn362376")))
  newEnt <- addObject(newEnt, bfrmRes)
  newEnt <- addObject(newEnt, evolveRes)
  newEnt <- storeEntity(newEnt)
}

