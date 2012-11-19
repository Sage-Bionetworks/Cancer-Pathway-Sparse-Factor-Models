## perturbBFRMParse.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Required libraries
require(synapseClient)
require(bfrm)
require(doMC)
require(rGithubClient)

registerDoMC() 

## Using rGithubClient to source in the 'bfrmHelperFunctions' script
githubRepo <- 
  getRepo('Sage-Bionetworks/Cancer-Pathway-Sparse-Factor-Models')

functionEnv <- sourceRepoFile(githubRepo, 'helperFunctions/bfrmHelperFunctions.R')
  ## Right now the rGithubClient sources into a local environment. Will change this when the new
  ## version sources into globalEnv

## Bring down the bfrm result objects
bfrmResultIDs <- 
  synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn1394611"')

bfrmResultIDList <- as.list(bfrmResultIDs[ , 2])

bfrmResultList <- mclapply(bfrmResultIDList, function(x){
  bfrmResultEnt <- loadEntity(x)
})

names(bfrmResultList) <- sapply(strsplit(bfrmResultIDs[ , 1], ' '), '[[', 1)

## Take the list of result objects and annotate them with gene symbols
# parsedBFRMResults <- mclapply(bfrmResultList, functionEnv$parseBFRM)
        ## For some reason there's unreliable performance with parallelization
        ## of the final BFRM result object annotation, hence a for loop

for(i in 1:length(bfrmResultList)){
  bfrmTarget <- bfrmResultList[[i]]
  bfrmName <- names(bfrmResultList)[i]
  parsedResult <- functionEnv$parseBFRM(bfrmTarget)
  varName <- paste(bfrmName, 'annotatedFactors', sep = '')
  assign(varName, parsedResult)
  newEnt <- Data(list(name = paste(varName),
                      parentId = 'syn1513398'))
  newEnt <- createEntity(newEnt)
  newEnt <- addObject(newEnt, get(varName), varName)
  newEnt$annotations$assayPlatform <- bfrmTarget$annotations$assayPlatform
  newEnt$annotations$derivedFrom <- bfrmTarget$properties$id
  newEnt <- storeEntity(newEnt)
}



