## bfrmStudy.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## Required libraries
require(synapseClient)
require(bfrm)
require(doMC)
require(rGithubClient)

# Using rGithubClient to source in the 'generateFacs' function
genFacsBlob <- 
  downloadRepoBlob('Sage-Bionetworks/Cancer-Pathway-Sparse-Factor-Models',
                   'factorGeneration/generateFacs.R')
sourceBlob(genFacsBlob, 'factorGeneration/generateFacs.R')

## Pull in the IDs of all normalized studies
normalizedEsets <- 
  synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn1091010"')

childIDs <- normalizedEsets[ , 2]

## Parallelized generation of BFRM objects
registerDoMC()
numProc <- getDoParWorkers()

bfrmList <- mclapply(childIDs, generateFacs)
