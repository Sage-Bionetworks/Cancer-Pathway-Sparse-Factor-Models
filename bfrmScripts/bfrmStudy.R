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

## Source in the 'generateFacs' function from Github
# githubRepo <- 
#   downloadRepo('Sage-Bionetworks/Cancer-Pathway-Sparse-Factor-Models')
# 
# source(paste(githubRepo@localPath, 
#              '/factorGeneration/generateFacs.R', sep = ''))

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
