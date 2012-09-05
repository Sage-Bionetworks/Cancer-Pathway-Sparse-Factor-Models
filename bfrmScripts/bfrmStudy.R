## bfrmStudy.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org


## Pull in the IDs of all normalized studies
normalizedEsets <- 
  synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn1091010"')

childIDs <- normalizedEsets[ , 2]

bfrmList <- lapply(testIDs, generateFacs)