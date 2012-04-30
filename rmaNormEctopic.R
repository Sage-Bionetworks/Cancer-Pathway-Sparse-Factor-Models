#####
## PROGRAM WRITTEN TO NORMALIZE DATA
#####

options(stringsAsFactors=F)

require(synapseClient)
require(affy)

theseData <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn162399"')

res <- lapply(as.list(theseData$entity.id), function(x){
  tmp <- downloadEntity(x)
  tmpAB <- ReadAffy( celfile.path = tmp$cacheDir )
  rmaEset <- rma(tmpAB, normalize=T, background=F)
  
  tmpEnt <- ExpressionData(list(name=paste(propertyValue(tmp, "name"), "- RMA normalized"),
                                parentId="syn134305"))
  tmpEnt <- createEntity(tmpEnt)
  tmpEnt <- addObject(tmpEnt, rmaEset)
  tmpEnt <- storeEntity(tmpEnt)
  tmpEnt
})



