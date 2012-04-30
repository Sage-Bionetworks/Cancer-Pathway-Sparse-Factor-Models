#####
## PROGRAM WRITTEN TO NORMALIZE DATA
#####

options(stringsAsFactors=F)

require(synapseClient)
require(affy)

theseGSEs <- c("GSE5851", "GSE7146", "GSE8970", "GSE9633")

res <- lapply(as.list(theseGSEs), function(x){
  myQuery <- synapseQuery(paste("SELECT id, name FROM study WHERE name == '", x, "'", sep=""))
  myQuery2 <- synapseQuery(paste("SELECT id, name FROM entity WHERE entity.parentId=='", myQuery$study.id, "'", sep=""))
  
  tmp <- downloadEntity(myQuery2$entity.id)
  theseFiles <- file.path(tmp$cacheDir, tmp$files)[ grep("U133", sapply(as.list(file.path(tmp$cacheDir, tmp$files)), whatcdf)) ]
  tmpAB <- ReadAffy( filenames=theseFiles )
  rmaEset <- rma(tmpAB, normalize=T, background=F)
  
  tmpEnt <- ExpressionData(list(name=paste(x, "- RMA normalized"),
                                parentId="syn134305"))
  tmpEnt <- createEntity(tmpEnt)
  tmpEnt <- addObject(tmpEnt, rmaEset)
  tmpEnt <- storeEntity(tmpEnt)
  tmpEnt
})

