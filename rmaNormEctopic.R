#####
## PROGRAM WRITTEN TO NORMALIZE DATA
#####

options(stringsAsFactors=F)

require(synapseClient)
require(affy)
require(snm)
require(metaGEO)

theseData <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn162399"')

res <- lapply(as.list(theseData$entity.id), function(x){
  tmp <- downloadEntity(x)
  tmpAB <- ReadAffy( celfile.path = tmp$cacheDir )
  rmaEset <- rma(tmpAB, normalize=T, background=F)
  treatment <- ifelse(grepl(annotValue(tmp, "treatmentString"), tolower(sampleNames(rmaEset))), 1, 0)
  
  ## PLOTS CREATED AND STORED SUCH THAT CAN BE ATTACHED VIA WEB UI AS ATTACHMENTS ON YET SUPPORTED IN R CLIENT
  myDir <- tempfile(pattern="dir", tmpdir = path.expand("~/"), fileext="")
  dir.create(myDir)
  
  tmpFit <- snm(exprs(rmaEset), adj.var=model.matrix(~factor(treatment)), rm.adj=T)
  tmpMat <- tmpFit$norm.dat
  
  u <- fs(tmpMat)
  f1 <- file.path(myDir, "pctVar.png")
  png(f1)
  plot(u$d, main="percent variance explained after removing treatment variable", ylab="% variance explained")
  dev.off()
  
  f2 <- file.path(myDir, "svd1v2.png")
  png(f2)
  print(xyplot(u$v[, 1] ~ u$v[, 2], groups = treatment, main="svd after removing treatment variable", ylab="1st svd", xlab="2nd svd"))
  dev.off()
  
  tmpEnt <- ExpressionData(list(name=paste(propertyValue(tmp, "name"), "- RMA normalized"),
                                parentId="syn299160"))
  tmpEnt <- createEntity(tmpEnt)
  tmpEnt <- addObject(tmpEnt, rmaEset)
  tmpEnt <- addObject(tmpEnt, treatment)
  tmpEnt <- storeEntity(tmpEnt)
  tmpEnt
})



