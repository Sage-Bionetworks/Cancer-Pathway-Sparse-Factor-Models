#####
## PROGRAM WRITTEN TO CHECK UNSUPERVISED NORMALIZATION OF DATASETS
#####

require(synapseClient)
require(affy)
require(snm)
require(metaGEO)


synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn134305"')

## FUNCTION TO CREATE DIAGNOSTIC PLOTS
function(dataEnt, grepString){
  tmpEset <- dataEnt$rmaEset
  treatment <- ifelse(grepl(grepString, sampleNames(tmpEset)), 1, 0)
  dataEnt <- addObject(dataEnt, treatment)
  
  tmpFit <- snm(exprs(rmaEset), adj.var=model.matrix(~factor(treatment)), rm.adj=T)
  tmpMat <- tmpFit$norm.mat
  u <- fs(tmpMat)
  f <- tempfile(pattern="pctVar", fileext=".png")
  png(f)
  plot(u$d, main="percent variance explained after removing treatment variable", ylab="% variance explained")
  dev.off()
  dataEnt <- addFile(dataEnt, f)
  
  f <- tempfile(pattern="svd1v2", fileext=".png")
  png(f)
  print(xyplot(u$v[, 1] ~ u$v[, 2], groups = grepl("bcat", colnames(rawMat)), main="svd after removing treatment variable", ylab="1st svd", xlab="2nd svd"))
  dev.off()
  dataEnt <- addFile(dataEnt, f)
}




