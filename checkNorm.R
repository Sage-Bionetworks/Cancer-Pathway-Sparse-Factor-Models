#####
## PROGRAM WRITTEN TO CHECK UNSUPERVISED NORMALIZATION OF DATASETS
#####

options(stringsAsFactors=F)

require(synapseClient)
require(affy)
require(snm)
require(metaGEO)


synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn134305"')

## FUNCTION TO CREATE DIAGNOSTIC PLOTS
normEval <- function(dataEnt, trtString){
  tmpEset <- dataEnt$objects$rmaEset
  treatment <- ifelse(grepl(trtString, tolower(sampleNames(tmpEset))), 1, 0)
#  dataEnt <- addObject(dataEnt, treatment)
  
  tmpFit <- snm(exprs(tmpEset), adj.var=model.matrix(~factor(treatment)), rm.adj=T)
  tmpMat <- tmpFit$norm.dat
  u <- fs(tmpMat)
#  f <- tempfile(pattern="pctVar", fileext=".png")
#  png(f)
  plot(u$d, main="percent variance explained after removing treatment variable", ylab="% variance explained")
#  dev.off()
#  dataEnt <- addFile(dataEnt, f)
  
#  f <- tempfile(pattern="svd1v2", fileext=".png")
#  png(f)
  print(xyplot(u$v[, 1] ~ u$v[, 2], groups = treatment, main="svd after removing treatment variable", ylab="1st svd", xlab="2nd svd"))
#  dev.off()
#  dataEnt <- addFile(dataEnt, f)
#  dataEnt <- storeEntity(dataEnt)
}

## BETA CATENIN (one outlier)
checkEnt <- loadEntity("syn299099")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "bcat")

## E2F3
checkEnt <- loadEntity("syn299101")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "e2f3")

## EGFR (two outliers)
checkEnt <- loadEntity("syn299103")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "wtegf")

## Myc (3 outliers)
checkEnt <- loadEntity("syn299105")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "m")

## PI3k
checkEnt <- loadEntity("syn299107")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "p11")

## RAS
checkEnt <- loadEntity("syn299109")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "ras")

## SRC
checkEnt <- loadEntity("syn299111")
sampleNames(checkEnt$objects$rmaEset)
res <- normEval(checkEnt, "src")





# ## GEO SPECIFIC
# ## FUNCTION TO CREATE DIAGNOSTIC PLOTS
# normEvalGEO <- function(dataEnt, treatment){
#   tmpEset <- dataEnt$objects$rmaEset
#   dataEnt <- addObject(dataEnt, treatment)
#   
#   tmpFit <- snm(exprs(tmpEset), adj.var=model.matrix(~factor(treatment)), rm.adj=T)
#   tmpMat <- tmpFit$norm.dat
#   u <- fs(tmpMat)
#   f <- tempfile(pattern="pctVar", fileext=".png")
#   png(f)
#   plot(u$d, main="percent variance explained after removing treatment variable", ylab="% variance explained")
#   dev.off()
#   dataEnt <- addFile(dataEnt, f)
#   
#   f <- tempfile(pattern="svd1v2", fileext=".png")
#   png(f)
#   print(xyplot(u$v[, 1] ~ u$v[, 2], groups = treatment, main="svd after removing treatment variable", ylab="1st svd", xlab="2nd svd"))
#   dev.off()
#   dataEnt <- addFile(dataEnt, f)
# }
# 
# 
# ## GSE5851
# checkEnt <- loadEntity("syn299131")
# sampleNames(checkEnt$objects$rmaEset)
# res <- normEvalGEO(checkEnt, "src")
# 
