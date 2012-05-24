#####
## A SCRIPT FOR NORMALIZING THE RAW ONCOGENIC PERTURBATION DATA
#####

#####
## SET OPTIONS
#####
options(stringsAsFactors = F)

#####
## REQUIRE LIBRARIES
#####
require(synapseClient)
require(affy)
require(corpcor)
require(car)
require(ggplot2)

#####
## BRING IN ANY NECESSARY SYNAPSE CODE ENTITIES
####
multiPlotEnt <- loadEntity('syn274067')
attach(multiPlotEnt)

#####
## DEFINE FUNCTIONS
#####
# function [1] standardized fast svd
fs <-function(x){
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  u
}

# function [2] rma normalization and creating esets
makeEsets <- function(x){
  tmp <- downloadEntity(x)
  tmpAB <- ReadAffy(celfile.path = tmp$cacheDir)
  rmaEset <- rma(tmpAB, normalize=T, background=F)
  treatment <- ifelse(grepl(annotValue(tmp, "treatmentString"), 
                            tolower(sampleNames(rmaEset))), 1, 0)
  tmpPhen <- pData(rmaEset)
  tmpPhen$treatment <- treatment
  pData(rmaEset) <- tmpPhen
  return(rmaEset)
}

# function [3] removing treatment effects (to identify outliers)
removeTx <- function(x){
  tmpExpress <- exprs(x)
  tmpTreatment <- x$treatment
  treatMM <- model.matrix(~ factor(tmpTreatment))
  tmpMat <- tmpExpress - t(treatMM %*% solve(t(treatMM) %*% treatMM) %*% 
    t(treatMM) %*% t(tmpExpress))
  svdObj <- fs(tmpMat)
  return(list('svdObj' = svdObj, "tmpMat" <- tmpMat))
}

# function [4] for visualizing the data before adjustment
vizData <- function(rmaEset, pc = 1){
  rmaSVD <- fs(exprs(rmaEset))
  rmaDF <- as.data.frame(cbind(1:dim(rmaEset)[2], rmaSVD$v[ , 1:3]))
  rmaTx <- rmaEset$treatment
  colnames(rmaDF)[2:4] <- paste('Eigen', 1:3, sep = '')
  colnames(rmaDF)[1] <- 'Samples'
  rmaFig <- ggplot(rmaDF, aes(x = rmaDF[ , 1], rmaDF[ , pc + 1])) + 
    geom_point(aes(colour = factor(rmaTx))) +
    opts(title = 'SVD Plot of Unadjusted Data') +
    xlab('Samples') + ylab(paste('EigenGene', pc))
  return(rmaFig)
}

# function [5] for visualizing the data after adjustment
vizAdjData <- function(bleachList, pc){
  
}

#####
## THE ACTUAL SCRIPT
#####

# pull down the synapse entity references from the synapse project
theseData <- 
  synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn162399"')

# rma normalize the constituent raw datasets and wrap them up as esets
rmaList <- lapply(as.list(theseData$entity.id), makeEsets)
names(rmaList) <- theseData[ , 1]

# remove the treatment effects so we can identify outlier samples
bleachList <- lapply(rmaList, removeTx)
names(bleachList) <- theseData[ , 1]

# visualize the unadjusted data in principle component space
beforeViz <- lapply(rmaList, vizData)
names(beforeViz) <- theseData[ , 1]



# res <- lapply(as.list(theseData$entity.id), function(x){
#   tmp <- downloadEntity(x)
#   tmpAB <- ReadAffy( celfile.path = tmp$cacheDir )
#   rmaEset <- rma(tmpAB, normalize=T, background=F)
#   treatment <- ifelse(grepl(annotValue(tmp, "treatmentString"), 
#                             tolower(sampleNames(rmaEset))), 1, 0)
#   tmpPhen <- pData(rmaEset)
#   tmpPhen$treatment <- treatment
#   pData(rmaEset) <- tmpPhen
#   
#   ## PLOTS CREATED AND STORED SUCH THAT CAN BE ATTACHED VIA WEB UI AS 
#   ## ATTACHMENTS ON YET SUPPORTED IN R CLIENT
#   myDir <- tempfile(pattern=strsplit(propertyValue(tmp, "name"), " ")[[1]][2], 
#                     tmpdir = path.expand("~/"), fileext="")
#   dir.create(myDir)
#   
#   ## TO IDENTIFY OUTLIERS, FIRST MODEL OUT TREATMENT EFFECTS
#   treatMM <- model.matrix(~ factor(treatment))
#   expressMat <- exprs(rmaEset)
#   tmpMat <- expressMat - t(treatMM %*% solve(t(treatMM) %*% treatMM) %*% 
#     t(treatMM) %*% t(expressMat))
#   svdObj <- fs(tmpMat)
#   
#   ## THEN USE THE CAR PACKAGE outlierTest() TO ASSESS THE THREE FIRST PRINCIPAL
#   ## AXES AND RETURN THE SAMPLES THAT ARE OUTLIERS
#   top3Axes <- svdObj$v[ , 1:3]
#   outlierSamples <- apply(top3Axes, 2, function(column){
#     outlierObj <- outlierTest(lm(column ~ 1))
#     if (outlierObj$signif == 'TRUE'){
#       outSamp <- as.numeric(names(outlierObj$rstudent))
#     } else {
#       outSamp <- 0
#     }
#   })
#   
#   rmaEset[ -na.omit(as.numeric(outlierSamples)) ]
#   
#   
# 
#   
  
  
  
  
  
  
  
  
  
  
  
  
  
#   tmpFit <- snm(exprs(rmaEset), adj.var=model.matrix(~factor(treatment)), rm.adj=T)
#   tmpMat <- tmpFit$norm.dat
#   
#   u <- fs(tmpMat)
#   f1 <- file.path(myDir, "pctVar.png")
#   png(f1)
#   plot(u$d, main="percent variance explained after removing treatment variable", ylab="% variance explained")
#   dev.off()
#   
#   f2 <- file.path(myDir, "svd1v2.png")
#   png(f2)
#   print(xyplot(u$v[, 1] ~ u$v[, 2], groups = treatment, main="svd after removing treatment variable", ylab="1st svd", xlab="2nd svd"))
#   dev.off()
#   
#   tmpEnt <- ExpressionData(list(name=paste(propertyValue(tmp, "name"), "- RMA normalized"),
#                                 parentId="syn299160"))
#   tmpEnt <- createEntity(tmpEnt)
#   tmpEnt <- addObject(tmpEnt, rmaEset)
#   tmpEnt <- addObject(tmpEnt, treatment)
#   tmpEnt <- storeEntity(tmpEnt)
#   tmpEnt
})



