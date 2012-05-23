#####
## PROGRAM WRITTEN TO NORMALIZE DATA
#####

options(stringsAsFactors = F)

require(synapseClient)
require(affy)
require(corpcor)
require(car)

fs <-function (x){
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  u
}

theseData <- synapseQuery('SELECT id, name FROM entity WHERE entity.parentId == "syn162399"')

res <- lapply(as.list(theseData$entity.id), function(x){
  tmp <- downloadEntity(x)
  tmpAB <- ReadAffy( celfile.path = tmp$cacheDir )
  rmaEset <- rma(tmpAB, normalize=T, background=F)
  treatment <- ifelse(grepl(annotValue(tmp, "treatmentString"), 
                            tolower(sampleNames(rmaEset))), 1, 0)
  tmpPhen <- pData(rmaEset)
  tmpPhen$treatment <- treatment
  pData(rmaEset) <- tmpPhen
  
  ## PLOTS CREATED AND STORED SUCH THAT CAN BE ATTACHED VIA WEB UI AS 
  ## ATTACHMENTS ON YET SUPPORTED IN R CLIENT
  myDir <- tempfile(pattern=strsplit(propertyValue(tmp, "name"), " ")[[1]][2], 
                    tmpdir = path.expand("~/"), fileext="")
  dir.create(myDir)
  
  ## TO IDENTIFY OUTLIERS, FIRST MODEL OUT TREATMENT EFFECTS
  treatMM <- model.matrix(~ factor(treatment))
  expressMat <- exprs(rmaEset)
  tmpMat <- expressMat - t(treatMM %*% solve(t(treatMM) %*% treatMM) %*% 
    t(treatMM) %*% t(expressMat))
  svdObj <- fs(tmpMat)
  
  ## THEN USE THE CAR PACKAGE outlierTest() TO ASSESS THE THREE FIRST PRINCIPAL
  ## AXES AND RETURN THE SAMPLES THAT ARE OUTLIERS
  top3Axes <- svdObj$v[ , 1:3]
  outlierSamples <- apply(top3Axes, 2, function(column){
    outlierObj <- outlierTest(lm(column ~ 1))
    if (outlierObj$signif == 'TRUE'){
      outSamp <- names(outlierObj$rstudent)
    } else {
      outSamp <- NA
    }
  })
  
  expressionset[ -na.omit(as.numeric(outlierSamples)) ]
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
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



