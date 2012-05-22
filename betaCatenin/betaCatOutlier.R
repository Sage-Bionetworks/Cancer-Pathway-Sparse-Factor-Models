## BETA CATENIN Outlier Identification
# betaCatOutlier.R

require(car)
require(synapseClient)
require(Biobase)
require(corpcor)

## BRING DOWN THE RMA-NORMALIZED BETA-CATENIN PERTURBATION DATA
bCatEnt <- loadEntity('syn299161')
bCatMat <- exprs(bCatEnt$objects$rmaEset)
bCatTx <- bCatEnt$objects$treatment

## LOOK AT THE DATA IN PC SPACE
svdObj <- fast.svd(bCatMat)

## FIND OUTLIERS
outlierObj <- outlierTest(lm(svdObj$v[ , 1] ~ 1), n.max = 19)
