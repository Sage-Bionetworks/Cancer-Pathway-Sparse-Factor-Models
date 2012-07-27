## BETA CATENIN DATA PREPARATION
# betaCatDataPrep.R

require(car)
require(synapseClient)
require(Biobase)
require(corpcor)
require(snm)

## BRING DOWN THE RMA-NORMALIZED BETA-CATENIN PERTURBATION DATA
bCatEnt <- loadEntity('syn299161')

## PULL IN THE CODE ENTITY WITH HELPER FUNCTIONS
pathFacFuncEnt <- loadEntity('syn308412')

## USE THE 'removeTreatment()' FUNCTION TO TAKE THE DATA, REMOVE THE TREATMENT
## EFFECT AND IDENTIFY OUTLIERS
svdObj <- removeTreatment(bCatEnt)

## INSPECT THE TOP PRINCIPAL COMPONENTS FOR OUTLIER SAMPLES
outlierPC1 <- findOutlier(svdObj, 1)
outlierPC2 <- findOutlier(svdObj, 2)
outlierPC3 <- findOutlier(svdObj, 3)

  # Looking at the first three principal axes, only sample 15 appears to be an
  # outlier.

## USE SNM TO RENORMALIZE THE DATA
adjVar <- model.matrix( ~ ifelse(svdObj$v[ , 1] < -0.8, 1, 0))
  # The outlier sample 15 is < -0.8 on principal axis 1

bioVar <- model.matrix( ~ factor(bCatEnt$objects$treatment))
bCatFit <- snm(exprs(bCatEnt$objects$rmaEset),
               bio.var = bioVar,
               adj.var = adjVar,
               rm.adj = T)

## FOR A REALITY CHECK, LET'S LOOK AT THE OUTLIER ADJUSTED DATA, MODEL OUT THE
## TREATMENT EFFECT AGAIN, AND INSPECT FOR OUTLIERS
Z <- model.matrix(~ factor(bCatEnt$objects$treatment))
Y <- bCatFit$norm.dat
tmpMat <- Y - t(Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% t(Y))
svdObjNew <- fs(tmpMat)

## AGAIN, INSPECT THE TOP THREE PRINCIPAL COMPONENTS FOR OUTLIERS
outlierPC1N <- findOutlier(svdObjNew, 1)
outlierPC2N <- findOutlier(svdObjNew, 2)
outlierPC3N <- findOutlier(svdObjNew, 3)

  # Looks good.

snmBCatEnt <- loadEntity('syn308414')
snmBCatEnt <- addObject(snmBCatEnt, bCatFit$norm.dat)

