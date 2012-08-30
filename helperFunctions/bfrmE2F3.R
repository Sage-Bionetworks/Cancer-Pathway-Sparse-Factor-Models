## bfrmE2F3.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

##########
# RUN BFRM IN THE SPARSE ANOVA MODE
##########
require(bfrm)
require(synapseClient)



e2f3Anova <- bfrm(dat2, design = ifelse(treatment == 'E2F3', 1, 0))
mPPib <- e2f3Anova@results$mPostPib
topProbeLogical <- mPPib[ , 2] >= 0.99
topProbeInd <- grep("TRUE", topProbeLogical)

##########
# RUN BFRM IN THE FACTOR DISCOVERY MODE
##########

bCatEvolveFactor <- evolve(dat2, 
                           init = as.numeric(topProbeInd),
                           priorpsia = 2,
                           priorpsib = 0.005,
                           varThreshold = 0.85,
                           facThreshold = 0.95,
                           maxVarIter = 30,
                           minFacVars = 10,
                           maxFacVars = length(topProbeInd),
                           maxFacs = 50,
                           maxVars = length(topProbeInd)
                           )