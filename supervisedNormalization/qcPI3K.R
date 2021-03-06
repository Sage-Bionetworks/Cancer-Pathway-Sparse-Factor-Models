# qcPI3K.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

##########
# LOAD IN REQUIRED LIBRARIES
##########

require(mGenomics)
require(snm)
require(ggplot2)
require(Biobase)
require(synapseClient)

# synapseLogin('username', 'password')

##########
# SOURCE IN HELPER FUNCTIONS FROM GITHUB
##########

githubRepo <- 
  downloadRepo('Sage-Bionetworks/Cancer-Pathway-Sparse-Factor-Models')

source(paste(githubRepo@localPath, 
             '/helperFunctions/generateQCFigures.R', sep = ''))

source(paste(githubRepo@localPath, 
             '/helperFunctions/qcFunctions.R', sep = ''))

##########
# LOAD IN RAW DATA ENTITY FROM SCR
##########

pi3kEnt <- loadEntity('syn138520')
fits <- runWorkflow(pi3kEnt$cacheDir, workflow = 'snm')

##########
# PULL OUT THE EXPRESSION DATA
##########

exprDat <- exprs(fits[[1]][[1]])

# Create a treatment model matrix  (using the filename annotations)
treatment <- ifelse(grepl('P11', list.files(pi3kEnt$cacheDir)), "PI3K", "GFP")
X <- model.matrix(~ factor(treatment))
sigObj <- calcSig(exprDat, X)

##########
# GENERATE FIGURES ON DATA PRIOR TO SUPERVISED NORMALIZATION
##########

# custom function from generateFigures.R
varBarPlot <- qcFigureA(sigObj)

svdObj <- fs(exprDat)

# custom function from generateFigures.R
initPcPlots <- qcFigureB(svdObj)

##########
# REGRESS OUT THE KNOWN EXPERIMENTAL PERTURBATION EFFECT (TREATMENT)
##########

# In order to better understand perturbation-independent latent structure in
# the data

propSSQ <- removeExpEffect(exprDat, X)

# Generate a figure
propSSQBarPlot <- qcFigureC(propSSQ)

# Visual inspection of this figure suggests that the dependence kernel 
# should be set as rank 3
svaFit <- sva(exprDat, bio.var = X, n.sv = 3, num.iter = 30, diagnose = FALSE)

# Now, we'll visualize the estimated basis vectors 
subtractedPCPlot <- qcFigureD(svaFit)

# The plots suggest that even with the experimental perturbation effect
# removed there is still latent structure that correlates with the 
# experimental treatment.

# One way to manages this is to estimate a dependence kernel from transcripts
# that are non-varying, as determined using the pi0 statistic to define a 
# significance cutoff. Basically, a dependence kernel is generated by taking
# an unweighted SVD of the data restricted to the first pi0 largest ranked 
# (by p values) transcripts.
nullNorm <- nullProbeNorm(sigObj, 3, pi3kEnt)

# Visualize the normalized data
nullNormFig <- qcFigureE(nullNorm$uMatrix)

##########
# MAKE AN ESET OUT OF THE NORMALIZED DATA
##########

nPI3KEset <- nullNorm$newFit[[1]][[1]]
tempPhen <- pData(nPI3KEset)
tempPhen$treatment <- treatment
pData(nPI3KEset) <- tempPhen

##########
# UPLOAD TO SYNAPSE
##########

sNormProjEnt <- getEntity('syn1091010')
sNormPI3KEnt <- Data(list(name = "Supervised Normalized PI3K Perturbation Data",
                          parentId = propertyValue(sNormProjEnt, "id")))
annotValue(sNormPI3KEnt, 'treatmentString') <- 'PI3K'
annotValue(sNormPI3KEnt, 'assayPlatform') <- 
  fits[[1]]$eset@annotation
sNormPI3KEnt <- addObject(sNormPI3KEnt, nPI3KEset)
sNormPI3KEnt <- storeEntity(sNormPI3KEnt)

