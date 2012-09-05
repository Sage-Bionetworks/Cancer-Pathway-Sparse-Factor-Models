## generateQCFigures.R

## Erich S. Huang
## Sage Bionetworks
## erich.huang@sagebase.org


## First figure function
## p-value distribution of transcripts by perturbation
qcFigureA <- function(sigObj){
  require(ggplot2)
  histFig <- qplot(sigObj$pval, geom = 'histogram') + 
    opts(title = 'p-value Distribution of Transcripts\n') +
    xlab('\np-values') +
    ylab('Experiment Effect on Expression Variation\n')
  return(histFig)
}

## Second figure function
## Eigengene plots on data prior to supervised normalization
qcFigureB <- function(svdObj){
  require(ggplot2)
  require(reshape)
  
  dMat <- svdObj[[1]]
  vMat <- svdObj[[3]]
  pcDim <- dim(vMat)[1]
  svdDF <- as.data.frame(cbind(1:pcDim, svdObj$d,
                               svdObj$v))
  colnames(svdDF) <- c('sample', 'percentVariance',
                      paste('PC', 1:pcDim, sep = ''))
  barPlot <- ggplot(svdDF, aes(sample, percentVariance)) +
    geom_bar(stat = 'identity') +
    opts(title = 'Percent Variance Explained by Each Experimental Eigengene\n') +
    xlab('\nEigengene') +
    ylab('Percent Variance Explained\n')
  
  svdDF2 <- cbind(svdDF[ , 1], treatment, svdDF[ , 3:dim(svdDF)[2]])
  colnames(svdDF2)[1] <- 'sample'
  
  meltDF <- melt(svdDF2, id = c('sample', 'treatment'))
  colnames(meltDF)[3] <- 'pc'
  smeltDF <- subset(meltDF, pc %in% c('PC1', 'PC2', 'PC3', 'PC4'))

  facetPcPlot <- qplot(sample, value, data = smeltDF) +
    facet_grid(. ~ pc) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = 'Eigenvalues by Sample\n') +
    xlab('\nSample') +
    ylab('Eigenvalue\n')
  
  allFigs <- list('barPlot' = barPlot, 
                  'facetPlot' = facetPcPlot)
}

## Third figure function
## Visualizing the proportion of SSQ explained by each perturbation eigengene

qcFigureC <- function(propSSQ){
  propSSQDF <- as.data.frame(cbind(1:length(treatment), propSSQ))
  colnames(propSSQDF) <- c('eigenGene', 'propSSQ')
  
  propSSQPlot <- ggplot(propSSQDF, aes(eigenGene, propSSQ)) +
    geom_bar(stat = 'identity') +
    opts(title = 'Proportion of total SSQ\n') +
    xlab('\nEigengene') +
    ylab('Proportion of total Sum of Squares Explained\n')
}


## Fourth figure function
## Subtracting treatment effect from the data to understand whether there are
## other latent variables influencing the data

qcFigureD <- function(svaFit){
  svaDF <- data.frame(1:length(treatment), treatment,
                               svaFit$svd[[svaFit$num.iter]]$v)
  colnames(svaDF) <- c('sample', 'treatment', paste('subtractedPC', 
                                       1:svaFit$n.sv, sep = ''))
  transform(svaDF, sample = as.numeric(sample), 
            subtractedPC1 = as.numeric(subtractedPC1),
            subtractedPC2 = as.numeric(subtractedPC2))
  meltDF <- melt(svaDF, id = c('sample', 'treatment'))
  colnames(meltDF)[3] <- 'pc'
  meltDF <- meltDF[order(meltDF[ , 2]), ]
  subFacetPlot <- qplot(sample, value, data = meltDF) +
    facet_grid(. ~ pc) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = 'Subtraction of Biological Effect\n') +
    xlab('\nSample') +
    ylab('Eigenvalue\n')
}


## Fifth figure function
## Visualizing the data after normalizing out null probe latent structure

qcFigureE <- function(u2){
  nullNormDF <- data.frame(1:length(treatment), treatment, u2$v)
  colnames(nullNormDF) <- c('sample', 
                            'treatment',
                            paste('nullNormPC', 1:length(treatment), sep = ''))
  meltDF <- melt(nullNormDF, id = c('sample', 'treatment'))
  colnames(meltDF)[2:3] <- c('treatment', 'pc')
  smeltDF <- subset(meltDF, pc %in% c('nullNormPC1', 'nullNormPC2',
                                      'nullNormPC3', 'nullNormPC4'))
  facetPcPlot <- qplot(sample, value, data = smeltDF) +
    facet_grid(. ~ pc) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = 'Eigenvalues by Sample\n') +
    xlab('\nSample') +
    ylab('Eigenvalue\n')
}


