## generateFigures.R

## Erich S. Huang
## Sage Bionetworks
## erich.huang@sagebase.org


## Multiplot helper function
multiplot <- function(..., plotlist = NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols     # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols)      # Number of rows needed
                                             # calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

## First figure function
## p-value distribution of transcripts by perturbation
pvalHistFig1 <- function(sigObj){
  require(ggplot)
  histFig <- qplot(sigObj$pval, geom = 'histogram') + 
    opts(title = 'p-value Distribution of Transcripts\n') +
    xlab('\np-values') +
    ylab('Experiment Effect on Expression Variation\n')
  return(histFig)
}

## Second figure function
## Eigengene plots on data prior to supervised normalization
pcPlotsFig2 <- function(svdObj){
  require(ggplot)
  dMat <- svdObj[[1]]
  vMat <- svdObj[[3]]
  pcDim <- dim(vMat)[1]
  svdDF <- as.data.frame(cbind(1:pcDim, svdObj$d,
                               svdObj$v))
  colnames(svdDF) <- c('sample', 'percentVariance',
                      paste('PC', 1:pcDim, sep = ''))
  barPlot <- ggplot(svdDF, aes(eigenGene, percentVariance)) +
    geom_bar(stat = 'identity') +
    opts(title = 'Percent Variance Explained by Each Experimental Eigengene\n') +
    xlab('\nEigengene') +
    ylab('Percent Variance Explained\n')
  
  eigenPlot1 <- ggplot(svdDF, aes(eigenGene, PC1)) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = "Loadings (before Supervised Norm.)\n") +
    xlab('\nSamples') +
    ylab('Eigengene 1 Loadings\n')
  
  eigenPlot2 <- ggplot(svdDF, aes(eigenGene, PC2)) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = "Loadings (before Supervised Norm.)\n") +
    xlab('\nSamples') +
    ylab('Eigengene 2 Loadings\n')
  
  eigenPlot3 <- ggplot(svdDF, aes(eigenGene, PC3)) +
    geom_point(aes(colour = factor(treatment))) +
    opts(title = "Loadings (before Supervised Norm.)\n") +
    xlab('\nSamples') +
    ylab('Eigengene 3 Loadings\n')
  
  combinedPlot <- multiplot(barPlot, eigenPlot1, eigenPlot2,
                            eigenPlot3, cols = 2)
  
  allFigs <- list('barPlot' = barPlot, 
                  'eigenPlot1' = eigenPlot1,
                  'eigenPlot2' = eigenPlot2,
                  'eigenPlot3' = eigenPlot3,
                  'combinedPlot' = combinedPlot)
}


# 
# multiplot(barFig, eigenFig1, eigenFig2, cols = 2)