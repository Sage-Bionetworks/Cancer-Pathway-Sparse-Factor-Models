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
  require(ggplot2)
  histFig <- qplot(sigObj$pval, geom = 'histogram') + 
    opts(title = 'p-value Distribution of Transcripts\n') +
    xlab('\np-values') +
    ylab('Experiment Effect on Expression Variation\n')
  return(histFig)
}

## Second figure function
## Eigengene plots on data prior to supervised normalization
pcPlotsFig2 <- function(svdObj){
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


# 
# multiplot(barFig, eigenFig1, eigenFig2, cols = 2)