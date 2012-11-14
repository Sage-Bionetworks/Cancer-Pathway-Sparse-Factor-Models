## qcFunctions.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washingtonb
## erich.huang@sagebase.org

## Calculate the total sum of squares of the data and remove the biological
## treatment effect
removeExpEffect <- function(exprDat, X){
  datC <- sweep(exprDat, 1, rowMeans(exprDat))
  tSSQdat <- sum(datC^2)
  residuals <- exprDat - t(X %*% solve(t(X) %*% X) %*% t(X) %*% t(exprDat))
  rSSQdat <- sum(residuals^2)
  u <- fs(residuals)
  propSSQ <- round(rSSQdat * u$d, 3) / tSSQdat
}

## Use 'null probes' to build a dependence kernel and renormalize the data
nullProbeNorm <- function(sigObj, pcs, expEntity){
  nullProbes <- which(rank(1 - sigObj$pval) < 
    (length(sigObj$pval) * sigObj$pi0))
  u <- fs(exprDat[nullProbes, ])
  Z <- model.matrix(~ u$v[ , 1:pcs])
  print('Normalizing out remaining latent structure')
  fits2 <- runWorkflow(expEntity$cacheDir,
                       workflow = "snm", 
                       bio.var = X, 
                       adj.var = Z, 
                       rm.adj = TRUE)
  dat2 <- exprs(fits2[[1]][[1]])
  u2 <- fs(dat2)
  nullNorm <- list('uMatrix' = u2, 'newFit' = fits2)
}
