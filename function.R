# R function to return P values corresponding to pairwise t test after correction for multiple comparison by Benjamini-Hochberg method, based on a lmer fitted object

pairwise_ttest <- function(dat_lmer){
  estim <- dat_lmer$coefficients[,1]
  df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )
  vcov <- as.matrix(dat_lmer$vcov)

  idx <- order(estim)
  estim <- estim[idx]
  vcov <- vcov[idx,idx]

  group <- names(estim)

  group_comp <- c()    						
  for (i in 1:(length(group)-1)){							
    for (j in (i+1):length(group)){
      group_comp <- c(group_comp, paste(group[i], group[j], sep="-"))
    }
  }

  id.mat <- matrix(0, ncol=1, nrow = length(group) )
  p.val <- c()
  for ( i in 1:(length(estim)-1)) {
    for (j in (i+1):length(estim)){
      id.mat.x <- id.mat
      id.mat.x[ i, 1 ] <- 1
      id.mat.x[ j, 1 ] <- -1
      stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
      t.val <- abs( estim[i]-estim[j]) / stder
      p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
    }
  }

  names(p.val) <- group_comp

  adj_p.val <- p.adjust(p.val, "fdr")

  return(adj_p.val)
}
