# Ryosuke Sugiyama et al., PNAS
# Fig 1c
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(ggplot2)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5g/"

dat_all <- read.table(paste(data_dir, "R_Fig 5g.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$genotype <- str_replace(dat_all$genotype, " ", "")
dat_all$batch <- factor(dat_all$batch)

genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

for(x in unique(dat_all$batch)){
  idx <- dat_all$batch == x
  dat <- dat_all[idx,]

  # total
  pdf(paste(data_dir, x, "-data_nomral.pdf", sep=""))
  	hist(dat$value, breaks=20)
  	qqnorm(dat$value)
  	qqline(dat$value, col='red')

  	hist(log10(dat$value), breaks=20)
  	qqnorm(log10(dat$value))
  	qqline(log10(dat$value), col='red')

    hist(sqrt(dat$value), breaks=20)
    qqnorm(sqrt(dat$value))
    qqline(sqrt(dat$value), col='red')
  dev.off()

  lmer_fit <- lmer(log10(value) ~ sulfur:genotype  + sample- 1 +  (1|pot), dat)
  dat_lmer <- summary(lmer_fit)

  sink(paste(data_dir, x, "-lmer_summary.txt", sep=""))
  print(dat_lmer)
  sink()

  pdf(paste(data_dir, x, "-diagnosis.pdf", sep=""))
  	plot(fitted(lmer_fit), resid(lmer_fit))
  	abline(0, 0, col="red")

  	qqnorm(dat_lmer$resid)
  	qqline(dat_lmer$resid, col="red")
  dev.off()


  estim <- dat_lmer$coefficients[,1]
  df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )
  vcov <- as.matrix(dat_lmer$vcov)


  idx <- names(estim) != "sample"
  estim <- estim[idx]

  idx <- order(estim)
  estim <- estim[idx]

  idx <- match(names(estim), rownames(vcov))
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
  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

  sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))

  write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)
}
   