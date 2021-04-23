# Ryosuke Sugiyama et al., PNAS
# Fig 2b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(psych)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 2b/"

dat <- read.table(paste(data_dir, "R_Fig 2b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$genotype <- str_replace(dat$genotype, " ", "")
dat$batch <- factor(dat$batch)

dat$sulfur <- str_replace(dat$sulfur, "4MSB-34S", "4MSB_34S")

for(cpd in unique(dat$cpd)){
  idx <- dat$cpd == cpd
  temp <- dat[idx,]

  pdf(paste(data_dir, cpd, "-data_nomral.pdf", sep=""))
    hist(dat$value, breaks=20)
    qqnorm(dat$value)
    qqline(dat$value, col='red')

    hist(log10(dat$value), breaks=20)
    qqnorm(log10(dat$value))
    qqline(log10(dat$value), col='red')
  dev.off()

  # lmer_fit <- lmer(log10(value) ~ sulfur - 1 + (1|batch), dat)
  lmer_fit <- lmer(value ~ sulfur - 1 + (1|batch), dat)
  dat_lmer <- summary(lmer_fit)

  sink(paste(data_dir, cpd, "-lmer_summary.txt", sep=""))
  print(dat_lmer)
  sink()

  pdf(paste(data_dir, cpd, "-diagnosis.pdf", sep=""))
    plot(fitted(lmer_fit), resid(lmer_fit))
    abline(0, 0, col="red")

    qqnorm(dat_lmer$resid)
    qqline(dat_lmer$resid, col="red")
  dev.off()


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
  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
  # names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

  # sort <- apply(expand.grid(c("wt", "bglu28bglu30"), c("S1500", "S150", "S15", "S0"))[, 2:1], 1, function(x) paste(x, collapse=":"))
  # sort <- c("S1500", "S150", "S15", "S0", "4MSB")
  sort <- c("S1500", "S150", "4MSB", "4MSB_34S")

  write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}


