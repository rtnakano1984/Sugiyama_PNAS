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

dat <- read.table(paste(data_dir, "R_Fig 5g.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$genotype <- str_replace(dat$genotype, " ", "")
dat$batch <- factor(dat$batch)

genotype <- c("wt", "bglu28bglu30")
dat$genotype <- factor(dat$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat$sulfur <- factor(dat$sulfur, levels=sulfur)



# normalize
for(x in unique(dat$batch)){
  idx <- dat$batch == x
  norm_fac <- mean(dat$value[idx & dat$sulfur == "S1500" & dat$genotype == "wt"])
  dat$rel[idx] <- dat$value[idx]/norm_fac
}


# total
pdf(paste(data_dir, "data_nomral_rel.pdf", sep=""))
	hist(dat$rel, breaks=20)
	qqnorm(dat$rel)
	qqline(dat$rel, col='red')

	hist(log10(dat$rel), breaks=20)
	qqnorm(log10(dat$rel))
	qqline(log10(dat$rel), col='red')

  hist(sqrt(dat$rel), breaks=20)
  qqnorm(sqrt(dat$rel))
  qqline(sqrt(dat$rel), col='red')
dev.off()

lmer_fit <- lmer(sqrt(rel) ~ sulfur:genotype - 1 + (1|pot:sample) + (1|batch), dat)
dat_lmer <- summary(lmer_fit)

sink(paste(data_dir, "lmer_summary_rel.txt", sep=""))
print(dat_lmer)
sink()

pdf(paste(data_dir, "diagnosis_rel.pdf", sep=""))
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
names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))

write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters_rel.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

 
#
p <- ggplot(dat, aes(x=genotype, y=rel, group=genotype, shape=factor(pot), color=factor(batch))) +
  geom_boxplot(outlier.shape=NA, position=position_dodge()) +
  geom_point(, position=position_jitterdodge()) +
  facet_grid(. ~ sulfur)
ggsave(p, file=paste(data_dir, "boxplot_rel.pdf", sep=""), width=5, height=4, bg="transparent")
