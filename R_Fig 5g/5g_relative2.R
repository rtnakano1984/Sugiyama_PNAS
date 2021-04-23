# Ryosuke Sugiyama et al., PNAS
# Fig 1c
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(ggplot2)
library(multcompView)
library(reshape2)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5g/"

dat_all <- read.table(paste(data_dir, "R_Fig 5g.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

# normalize
dat <- lapply(unique(dat_all$batch), function(x){
  idx <- dat_all$batch == x
  temp <- dat_all[idx,] %>% group_by(sulfur, genotype, pot, batch) %>% summarise(value=median(value)) %>% data.frame(., stringsAsFactors=F)

  dcast <- dcast(temp, sulfur + pot + batch ~ genotype, value.var="value")
  dcast$value <- dcast$bglu28bglu30 / dcast$wt

  return(dcast)
}) %>% do.call(rbind, .)



# total
pdf(paste(data_dir, "data_nomral_rel2.pdf", sep=""))
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

lmer_fit <- lmer(log10(value) ~ sulfur - 1 + (1|pot:batch), dat)
dat_lmer <- summary(lmer_fit)

sink(paste(data_dir, "lmer_summary_rel2.txt", sep=""))
print(dat_lmer)
sink()

pdf(paste(data_dir, "diagnosis_rel2.pdf", sep=""))
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

write.table(as.data.frame(p.letters$Letters[sulfur]), file=paste(data_dir, "FDR_letters_rel2.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

 
#
p <- ggplot(dat, aes(x=sulfur, y=value, group=sulfur, color=factor(batch))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter())
ggsave(p, file=paste(data_dir, "boxplot_rel2.pdf", sep=""), width=5, height=4, bg="transparent")
