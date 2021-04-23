# Ryosuke Sugiyama et al., PNAS
# Fig S7a
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig5f/"
# data_dir <- "~/Desktop/ryosuke_temp/R_Fig4fgh_new/"

dat <- read.table(paste(data_dir, "R_Fig5f.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "bglu28 bglu30")
sulfur <- c("S200", "S100", "4MTB", "4MSB", "PhE")

idx <- dat$sulfur %in% sulfur
dat <- dat[idx,]

dat$genotype <- factor(dat$genotype, genotype)
dat$sulfur   <- factor(dat$sulfur, sulfur)

# pairwise_test <- function(dat_lmer){

# 	estim <- dat_lmer$coefficients[,1]
# 	df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )
# 	vcov <- as.matrix(dat_lmer$vcov)

# 	idx <- order(estim)
# 	estim <- estim[idx]
# 	vcov <- vcov[idx,idx]


# 	group <- names(estim)

# 	group_comp <- c()    						
# 	for (i in 1:(length(group)-1)){							
# 	  for (j in (i+1):length(group)){
# 	    group_comp <- c(group_comp, paste(group[i], group[j], sep="-"))
# 	  }
# 	}

# 	id.mat <- matrix(0, ncol=1, nrow = length(group) )
# 	p.val <- c()
# 	for ( i in 1:(length(estim)-1)) {
# 	  for (j in (i+1):length(estim)){
# 	    id.mat.x <- id.mat
# 	    id.mat.x[ i, 1 ] <- 1
# 	    id.mat.x[ j, 1 ] <- -1
# 	    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
# 	    t.val <- abs( estim[i]-estim[j]) / stder
# 	    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=f ) )
# 	  }
# 	}

# 	names(p.val) <- group_comp

# 	adj_p.val <- p.adjust(p.val, "fdr")
# 	p.letters <- multcompletters(adj_p.val)

# 	return(p.letters)
# }



pdf(paste(data_dir, "data_nomral.pdf", sep=""))
	hist(dat$value, breaks=20)
	qqnorm(dat$value)
	qqline(dat$value, col='red')

	hist(sqrt(dat$value), breaks=20)
	qqnorm(sqrt(dat$value))
	qqline(sqrt(dat$value), col='red')

	hist(log10(dat$value), breaks=20)
	qqnorm(log10(dat$value))
	qqline(log10(dat$value), col='red')
dev.off()

fit <- lm(value ~ 0 + genotype:sulfur + batch, dat)
anova <- anova(aov(fit))
# lmer_fit <- lmer(sqrt(value) ~ genotype:sulfur - 1 + (1|batch), dat)
# dat_lmer <- summary(lmer_fit)

sink(paste(data_dir, "lmer_summary.txt", sep=""))
	print(anova)
sink()

pdf(paste(data_dir, "diagnosis.pdf", sep=""))
	plot(fitted(fit), resid(fit))
	abline(0, 0, col="red")

	qqnorm(summary(fit)$resid)
	qqline(summary(fit)$resid, col="red")
dev.off()

tukey <- TukeyHSD(aov(fit))
idx <- order(-tukey$'genotype:sulfur'[, 'diff'])
p.letters <- multcompLetters(tukey$'genotype:sulfur'[idx, 'p adj'])

sort <- expand.grid(genotype, sulfur) %>% apply(1, paste, collapse=":")
write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

p <- ggplot(dat, aes(x=sulfur, y=value, colour=genotype)) +
	geom_point(position=position_jitterdodge(jitter.width=.2)) +
	geom_boxplot(fill=NA, outlier.shape=NA) +
	theme_bw()
ggsave(p, file=paste(data_dir, "boxplot.pdf", sep=""), width=6, height=4, bg="transparent")

# p.letters <- pairwise_test(dat_lmer)
# sort <- expand.grid(paste("genotype", genotype, sep=""), paste("sulfur", sulfur, sep="")) %>% apply(1, paste, collapse=":")
# write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)


