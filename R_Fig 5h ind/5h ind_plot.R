# Ryosuke Sugiyama et al., PNAS
# Fig S10b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)
library(ggplot2)
library(patchwork)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")


data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5h ind/"
GL_names <- read.table(paste(data_dir, "../GL names.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
mean     <- read.table(paste(data_dir, "logFC.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

# some stylish adjustments for plotting
split <- as.data.frame(matrix(unlist(str_split(mean$group, ":")), ncol=2, byrow=T))
names(split) <- c("sulfur", "genotype")
mean <- cbind(mean, split)

mean$batch <- factor(mean$batch)

genotype <- c("wt", "bglu28bglu30")
geno_label <- c(
	expression(`Col-0`),
	expression(italic("bglu28"~"bglu30")))
mean$genotype <- factor(mean$genotype, levels=rev(genotype))

sulfur <- c("S1500", "S30", "S3")
mean$sulfur <- factor(mean$sulfur, levels=sulfur)

idx <- match(mean$cpd, GL_names$code)
mean$cpd   <- factor(GL_names$abbr[idx],  levels=GL_names$abbr)
mean$class <- factor(GL_names$class[idx], levels=unique(GL_names$class))

min_q <- round(quantile(mean$logFC, .01), digits=1)
mean$logFC_fill <- mean$logFC
mean$logFC_fill[mean$logFC < min_q] <- min_q

# significance 
for(i in 1:(nrow(mean))){
	if(i %% 2 == 1){
		mean$sig[i] <- 0
	} else {
		mean$sig[i] <- sum(unlist(str_split(mean$letters[i], "")) %in% unlist(str_split(mean$letters[i-1], ""))) == 0
	}
}

# plotting
p1 <- ggplot(mean, aes(x=cpd:batch, y=genotype, fill=logFC_fill)) +
	geom_tile(aes(size=factor(sig)), colour="black") +
	geom_text(aes(label=letters), size=2, colour=c_black) +
	facet_grid(sulfur ~ class, switch="both", space="free", drop=T, scales="free") +
	scale_size_manual(values=c("0"=0.25, "1"=1), guide=F) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, breaks=c(2, 0, min_q/2, min_q), labels=c(2, 0, min_q/2, paste("<", min_q, sep=""))) +
	scale_y_discrete(breaks=rev(genotype), labels=rev(geno_label)) +
	theme_RTN +
	theme(
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		strip.text.x=element_blank()) +
	labs(x="Glucosinolate species", y="")
# ggsave(p1, file=paste(data_dir, "heatmap.pdf", sep=""), width=10, height=8, bg="transparent")

p2 <- ggplot(mean, aes(x=genotype, y=logFC, group=genotype, shape=batch)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=rev(genotype), labels=rev(geno_label)) +
	facet_grid(sulfur ~ .) +
	theme_RTN +
	theme(axis.text.y=element_blank(),
		strip.background=element_blank(),
		strip.text.y=element_blank(),
		legend.position="none") +
	labs(y="mean logFC", x="") +
	coord_flip()
# ggsave(p2, file=paste(data_dir, "FC_boxplot.pdf", sep=""), width=4, height=8, bg="transparent")

p <- ( p1 | p2 ) + plot_layout(widths = c(6, 1))
ggsave(p, file=paste(data_dir, "logFC_composite.pdf", sep=""), width=12, height=3.5, bg="transparent")





# transposed
mean$genotype <- factor(mean$genotype, levels=genotype)
p1 <- ggplot(mean, aes(x=genotype, y=cpd:batch, fill=logFC_fill)) +
	geom_tile(aes(size=factor(sig)), colour="black") +
	geom_text(aes(label=letters), size=2, colour=c_black) +
	facet_grid(class ~ sulfur, switch="both", space="free", drop=T, scales="free") +
	scale_size_manual(values=c("0"=0.25, "1"=1), guide=F) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, breaks=c(2, 0, min_q/2, min_q), labels=c(2, 0, min_q/2, paste("<", min_q, sep=""))) +
	scale_x_discrete(breaks=genotype, labels=geno_label) +
	theme_RTN +
	theme(
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		strip.text=element_blank(),
		legend.position="bottom") +
	labs(y="Glucosinolate species", x="")
# ggsave(p1, file=paste(data_dir, "heatmap.pdf", sep=""), width=10, height=8, bg="transparent")


p2 <- ggplot(mean, aes(x=genotype, y=logFC, group=genotype, shape=batch)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=genotype, labels=geno_label) +
	facet_grid(. ~ sulfur) +
	theme_RTN +
	theme(axis.text.x=element_blank(),
		strip.background=element_blank(),
		legend.position="none") +
	labs(y="mean logFC", x="")
# ggsave(p2, file=paste(data_dir, "FC_boxplot.pdf", sep=""), width=4, height=8, bg="transparent")

p <- ( p2 / p1 ) + plot_layout(height = c(1, 6))
ggsave(p, file=paste(data_dir, "logFC_composite.t.pdf", sep=""), width=3, height=12, bg="transparent")







# stats
mean$FC <- 2^mean$logFC
pdf(paste(data_dir, "logFC-data_nomral.pdf", sep=""))
	hist(mean$FC, breaks=20)
	qqnorm(mean$FC)
	qqline(mean$FC, col='red')

    hist(sqrt(mean$FC), breaks=20)
    qqnorm(sqrt(mean$FC))
    qqline(sqrt(mean$FC), col='red')

    hist(log10(mean$FC), breaks=20)
    qqnorm(log10(mean$FC))
    qqline(log10(mean$FC), col='red')
dev.off()

lmer_fit <- lmer(FC ~ sulfur:genotype - 1 + (1|batch), mean)
dat_lmer <- summary(lmer_fit)

pdf(paste(data_dir, "logFC-diagnosis.pdf", sep=""))
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
p.letters_df <- as.data.frame(p.letters$Letters[sort])
write.table(p.letters_df, file=paste(data_dir,"FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)


