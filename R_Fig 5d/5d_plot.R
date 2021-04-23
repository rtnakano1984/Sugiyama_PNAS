# Ryosuke Sugiyama et al., PNAS
# Fig 5d (for plotting)
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)
library(ggplot2)
library(patchwork)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")


data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5d/"
GL_names <- read.table(paste(data_dir, "../GL names.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
mean     <- read.table(paste(data_dir, "logFC.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)


# some data rearrangement for plotting
split <- as.data.frame(matrix(unlist(str_split(mean$group, ":")), ncol=2, byrow=T))
names(split) <- c("time", "genotype")
mean <- cbind(mean, split)


genotype <- c("wt", "bglu28", "bglu30", "bglu28bglu30")
geno_label <- c(
	expression(`Col-0`),
	expression(italic("bglu28")),
	expression(italic("bglu30")),
	expression(italic("bglu28"~"bglu30")))
mean$genotype <- factor(mean$genotype, levels=rev(genotype))

time <- c("S1500_0d", "S1500_5d", "S0_5d")
mean$time <- factor(mean$time, levels=time)


idx <- match(mean$cpd, GL_names$code)
mean$cpd   <- factor(GL_names$abbr[idx],  levels=GL_names$abbr)
mean$class <- factor(GL_names$class[idx], levels=unique(GL_names$class))



min_q <- round(quantile(mean$logFC, .10), digits=1)
mean$logFC_fill <- mean$logFC
mean$logFC_fill[mean$logFC < min_q] <- min_q

# significance 
for(x in unique(mean$cpd)){
	for(y in unique(mean$time)){
		idx_0 <- mean$cpd == x & mean$time == y
		temp <- mean[idx_0,]
		sig <- 0
		for(i in 2:nrow(temp)){
			idx <- temp$genotype == "wt"
			sig <- c(sig, sum(unlist(str_split(temp$letters[i], "")) %in% unlist(str_split(temp$letters[idx], ""))) == 0)
		}
		mean$sig[idx_0] <- sig
	}
}

levels(mean$time) <- c("Day 0", "Day 5\nin S1500", "Day 5\nin S0")

# plotting
p1 <- ggplot(mean, aes(x=cpd, y=genotype, fill=logFC_fill)) +
	geom_tile(aes(size=factor(sig)), colour="black") +
	geom_text(aes(label=letters), size=2, colour=c_black) +
	facet_grid(time ~ class, switch="both", space="free", drop=T, scales="free") +
	scale_size_manual(values=c("0"=0.25, "1"=1), guide=F) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, breaks=c(2, 0, min_q/2, min_q), labels=c(2, 0, min_q/2, paste("<", min_q, sep=""))) +
	scale_y_discrete(breaks=rev(genotype), labels=rev(geno_label)) +
	theme_RTN +
	theme(
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		strip.text.x=element_blank()) +
	labs(x="Glucosinolate species", y="", fill="logFC (vs Col-0, Day 0 in S1500)")
# ggsave(p1, file=paste(data_dir, "heatmap.pdf", sep=""), width=10, height=8, bg="transparent")


p2 <- ggplot(mean, aes(x=genotype, y=logFC)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=rev(genotype), labels=rev(geno_label)) +
	scale_y_continuous(breaks=c(2, 0, -4, -8)) +
	facet_grid(time ~ .) +
	theme_RTN +
	theme(axis.text.y=element_blank(),
		strip.background=element_blank(),
		strip.text.y=element_blank()) +
	labs(y="mean logFC", x="") +
	coord_flip()
# ggsave(p2, file=paste(data_dir, "FC_boxplot.pdf", sep=""), width=4, height=8, bg="transparent")



p <- ( p1 | p2 ) + plot_layout(widths = c(6, 1))
ggsave(p, file=paste(data_dir, "logFC_composite.pdf", sep=""), width=9, height=5, bg="transparent")




# transposed
mean$genotype <- factor(mean$genotype, levels=genotype)
p1 <- ggplot(mean, aes(x=genotype, y=cpd, fill=logFC_fill)) +
	geom_tile(aes(size=factor(sig)), colour="black") +
	geom_text(aes(label=letters), size=2, colour=c_black) +
	facet_grid(class ~ time, switch="both", space="free", drop=T, scales="free") +
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


p2 <- ggplot(mean, aes(x=genotype, y=logFC)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=genotype, labels=geno_label) +
	facet_grid(. ~ time) +
	theme_RTN +
	theme(axis.text.x=element_blank(),
		strip.background=element_blank()) +
	labs(y="mean logFC", x="")
# ggsave(p2, file=paste(data_dir, "FC_boxplot.pdf", sep=""), width=4, height=8, bg="transparent")

p <- ( p2 / p1 ) + plot_layout(height = c(1, 6))
ggsave(p, file=paste(data_dir, "logFC_composite.t.pdf", sep=""), width=4.5, height=9, bg="transparent")





# some stats using mean logFC
sort <- apply(expand.grid(genotype, time)[, 2:1], 1, function(x) paste(x, collapse=":"))
mean$group <- factor(mean$group, levels=sort)
kruskal.test(logFC ~ group, mean)

mean_summary <- mean %>% group_by(group) %>% summarise(mean=mean(logFC))

idx <- order(mean_summary$mean)
group_ordered <- sort[idx]

mean$group <- factor(mean$group, levels=group_ordered)
wilcox <- pairwise.wilcox.test(mean$logFC, mean$group, "fdr")

p_val <- matrix(NA, ncol=length(sort), nrow=length(sort))
colnames(p_val) <- group_ordered
rownames(p_val) <- group_ordered

for(i in 1:(length(group_ordered)-1)){
	for(j in (i+1):length(group_ordered)){
		p_val[j,i] <- wilcox$p.value[j-1, i]
		p_val[i,j] <- wilcox$p.value[j-1, i]		
	}
}
p_letters <- multcompLetters(p_val)$Letters[sort]
write.table(as.data.frame(p_letters), file=paste(data_dir, "pairwise_wilcox-Letters.txt", sep=""), col.names=F, row.names=T, quote=F, sep="\t")


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

aov <- aov(logFC ~ time:genotype - 1, mean)
pdf(paste(data_dir, "logFC-diagnosis.pdf", sep=""))
	plot(fitted(aov), resid(aov))
	abline(0, 0, col="red")

	qqnorm(aov$resid)
	qqline(aov$resid, col="red")
dev.off()













