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


data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5h/"
GL_names <- read.table(paste(data_dir, "../GL names.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
mean     <- read.table(paste(data_dir, "logFC.txt", sep=""), sep="\t", header=T, stringsAsFactors=F)

# load stat data
stats_letters <- lapply(unique(mean$cpd), function(x){
	temp <- read.table(paste(data_dir, x, "-FDR_letters.txt", sep=""), header=T, stringsAsFactors=F, row.names=NULL)
	names(temp) <- c("group", "letters")
	temp$cpd <- x
	return(temp)
}) %>% do.call(rbind, .)

stats <- lapply(unique(mean$cpd), function(x){
	temp <- read.table(paste(data_dir, x, "-P_values.txt", sep=""), header=T, stringsAsFactors=F, row.names=1)
	rownames(temp) <- str_replace_all(rownames(temp), "sulfur", "")
	rownames(temp) <- str_replace_all(rownames(temp), "genotype", "")
	group_comp <- matrix(unlist(str_split(rownames(temp), ":|-")), ncol=4, byrow=T)
	colnames(group_comp) <- c("sulfur1", "genotype1", "sulfur2", "genotype2")
	temp <- data.frame(temp, group_comp, cpd=x)
	return(temp)
}) %>% do.call(rbind, .)

# some data rearrangement for plotting
split <- as.data.frame(matrix(unlist(str_split(mean$group, ":")), ncol=2, byrow=T))
names(split) <- c("sulfur", "genotype")
mean <- cbind(mean, split)

genotype <- c("wt", "bglu28bglu30")
geno_label <- c(
	expression(`Col-0`),
	expression(italic("bglu28"~"bglu30")))
mean$genotype <- factor(mean$genotype, levels=rev(genotype))

sulfur <- c("S1500", "S30", "S3")
mean$sulfur <- factor(mean$sulfur, levels=sulfur)

pdf(paste(data_dir, "hist_logFC.pdf", sep=""))
	hist(mean$logFC, breaks=20)
dev.off()

min_q <- round(quantile(mean$logFC, .01), digits=1)
mean$logFC_fill <- mean$logFC
mean$logFC_fill[mean$logFC < min_q] <- min_q

# merge with stats data (letters)
idx <- match(paste(mean$sulfur, mean$genotype, mean$cpd, sep=":"), paste(stats_letters$group, stats_letters$cpd, sep=":"))
mean <- data.frame(mean, letters=stats_letters$letters[idx])

# merge with stats data (p values)
idx <- stats$sulfur1 == stats$sulfur2
stats <- stats[idx,]

idx <- stats$genotype1 == "wt" | stats$genotype2 == "wt"
stats <- stats[idx,]

for(i in 1:nrow(mean)){
	if(mean$genotype[i] == "wt"){
		mean$sig[i] <- 0
	} else {
		idx <- stats$sulfur1 == mean$sulfur[i] & ( stats$genotype1 == mean$genotype[i] | stats$genotype2 == mean$genotype[i] ) & stats$cpd == mean$cpd[i]
		mean$sig[i] <- as.numeric(stats$adj_p[idx] < .05)
	}
}


# some stylish adjustments for plotting
idx <- match(mean$cpd, GL_names$code)
mean$cpd   <- factor(GL_names$abbr[idx],  levels=GL_names$abbr)
mean$class <- factor(GL_names$class[idx], levels=unique(GL_names$class))

# plotting
p1 <- ggplot(mean, aes(x=cpd, y=genotype, fill=logFC_fill)) +
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


p2 <- ggplot(mean, aes(x=genotype, y=logFC)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=rev(genotype), labels=rev(geno_label)) +
	facet_grid(sulfur ~ .) +
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


p2 <- ggplot(mean, aes(x=genotype, y=logFC)) +
	geom_hline(yintercept=0) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(colour=c_black, position=position_jitter(width=.2), size=.5) +
	scale_x_discrete(breaks=genotype, labels=geno_label) +
	facet_grid(. ~ sulfur) +
	theme_RTN +
	theme(axis.text.x=element_blank(),
		strip.background=element_blank()) +
	labs(y="mean logFC", x="")
# ggsave(p2, file=paste(data_dir, "FC_boxplot.pdf", sep=""), width=4, height=8, bg="transparent")

p <- ( p2 / p1 ) + plot_layout(height = c(1, 6))
ggsave(p, file=paste(data_dir, "logFC_composite.t.pdf", sep=""), width=3, height=9, bg="transparent")







# some stats using mean logFC
sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))
mean$group <- factor(paste(mean$sulfur, mean$genotype, sep=":"), levels=sort)
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

aov <- aov(logFC ~ sulfur:genotype - 1, mean)
pdf(paste(data_dir, "logFC-diagnosis.pdf", sep=""))
	plot(fitted(aov), resid(aov))
	abline(0, 0, col="red")

	qqnorm(aov$resid)
	qqline(aov$resid, col="red")
dev.off()





