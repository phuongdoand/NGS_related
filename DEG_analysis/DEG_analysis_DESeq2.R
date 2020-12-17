setwd('Documents/translational genomics/cryocord/')
library(RColorBrewer)
library(Rtsne)
library(DESeq2)
library(gplots)
if(!require("affy")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")  
}
library(affy)
library(sva)

# Load data and rename the column
data = read.delim("mirna_ordered.txt")
row.names(data) <- make.unique(as.character(data[,1]))
data = data[,-1]
group = c('A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B', 'B', 
           'C', 'C', 'C', 'C', 'C', 'C', 'D', 'D', 'D', 'D', 'D', 'D')
colors = colorRampPalette(c("red", "green", "blue", "magenta"))(length(unique(as.factor(group))))[factor(as.factor(group))]
colnames(data) = c('0619_A_1', '0619_A', '3008_A', '3008_A_1', '3678_A', '3678_A_1', '0619_B_1', '0619_B',
                   '3008_B_1', '3008_B', '3678_B', '3678_B_1', '0619_C', '0619_C_1', '3008_C', '3008_C_1',
                   '3678_C_1', '3678_C', '0619_D_1', '0619_D', '3008_D', '3008_D_1', '3678_D_1', '3678_D')

# Clear the one with no useful data, i.e null data
sum(rowSums(data) > 0) # 796
sum(rowSums(data) == 0) # 2087
new_data = data[rowSums(data) > 0,]

# Draw some histograms, boxplot, and density plot to view the data
par(mfrow=c(2,1))
hist(as.matrix(new_data), col="blue", border="white", breaks=100)
hist(as.matrix(log2(new_data + 1)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

par(mfrow=c(1,1))
boxplot(log2(new_data + 1), col=colors, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(Counts +1)")

plotDensity(log2(new_data + 1), lty=1, col=colors, lwd=2)
grid()
legend("topright", legend=unique(group), col=unique(colors), lwd=2)

plotFun <- function(x,y){ 
  dns <- densCols(x,y); 
  points(x,y, col=dns, pch=".", panel.first=grid());  
  #  abline(a=0, b=1, col="brown")
}
pairs(log2(new_data[,sample(ncol(new_data), 6)] + 1), 
      panel=plotFun, lower.panel = NULL)

# Draw null percentage in all the samples
prop.null <- apply(data, 2, function(x) 100*mean(x==0))
tiff(filename = 'null_percentage.tiff')
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=colors, ylab='Samples', xlab='% of null counts')
dev.off()

# Draw MDS plot of the raw data log transformed
fit <- cmdscale(dist(t(log2(new_data+1))))
tiff(filename = "MDS_raw_log.tiff", width = 480, height = 480)
plot(fit, main = 'MDS raw data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-35, 35), ylim = c(-20,30), type = 'n')
for(i in 1:24){
  text(fit[i,1], fit[i,2], labels = colnames(new_data)[i], col = colors[i])
}
dev.off()

# Drawing t-SNE plot of the raw data log transformed
labels = as.factor(group)
tnse_raw = Rtsne(t(log2(new_data+1)), dims = 2, perplexity = 5, verbose = TRUE,
                 max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_raw_log.tiff", width = 480, height = 480)
plot(tnse_raw$Y, main="tSNE raw data", col=colors,  pch=18, xlab = 't-SNE1', 
     ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_raw$Y[i,1], tnse_raw$Y[i,2], labels = colnames(new_data)[i], col = colors[i])
}
#dev.off()

# Load data to DESeq2 and estimate the normalize factors
dds = DESeqDataSetFromMatrix(as.matrix(new_data), data.frame(group), 
                             design = ~ group)
dds = estimateSizeFactors(dds)
# Get the normalized counts using the normalize method of DESeq2
counts = as.data.frame(counts(dds, normalized = TRUE))

# Draw MDS plot of normalized data log transformation
normed_fit = cmdscale(dist(t(log2(x = counts + 1))))
tiff(filename = "MDS_norm_log.tiff", width = 480, height = 480)
plot(normed_fit, main = 'MDS norm data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-35, 35), ylim = c(-30, 30), type = 'n')
for(i in 1:24){
  text(normed_fit[i,1], normed_fit[i,2], labels = colnames(new_data)[i], col = colors[i])
}
dev.off()

# Draw t-SNE plot of normalized data log transform
tnse_norm = Rtsne(t(log2(counts+1)), dims = 2, perplexity = 5, verbose = TRUE,
                  max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_norm_log.tiff", width = 480, height = 480)
plot(tnse_norm$Y, main="tSNE norm data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-35, 45), ylim = c(-35,45), ylab='t-SNE2', type='n')
for(i in 1 : 24) {
  text(tnse_norm$Y[i,1], tnse_norm$Y[i,2], labels = colnames(new_data)[i], col = colors[i])
}
#dev.off()

# Draw boxplot of raw and normed data
tiff(filename = 'boxplot_raw.tiff')
boxplot(log2(new_data+1),  col=colors, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
dev.off()
tiff(filename = 'boxplot_norm.tiff')
boxplot(log2(counts(dds, normalized=TRUE)+1),  col=colors, cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
dev.off()

# Draw density distribution before and after normalization
tiff(filename = 'density_raw.tiff')
plotDensity(log2(counts(dds)+1),  col=colors, lty=1, lwd=2, xlab="log2(counts)", 
            main="Raw read count distribution", cex.lab=0.7, panel.first=grid()) 
legend("topright", legend=unique(group), col=unique(colors), lwd=2)
dev.off()
tiff(filename = 'density_norm.tiff')
plotDensity(log2(counts(dds, normalized=TRUE)+1), col=colors, lty=1, lwd=2,
            main="Norm read count distribution",
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
legend("topright", legend=unique(group), col=unique(colors), lwd=2)
dev.off()

# log transform data 
rld = rlogTransformation(dds, blind = TRUE)

# Calculate some statistics
norm.counts = counts(dds, normalized=TRUE)
mean.counts = rowMeans(norm.counts)
variance.counts = apply(norm.counts, 1, var)

norm.counts.stats <- data.frame(
  min=apply(norm.counts, 2, min),
  mean=apply(norm.counts, 2, mean),
  median=apply(norm.counts, 2, median),
  max=apply(norm.counts, 2, max),
  zeros=apply(norm.counts==0, 2, sum),
  percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
  perc05=apply(norm.counts, 2, quantile, 0.05),
  perc10=apply(norm.counts, 2, quantile, 0.10),
  perc90=apply(norm.counts, 2, quantile, 0.90),
  perc95=apply(norm.counts, 2, quantile, 0.95)
)

# Draw the mean and variance distribution of the data
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
tiff(filename = 'mean_variance_cor.tiff')
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")
dev.off()

# Estimate the dispersions for statistic test of DEG 
dds = estimateDispersions(dds)
tiff(filename = "dispersion.tiff", width = 480, height = 480)
plotDispEsts(dds)
dev.off()

# Perform the statistic test
dds = nbinomWaldTest(dds)
dds = DESeqDataSetFromMatrix(as.matrix(new_data), data.frame(group), 
                             design = ~ group)
dds2 = DESeq(dds)

# Extract result for each groups with a = 0.05
res1 = results(dds2, contrast = c("group","A","B"), alpha = 0.05)
res2 = results(dds2, contrast = c("group","A","C"), alpha = 0.05)
res3 = results(dds2, contrast = c("group","A","D"), alpha = 0.05)
res4 = results(dds2, contrast = c("group","B","C"), alpha = 0.05)
res5 = results(dds2, contrast = c("group","B","D"), alpha = 0.05)
res6 = results(dds2, contrast = c("group","C","D"), alpha = 0.05)

## Extract results for each group with a = 0.1 (default)
# res1 = results(dds2, contrast = c("group","A","B"))
# res2 = results(dds2, contrast = c("group","A","C"))
# res3 = results(dds2, contrast = c("group","A","D"))
# res4 = results(dds2, contrast = c("group","B","C"))
# res5 = results(dds2, contrast = c("group","B","D"))
# res6 = results(dds2, contrast = c("group","C","D"))

# Draw p-value histogram 
hist(res1$pvalue[res1$baseMean > 1], col = "grey50",main = "A vs B", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res2$pvalue[res2$baseMean > 1], col = "grey50",main = "A vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res3$pvalue[res3$baseMean > 1], col = "grey50",main = "A vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res4$pvalue[res4$baseMean > 1], col = "grey50",main = "B vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res5$pvalue[res5$baseMean > 1], col = "grey50",main = "B vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res6$pvalue[res6$baseMean > 1], col = "grey50",main = "C vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")

# Show the table to see how many miRNA passed the test
table(res1$padj<0.05)
table(res2$padj<0.05)
table(res3$padj<0.05)
table(res4$padj<0.05)
table(res5$padj<0.05)
table(res6$padj<0.05)

# table(res1$padj<0.1)
# table(res2$padj<0.1)
# table(res3$padj<0.1)
# table(res4$padj<0.1)
# table(res5$padj<0.1)
# table(res6$padj<0.1)

# Order the result by sorting the p-adjusted value
res1 <- res1[order(res1$padj),]
res2 <- res2[order(res2$padj),]
res3 <- res3[order(res3$padj),]
res4 <- res4[order(res4$padj),]
res5 <- res5[order(res5$padj),]
res6 <- res6[order(res6$padj),]

# Create the subset with the miRNA that passed the test (aka significantly different)
deg1 <- subset(res1,padj < 0.05)
deg2 <- subset(res2,padj < 0.05)
deg3 <- subset(res3,padj < 0.05)
deg4 <- subset(res4,padj < 0.05)
deg5 <- subset(res5,padj < 0.05)
deg6 <- subset(res6,padj < 0.05)

# deg1 <- subset(res1,padj < 0.1)
# deg2 <- subset(res2,padj < 0.1)
# deg3 <- subset(res3,padj < 0.1)
# deg4 <- subset(res4,padj < 0.1)
# deg5 <- subset(res5,padj < 0.1)
# deg6 <- subset(res6,padj < 0.1)

table(deg1$log2FoldChange>1)
table(deg1$log2FoldChange<(-1))

## Write data to csv files for further analysis
# write.csv(as.data.frame(deg1[deg1$log2FoldChange>1,]), file = "Aup_Bdown.csv")
# write.csv(as.data.frame(deg1[deg1$log2FoldChange<(-1),]), file = "Adown_Bup.csv")
# write.csv(as.data.frame(deg2[deg2$log2FoldChange>1,]), file = "Aup_Cdown.csv")
# write.csv(as.data.frame(deg2[deg2$log2FoldChange<(-1),]), file = "Adown_Cup.csv")
# write.csv(as.data.frame(deg3[deg3$log2FoldChange>1,]), file = "Aup_Ddown.csv")
# write.csv(as.data.frame(deg3[deg3$log2FoldChange<(-1),]), file = "Adown_Dup.csv")
# write.csv(as.data.frame(deg4[deg4$log2FoldChange>1,]), file = "Bup_Cdown.csv")
# write.csv(as.data.frame(deg4[deg4$log2FoldChange<(-1),]), file = "Bdown_Cup.csv")
# write.csv(as.data.frame(deg5[deg5$log2FoldChange>1,]), file = "Bup_Ddown.csv")
# write.csv(as.data.frame(deg5[deg5$log2FoldChange<(-1),]), file = "Bdown_Dup.csv")
# write.csv(as.data.frame(deg6[deg6$log2FoldChange>1,]), file = "Cup_Ddown.csv")
# write.csv(as.data.frame(deg6[deg6$log2FoldChange<(-1),]), file = "Cdown_Dup.csv")

write.csv(as.data.frame(deg1), file = "AvsB.csv")
write.csv(as.data.frame(deg2), file = "AvsC.csv")
write.csv(as.data.frame(deg3), file = "AvsD.csv")
write.csv(as.data.frame(deg4), file = "BvsC.csv")
write.csv(as.data.frame(deg5), file = "BvsD.csv")
write.csv(as.data.frame(deg6), file = "CvsD.csv")


## Trying ComBat
pheno = data.frame("sample" = colnames(data), 
                   "group"=c(rep("A", 6), rep("B", 6), rep("C", 6), rep("D", 6)), 
                   "batch"=c(2,1,1,2,1,2,2,1,2,1,1,2,1,2,1,2,2,1,2,1,1,2,2,1))
batch = pheno$batch
modcombat = model.matrix(~as.factor(pheno$group), data=pheno)
combat_data = ComBat(dat = as.matrix(new_data), batch = batch, mod = modcombat, 
                     par.prior = FALSE)

combat_shift = combat_data - min(combat_data) # Shift the data so no more negative value

# Draw MDS plot for combat data
fit_combat = cmdscale(dist(t(log2(combat_shift+1))))
tiff(filename = "MDS_combat_raw_log.tiff", width = 480, height = 480)
plot(fit_combat, main = 'MDS Combat raw data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-15, 20), ylim = c(-15,15), type = 'n')
for(i in 1:24){
  text(fit_combat[i,1], fit_combat[i,2], labels = colnames(combat_shift)[i], col = colors[i])
}
dev.off()

# Draw t-SNE plot for the combat data
tnse_raw_combat = Rtsne(t(log2(combat_shift+1)), dims = 2, perplexity = 5, verbose = TRUE,
                        max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_raw_log.tiff", width = 480, height = 480)
plot(tnse_raw_combat$Y, main="tSNE Combat raw data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-40, 40), ylim = c(-55,55), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_raw_combat$Y[i,1], tnse_raw_combat$Y[i,2], labels = colnames(combat_data)[i], col = colors[i])
}
#dev.off()

# Use the combat result for DEG
dds_combat = DESeqDataSetFromMatrix(as.matrix(round(combat_shift)), 
                                    data.frame(group), design = ~ group)
dds_combat = estimateSizeFactors(dds_combat)

counts_norm_combat = as.data.frame(counts(dds_combat, normalized = TRUE))

plotDensity(log2(counts(dds_combat)+1),  col=colors, lty=1, lwd=2, xlab="log2(counts)", 
            main="Combat raw read count distribution", cex.lab=0.7, panel.first=grid()) 
legend("topright", legend=unique(group), col=unique(colors), lwd=2)
plotDensity(log2(counts(dds_combat, normalized=TRUE)+1), col=colors, lty=1, lwd=2,
            main="Combat norm read count distribution",
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
legend("topright", legend=unique(group), col=unique(colors), lwd=2)

combat_data_normed = ComBat(dat = as.matrix(counts(dds, normalized = TRUE)), batch = batch, mod = modcombat, 
                     par.prior = FALSE)
combat_shift_normed = combat_data_normed - min(combat_data_normed) 

fit_norm_combat = cmdscale(dist(t(log2(x = combat_shift_normed + 1))))
tiff(filename = "t-MDS_combat_norm_log.tiff", width = 480, height = 480)
plot(fit_norm_combat, main = 'MDS Norm-Combat data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-20, 20), ylim = c(-15,15), type = 'n')
for(i in 1:24){
  text(fit_norm_combat[i,1], fit_norm_combat[i,2], labels = colnames(counts_norm_combat)[i], col = colors[i])
}
dev.off()

#tnse_norm_combat = Rtsne(t(log2(combat_shift_normed+1)), dims = 2, perplexity = 5, verbose = TRUE,
#                        max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_norm_combat$Y, main="tSNE Norm-Combat data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-65, 55), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_norm_combat$Y[i,1], tnse_norm_combat$Y[i,2], labels = colnames(counts_norm_combat)[i], col = colors[i])
}
#dev.off()

dds_combat = estimateDispersions(dds_combat)
dds_combat = nbinomWaldTest(dds_combat)

res1_com = results(dds_combat, contrast = c("group","A","B"), alpha = 0.05)
res2_com = results(dds_combat, contrast = c("group","A","C"), alpha = 0.05)
res3_com = results(dds_combat, contrast = c("group","A","D"), alpha = 0.05)
res4_com = results(dds_combat, contrast = c("group","B","C"), alpha = 0.05)
res5_com = results(dds_combat, contrast = c("group","B","D"), alpha = 0.05)
res6_com = results(dds_combat, contrast = c("group","C","D"), alpha = 0.05)

hist(res1_com$pvalue[res1_com$baseMean > 1], col = "grey50",main = "A vs B", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res2_com$pvalue[res2_com$baseMean > 1], col = "grey50",main = "A vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res3_com$pvalue[res3_com$baseMean > 1], col = "grey50",main = "A vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res4_com$pvalue[res4_com$baseMean > 1], col = "grey50",main = "B vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res5_com$pvalue[res5_com$baseMean > 1], col = "grey50",main = "B vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res6_com$pvalue[res6_com$baseMean > 1], col = "grey50",main = "C vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")

table(res1_com$padj<0.05)
table(res2_com$padj<0.05)
table(res3_com$padj<0.05)
table(res4_com$padj<0.05)
table(res5_com$padj<0.05)
table(res6_com$padj<0.05)

###############################################
# Testing RPM
data_rpm = as.data.frame(apply(new_data, 2, function(x) x/(sum(x)/1000000)))
data_norm_rpm = as.data.frame(apply(counts, 2, function(x) x/(sum(x)/1000000)))

distance_rpm = cmdscale(dist(t(log2(data_rpm+1))))

#tiff(filename = "MDS_raw_log.tiff", width = 480, height = 480)
plot(distance_rpm, main = 'MDS RPM data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-35, 35), ylim = c(-30,30), type = 'n')
for(i in 1:24){
  text(distance_rpm[i,1], distance_rpm[i,2], labels = colnames(data_rpm)[i], col = colors[i])
}
#dev.off()

#tnse_rpm = Rtsne(t(log2(data_rpm+1)), dims = 2, perplexity = 5, verbose = TRUE,
#                         max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_rpm$Y, main="tSNE RPM data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-55, 50), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_rpm$Y[i,1], tnse_rpm$Y[i,2], labels = colnames(data_rpm)[i], col = colors[i])
}
#dev.off()

# Using RPM data for DeSeq2
# Raw => RPM => DeSeq2
dds_rpm = DESeqDataSetFromMatrix(as.matrix(round(data_rpm)), 
                                 as.data.frame(group), design = ~group)
dds_rpm = DESeq(dds_rpm)
norm_rpm = counts(dds_rpm, normalized=TRUE)

plotDensity(log2(counts(dds_rpm)+1),  col=colors, lty=1, lwd=2, xlab="log2(counts)", 
            main="RPM read count distribution", cex.lab=0.7, panel.first=grid()) 
legend("topright", legend=unique(group), col=unique(colors), lwd=2)
plotDensity(log2(counts(dds_rpm, normalized=TRUE)+1), col=colors, lty=1, lwd=2,
            main="RPM norm read count distribution",
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
legend("topright", legend=unique(group), col=unique(colors), lwd=2)

distance_norm_rpm = cmdscale(dist(t(log2(norm_rpm+1))))
plot(distance_norm_rpm, main = 'MDS RPM norm data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-35, 35), ylim = c(-30,30), type = 'n')
for(i in 1:24){
  text(distance_norm_rpm[i,1], distance_norm_rpm[i,2], labels = colnames(data_norm_rpm)[i], col = colors[i])
}

#tnse_rpm_norm = Rtsne(t(log2(norm_rpm+1)), dims = 2, perplexity = 5, verbose = TRUE,
#                         max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_rpm_norm$Y, main="tSNE RPM norm data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-55, 50), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_rpm_norm$Y[i,1], tnse_rpm_norm$Y[i,2], labels = colnames(norm_rpm)[i], col = colors[i])
}
#dev.off()

res1_rpm = results(dds_rpm, contrast = c("group","A","B"), alpha = 0.05)
res2_rpm = results(dds_rpm, contrast = c("group","A","C"), alpha = 0.05)
res3_rpm = results(dds_rpm, contrast = c("group","A","D"), alpha = 0.05)
res4_rpm = results(dds_rpm, contrast = c("group","B","C"), alpha = 0.05)
res5_rpm = results(dds_rpm, contrast = c("group","B","D"), alpha = 0.05)
res6_rpm = results(dds_rpm, contrast = c("group","C","D"), alpha = 0.05)

hist(res1_rpm$pvalue[res1_rpm$baseMean > 1], col = "grey50",main = "A vs B", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res2_rpm$pvalue[res2_rpm$baseMean > 1], col = "grey50",main = "A vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res3_rpm$pvalue[res3_rpm$baseMean > 1], col = "grey50",main = "A vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res4_rpm$pvalue[res4_rpm$baseMean > 1], col = "grey50",main = "B vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res5_rpm$pvalue[res5_rpm$baseMean > 1], col = "grey50",main = "B vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res6_rpm$pvalue[res6_rpm$baseMean > 1], col = "grey50",main = "C vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")

table(res1_rpm$padj<0.05)
table(res2_rpm$padj<0.05)
table(res3_rpm$padj<0.05)
table(res4_rpm$padj<0.05)
table(res5_rpm$padj<0.05)
table(res6_rpm$padj<0.05)

####### 
pheno = data.frame("sample" = colnames(data), 
                   "group"=c(rep("A", 6), rep("B", 6), rep("C", 6), rep("D", 6)), 
                   "batch"=c(2,1,1,2,1,2,2,1,2,1,1,2,1,2,1,2,2,1,2,1,1,2,2,1))
batch = pheno$batch
modcombat = model.matrix(~as.factor(pheno$group), data=pheno)
combat_data_rpm = ComBat(dat = as.matrix(data_rpm), batch = batch, mod = modcombat, 
                     par.prior = FALSE)

combat_rpm_shift = combat_data_rpm - min(combat_data_rpm)

distance_combat_rpm = cmdscale(dist(t(log2(combat_rpm_shift+1))))
#tiff(filename = "MDS_raw_log.tiff", width = 480, height = 480)
plot(distance_combat_rpm, main = 'MDS RPM-Combat data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-20, 20), ylim = c(-15,15), type = 'n')
for(i in 1:24){
  text(distance_combat_rpm[i,1], distance_combat_rpm[i,2], labels = colnames(combat_data_rpm)[i], col = colors[i])
}
#dev.off()

#tnse_rpm_combat = Rtsne(t(log2(combat_rpm_shift+1)), dims = 2, perplexity = 5, verbose = TRUE,
#                         max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_rpm_combat$Y, main="tSNE RPM-Combat data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-55, 40), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_rpm_combat$Y[i,1], tnse_rpm_combat$Y[i,2], labels = colnames(combat_rpm_shift)[i], col = colors[i])
}
#dev.off()

## RPM-Norm-Combat
combat_rpm_norm = ComBat(dat = as.matrix(round(norm_rpm[rowSums(norm_rpm)>1,])), batch = batch, mod = modcombat, 
                         par.prior = FALSE)
combat_rpm_norm_shift = combat_rpm_norm - min(combat_rpm_norm)

distance_combat_rpm_norm = cmdscale(dist(t(log2(combat_rpm_norm_shift+1))))
#tiff(filename = "MDS_raw_log.tiff", width = 480, height = 480)
plot(distance_combat_rpm_norm, main = 'MDS RPM-Norm-Combat data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-20, 20), ylim = c(-15,15), type = 'n')
for(i in 1:24){
  text(distance_combat_rpm_norm[i,1], distance_combat_rpm_norm[i,2], labels = colnames(combat_rpm_norm_shift)[i], col = colors[i])
}
#dev.off()

tnse_rpm_combat_norm = Rtsne(t(log2(combat_rpm_norm_shift+1)), dims = 2, perplexity = 5, verbose = TRUE,
                         max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_rpm_combat_norm$Y, main="tSNE RPM-Norm-Combat data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-55, 50), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {ro
  text(tnse_rpm_combat_norm$Y[i,1], tnse_rpm_combat_norm$Y[i,2], labels = colnames(combat_rpm_norm_shift)[i], col = colors[i])
}
#dev.off()

#### DeSeq2 Wrong !!! 
dds_rpm_combat = DESeqDataSetFromMatrix(as.matrix(round(combat_rpm_shift)), 
                                        as.data.frame(group), design = ~group)
dds_rpm_combat = DESeq(dds_rpm_combat)
rpm_combat_norm = counts(dds_rpm_combat, normalized=TRUE)

plotDensity(log2(counts(dds_rpm_combat)+1),  col=colors, lty=1, lwd=2, xlab="log2(counts)", 
            main="RPM-Combat read count distribution", cex.lab=0.7, panel.first=grid()) 
legend("topright", legend=unique(group), col=unique(colors), lwd=2)
plotDensity(log2(counts(dds_rpm_combat, normalized=TRUE)+1), col=colors, lty=1, lwd=2,
            main="RPM-Combat norm read count distribution",
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
legend("topright", legend=unique(group), col=unique(colors), lwd=2)

distance_combat_rpm_norm = cmdscale(dist(t(log2(rpm_combat_norm+1))))
#tiff(filename = "MDS_raw_log.tiff", width = 480, height = 480)
plot(distance_combat_rpm_norm, main = 'MDS RPM-Combat norm data', xlab = 'MDS1', ylab = 'MDS2', 
     xlim = c(-20, 20), ylim = c(-15,15), type = 'n')
for(i in 1:24){
  text(distance_combat_rpm_norm[i,1], distance_combat_rpm_norm[i,2], labels = colnames(rpm_combat_norm)[i], col = colors[i])
}
#dev.off()

#tnse_rpm_combat_norm = Rtsne(t(log2(rpm_combat_norm+1)), dims = 2, perplexity = 5, verbose = TRUE,
#                         max_iter = 2000, check_duplicates = FALSE, theta=0.5)
#tiff(filename = "t-SNE_combat_norm_log.tiff", width = 480, height = 480)
plot(tnse_rpm_combat_norm$Y, main="tSNE RPM-Combat norm data", col=colors,  pch=18, xlab = 't-SNE1', 
     xlim = c(-55, 50), ylim = c(-55,65), ylab='t-SNE2', type = "n")
for(i in 1 : 24) {
  text(tnse_rpm_combat_norm$Y[i,1], tnse_rpm_combat_norm$Y[i,2], labels = colnames(rpm_combat_norm)[i], col = colors[i])
}
#dev.off()

res1_rpm_combat = results(dds_rpm_combat, contrast = c("group","A","B"), alpha = 0.05)
res2_rpm_combat = results(dds_rpm_combat, contrast = c("group","A","C"), alpha = 0.05)
res3_rpm_combat = results(dds_rpm_combat, contrast = c("group","A","D"), alpha = 0.05)
res4_rpm_combat = results(dds_rpm_combat, contrast = c("group","B","C"), alpha = 0.05)
res5_rpm_combat = results(dds_rpm_combat, contrast = c("group","B","D"), alpha = 0.05)
res6_rpm_combat = results(dds_rpm_combat, contrast = c("group","C","D"), alpha = 0.05)

hist(res1_rpm_combat$pvalue[res1_rpm_combat$baseMean > 1], col = "grey50",main = "A vs B", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res2_rpm_combat$pvalue[res2_rpm_combat$baseMean > 1], col = "grey50",main = "A vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res3_rpm_combat$pvalue[res3_rpm_combat$baseMean > 1], col = "grey50",main = "A vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res4_rpm_combat$pvalue[res4_rpm_combat$baseMean > 1], col = "grey50",main = "B vs C", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res5_rpm_combat$pvalue[res5_rpm_combat$baseMean > 1], col = "grey50",main = "B vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")
hist(res6_rpm_combat$pvalue[res6_rpm_combat$baseMean > 1], col = "grey50",main = "C vs D", 
     xlab = "p-values", breaks = 0:20/20, border = "white")

table(res1_rpm_combat$padj<0.05)
table(res2_rpm_combat$padj<0.05)
table(res3_rpm_combat$padj<0.05)
table(res4_rpm_combat$padj<0.05)
table(res5_rpm_combat$padj<0.05)
table(res6_rpm_combat$padj<0.05)



