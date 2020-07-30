#10pc data
f_raw <- read.table("permutation_featureCounts.all.txt.gz")
#f_raw <- read.table("nominal_featureCounts_full.txt.gz")
k_raw <- read.table("permutation_kallisto_scaled_tpm.all.txt.gz")
#k_raw <- read.table("nominal_kallisto_full.txt.gz")
f_raw <- na.omit(f_raw)
k_raw <- na.omit(k_raw)
#making histograms of p-values for both methods
# pdf("kallisto_p_val.pdf")
# hist(k_raw$V19,breaks=100,xlab="raw p-values",main="raw p-values(Permutation Pass for kallisto-10pcs)")

#pdf("featureCounts_p_val.pdf")
#hist(f_raw$V19,breaks=100,xlab="raw p-values",main="raw p-values(Permutation Pass for featureCounts-10pcs)")


f <- f_raw[,c(1,17,19)]  #15325 genes
k <- k_raw[,c(1,17,19)]  #16990 genes
common_gene <- intersect(f$V1,k$V1) #14816 genes

k_common<-k[which(k$V1%in%common_gene),]
f_common<-f[which(f$V1%in%common_gene),]
new <- cbind(k_common$V1, k_common$V17, f_common$V17,k_common$V19,f_common$V19)
new <- as.data.frame(new)
library(ggplot2)
#making a density scatterplot of the p-vals
pdf("p-val density scatterplot.pdf")
sp <- ggplot(new, aes(f_common$V19, y=k_common$V19)) +
  geom_point()
#sp + geom_density_2d()
sp+stat_density_2d(aes(fill = ..level..), geom="polygon")+
  xlab("featureCounts")+ylab("kallisto")+ggtitle("p-value Density Scatterplot")
