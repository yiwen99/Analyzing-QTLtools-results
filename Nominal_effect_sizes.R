f_raw <- read.table("nominal_featureCounts_full.txt.gz")
#f_raw <- read.table("/u/project/zarlab/yiwen99/featureCounts_qtl/nominal_10/results/nominal_featureCounts_full.txt.gz")
k_raw <- read.table("nominal_kallisto_full.txt.gz")
#k_raw <- read.table("/u/project/zarlab/yiwen99/QTLtools_pipeline/nominal_10/results/nominal_kallisto_full.txt.gz")
f_raw <- na.omit(f_raw)
k_raw <- na.omit(k_raw)
#add a column of FDR values

#genename, 12pvalue, 13effects
f <- f_raw[,c(1,12,13)]
k <- k_raw[,c(1,12,13)]

common_gene <- intersect(f$V1,k$V1)

k_common<-k[which(k$V1%in%common_gene),]
f_common<-f[which(f$V1%in%common_gene),]
new <- cbind(k_common$V1, k_common$V12, f_common$V12,k_common$V13,f_common$V13)
#new$V1 common gene names
#new$V2 kallisto p-val
#new$V3 featureCounts p-val
#new$V4 kallisto effect size
#new$V5 featureCounts effect size
new <- as.data.frame(new)

#create a column of kallisto FDR values
new$V6 <- p.adjust(new$V2,method="BH")
#create a column of featureCounts FDR values
new$V7 <- p.adjust(new$V3,method="BH")

axis_max <- max(max(new$V4),max(new$V5))
axis_min <- min(min(new$V4),min(new$V5))

library(ggplot2)
pdf("Nominal Effect Sizes_k_FDR.pdf")
FDR_kallisto <- as.numeric(new$V6)
ggplot(data = new,aes(x=as.numeric(V4), y=as.numeric(V5)))+geom_point(aes(color=FDR_kallisto))+
  ggtitle("Nominal Effect Sizes")+
  xlab("featureCounts")+ylab("kallisto")+
  scale_x_continuous(breaks=seq(axis_min, axis_max, 0.5))+
  scale_y_continuous(breaks=seq(axis_min, axis_max, 0.5))+
  theme(axis.text.x = element_text(angle=90, hjust=1))

dev.off()

pdf("Nominal Effect Sizes_f_FDR.pdf")
FDR_featureCounts <- as.numeric(new$V7)
ggplot(data = new,aes(x=as.numeric(V4), y=as.numeric(V5)))+geom_point(aes(color=FDR_featureCounts))+
  ggtitle("Nominal Effect Sizes")+
  xlab("featureCounts")+ylab("kallisto")+
  scale_x_continuous(breaks=seq(axis_min, axis_max, 0.5))+
  scale_y_continuous(breaks=seq(axis_min, axis_max, 0.5))+
  theme(axis.text.x = element_text(angle=90, hjust=1))

dev.off()


