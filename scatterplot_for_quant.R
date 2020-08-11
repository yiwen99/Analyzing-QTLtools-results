#randomly choose 5 samples out of 87
#compare quantification data between 2 methods
#chromosome 1 chromosome 6 chromosome 19
intersection <- read.table("yri_sample_intersection.txt")
samples <- intersection$V1
chosen<-sample(samples,5)
#"NA19210" "NA18917" "NA19131" "NA19099" "NA18933"

bed_k_1 <- read.delim(paste0('kallisto_qqnorm','/qqnorm_',1,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
bed_f_1 <- read.delim(paste0('featureCounts_qqnorm','/qqnorm_',1,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
library(dplyr)
k_1 <- bed_k_1 %>% select("X.Chr", "gid", chosen)
f_1 <- bed_f_1 %>% select("X.Chr", "gid", chosen)
chrom_1 <- inner_join(f_1,k_1,by=c("X.Chr","gid"))

library(ggplot2)
pdf("NA19210 chrom1.pdf")
ggplot(data=chrom_1, aes(x=NA19210.x, y= NA19210.y))+
  geom_point(size=0.8)+
  ggtitle("Gene Quantification(scaled TPM) \nNA19210 chromosome1")+
  xlab("featureCounts")+ylab("kallisto")+
  geom_text(data=chrom_1[which(abs(chrom_1$NA19210.x-chrom_1$NA19210.y)>=2),],
            label=chrom_1[which(abs(chrom_1$NA19210.x-chrom_1$NA19210.y)>=2),]$gid,
            size=2,color="red3")
dev.off()

pdf("NA18917 chrom1.pdf")
ggplot(data=chrom_1, aes(x=NA18917.x, y= NA18917.y))+
  geom_point(size=0.8)+
  ggtitle("Gene Quantification(scaled TPM) \nNA18917 chromosome1")+
  xlab("featureCounts")+ylab("kallisto")+
  geom_text(data=chrom_1[which(abs(chrom_1$NA18917.x-chrom_1$NA18917.y)>=2),],
            label=chrom_1[which(abs(chrom_1$NA18917.x-chrom_1$NA18917.y)>=2),]$gid,
            size=2,color="red3")
dev.off()

pdf("NA19131 chrom1.pdf")
ggplot(data=chrom_1, aes(x=NA19131.x, y= NA19131.y))+
  geom_point(size=0.8)+
  ggtitle("Gene Quantification(scaled TPM) \nNA19131 chromosome1")+
  xlab("featureCounts")+ylab("kallisto")+
  geom_text(data=chrom_1[which(abs(chrom_1$NA19131.x-chrom_1$NA19131.y)>=2),],
            label=chrom_1[which(abs(chrom_1$NA19131.x-chrom_1$NA19131.y)>=2),]$gid,
            size=2,color="red3")
dev.off()

pdf("NA19099 chrom1.pdf")
ggplot(data=chrom_1, aes(x=NA19099.x, y= NA19099.y))+
  geom_point(size=0.8)+
  ggtitle("Gene Quantification(scaled TPM) \nNA19099 chromosome1")+
  xlab("featureCounts")+ylab("kallisto")+
  geom_text(data=chrom_1[which(abs(chrom_1$NA19099.x-chrom_1$NA19099.y)>=2),],
            label=chrom_1[which(abs(chrom_1$NA19099.x-chrom_1$NA19099.y)>=2),]$gid,
            size=2,color="red3")
dev.off()

pdf("NA18933 chrom1.pdf")
ggplot(data=chrom_1, aes(x=NA18933.x, y= NA18933.y))+
  geom_point(size=0.8)+
  ggtitle("Gene Quantification(scaled TPM) \nNA18933 chromosome1")+
  xlab("featureCounts")+ylab("kallisto")+
  geom_text(data=chrom_1[which(abs(chrom_1$NA18933.x-chrom_1$NA18933.y)>=2),],
            label=chrom_1[which(abs(chrom_1$NA18933.x-chrom_1$NA18933.y)>=2),]$gid,
            size=2,color="red3")
dev.off()
