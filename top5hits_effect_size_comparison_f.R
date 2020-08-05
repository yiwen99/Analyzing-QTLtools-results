quant_method <- "featureCounts_qqnorm"
bed <- list()
for (i in 1:22){
  #bed[[i]] <- read.delim(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
  bed[[i]] <- read.delim(paste0(quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
}

library(data.table)
vcf <- fread('genotypes_yri.vcf.gz') 

intersection <- read.table("yri_sample_intersection.txt")
samples <- intersection$V1

k_raw <- read.table("permutation_kallisto_scaled_tpm.all.txt.gz")
f_raw <- read.table("permutation_featureCounts.all.txt.gz")
k <- na.omit(k_raw)
f <- na.omit(f_raw)

#find the top 5 hits (5 least FDR-values)
#build a column with FDR values

f$V20 <- p.adjust(f$V19,method="BH")

top5_f <- f[order(f$V20,decreasing=F)[1:5],]
top5_f <- top5_f[,c(1,8,17)]
effect_sizes_featureCounts <- top5_f$V17
sig_snp_f <- top5_f$V8
sig_gene_f <- top5_f$V1
#"ENSG00000164308" "ENSG00000075234" "ENSG00000051620" "ENSG00000198468" "ENSG00000254184"

top5_f_in_k <- k[which(k$V1 %in% sig_gene_f),]
top5_f_in_k

#effect_sizes_kallisto <- top5_f_in_k$V17
effect_sizes_kallisto <- c(-1.103960,0.968221,-1.147010,1.381500,1.330890)

for_plot <- data.frame(cbind(sig_gene_f,effect_sizes_featureCounts,effect_sizes_kallisto))

library(ggplot2)
pdf("top 5 hits from featureCounts permutation pass.pdf")
ggplot(for_plot, aes(x=as.numeric(effect_sizes_featureCounts), y=as.numeric(effect_sizes_kallisto))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_plot$sig_gene_f, 
    nudge_x = -0.25, nudge_y = 0.05, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits from featureCounts")

dev.off()
