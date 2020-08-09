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
top5_f <- top5_f[,c(1,2,8,17)]
effect_sizes_featureCounts <- top5_f$V17
sig_snp_f <- top5_f$V8 
#"rs2548531"         "rs57866552"       "rs6928273"        "rs12023052"   "rs56176828"
#chrom 5              22x                 6x                  1              7x
sig_gene_f <- top5_f$V1
#"ENSG00000164308" "ENSG00000075234" "ENSG00000051620" "ENSG00000198468" "ENSG00000254184"

#top5_f_in_k_test <- k[which(k$V8 %in% sig_gene_f)]
top5_f_in_k_5 <- k[which(k$V8 =="rs56176828"),]
top5_f_in_k_5
top5_f_in_k_2 <- k[which(k$V8 == "rs57866552"),]
top5_f_in_k_2
top5_f_in_k_3 <- k[which(k$V8 == "rs6928273"),]
top5_f_in_k_3
top5_f_in_k_4 <- k[which(k$V8 == "rs12023052"),]
top5_f_in_k_4



#on putty kallisto nominal_10
#fn5 <- read.table("nominal_5.txt")
#top5_k_in_f_1 <- fn5[which(fn5$V8=="rs2548531"),]
#top5_k_in_f_1[which(top5_k_in_f_1$V1 =="ENSG00000164308"),]$V13
#-1.08206




#effect_sizes_kallisto <- top5_f_in_k$V17
#effect_sizes_kallisto <- c(-1.103960,0.968221,-1.147010,1.381500,1.330890)
effect_sizes_kallisto <- c(-1.08206, 0.968221, -1.14701, 1.3815, 1.33089)

for_plot_fsig <- data.frame(cbind(sig_snp_f,sig_gene_f,effect_sizes_featureCounts,effect_sizes_kallisto))

library(ggplot2)
pdf("top 5 hits from featureCounts permutation pass.pdf")
ggplot(for_plot_fsig, aes(x=as.numeric(effect_sizes_featureCounts), y=as.numeric(effect_sizes_kallisto))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_plot_fsig$sig_snp_f, 
    nudge_x = -0.25, nudge_y = 0.05, 
    check_overlap = F
  )+
  geom_text(
    size=2,
    label=for_plot_fsig$sig_gene_f,
    nudge_x = -0.75, nudge_y = 0.05, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits from featureCounts")

dev.off()
