library(dplyr)
quant_method <- "kallisto_qqnorm"
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

k$V20 <- p.adjust(k$V19,method="BH")

top5_k <- k[order(k$V20,decreasing=F)[1:5],]
top5_k <- top5_k[,c(1,2,8,17)]
effect_sizes_kallisto <- top5_k$V17
sig_snp_k <- top5_k$V8
#"rs2608824"  "rs6928273"  "rs12023052" "rs2247650"  "rs1739652"
#chromosome 4 6x 1x 5 20x
#need those files from nominal pass

sig_gene_k <- top5_k$V1
#"ENSG00000163682" "ENSG00000051620" "ENSG00000198468" "ENSG00000164308" "ENSG00000196756"

top5_k_in_f_2 <- f[which(f$V8 == "rs6928273"),] #effect size -1.16801
top5_k_in_f_3 <- f[which(f$V8 == "rs12023052"),] #effect size 1.39044
top5_k_in_f_5 <- f[which(f$V8 == "rs1739652"),] #effect size -0.942852

##########on putty
#fn4 <- read.table("nominal_4.txt")
#top5_k_in_f_1 <- fn4[which(fn4$V8=="rs2608824"),]
#top5_k_in_f_1[which(top5_k_in_f_1$V1 =="ENSG00000163682"),]$V13
#-0.656625


#fn5 <- read.table("nominal_5.txt")
#top5_k_in_f_4 <- fn5[which(fn5$V8=="rs2247650"),]
#top5_k_in_f_4[which(top5_k_in_f_4$V1 =="ENSG00000164308"),]$V13
#-1.15049


#top5_k_in_f <- f[which(f$V1 %in% sig_gene_k),]
#top5_k_in_f
#"ENSG00000198468" "ENSG00000163682" "ENSG00000164308" "ENSG00000051620" "ENSG00000196756"
#"rs12023052"         "rs1012352"         "rs2548531"    "rs6928273"      "rs1739652" 

#effect_sizes_featureCounts <- top5_k_in_f$V17

#effect_sizes_featureCounts <- c(-0.652825,-1.168010,1.390440,-1.133510,-0.942852)
effect_sizes_featureCounts <- c( -0.656625, -1.168010, 1.390440, -1.15049, -0.942852)

for_plot_ksig <- data.frame(cbind(sig_snp_k,sig_gene_k,effect_sizes_kallisto,effect_sizes_featureCounts))


library(ggplot2)
pdf("top 5 hits from kallisto permutation pass.pdf")
ggplot(for_plot_ksig, aes(x=as.numeric(effect_sizes_featureCounts), y=as.numeric(effect_sizes_kallisto))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_plot_ksig$sig_snp_k, 
    nudge_x = -0.25, nudge_y = 0.05, 
    check_overlap = F
  )+
  geom_text(
    size=2,
    label=for_plot_ksig$sig_gene_k,
    nudge_x = -0.75, nudge_y = 0.05, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits from kallisto")

dev.off()
