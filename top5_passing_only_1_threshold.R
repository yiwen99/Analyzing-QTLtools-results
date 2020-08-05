k_raw <- read.table("permutation_kallisto_scaled_tpm.all.txt.gz")
f_raw <- read.table("permutation_featureCounts.all.txt.gz")
k <- na.omit(k_raw)
f <- na.omit(f_raw)

k$V20 <- p.adjust(k$V19,method="BH")
f$V20 <- p.adjust(f$V19,method="BH")

pass_k <- k[which(k$V20<=0.1),]
not_pass_k <- k[which(k$V20>0.1),]
pass_f <- f[which(f$V20<=0.1),]
not_pass_f <- f[which(f$V20>0.1),]  

intersect_knotf_gene <- intersect(pass_k$V1,not_pass_f$V1) #  ~400
intersect_fnotk_gene <- intersect(pass_f$V1,not_pass_k$V1) #  ~1000

#the part if passing featuerCounts FDR threshold but not passing kallisto
f_not_k <- pass_f[which(pass_f$V1 %in% intersect_fnotk_gene),]
top5_f_not_k <- f_not_k[order(f_not_k$V20,decreasing=F)[1:5],]
f_sig_only_gene <- top5_f_not_k$V1
#"ENSG00000181126" "ENSG00000170632" "ENSG00000176681" "ENSG00000269501" "ENSG00000251287"
feffect_f_only <- top5_f_not_k$V17
#1.466150 -0.825901 -0.897273 -1.106160 -1.036260
f_not_k_in_k <- not_pass_k[which(not_pass_k$V1 %in% f_sig_only_gene),]
f_not_k_in_k$V17
#-0.869195 -0.927947  0.411740 -0.373029 -0.341884
keffect_f_only <- c(0.411740, -0.373029, -0.341884, -0.869195, -0.927947)
for_f_only <- data.frame(cbind(f_sig_only_gene, keffect_f_only, feffect_f_only))
library(ggplot2)
pdf("top 5 hits only passing featureCounts threshold.pdf")
ggplot(for_f_only, aes(x=as.numeric(feffect_f_only), y=as.numeric(feffect_k_only))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_f_only$f_sig_only_gene, 
    nudge_x = 0.25, nudge_y = 0.1, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits only passing featureCounts threshold")

dev.off()





#the part if passing kallisto FDR threshold but not passing featureCounts
k_not_f <- pass_k[which(pass_k$V1 %in% intersect_knotf_gene),]
top5_k_not_f <- k_not_f[order(k_not_f$V20,decreasing=F)[1:5],]
k_sig_only_gene <- top5_k_not_f$V1
#"ENSG00000233927" "ENSG00000168894" "ENSG00000176261" "ENSG00000266967" "ENSG00000170291"
keffect_k_only <- top5_k_not_f$V17
#-1.475170 -0.604513 -0.656313  0.664288 -1.020830
k_not_f_in_f <- not_pass_f[which(not_pass_f$V1 %in% k_sig_only_gene),]
k_not_f_in_f$V17
#0.384086  0.237240 -0.168156 -0.575646 -0.302949
feffect_k_only <- c(-0.302949, 0.237240, 0.384086, -0.575646, -0.168156)
for_k_only <-data.frame(cbind(k_sig_only_gene, keffect_k_only,feffect_k_only))
pdf("top 5 hits only passing kallisto threshold.pdf")
ggplot(for_k_only, aes(x=as.numeric(feffect_k_only), y=as.numeric(keffect_k_only))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_k_only$k_sig_only_gene, 
    nudge_x = 0.25, nudge_y = 0.1, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits only passing kallisto threshold")

dev.off()



