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
f_sig_only_snp <- top5_f_not_k$V8
#"rs112731253" "rs1404213"   "rs182191737" "rs2814294"   "rs4610271"  
#chrom 6        7             17            1             3
feffect_f_only <- top5_f_not_k$V17
#1.466150 -0.825901 -0.897273 -1.106160 -1.036260

##########on putty
#fn6 <- read.table("nominal_6.txt")
#k_effect_1 <- fn6[which(fn6$V8== "rs112731253"),]
#k_effect_1[which(k_effect_1$V1 =="ENSG00000181126"),]$V13
# 0.335983


#fn7 <- read.table("nominal_7.txt")
#k_effect_2 <- fn7[which(fn7$V8=="rs1404213"),]
#k_effect_2[which(k_effect_2$V1 =="ENSG00000170632"),]$V13
#-0.113692

#fn17 <- read.table("nominal_17.txt")
#k_effect_3 <- fn17[which(fn17$V8=="rs182191737"),]
#k_effect_3[which(k_effect_3$V1 =="ENSG00000176681"),]$V13
#-0.333166

#fn1 <- read.table("nominal_1.txt")
#k_effect_4 <- fn1[which(fn1$V8=="rs2814294"),]
#k_effect_4[which(k_effect_4$V1 =="ENSG00000269501"),]$V13
# -0.481215

#fn3 <- read.table("nominal_3.txt")
#k_effect_5 <- fn3[which(fn3$V8=="rs4610271"),]
#k_effect_5[which(k_effect_5$V1 =="ENSG00000251287"),]$V13
#-0.55163


keffect_f_only <- c(0.335983,-0.113692,-0.333166, -0.481215,-0.55163)

# f_not_k_in_k <- not_pass_k[which(not_pass_k$V1 %in% f_sig_only_gene),]
# f_not_k_in_k$V17
# #-0.869195 -0.927947  0.411740 -0.373029 -0.341884
# keffect_f_only <- c(0.411740, -0.373029, -0.341884, -0.869195, -0.927947)
for_f_only <- data.frame(cbind(f_sig_only_snp,f_sig_only_gene, keffect_f_only, feffect_f_only))
library(ggplot2)
pdf("top 5 hits only passing featureCounts threshold.pdf")
ggplot(for_f_only, aes(x=as.numeric(feffect_f_only), y=as.numeric(keffect_f_only))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_f_only$f_sig_only_snp, 
    nudge_x = -0.25, nudge_y = 0, 
    check_overlap = F
  )+
  geom_text(
    size=2,
    label=for_f_only$f_sig_only_gene, 
    nudge_x = -0.75, nudge_y = 0, 
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
k_sig_only_snp <- top5_k_not_f$V8
#"rs732217"  "rs2232747" "rs704878"  "rs7219730" "rs4562"  
#chrom 19      2          1           17          17
keffect_k_only <- top5_k_not_f$V17
#-1.475170 -0.604513 -0.656313  0.664288 -1.020830

#k_not_f_in_f <- not_pass_f[which(not_pass_f$V1 %in% k_sig_only_gene),]
#k_not_f_in_f$V17

##########on putty
#fn19 <- read.table("nominal_19.txt")
#f_effect_1 <- fn19[which(fn19$V8== "rs732217"),]
#f_effect_1[which(f_effect_1$V1 =="ENSG00000233927"),]$V13
#-0.301818

#fn2 <- read.table("nominal_2.txt")
#f_effect_2 <- fn2[which(fn2$V8=="rs2232747"),]
#f_effect_2[which(f_effect_2$V1 =="ENSG00000168894"),]$V13
#0.0356306

#fn1 <- read.table("nominal_1.txt")
#f_effect_3 <- fn1[which(fn1$V8=="rs704878"),]
#f_effect_3[which(f_effect_3$V1 =="ENSG00000176261"),]$V13
# -0.0995719

#fn17 <- read.table("nominal_17.txt")
#f_effect_4 <- fn17[which(fn17$V8=="rs7219730"),]
#f_effect_4[which(f_effect_4$V1 =="ENSG00000266967"),]$V13
# 0.126137


#f_effect_5 <- fn17[which(fn17$V8=="rs4562"),]
#f_effect_5[which(f_effect_5$V17 =="ENSG00000170291"),]$V13
#NA ?????



feffect_k_only <- c(-0.301818, 0.0356306, -0.0995719, 0.126137, NA)
#feffect_k_only <- c(-0.302949, 0.237240, 0.384086, -0.575646, -0.168156)

for_k_only <-data.frame(cbind(k_sig_only_gene,k_sig_only_snp, keffect_k_only,feffect_k_only))
library(ggplot2)
pdf("top 5 hits only passing kallisto threshold(one pair not found).pdf")
ggplot(for_k_only, aes(x=as.numeric(feffect_k_only), y=as.numeric(keffect_k_only))) +
  geom_point()+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_text(
    size=2,
    label=for_k_only$k_sig_only_snp, 
    nudge_x = -0.25, nudge_y = 0.05, 
    check_overlap = F
  )+
  geom_text(
    size=2,
    label=for_k_only$k_sig_only_gene, 
    nudge_x = -0.75, nudge_y = 0.05, 
    check_overlap = F
  )+
  labs(y= "effect sizes from kallisto",x="effect sizes from featureCounts")+
  ggtitle("Effect Sizes for Top 5 Hits only passing kallisto threshold\n rs4562-ENSG00000170291 not fount in featureCounts")

dev.off()
