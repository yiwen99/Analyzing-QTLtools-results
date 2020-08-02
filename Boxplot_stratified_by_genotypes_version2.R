library(dplyr)
library(ggplot2)


quant_method = 'kallisto_qqnorm'
#quant_method = "featureCounts"

bed <- list()
for (i in 1:22){
  #bed[[i]] <- read.delim(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
  bed[[i]] <- read.delim(paste0(quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
}


library(data.table)
vcf <- fread('genotypes_yri.vcf.gz')

intersection <- read.table("yri_sample_intersection.txt")
samples <- intersection$V1


box_genotype_full <- function(snp){
  vcf_row <- which(vcf$ID==snp)
  chrom_num <- vcf$`#CHROM`[vcf_row]
  
  nominal_kallisto <- read.table(paste0("nominal_kallisto_",chrom_num,".txt"),sep=" ",stringsAsFactors=FALSE, quote="")

  permutation_kallisto <- read.table(paste0("permutation_kallisto_",chrom_num,".txt"),sep=" ",stringsAsFactors=FALSE, quote="")
  
  nominal_kallisto <- na.omit(nominal_kallisto)
  permutation_kallisto <-na.omit(permutation_kallisto)
  
  output <- left_join(nominal_kallisto, permutation_kallisto,by="V1") 
  output_needed <- output[,c(1,2,8)]
  
  
  bed_file <- bed[[chrom_num]] #bed[[chromosome#]]
  snp_row <- which(output_needed$V8.x==snp)
  gene <- output_needed[snp_row,1] #the gene whose expression level we want to look for
  
  #dosage
  dosage_row <- vcf_row
  dosage_samples <- vcf %>% select(samples)
  dosage <- round(dosage_samples[dosage_row,])
  
  #expression level
  level_row <- which(bed_file$gid==gene) 
  levels <- bed_file[level_row,c(7:93)]
  
  for_boxplot <- as.data.frame(t(rbind(dosage,levels)))
  
  for_boxplot$V1 <- as.factor(for_boxplot$V1)
  my.bp <- ggplot(data=for_boxplot, aes(y= V2, x=V1 ) )+geom_boxplot()+
    ggtitle(snp)+ 
    ylab("gene expression level") + xlab("snp dosage")+
    geom_smooth(aes(group=1),method="lm")
  my.bp # displays the boxplots
  
}
