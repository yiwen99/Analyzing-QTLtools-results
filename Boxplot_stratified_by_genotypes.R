#a function that will take a specific snp 
#and make a figure with 3 box plots varying the dosage
#plotting the standardize expression that is used for QTLtools (qqnorm_*.*)

#.bed.gz file includes the phenotype files-expression levels
#.vcf.gz file includes the dosage
#find the pair-up information between snp and gene expression level

#snp is the ID column of .vcf.gz file
#dosage of the snp are in the .vcf.gz file, under the sample columns
#expression level are 

#BASE = '/u/project/zarlab/hjp/geuvadis_data'
quant_method = 'kallisto_qqnorm'
#quant_method = "featureCounts"

#bednames = Sys.glob(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'))

bed <- list()
for (i in 1:22){
  #bed[[i]] <- read.delim(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
  bed[[i]] <- read.delim(paste0(quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
}

library(data.table)
vcf <- fread('genotypes_yri.vcf.gz') #128 cols


#results file left join output from nominal pass and permutation pass
#results file
# 1. The phenotype ID                     *** gene id
# 2. The chromosome ID of the phenotype   ***chromosome id
# 3. The start position of the phenotype
# 4. The end position of the phenotype
# 5. The strand orientation of the phenotype
# 6. The total number of variants tested in cis
# 7. The distance between the phenotype and the tested variant (accounting for strand orientation)
# 8. The ID of the top variant            *** snp id
# 9. The chromosome ID of the top variant 

nominal21_kallisto <- read.table("nominal_kallisto_21.txt",sep=" ",stringsAsFactors=FALSE, quote="")
#nominal21_featureCounts <- read.table("nominal_featureCounts_21.txt",sep=" ",stringsAsFactors=FALSE, quote="")
permutation21_kallisto <- read.table("permutation_kallisto_21.txt",sep=" ",stringsAsFactors=FALSE, quote="")

nominal21_kallisto <- na.omit(nominal21_kallisto)
permutation21_kallisto <-na.omit(permutation21_kallisto)

library(dplyr)
#nominal21_k <- nominal21_kallisto[,c(1,2,8)] #1048915 rows
#permutation21_k <-permutation21_kallisto[,c(1,2,8)] #176 rows
output <- left_join(nominal21_kallisto, permutation21_kallisto,by="V1") #1048915 rows
output_needed <- output[,c(1,2,8)]
#col1 gene id           col2 chromosome        col3 snp id

#snp <- "rs2821672"
#which(output_needed$V8.x==sp) #row

intersection <- read.table("yri_sample_intersection.txt")
samples <- intersection$V1

library(ggplot2)
box_genotype <- function(snp){
  bed <- bed[[21]] #bed[[chromosome#]]
  #pair the snp ID with gene ID using the results file, also get chromosome #
  #find the dosage(s) of that snp in vcf file of each sample, round dosage to 0,1,2
  #find the gene expression level associated with that snp in bed[[#chrom]] file of each sample
  snp_row <- which(output_needed$V8.x==snp)
  gene <- output_needed[snp_row,1] #the gene whose expression level we want to look for
  
  #dosage
  dosage_row <- which(vcf$ID==snp)
  dosage_samples <- vcf %>% select(samples)
  dosage <- round(dosage_samples[dosage_row,])
  
  #expression level
  level_row <- which(bed$gid==gene) 
  levels <- bed[level_row,c(7:93)]
  
  for_boxplot <- as.data.frame(t(rbind(dosage,levels)))
  
  
  #boxplot(for_boxplot$V2~ for_boxplot$V1, main=sp, xlab="snp dosage", ylab="expression level", horizontal=FALSE) 
  
  for_boxplot$V1 <- as.factor(for_boxplot$V1)
  my.bp <- ggplot(data=for_boxplot, aes(y= V2, x=V1 ) )+geom_boxplot()+
    ggtitle(snp)+ 
    ylab("gene expression level") + xlab("snp dosage")+
    geom_smooth(aes(group=1),method="lm")
  my.bp # displays the boxplots
  
  #input a snp from chromosome 21
  #find the corresponding gene in the result file, the line with that snp, find geneid
  #find dosage in vcf: that line with chromosome# 21 and snpid==snp (87)
  #find gene expression level in bed[[21]]  (87)
  #build dataframe:1col of dosage, 1col of expression, rownames being the sample names
}
