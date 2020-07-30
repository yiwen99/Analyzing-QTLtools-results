#a function that will take a specific snp 
#and make a figure with 3 box plots varying the dosage
#plotting the standardize expression that is used for QTLtools (qqnorm_*.*)

#.bed.gz file includes the phenotype files-expression levels
#.vcf.gz file includes the dosage

#snp is the ID column of .vcf.gz file
#dosage of the snp are in the .vcf.gz file, under the sample columns
#expression level are 

#snp ------ sample dosage(need to round) --(find that sample in qqnorm file)
#-----expression level???
BASE = '/u/project/zarlab/hjp/geuvadis_data'
quant_method = 'kallisto_scaled_tpm'
#quant_method = "featureCounts"
bednames = Sys.glob(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'))

bed <- list()
for (i in 1:22){
  bed[[i]] <- read.delim(paste0(BASE,'/qtltools/',quant_method,'/qqnorm_',i,'.bed.gz'), sep="\t",stringsAsFactors=FALSE, quote="")
}

vcf <- fread(BASE,'/annotation/genotypes_yri.vcf.gz')

#results file left join output from nominal pass and permutation pass

#bed <- as.data.frame(read.delim("qqnorm_1.bed.gz", sep="\t",stringsAsFactors=FALSE, quote=""))

# 1. The phenotype ID *** gene id
# 2. The chromosome ID of the phenotype
# 3. The start position of the phenotype
# 4. The end position of the phenotype
# 5. The strand orientation of the phenotype
# 6. The total number of variants tested in cis
# 7. The distance between the phenotype and the tested variant (accounting for strand orientation)
# 8. The ID of the top variant *** snp id
# 9. The chromosome ID of the top variant ***chromosome id

box_genotype <- function(snp){
  #pair the snp ID with gene ID using the results file, also get chromosome #
  #find the dosage(s) of that snp in vcf file of each sample
  #find the gene expression level associated with that snp in bed[[#chrom]] file of each sample
  #save the information in a data frame with 1 col of snp dosage, 1 col of gene expresison levels
  #plot each dot for a sample with x axis as the dosage and y axis as the gene expression level
}
