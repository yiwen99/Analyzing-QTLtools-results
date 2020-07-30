# library('dplyr')
# 
# BASE = '/u/home/h/hjp/project-zarlab/qtltest/results/qtltools'
# 
# quant_method = 'kallisto_scaled_tpm'
# # quant_method = 'featureCounts'
# fnames = Sys.glob(paste0(BASE, '/epc_*/', quant_method, '/permutation_*.all.txt.gz'))
# npc = sub('/.*$', '', sub('^.*epc_', '', fnames))
# 
# # 1. The phenotype ID
# # 2. The chromosome ID of the phenotype
# # 3. The start position of the phenotype
# # 4. The end position of the phenotype
# # 5. The strand orientation of the phenotype
# # 6. The total number of variants tested in cis
# # 7. The distance between the phenotype and the tested variant (accounting for strand orientation)
# # 8. The ID of the top variant
# # 9. The chromosome ID of the top variant
# # 10. The start position of the top variant
# # 11. The end position of the top variant
# # 12. The number of degrees of freedom used to compute the P-values
# # 13. Dummy
# # 14. The first parameter value of the fitted beta distribution
# # 15. The second parameter value of the fitted beta distribution (it also gives the effective number of independent tests in the region)
# # 16. The nominal P-value of association between the phenotype and the top variant in cis
# # 17. The corresponding regression slope
# # 18. The P-value of association adjusted for the number of variants tested in cis given by the direct method (i.e. empirircal P-value)
# # 19. The P-value of association adjusted for the number of variants tested in cis given by the fitted beta distribution. We strongly recommend to use this adjusted P-value in any downstream analysis
# 
# all_permutations = lapply(seq_along(fnames),
#                           function(i) {
#                             ret = read.table(fnames[i], stringsAsFactors = FALSE, header = FALSE)
#                             colnames(ret) = c('pid', 'chr', 'start', 'end', 'strand', 'n_snps',
#                                               'distance', 'best_snp', 'snp_chr', 'snp_start', 'snp_end', 'df', 'dummy',
#                                               'beta1', 'beta-2', 'p_nominal', 'effect_size', 'p_empirical', 'p_beta')
#                             dplyr::mutate(ret, n_epc = npc[i], fdr = p.adjust(p_beta, method = 'BH'))
#                           })
# all_permutations = bind_rows(all_permutations)
# 
# summarize(group_by(all_permutations, n_epc), n_fdr10 = sum(fdr < 0.10, na.rm = TRUE))
##The above code are provided by Professor Pimentel. It gives the number of eQTLs discovered with each pc's used.

kallisto <- read.delim("kallisto_scaled_tpm.txt")
featureCounts <- read.delim("featureCounts.txt")

#colnames(kallisto) <- c("k_n_epc","k_n_fdr10")
method <- c("kallisto","featureCounts")
k <- cbind(kallisto,method=method[1])
f <- cbind(featureCounts,method=method[2])
all <- rbind(k,f)

library(ggplot2)

#plot for kallisto
pdf("kallisto-pc-qtl.pdf")
ggplot(data = kallisto,aes(x=n_epc, y=n_fdr10))+geom_point()+ggtitle("kallisto")
dev.off()

#plot for featureCounts
pdf("featureCounts-pc-qtl.pdf")
ggplot(data = featureCounts,aes(x=n_epc, y=n_fdr10))+geom_point()+ggtitle("featureCounts")
dev.off()

pdf("bothMethod.pdf")
ggplot(data=all,aes(x=n_epc, y=n_fdr10))+geom_point(aes(color=method))+
  scale_color_manual(values=c("orange","blue"))

dev.off() #detach the pdf
