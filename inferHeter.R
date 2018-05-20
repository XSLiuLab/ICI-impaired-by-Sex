setwd("~/wangshx/projects/")
library(maftools)
load("data/TCGA_LUAD_Maf.RData")
seg <- "data/segmented_scna_minus_germline_cnv_hg19_name_chopoff.seg.txt"
luad_het <- inferHeterogeneity(maf=luad_maf, tsb = getSampleSummary(luad_maf)$Tumor_Sample_Barcode, segFile = seg)

luad_het_noCNV <- inferHeterogeneity(maf=luad_maf, tsb = getSampleSummary(luad_maf)$Tumor_Sample_Barcode)
#plotClusters(clusters = luad_het)

save(luad_het, file="data/inferHeter_remove_CNV.RData")
save(luad_het_noCNV, file="data/inferHeter.RData")
