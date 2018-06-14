library(maftools)
library(tidyverse)
library(data.table)

load("../TCGA_LUAD_Maf.RData")
#load("../TCGA_LUAD_clin.RData")

# plot MAF summary
plotmafSummary(maf = luad_maf, rmOutlier = FALSE, addStat = "median", dashboard = TRUE, titvRaw = FALSE,
               top = 20, fs=.9)

# oncoplot
oncoplot(maf = luad_maf, top = 10, fontSize = 12)
# luad_mutsig <- "../mutsig/my_LUAD_result.sig_genes.txt"
luad_mutsig <- "../mutsig/new_LUAD_sig_genes_top.txt"
pdf(file = "onc2.pdf", width = 15, height = 15, paper = "special", bg = "white", pointsize = 9)
oncoplot(maf = luad_maf, top = 20,  mutsig = luad_mutsig, mutsigQval = 1e-12 )
dev.off()

# Transition and transversions
luad_titv <- titv(maf = luad_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = luad_titv)

# plot VAF
plotVaf(maf = luad_maf)

# somatic analyses
somaticInteractions(maf = luad_maf, top = 30, pvalue = c(0.05, 0.1))

# Detecting cancer driver genes based on positional clustering
luad.sig <- oncodrive(maf = luad_maf, minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = luad.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 5)

# Pan-cancer comparison about mutsigCV result
luad_mutsig <- "../mutsig/my_LUAD_result.sig_genes.txt"
pancanComparison(mutsigResults = luad_mutsig, qval = 0.01, cohortName = "TCGA-LUAD")

##> Mutation Signature Analysis
require(NMF)
luad_tnm <- trinucleotideMatrix(maf=luad_maf, ref_genome = "G:/biodata/reference/hg38.fa", ignoreChr = "chrM")
plotApobecDiff(tnm=luad_tnm, maf=luad_maf)
luad_signature <- extractSignatures(luad_tnm, nTry = 5, plotBestFitRes = TRUE)
plotSignatures(luad_signature, title_size = 1.5)
luad_se <- signatureEnrichment(maf=luad_maf, sig_res = luad_signature)
