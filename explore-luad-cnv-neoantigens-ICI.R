# explore the influence of CNV/TMB/Neoantigens/heter.. on patient survival and ICI response in LUAD

## mark data:
# JCO and Science Rizvi all got CNV segmentation files, Forde dataset has gene level cnv
  
setwd("G:/biodata/immunotherapyDatasets")
library(tidyverse)
# load cache MAF and sampleInfo
load("Rdata/cache_data_code.RData")
load("Rdata/cache_Unify_sampleinfo.RData")

#----------- keep only LUAD ---------------#

# check nsclc data subtype first
Forde_tsb <- Forde_maf@clinical.data %>% as.tibble() %>% select(Histology, Tumor_Sample_Barcode) %>% 
    filter(Histology == "Adenocarcinoma") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character

# Hellmann_maf has no label of luad, we treat non-squamous as luad (NSCLC include LUAD and LUSC)
Hellmann_tsb <- Hellmann_maf@clinical.data %>% as.tibble() %>% select(Histology, Tumor_Sample_Barcode) %>% 
    filter(Histology == "non-squamous") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character()

# It is a little strange that there is a lable about NSCLC, while LUAD, LUSC exist
JCO_Rizvi_tsb <- JCO_Rizvi_maf@clinical.data %>% as.tibble() %>% select(Histology, Tumor_Sample_Barcode) %>% 
    filter(Histology == "LUAD") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character()

Science_Rizvi_tsb <- Science_Rizvi_maf@clinical.data %>% as.tibble() %>% select(Histology, Tumor_Sample_Barcode) %>% 
    filter(Histology == "Adeno") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character()

# subset LUAD maf
luad_forde <- subsetMaf(Forde_maf, tsb = Forde_tsb, mafObj = TRUE)
luad_hellmann <- subsetMaf(Hellmann_maf, tsb = Hellmann_tsb, mafObj = TRUE)
luad_jcoRizvi <- subsetMaf(JCO_Rizvi_maf, tsb = JCO_Rizvi_tsb, mafObj = TRUE)
luad_sciRizvi <- subsetMaf(Science_Rizvi_maf, tsb = Science_Rizvi_tsb, mafObj = TRUE)

# save
save(luad_forde, luad_hellmann, luad_jcoRizvi, luad_sciRizvi, file = "Rdata/ICI_LUAD_MAF.RData")

#------------ Signature Analysis and Assignment ----------#
require(NMF)
tnm <- trinucleotideMatrix(maf=luad_forde, ref_genome = ref_19genome, ignoreChr = "chrM",
                           prefix = "chr", add = TRUE)
plotApobecDiff(tnm=tnm, maf=luad_forde)
signature <- extractSignatures(tnm, plotBestFitRes = TRUE, nTry = 5)
plotSignatures(signature, title_size = 1.5)
se <- signatureEnrichment(maf=luad_forde, sig_res = signature)

signatureFlow <- function(maf, ref_genome, filename=NULL, prefix = "chr", add = TRUE, ignoreChr = "chrM", 
                          useSyn = TRUE, nTry = 5, title_size = 1.2, parallel = "P4"){
    require(NMF)
    tnm <- trinucleotideMatrix(maf=maf, ref_genome = ref_genome, ignoreChr = ignoreChr,
                               prefix = prefix, add = add, useSyn = useSyn)
    pdf(file = paste0(filename, "APOBECcompare.pdf"), width = 6, height = 3)
    plotApobecDiff(tnm=tnm, maf=maf)
    dev.off()
    signature <- extractSignatures(tnm, nTry = nTry, plotBestFitRes = FALSE, parallel = parallel)
    pdf(file = paste0(filename, "MutationSignature.pdf"), width = 6, height = 6)
    plotSignatures(signature, title_size = title_size)
    dev.off()
    pdf(file = paste0(filename, "SignatureCluster.pdf"), width = 5, height = 3)
    se <- signatureEnrichment(maf=maf, sig_res = signature)
    dev.off()
    res <- list(signature=signature, se=se, tnm=tnm)
    return(res)
}
# signatureFlow(maf = luad_forde, ref_genome = ref_19genome)
sig_hellmann <- signatureFlow(maf = luad_hellmann, ref_genome = ref_19genome, filename = "LUAD_Hellmann_Signature")
sig_sciRizvi <- signatureFlow(maf = luad_sciRizvi, ref_genome = ref_19genome, filename = "LUAD_sciRizvi_Signature")
sig_jcoRizvi <- signatureFlow(maf = luad_jcoRizvi, ref_genome = ref_19genome, filename = "LUAD_jcoRizvi_Signature")

# check NSCLC
signatureFlow(maf = Hellmann_maf, ref_genome = ref_19genome, filename = "NSCLC_Hellmann_Signature")
signatureFlow(maf = Science_Rizvi_maf, ref_genome = ref_19genome, filename = "NSCLC_sciRizvi_Signature")
signatureFlow(maf = JCO_Rizvi_maf, ref_genome = ref_19genome, filename = "NSCLC_jcoRizvi_Signature")

