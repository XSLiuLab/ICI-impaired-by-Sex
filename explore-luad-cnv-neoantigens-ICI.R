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

# # check NSCLC
# signatureFlow(maf = Hellmann_maf, ref_genome = ref_19genome, filename = "NSCLC_Hellmann_Signature")
# signatureFlow(maf = Science_Rizvi_maf, ref_genome = ref_19genome, filename = "NSCLC_sciRizvi_Signature")
# signatureFlow(maf = JCO_Rizvi_maf, ref_genome = ref_19genome, filename = "NSCLC_jcoRizvi_Signature")
# # validata result shown in APOBEC3B paper
# tsb_DCB <- list()
# tsb_NDB <- list()
# 
# tsb_DCB$hellmann <- sampleInfo_Hellmann %>% filter(Clinical_Benefit == "DCB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()
# tsb_NDB$hellmann <- sampleInfo_Hellmann %>% filter(Clinical_Benefit == "NDB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()    
# tsb_DCB$jco <- sampleInfo_JCO_Rizvi %>% filter(Clinical_Benefit == "DCB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()
# tsb_NDB$jco <- sampleInfo_JCO_Rizvi %>% filter(Clinical_Benefit == "NDB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()   
# tsb_DCB$sci <- sampleInfo_Sci_Rizvi %>% filter(Clinical_Benefit == "DCB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()
# tsb_NDB$sci <- sampleInfo_Sci_Rizvi %>% filter(Clinical_Benefit == "NDB") %>% 
#     select(Tumor_Sample_Barcode) %>% unlist %>% as.character()   
# 
# signatureFlow(maf = subsetMaf(Hellmann_maf, tsb = tsb_DCB$hellmann, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_Hellmann_DCB")
# signatureFlow(maf = subsetMaf(Hellmann_maf, tsb = tsb_NDB$hellmann, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_Hellmann_NDB")
# signatureFlow(maf = subsetMaf(JCO_Rizvi_maf, tsb = tsb_DCB$jco, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_JCO_DCB")
# signatureFlow(maf = subsetMaf(JCO_Rizvi_maf, tsb = tsb_NDB$jco, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_JCO_NDB")
# signatureFlow(maf = subsetMaf(Science_Rizvi_maf, tsb = tsb_DCB$sci, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_Sci_DCB")
# signatureFlow(maf = subsetMaf(Science_Rizvi_maf, tsb = tsb_NDB$sci, mafObj = TRUE),
#               ref_genome = ref_19genome, filename = "vAPOBEC/NSCLC_Sci_NDB")

# further test TCW and total mutation correlation
apobec_jcoRizvi <- sig_jcoRizvi$tnm$APOBEC_scores %>% select(Tumor_Sample_Barcode, APOBEC_Enrichment, APOBEC_Enriched) %>% 
    dplyr::left_join(sampleInfo_JCO_Rizvi %>% select(Tumor_Sample_Barcode, sTMB), 
                     by='Tumor_Sample_Barcode')
apobec_sciRizvi <- sig_sciRizvi$tnm$APOBEC_scores %>% select(Tumor_Sample_Barcode, APOBEC_Enrichment, APOBEC_Enriched) %>% 
    dplyr::left_join(sampleInfo_Sci_Rizvi %>% select(Tumor_Sample_Barcode, sTMB), 
                     by='Tumor_Sample_Barcode')
apobec_hellmann <- sig_hellmann$tnm$APOBEC_scores %>% select(Tumor_Sample_Barcode, APOBEC_Enrichment, APOBEC_Enriched) %>% 
    dplyr::left_join(sampleInfo_Hellmann %>% select(Tumor_Sample_Barcode, sTMB), 
                     by='Tumor_Sample_Barcode')

apobec_jcoRizvi %>% filter(!is.na(APOBEC_Enriched)) %>% 
    summarise(CR = cor(APOBEC_Enrichment, sTMB))
apobec_sciRizvi %>% filter(!is.na(APOBEC_Enriched)) %>% 
    summarise(CR = cor(APOBEC_Enrichment, sTMB))
apobec_hellmann %>% filter(!is.na(APOBEC_Enriched)) %>% 
    summarise(CR = cor(APOBEC_Enrichment, sTMB))
##
## APOBEC enrichment is not correlated with TMB


#----------- handly assign the signature to factors ------#
sig_Assignment <- list()
sig_Assignment$hellmann <- sig_hellmann$se$Signature_Assignment %>% 
    filter(Signature != "Signature_1") %>%  # filter UV signature
    mutate(Aetiology = case_when(
        Signature == "Signature_2" ~ "dMMR",
        Signature %in% c("Signature_3", "Signature_5") ~ "APOBEC-C>T",
        Signature == "Signature_4" ~ "Smoking"
    ))
sig_Assignment$jcoRizvi <- sig_jcoRizvi$se$Signature_Assignment %>% 
    mutate(Aetiology = case_when(
        Signature == "Signature_1" ~ "APOBEC-C>G",
        Signature == "Signature_2" ~ "APOBEC-C>T",
        Signature == "Signature_3" ~ "Smoking"
    ))
sig_Assignment$sciRizvi <- sig_sciRizvi$se$Signature_Assignment %>% 
    mutate(Aetiology = case_when(
        Signature == "Signature_1" ~ "Smoking",
        Signature == "Signature_2" ~ "dMMR",
        Signature == "Signature_3" ~ "APOBEC-C>G"
    ))

luad_hellmann_info <- sig_Assignment$hellmann %>% dplyr::inner_join(sampleInfo_Hellmann,
                                                                    by = "Tumor_Sample_Barcode")
luad_jcoRizvi_info <- sig_Assignment$jcoRizvi %>% dplyr::inner_join(sampleInfo_JCO_Rizvi,
                                                                    by = "Tumor_Sample_Barcode")
luad_sciRizvi_info <- sig_Assignment$sciRizvi %>% dplyr::inner_join(sampleInfo_Sci_Rizvi,
                                                                    by = "Tumor_Sample_Barcode")

# save info
save(luad_hellmann_info, luad_jcoRizvi_info, luad_sciRizvi_info, file = "Rdata/ICI_LUAD_Info.RData")

# merge three dataset into one 
luad <- bind_rows(luad_hellmann_info %>% select(Tumor_Sample_Barcode:Aetiology, Age:PFS_Event, Clinical_Benefit, 
                                        TMB_NonsynVariants, NQuality:sNeo) %>% mutate(PDL1 = PDL1_Expression)
                                        %>% select(-PDL1_Expression) %>% mutate(PDL1 = ifelse(PDL1 == 'Unknown', NA, as.numeric(PDL1))) ,
          luad_jcoRizvi_info %>% select(Tumor_Sample_Barcode:Aetiology, Gender:Smoking_History, 
                                        PFS_Months:Clinical_Benefit, PDL1_Score, Histology, 
                                        TMB_NonsynVariants, sTMB) %>% mutate(PDL1 = PDL1_Score) %>% select(-PDL1_Score), 
          luad_sciRizvi_info %>% select(Tumor_Sample_Barcode:Smoking_History, `PD-L1`, PFS_Months, PFS_Event,
                                         Clinical_Benefit,TMB_NonsynVariants, NQuality:sNeo)  %>% mutate(PDL1_Status = `PD-L1`) %>% select(-`PD-L1`),
          .id = "Dataset")

luad <- luad %>% mutate(Dataset = case_when(
    Dataset == 1 ~ "Hellmann et al",
    Dataset == 2 ~ "Rizvi (2018) et al",
    Dataset == 3 ~ "Rizvi (2015) et al"
))

#save
save(luad, file="Rdata/ICI_luad_merge.RData")









#---------- Analyses -----------#
## ROC
plotROC <- function(.data, predict_col, target, group, positive="DCB", all=TRUE){
    if(!(require(tidyverse) & require(plotROC))){
        stop("--> tidyverse and plotROC packages are required..")
    } 
    
    predict_col <- enquo(predict_col)
    target <- enquo(target)
    group  <- enquo(group)
    
    predictN <- quo_name(predict_col)
    groupN   <- quo_name(group)
    
    df <- .data %>% dplyr::select(!! predict_col, !! target, !! group) %>%
        mutate(targetN = ifelse(!! target == positive, 1, 0)) %>% as.data.frame()
    if (all){
        df2 <- df 
        df2[, groupN] <- "ALL"
        
        df <- rbind(df, df2)
    }
    p  <- df %>%  ggplot(aes_string(m = predictN, 
                                    d = "targetN",
                                    color = groupN)) + geom_roc(show.legend = TRUE, labels=FALSE)
    p <- p + ggpubr::theme_classic2()
    
    ng <- levels(factor(df[, groupN]))
    if(length(ng) == 3){
        auc <- calc_auc(p)$AUC
        names(auc) <- ng
        auc <- base::sort(auc, decreasing = TRUE)
        p <- p + annotate("text", x = .75, y = .25, 
                          label = paste(names(auc)[1], " AUC =", round(auc[1], 3), "\n",
                                        names(auc)[2], " AUC =", round(auc[2], 3), "\n",
                                        names(auc)[3], " AUC =", round(auc[3], 3), "\n"),
                          size = 4)
    }
    
    p + xlab("1 - Specificity") + ylab("Sensitivity") + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
}

plotROC(luad_sciRizvi_info, sTMB, Clinical_Benefit, Aetiology, all = FALSE)
plotROC(luad_jcoRizvi_info, sTMB, Clinical_Benefit, Aetiology, all = FALSE)
plotROC(luad_hellmann_info, sTMB, Clinical_Benefit, Aetiology, all = FALSE)


## TMB
luad_sciRizvi_info %>% filter(Clinical_Benefit!="NR") %>% 
    ggplot(aes(x=Clinical_Benefit, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw() 
luad_jcoRizvi_info  %>% 
    ggplot(aes(x=Clinical_Benefit, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw()
luad_hellmann_info  %>% 
    ggplot(aes(x=Clinical_Benefit, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw()


luad_sciRizvi_info %>% filter(Clinical_Benefit!="NR") %>%
    ggplot(aes(x=Aetiology, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw() 
luad_jcoRizvi_info %>% ggplot(aes(x=Aetiology, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw()
luad_hellmann_info %>% ggplot(aes(x=Aetiology, y=sTMB, fill=Clinical_Benefit)) + 
    geom_boxplot() + theme_bw()

luad %>% ggplot(aes(x = Clinical_Benefit, y = sTMB)) +
    ggplot()
  

## survival
library(survival)
library(survminer)

ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ Aetiology, luad_sciRizvi_info), pval = TRUE)
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ Aetiology, luad_jcoRizvi_info), pval = TRUE)
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ Aetiology, luad_hellmann_info), pval = TRUE)


# check NQuality on survival
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          sampleInfo_Sci_Rizvi %>% 
                              mutate(group = ifelse(NQuality >= median(NQuality),
                                                                       "High", 
                                                                       "Low"))), pval = TRUE)
neo_fitness = data.table::fread("G:/biodata/immunotherapyDatasets/Neo_result/neoantigen_fitness_Rizvi.txt")
neo_fitness = neo_fitness %>% group_by(Sample) %>% 
    summarise(Neo_fitness = max(NeoantigenRecognitionPotential))
new_dat = dplyr::full_join(sampleInfo_Sci_Rizvi, neo_fitness, by=c("Tumor_Sample_Barcode"="Sample"))
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          new_dat %>% 
                              filter(!is.na(Neo_fitness)) %>%
                              mutate(group = ifelse(Neo_fitness >= median(Neo_fitness),
                                                                       "High", 
                                                                       "Low"))), pval = TRUE)
neo_fitness2 = data.table::fread("G:/biodata/immunotherapyDatasets/Neo_result/neoantigen_Quality_Sci_Rizvi.txt")
length(intersect(neo_fitness$Mutation, neo_fitness2$Mutation))
neo_fitness3 = data.table::fread("G:/biodata/immunotherapyDatasets/Neo_result/neoantigen_Quality2.txt")
length(intersect(neo_fitness$Mutation, neo_fitness3$Mutation))

neo_fitness3 = neo_fitness3 %>% group_by(Sample) %>% 
    summarise(Neo_fitness = max(NeoantigenRecognitionPotential))
new_dat = dplyr::full_join(sampleInfo_Sci_Rizvi, neo_fitness3, by=c("Tumor_Sample_Barcode"="Sample"))
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          new_dat %>% 
                              filter(!is.na(Neo_fitness)) %>%
                              mutate(group = ifelse(Neo_fitness >= median(Neo_fitness),
                                                    "High", 
                                                    "Low"))), pval = TRUE)
# neos = data.table::fread("C:/Users/wangshx/Desktop/source_code/SupplementaryDataFile7/InputData/neoantigens_Rizvi.txt")
# max(neos$MT.Score)

# check hellmann
neo_fitness4 = data.table::fread("G:/biodata/immunotherapyDatasets/Neo_result/neoantigen_Quality2_hellmann.txt")

neo_fitness4 = neo_fitness4 %>% group_by(Sample) %>% 
    summarise(Neo_fitness = max(NeoantigenRecognitionPotential))
new_dat = dplyr::full_join(sampleInfo_Hellmann, neo_fitness4, by=c("Tumor_Sample_Barcode"="Sample"))
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          new_dat %>% 
                              filter(!is.na(Neo_fitness)) %>%
                              mutate(group = ifelse(Neo_fitness >= median(Neo_fitness),
                                                    "High", 
                                                    "Low"))), pval = TRUE)
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          sampleInfo_Hellmann %>% 
                              filter(!is.na(NQuality)) %>% 
                              mutate(group = ifelse(NQuality >= median(NQuality),
                                                    "High", 
                                                    "Low"))), pval = TRUE)
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          sampleInfo_Hellmann %>% 
                              filter(!is.na(sTMB)) %>% 
                              mutate(group = ifelse(sTMB >= median(sTMB),
                                                    "High", 
                                                    "Low"))), pval = TRUE)
ggsurvplot(fit = surv_fit(Surv(PFS_Months, PFS_Event) ~ group,
                          sampleInfo_Sci_Rizvi %>% 
                              filter(!is.na(sTMB)) %>% 
                              mutate(group = ifelse(sTMB >= median(sTMB),
                                                    "High", 
                                                    "Low"))), pval = TRUE)
