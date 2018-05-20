setwd("C:/Users/wangshx/Desktop")

#>> Load pkgs and configuration
#options(stringsAsFactors=FALSE)
library(tidyverse)
library(data.table)
library(maftools)


source("data/neoQ/Allfunctions.R")
#>> prepare Neoantigen data
NQ_obj   <- summaryNQres("data/luad/neoResults/neoantigen_Quality_TCGA-LUAD.txt")
Neo_obj  <- summaryMergedNeos("data/luad/neoResults/merged_neoantigens.tsv")
neo <- dplyr::full_join(NQ_obj$NQsummary, Neo_obj$Neosummary, by=c("Sample"="Sample.Name")) 

#>> prepare Tumor-infiltrating lymphocytes data from three method:
#       Timer, Cibersort, QuanTIseq
TILtimer <- fread(input = "data/TIL/TCGA_sixcell_population_TIMER.txt")
TILtcia  <- fread(input = "data/TIL/TCGA-all-cellTypeFractionsAll.tsv")
TILtcia  <- TILtcia[!duplicated(TILtcia),]
TILtimer <- TILtimer[!duplicated(TILtimer),]

TILtcia_luad <- subset(TILtcia, disease=="LUAD")
TILtcia_wide_quanTIseq <- TILtcia_luad %>% data.frame %>%
                            select(patientBarcode:quanTIseq_lsei_TIL10) %>%  
                            spread(cell_type, quanTIseq_lsei_TIL10) %>% 
                            rename(Tumor_Sample_Barcode=patientBarcode)
TILtcia_wide_cibersort <- TILtcia_luad %>% data.frame %>%
                            select(-quanTIseq_lsei_TIL10) %>%  
                            spread(cell_type,  cibersort_LM22) %>% 
                            rename(Tumor_Sample_Barcode=patientBarcode)

colnames(TILtcia_wide_quanTIseq) <- make.names(colnames(TILtcia_wide_quanTIseq))
colnames(TILtcia_wide_cibersort) <- make.names(colnames(TILtcia_wide_cibersort))

TILtimer <- TILtimer %>% data.frame %>% rename(Tumor_Sample_Barcode=sample)

#any(duplicated(TILtcia_wide_cibersort$Tumor_Sample_Barcode))
#any(duplicated(TILtimer$Tumor_Sample_Barcode))
dup_rows <- which(duplicated(TILtimer$Tumor_Sample_Barcode))
dup_ids  <- TILtimer[dup_rows, ]$Tumor_Sample_Barcode
#TILtimer[TILtimer$Tumor_Sample_Barcode%in%dup_ids, ]

# data of these duplacated samples are variable, remove them
TILtimer <- TILtimer %>% filter(!(Tumor_Sample_Barcode%in%dup_ids))

# unify tumor sample barcode
TILtcia_wide_cibersort$Tumor_Sample_Barcode <- paste0(TILtcia_wide_cibersort$Tumor_Sample_Barcode, "-01")
TILtcia_wide_quanTIseq$Tumor_Sample_Barcode <- paste0(TILtcia_wide_quanTIseq$Tumor_Sample_Barcode, "-01")
TILtimer$Tumor_Sample_Barcode <- substr(TILtimer$Tumor_Sample_Barcode, 1, 15)

# create TIL list for storing three datasets
TIL <- list(timer=TILtimer, cibersort=TILtcia_wide_cibersort, quanTIseq=TILtcia_wide_quanTIseq)

# remove objects
rm(TILtcia, TILtcia_luad, TILtcia_wide_cibersort, TILtcia_wide_quanTIseq, TILtimer, dup_ids, dup_rows); gc()


#> prepare Maf, Gistic2 and Clinical data
persCli  <- fread(input = "data/luad/LUAD_clinicalMatrix/data")

# update TIL timer
TIL$timer <- subset(TIL$timer, Tumor_Sample_Barcode%in%persCli$sampleID)

save(TIL,file = "data/TCGA_LUAD_TIL.RData")

# select and filter clinical data
clin <- persCli %>% as.tibble %>% 
    select(sampleID:STK11, `_EVENT`, `_OS`:`_OS_UNIT`, additional_pharmaceutical_therapy:anatomic_neoplasm_subdivision_other, days_to_additional_surgery_locoregional_procedure:days_to_new_tumor_event_after_initial_treatment, gender:icd_o_3_site, new_neoplasm_event_type:number_pack_years_smoked, pathologic_M:pathologic_stage, person_neoplasm_cancer_status, primary_therapy_outcome_success, radiation_therapy, residual_tumor,sample_type, stopped_smoking_year, targeted_molecular_therapy, tobacco_smoking_history, tobacco_smoking_history_indicator, vital_status) %>% rename(EVENT=`_EVENT`, OS=`_OS`, OS_IND=`_OS_IND`, OS_UNIT=`_OS_UNIT`, Tumor_Sample_Barcode=sampleID) %>% as.data.frame

save(clin, file="data/TCGA_LUAD_clin.RData")
# EVENT, OS_IND and vital_status column provide sample info and has been checked

# read maf file and gistic2 files
maf <- fread("data/TCGA.LUAD.mutect.hg38.somatic.maf", stringsAsFactors = FALSE)
#maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 15)

all.lesions <- "data/luad/gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_lesions.conf_99.txt"
amp.genes   <- "data/luad/gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt"
del.genes   <- "data/luad/gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt"
score.gis   <- "data/luad/gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/scores.gistic"

# retrieve primary tumor data
clin2 <- clin %>% filter(sample_type == "Primary Tumor") 
clin2$Tumor_Sample_Barcode <- substr(clin2$Tumor_Sample_Barcode, 1, 12)

luad_maf <- read.maf(maf = maf, clinicalData = clin2, gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = score.gis, isTCGA = TRUE)

save(luad_maf, file="data/TCGA_LUAD_Maf.RData")

luad_gistic <- readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = score.gis, isTCGA = TRUE)
save(luad_gistic, file="data/TCGA_LUAD_gistic.RData")

# calculate mutation load, nonsynonymous substitutions (which can be divided into missense or nonsense mutations) 
sampleSummary <- getSampleSummary(luad_maf)
sampleSummary[, Nonsynonmous_Mutation:=(Missense_Mutation+Nonsense_Mutation)]
sampleSummary$Tumor_Sample_Barcode <- paste0(sampleSummary$Tumor_Sample_Barcode, "-01")
# merge sampleSummary with neo
MutationSummary <- dplyr::full_join(sampleSummary, neo, by=c("Tumor_Sample_Barcode"="Sample")) 

save(MutationSummary, file = "data/summary_of_Mutations_and_Neoantigens.RData")

cnv      <- fread(input = "data/luad/Gistic2_CopyNumber_Gistic2_all_data_by_genes/data")
geneExpr <- fread(input = "data/luad/HiSeqV2/data")
geneMut  <- fread(input = "data/luad/mutation_broad_gene/data")
genePro  <- fread(input = "data/luad/RPPA/data")

cnv <- cnv %>% rename(GeneSymbol = `Gene Symbol`)
geneExpr <- geneExpr %>% rename(GeneSymbol = sample)
geneMut  <- geneMut  %>% rename(GeneSymbol = sample)
genePro  <- genePro  %>% rename(GeneSymbol = sample)

GeneSummary <- list(CNV=cnv, Expr=geneExpr, Mut=geneMut, Protein=genePro) 
save(GeneSummary, file="data/summary_of_genes.RData")

##> Mutation Signature Analysis
require(NMF)
luad_tnm <- trinucleotideMatrix(maf=luad_maf, ref_genome = "G:/biodata/reference/hg38.fa", ignoreChr = "chrM")
#plotApobecDiff(tnm=luad_tnm, maf=luad_maf)
luad_signature <- extractSignatures(luad_tnm, nTry = 5, plotBestFitRes = TRUE)
#plotSignatures(luad_signature, title_size = 1.5)
luad_se <- signatureEnrichment(maf=luad_maf, sig_res = luad_signature)
#plotEnrichmentResults(enrich_res = luad_se, pVal = 0.0000000001)

mut_signatures <- c(luad_tnm, luad_signature, luad_se)
save(mut_signatures, file="data/summary_of_mutation_signatures.RData")

rm(clin, clin2, cnv, geneExpr, geneMut, genePro, GeneSummary, luad_gistic, luad_se, luad_signature, maf, mut_signatures, MutationSummary, neo, Neo_obj, NQ_obj, persCli, sampleSummary, TIL, luad_tnm);gc()
##> info tumor heterogeneity | This run on HPC
# seg <- "data/luad/broadinstitute_seg/LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"
# rawSeg <- read_tsv(seg)
# rawSeg$Sample <- substr(rawSeg$Sample, 1, 12)
# write_tsv(rawSeg, path="data/luad/broadinstitute_seg/segmented_scna_minus_germline_cnv_hg19_name_chopoff.seg.txt")
# 
# seg <- "data/luad/broadinstitute_seg/segmented_scna_minus_germline_cnv_hg19_name_chopoff.seg.txt"
# luad_het <- inferHeterogeneity(maf=luad_maf, tsb = getSampleSummary(luad_maf)$Tumor_Sample_Barcode, segFile = seg)
# 
# luad_het_noCNV <- inferHeterogeneity(maf=luad_maf, tsb = getSampleSummary(luad_maf)$Tumor_Sample_Barcode)
# #plotClusters(clusters = luad_het)

# load result of heterogeneity
load("C:/Users/wangshx/Desktop/data/inferHeter_remove_CNV.RData")
load("C:/Users/wangshx/Desktop/data/inferHeter.RData")

math_rm_CNA <- luad_het$clusterData %>% select(Tumor_Sample_Barcode, MATH) %>% filter(!is.na(MATH)) %>% distinct() %>% data.table::setDT() %>% rename(MATH_rm_CNA=MATH)
math_con_CNA <- luad_het_noCNV$clusterData %>% select(Tumor_Sample_Barcode, MATH) %>% distinct()
MATH <- dplyr::full_join(math_rm_CNA, math_con_CNA, by="Tumor_Sample_Barcode")

heter <- list(Heter=luad_het, Heter_noCNV=luad_het_noCNV, MATH=MATH)
save(heter, file="data/TCGA_LUAD_heterogeneity.RData")


# #>> test maftools plots
# plotmafSummary(maf=luad_maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
# oncoplot(maf = luad_maf, top = 10, fontSize = 12)
# # changing colors, including MutSig Q-Values, adding annotations and sorting
# col <- RColorBrewer::brewer.pal(n=10, name="Paired")
# names(col) <- c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
#                 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del', "Amp", "Del")
# tnmcolors <- RColorBrewer::brewer.pal(n=9, name="Spectral")
# names(tnmcolors) <- names(table(getClinicalData(luad_maf)$pathologic_stage))[-1]
# tnmcolors <- list(pathlogic_stage = tnmcolors)
# 
# # mutsig result
# luad_mutsig <- "data/mutsig/my_LUAD_result.sig_genes.txt"
# 
# oncoplot(maf = luad_maf, colors = col, top=10, mutsig = luad_mutsig, mutsigQval = 0.01, clinicalFeatures = "pathologic_stage", sortByAnnotation = TRUE, annotationColor = tnmcolors )
# somaticInteractions(maf=luad_maf, top=25, pvalue = c(0.05, 0.1))
# mafSurvival(maf=luad_maf, samples=getSampleSummary(luad_maf)$Tumor_Sample_Barcode[1:250], time="OS", Status = "EVENT", isTCGA = TRUE)
#tcgaCompare(luad_maf, cohortName = "test")
# gisticChromPlot(luad_gistic, markBands = "all")
# gisticBubblePlot(luad_gistic)
