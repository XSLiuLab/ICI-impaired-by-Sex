# Explore analysis for 4 immunotherapy datasets
library(tidyverse)
library(data.table)
library(maftools)
library(NMF)
library(ggsci)
library(gridExtra)
library(survival)
library(survminer)

ref_38genome = "G:/biodata/reference/hg38.fa"
ref_19genome = "G:/biodata/reference/hg19.fa"

setwd("G:/biodata/immunotherapyDatasets/")
# load functions
source("C:/Users/wangshx/Desktop/data/neoQ/Allfunctions.R")

##> load immunotherapy datasets
# clinical datasets
load("Rdata/Science_Rizvi_sampleInfo.RData")
load("Rdata/JCO_Rizvi_sampleInfo.RData")
load("Rdata/Hellmann_sampleInfo.RData")
load("Rdata/Forde_sampleInfo.RData")


##> load and calculate tumor mutation burden from maf data
load("Rdata/JCO_Rizvi_maf.RData")
sampleInfo_JCO_Rizvi <- getSampleTMB(JCO_Rizvi_maf) %>% 
    left_join(x=JCO_Rizvi_sampleInfo, by="Tumor_Sample_Barcode")
load("Rdata/Hellmann_maf.RData")
sampleInfo_Hellmann <- getSampleTMB(Hellmann_maf) %>%
    left_join(x=Hellmann_sampleInfo, by="Tumor_Sample_Barcode")
load("Rdata/Forde_maf.RData")
sampleInfo_Forde <- getSampleTMB(Forde_maf) %>%
    left_join(x=Forde_sampleInfo, by="Tumor_Sample_Barcode")
load("Rdata/Science_Rizvi_maf.RData")
sampleInfo_Sci_Rizvi <-  getSampleTMB(Science_Rizvi_maf) %>% 
    left_join(x=Science_Rizvi_sampleInfo, by="Tumor_Sample_Barcode")
# Science Rizvi mutation only contatins non-syn mutation data
sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% select(-TMB_Total) %>% dplyr::rename(TMB_Total = Total_Mutation_Burden)


sampleInfo_Sci_Rizvi <- setDT(sampleInfo_Sci_Rizvi)
sampleInfo_Sci_Rizvi[, Gender:=ifelse(Gender=="F", "Female", "Male")]
p1_1 <- compareMutPlot(sampleInfo_Sci_Rizvi[Clinical_Benefit %in% c("DCB", "NDB")])
sampleInfo_JCO_Rizvi <- setDT(sampleInfo_JCO_Rizvi)
p1_2 <- compareMutPlot(sampleInfo_JCO_Rizvi[Clinical_Benefit %in% c("DCB", "NDB")])
sampleInfo_Hellmann <- setDT(sampleInfo_Hellmann)
sampleInfo_Hellmann[, Gender:=ifelse(Gender=="female", "Female", "Male")]
p1_3 <- compareMutPlot(sampleInfo_Hellmann[Clinical_Benefit %in% c("DCB", "NDB")])
sampleInfo_Forde <- setDT(sampleInfo_Forde)
sampleInfo_Forde[, Gender:=ifelse(Gender=="F", "Female", "Male")]
compareMutPlot(sampleInfo_Forde[Clinical_Benefit %in% c("DCB", "NDB")])

p2_1 <- compareMutPlot(subset(sampleInfo_Sci_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_2 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_3 <- compareMutPlot(subset(sampleInfo_Hellmann, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_4 <- compareMutPlot(subset(sampleInfo_Forde, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")

p3_1 <- compareMutPlot(sampleInfo_Sci_Rizvi[Clinical_Benefit %in% c("DCB", "NDB")], value = "TMB_NonsynVariants")
p3_2 <- compareMutPlot(sampleInfo_JCO_Rizvi[Clinical_Benefit %in% c("DCB", "NDB")], value = "TMB_NonsynVariants")
p3_3 <- compareMutPlot(sampleInfo_Hellmann[Clinical_Benefit %in% c("DCB", "NDB")], value = "TMB_NonsynVariants")


groupSummary(sampleInfo_Sci_Rizvi, summarise_var = TMB_NonsynVariants, Gender)
groupSummary(sampleInfo_JCO_Rizvi, summarise_var = TMB_NonsynVariants, Gender)
groupSummary(sampleInfo_Hellmann, summarise_var = TMB_NonsynVariants, Gender)
groupSummary(sampleInfo_Forde, summarise_var = TMB_NonsynVariants, Gender)

p1 <- grid.arrange(p1_1, p1_2, p1_3, nrow=1, ncol=3)
ggsave(filename = "Total_Mutation_asTMB_across_datasets.pdf", plot=p1, height = 3,
       width = 8)
p2 <- grid.arrange(p2_1, p2_2, p2_3, nrow=1, ncol=3)
ggsave(filename = "Nonsyn_SNVs_asTMB_across_datasets.pdf", plot=p2, height = 3,
       width = 8)

p3 <- grid.arrange(p3_1, p3_2, p3_3, nrow=1, ncol=3)
ggsave(filename = "Nonsyn_Variants_asTMB_across_datasets.pdf", plot=p3, height = 3,
       width = 8)


# JCO dataset has different gene panels
table(sampleInfo_JCO_Rizvi$Gene_Panel)
p4_1 <- compareMutPlot(sampleInfo_JCO_Rizvi[Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT341"], value = "TMB_NonsynSNP",label_name = "p.format")
p4_2 <- compareMutPlot(sampleInfo_JCO_Rizvi[Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT410"], value = "TMB_NonsynSNP", label_name = "p.format")
p4_3 <- compareMutPlot(sampleInfo_JCO_Rizvi[Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT468"], value = "TMB_NonsynSNP", label_name = "p.format")

p4 <- grid.arrange(p4_1, p4_2, p4_3, nrow=1, ncol=3)
ggsave(filename = "JCO_dataset_TMB_nonsynSNP_across_Gene_Panels.pdf", plot=p4, height = 3,
       width = 8)

# compare different ways for TMB representation
corr_dat <- rbind(sampleInfo_Hellmann[, .(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants)],
                  sampleInfo_JCO_Rizvi[, .(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants)],
                  sampleInfo_Sci_Rizvi[, .(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants)])

p5 <- ggstatsplot::ggcorrmat(
    data = corr_dat,
    corr.method = "spearman",                # correlation method
    sig.level = 0.005,                        # threshold of significance
    cor.vars = TMB_Total:TMB_NonsynVariants,     # a range of variables can be selected  
    #cor.vars.names = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Width"),
    title = "Correlalogram for Different TMB Calculation Ways",
    #subtitle = "Iris dataset by Anderson",
    caption = expression(
        paste(
            italic("Note"),
            ": X denotes correlation non-significant at ",
            italic("p "),
            "< 0.005; adjusted alpha"
        )
    )
)

# ggsave(filename = "compare_different_ways_for_TMB_representation.pdf", plot = p5,
#        width = 4, height = 4)

p6_1 <- compareBoxplot(sampleInfo_Forde, x="Gender", y="TMB_NonsynVariants")
p6_2 <- compareBoxplot(sampleInfo_Hellmann, x="Gender", y="TMB_NonsynVariants")
p6_3 <- compareBoxplot(sampleInfo_Sci_Rizvi, x="Gender", y="TMB_NonsynVariants")
p6_4 <- compareBoxplot(sampleInfo_JCO_Rizvi, x="Gender", y="TMB_NonsynVariants")

library(cowplot)
p6 <- plot_grid(p6_1, p6_2, p6_3, p6_4, nrow=2, ncol=2, align = "v")
ggsave("Compare-TMB-between-F-and-M.pdf", plot = p6, width = 5, height = 5)

rm(list = ls(pattern = "^p"))

generateHLAfile(sampleInfo_Forde, path="Rdata/Forde_HLA.tsv")
generateHLAfile(sampleInfo_Hellmann, path="Rdata/Hellmann_HLA.tsv")
generateHLAfile(sampleInfo_Sci_Rizvi, path="Rdata/Sci_Rizvi_HLA.tsv")

##> retrieve functional mutation for neoantigen load and neoantigen quality computation
write_tsv(Forde_maf@data, path="Rdata/NQ_Forde.maf")
write_tsv(Hellmann_maf@data, path="Rdata/NQ_Hellmann.maf")
write_tsv(Science_Rizvi_maf@data, path="Rdata/NQ_Sci_Rizvi.maf")


##> cache data and code for now
save(sampleInfo_Forde, sampleInfo_Hellmann, sampleInfo_JCO_Rizvi, sampleInfo_Sci_Rizvi,
     Forde_maf, Hellmann_maf, JCO_Rizvi_maf, Science_Rizvi_maf, ref_19genome, ref_38genome,
     autoMutSig, autoTumorHeter, file = "Rdata/cache_data_code.RData")
rm(list = ls()); gc()
        
