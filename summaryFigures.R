# summaryFigures


# mutation compare
p2_1 <- compareMutPlot(subset(sampleInfo_Sci_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP", label_name = "p.format")
p2_2 <- compareMutPlot(subset(sampleInfo_Hellmann, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP", label_name = "p.format")
p2_3 <- compareMutPlot(subset(sampleInfo_Forde, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP", label_name = "p.format")

p2_1 <- ggpar(p2_1, ylab="Tumor Mutation Burden",  xlab = "")
p2_2 <- ggpar(p2_3, ylab="Tumor Mutation Burden",  xlab = "") 
p2_3 <- ggpar(p2_3, ylab="Tumor Mutation Burden",  xlab = "")

ggsave("summaryFigures/突变比较-NSCLC-Science_Rizvi.tiff", plot = p2_1, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-NSCLC-CancerCell_Hellmann.tiff", plot = p2_2, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-NSCLC-NEJM-Forde.tiff", plot = p2_3, width = 5, height = 4, dpi=300)

p3_1 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT341"), value = "TMB_NonsynSNP",label_name = "p.format")
p3_2 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT410"), value = "TMB_NonsynSNP", label_name = "p.format")
p3_3 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB") & Gene_Panel == "IMPACT468"), value = "TMB_NonsynSNP", label_name = "p.format")

p3_1 <- ggpar(p3_1, ylab="Tumor Mutation Burden",  xlab = "")
p3_2 <- ggpar(p3_2, ylab="Tumor Mutation Burden",  xlab = "") 
p3_3 <- ggpar(p3_3, ylab="Tumor Mutation Burden",  xlab = "")

ggsave("summaryFigures/突变比较-NSCLC-JCO_Rizvi-341Gene.tiff", plot = p3_1, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-NSCLC-JCO_Rizvi-410Gene.tiff", plot = p3_2, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-NSCLC-JCO_Rizvi-468Gene.tiff", plot = p3_3, width = 5, height = 4, dpi=300)

p4_1 <- compareMutPlot(MELA_NEJM, group1 = "Gender", 
                                group2 = "Clinical_Benefit", 
                                value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
p4_2 <- compareMutPlot(MELA_Cell, group1 = "Gender", 
                                group2 = "Clinical_Benefit", 
                                value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
p4_3 <- compareMutPlot(MELA_Science, group1 = "Gender", group2 = "Clinical_Benefit",
                               value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")

ggsave("summaryFigures/突变比较-MELA-NEJM.tiff", plot = p4_1, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-MELA-Cell.tiff", plot = p4_2, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/突变比较-MELA-Science.tiff", plot = p4_3, width = 5, height = 4, dpi=300)

# ROC compare
p5_1 <- plotROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
p5_2 <- plotROC(sampleInfo_Hellmann, TMB_NonsynSNP, Clinical_Benefit, Gender)
p5_3 <- plotROC(sampleInfo_Forde, TMB_NonsynSNP, Clinical_Benefit, Gender)

p5_4 <- plotROC(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT341"), TMB_NonsynSNP,
                Clinical_Benefit, Gender)
p5_5 <- plotROC(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT410"), TMB_NonsynSNP,
                Clinical_Benefit, Gender)
p5_6 <- plotROC(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT468"), TMB_NonsynSNP,
                Clinical_Benefit, Gender)

p5_7 <- plotROC(MELA_NEJM, mutation_load, Clinical_Benefit, Gender, positive = "DCB")
p5_8 <- plotROC(MELA_Cell, mutation_load, Clinical_Benefit, Gender)
p5_9 <- plotROC(MELA_Science, mutation_load, Clinical_Benefit, Gender)

ggsave("summaryFigures/ROC-NSCLC-Science-Rizvi.tiff", plot = p5_1, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC-CancerCell-Hellmann.tiff", plot = p5_2, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC-NEJM-Forde.tiff", plot = p5_3, width = 5, height = 4, dpi=300)

ggsave("summaryFigures/ROC-NSCLC-JCO-Rizvi_341Gene.tiff", plot = p5_4, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC-JCO-Rizvi_410Gene.tiff", plot = p5_5, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC-JCO-Rizvi_468Gene.tiff", plot = p5_6, width = 5, height = 4, dpi=300)

ggsave("summaryFigures/ROC-MELA-NEJM.tiff", plot = p5_7, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-MELA-Cell.tiff", plot = p5_8, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-MELA-Science.tiff", plot = p5_9, width = 5, height = 4, dpi=300)

p6_1 <- sampleInfo_Hellmann %>% mutate(Smoke=ifelse(Smoking_History == "never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
p6_2 <- sampleInfo_Sci_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
p6_3 <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
p6_4 <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>%
    filter(Smoke == "Never", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")

ggsave("summaryFigures/ROC-NSCLC吸烟-CancerCell-Hellmann.tiff", plot = p6_1, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC吸烟-Science_Rizvi.tiff", plot = p6_2, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC吸烟-JCO-Rizvi-410Gene.tiff", plot = p6_3, width = 5, height = 4, dpi=300)
ggsave("summaryFigures/ROC-NSCLC从不吸烟-JCO-Rizvi-410Gene.tiff", plot = p6_4, width = 5, height = 4, dpi=300)


# summary data from Survival_HR_for... R

