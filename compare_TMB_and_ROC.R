setwd("G:/biodata/immunotherapyDatasets/")
library(ggpubr)
load("Rdata/sampleInfo_cache.RData")

p2_1 <- compareMutPlot(subset(sampleInfo_Sci_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_2 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_3 <- compareMutPlot(subset(sampleInfo_Hellmann, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")
p2_4 <- compareMutPlot(subset(sampleInfo_Forde, Clinical_Benefit %in% c("DCB", "NDB")), value = "TMB_NonsynSNP")

p2_1 <- ggpar(p2_1 , title="Science Rizvi etc. dataset", ylab="Nonsynonymous Mutation", font.main = 12, xlab = "")
p2_2 <- ggpar(p2_2 , title="JCO Rizvi etc. dataset", ylab="Nonsynonymous Mutation", font.main = 12, xlab = "")
p2_3 <- ggpar(p2_3 , title="Hellmann etc. dataset", ylab="Nonsynonymous Mutation", font.main = 12, xlab = "") 
p2_4 <- ggpar(p2_4 , title="Forde etc. dataset", ylab="Nonsynonymous Mutation", font.main = 12, xlab = "")


library(cowplot)
p2 <- plot_grid(p2_1, p2_2, p2_3, p2_4, nrow=2, ncol=2,  align = "v")
ggsave(filename = "NonsynonmousMutation_between_Gender.pdf", plot=p2, height = 5,
       width = 6)

p3_1 <- compareMutPlot(subset(sampleInfo_Sci_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "NeoCounts")
#p3_2 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "NeoCounts")
p3_3 <- compareMutPlot(subset(sampleInfo_Hellmann, Clinical_Benefit %in% c("DCB", "NDB")), value = "NeoCounts")
p3_4 <- compareMutPlot(subset(sampleInfo_Forde, Clinical_Benefit %in% c("DCB", "NDB")), value = "NeoCounts")

p3_1 <- ggpar(p3_1 , title="Science Rizvi etc. dataset", ylab="Neoantigen Load", font.main = 12, xlab = "")
#p3_2 <- ggpar(p3_2 , title="JCO Rizvi etc. dataset", ylab="Neoantigen Load", font.main = 12, xlab = "")
p3_3 <- ggpar(p3_3 , title="Hellmann etc. dataset", ylab="Neoantigen Load", font.main = 12, xlab = "") 
p3_4 <- ggpar(p3_4 , title="Forde etc. dataset", ylab="Neoantigen Load", font.main = 12, xlab = "")
p3 <- plot_grid(p3_1, p3_3, p3_4, nrow=1, ncol=3,  align = "v")
ggsave(filename = "NeoantigenLoad_between_Gender.pdf", plot=p3, height = 3,
       width = 8)

##> calculate coefficient of variation (CV)
calcCV <- function(.data, target, ...){
    require(tidyverse)
    target <- enquo(target)
    group_var <- quos(...)
    
    .data %>% group_by(!!! group_var) %>% 
        summarise(
            N = n(),
            CV = sd(!! target, na.rm = TRUE)/ mean(!! target, na.rm = TRUE),
            Mean = mean(!! target, na.rm = TRUE),
            SD = sd(!! target, na.rm = TRUE)
        )
    
}


load("C:/Users/wangshx/Desktop/data/TCGA_LUAD_Maf.RData")

LUAD_TMB <- getSampleTMB(luad_maf)
LUAD_clin <- luad_maf@clinical.data %>% select(Tumor_Sample_Barcode, gender, pathologic_stage) %>% 
    filter(grepl("Stage", pathologic_stage)) %>% 
    rename(Gender = gender) %>%  mutate(Gender = ifelse(Gender=="MALE", "Male", "Female"),
                                        TNM = ifelse(pathologic_stage %in% c("Stage_I", "Stage_IA", "Stage_IB"),
                                                     "early", ifelse(pathologic_stage %in% c("Stage_II", "Stage_IIA", "Stage_IIB", "Stage_IIIA"),
                                                                     "mid", "late")))
LUAD_TMB2 <-  dplyr::full_join(LUAD_TMB, LUAD_clin, by="Tumor_Sample_Barcode") %>% 
    filter(Gender %in% c("Male", "Female"))


CV.nonsyn <- list()
CV.nonsyn$Sci_Rizvi <- calcCV(sampleInfo_Sci_Rizvi, target = TMB_NonsynSNP, Gender)
CV.nonsyn$JCO_Rizvi <- calcCV(sampleInfo_JCO_Rizvi, target = TMB_NonsynSNP, Gender, Gene_Panel)
CV.nonsyn$Hellmann  <- calcCV(sampleInfo_Hellmann, target = TMB_NonsynSNP, Gender)
CV.nonsyn$Forde     <- calcCV(sampleInfo_Forde, target = TMB_NonsynSNP, Gender)

CV.nonsyn$LUAD      <- calcCV(LUAD_TMB2, target =TMB_NonsynSNP, Gender, TNM )

calcCV(sampleInfo_Sci_Rizvi, target = NeoCounts, Gender)
calcCV(sampleInfo_Hellmann, target = NeoCounts, Gender)
calcCV(sampleInfo_Forde, target = NeoCounts, Gender)


# plot ROC curve
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
    
    df2 <- df 
    df2[, groupN] <- "ALL"
    
    df <- rbind(df, df2)
    
    p  <- df %>%  ggplot(aes_string(m = predictN, 
                                    d = "targetN",
                                    color = groupN)) + geom_roc(show.legend = TRUE, labels=FALSE)
    p <- p + ggpubr::theme_classic2()
    
    # p <- direct_label(p) + ggpubr::theme_classic2()
    # annotate("text", x = .75, y = .25, 
    #          label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)))
    # pairplot <- ggplot(longtest, aes(d = D, m = M, color = name)) + 
    #     geom_roc(show.legend = FALSE) + style_roc()
    # direct_label(pairplot)
    
    ng <- levels(factor(df[, groupN]))
    if(length(ng) == 3){
        auc <- calc_auc(p)$AUC
        names(auc) <- ng
        auc <- base::sort(auc, decreasing = TRUE)
        p <- p + annotate("text", x = .75, y = .25, 
                          label = paste(names(auc)[1], " AUC =", round(auc[1], 3), "\n",
                                        names(auc)[2], " AUC =", round(auc[2], 3), "\n",
                                        names(auc)[3], " AUC =", round(auc[3], 3), "\n"))
    }
    
    p + xlab("1 - Specificity") + ylab("Sensitivity") + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
}

# tmb ROC
tmbROC_SciRizvi <- plotROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggsave("tmbROC_SciRizvi.pdf", plot=tmbROC_SciRizvi, width=5, height=4)
tmbROC_JCORizvi <- plotROC(sampleInfo_JCO_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggsave("tmbROC_JCORizvi.pdf", plot=tmbROC_JCORizvi, width=5, height=4)
tmbROC_Hellmann <- plotROC(sampleInfo_Hellmann, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggsave("tmbROC_Hellmann.pdf", plot=tmbROC_Hellmann, width=5, height=4)
tmbROC_Forde <- plotROC(sampleInfo_Forde, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggsave("tmbROC_Forde.pdf", plot=tmbROC_Forde, width=5, height=4)

tmbROC_JCORizvi_IMPACT341 <- plotROC(sampleInfo_JCO_Rizvi[Gene_Panel=="IMPACT341"], TMB_NonsynSNP,
                           Clinical_Benefit, Gender)
tmbROC_JCORizvi_IMPACT410 <- plotROC(sampleInfo_JCO_Rizvi[Gene_Panel=="IMPACT410"], TMB_NonsynSNP,
                                     Clinical_Benefit, Gender)
tmbROC_JCORizvi_IMPACT468 <- plotROC(sampleInfo_JCO_Rizvi[Gene_Panel=="IMPACT468"], TMB_NonsynSNP,
                                     Clinical_Benefit, Gender)

tmbROC_JCORizvi_com3 <- cowplot::plot_grid(tmbROC_JCORizvi_IMPACT341, tmbROC_JCORizvi_IMPACT410, tmbROC_JCORizvi_IMPACT468,
                   nrow = 3, align = 'v')
ggsave("tmbROC_JCORizvi_com3.pdf", plot=tmbROC_JCORizvi_com3, width=5, height=12)

tmbROC_JCORizvi_com2 <- cowplot::plot_grid(tmbROC_JCORizvi_IMPACT341, tmbROC_JCORizvi_IMPACT410, 
                                           nrow = 1, align = 'h')
ggsave("tmbROC_JCORizvi_com2.pdf", plot=tmbROC_JCORizvi_com2, width=10, height=4)

# neoantigen ROC
neoROC_SciRizvi <- plotROC(sampleInfo_Sci_Rizvi, NeoCounts, Clinical_Benefit, Gender)
ggsave("neoROC_SciRizvi.pdf", plot=neoROC_SciRizvi, width=5, height=4)

neoROC_Hellmann <- plotROC(sampleInfo_Hellmann, NeoCounts, Clinical_Benefit, Gender)
ggsave("neoROC_Hellmann.pdf", plot=neoROC_Hellmann, width=5, height=4)
neoROC_Forde <- plotROC(sampleInfo_Forde, NeoCounts, Clinical_Benefit, Gender)
ggsave("neoROC_Forde.pdf", plot=neoROC_Forde, width=5, height=4)


rm(list = ls(pattern = "^tmb|^roc|^neo"))


# compare mutation between Gender in TCGA dataset
load("C:/Users/wangshx/Desktop/data/TCGA_LUAD_Maf.RData")

LUAD_TMB <- getSampleTMB(luad_maf)
LUAD_clin <- luad_maf@clinical.data %>% select(Tumor_Sample_Barcode, gender) %>% 
    rename(Gender = gender) %>%  mutate(Gender = ifelse(Gender=="MALE", "Male", "Female"))
                                        
LUAD_TMB2 <-  dplyr::full_join(LUAD_TMB, LUAD_clin, by="Tumor_Sample_Barcode") %>% 
    filter(Gender %in% c("Male", "Female"))

LUAD_TMB3 <- LUAD_TMB2 %>% mutate(TMB_Status = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP)[4],
                                         "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP)[2],
                                                      "Low", "Mid"))) %>% filter(TMB_Status != "Mid")


# LUAD_TMB3 <- LUAD_TMB2 %>% mutate(TMB_Status = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP, 0.75),
#                                                       "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP, 0.25),
#                                                                      "Low", "Mid"))) %>% filter(TMB_Status != "Mid")

LUAD_TMB3 %>% group_by(Gender) %>% 
    summarise(N=n(), N_DOWN=length(which(TMB_Status=="Low")), N_UP=length(which(TMB_Status=="High")))

rm(luad_maf); gc()
# load gene levels
load("C:/Users/wangshx/Desktop/data/summary_of_genes.RData")
selt_genes <- GeneSummary$Expr[GeneSymbol %in% c("PDCD1", "CD274", "PDCD1LG2", "CTLA4")]

selt_genes <- selt_genes %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue) %>% 
    rename(PD1 = PDCD1, PDL1 = CD274, PDL2 = PDCD1LG2)

# 
# LUAD_TMB3 <- LUAD_TMB2
LUAD_TMB3$Tumor_Sample_Barcode <- paste0(LUAD_TMB3$Tumor_Sample_Barcode, "-01")

luad_merge <- dplyr::left_join(LUAD_TMB3, selt_genes, by="Tumor_Sample_Barcode")
luad_info <- LUAD_TMB2 %>% 
    mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode, "-01") )%>%
    dplyr::left_join(y= selt_genes, by="Tumor_Sample_Barcode")

ggstatsplot::ggcorrmat(
    data = luad_info,
    corr.method = "spearman",                # correlation method
    sig.level = 0.05,                        # threshold of significance
    cor.vars = c(TMB_NonsynSNP, PD1, PDL1, PDL2, CTLA4),     # a range of variables can be selected  
    cor.vars.names = c("TMB", "PD1 expression", "PDL1 expression", "PDL2 expression", "CTLA4"),
    #subtitle = "Iris dataset by Anderson",
    caption = expression(
        paste(
            italic("Note"),
            ": X denotes correlation non-significant at ",
            italic("p "),
            "< 0.05; adjusted alpha"
        )
    ))

luad_info2 <- luad_info %>% mutate(TMB_Status = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP)[4],
                                         "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP)[2],
                                                        "Low", "Mid"))) 

compareBoxplot(luad_info2 %>% filter(TMB_Status != "Mid"),
               x = "TMB_Status", y = "PD1")
compareBoxplot(luad_info2 %>% filter(TMB_Status != "Mid"),
               x = "TMB_Status", y = "PDL1")
compareBoxplot(luad_info2 %>% filter(TMB_Status != "Mid"),
               x = "TMB_Status", y = "PDL2")
compareBoxplot(luad_info2 %>% filter(TMB_Status != "Mid"),
               x = "TMB_Status", y = "CTLA4")


luad_info2 %>% filter(!is.na(PD1)) %>% group_by(Gender) %>% summarise(N=n())

luad_info2 %>% filter(!is.na(PD1)) %>% tidyr::gather(Genes, mValue, starts_with("PD"), CTLA4) %>% 
    ggbarplot(x="Genes", y="mValue", add="mean_se", color="Gender", palette = "jco",
              position = position_dodge(0.8)) + 
    stat_compare_means(aes(group=Gender), label = "p.format") +
    xlab("") + ylab("mRNA expression (log2)")

# , label = format.pval(..p.adj.., digits = 3)


compareBoxplot(luad_merge,
               x = "Gender", y = "CTLA4", method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PD1", method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PDL1",  method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PDL2",  method = "wilcox.test")




# # test top 10 percent 
# top10_male <- luad_info2 %>% filter(Gender=="Male") %>% 
#     mutate(MS = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP, probs = 0.75) ,"Male-Top10%", "NO")) %>%
#                filter(MS != "NO")
# top10_female <- luad_info2 %>% filter(Gender=="Female") %>% 
#     mutate(MS = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP, probs = 0.75) ,"Female-Top10%", "NO")) %>%
#     filter(MS != "NO")
# 
# top10 <- rbind(top10_female, top10_male)
# compareMutPlot(top10, group1="MS", group2 = "Gender", value = "PD1",  label_name = "p.format", method = "t.test")
# compareMutPlot(top10, group1="MS", group2 = "Gender", value = "PDL1",  label_name = "p.format", method = "t.test")
# compareMutPlot(top10, group1="MS", group2 = "Gender", value = "PDL2",  label_name = "p.format", method = "t.test")
# compareMutPlot(top10, group1="MS", group2 = "Gender", value = "CTLA4",  label_name = "p.format", method = "t.test")


compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PD1",  label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PDL2", label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "CTLA4", label_name = "p.format", method = "t.test")


# compareMutPlot(luad_merge, group2="TMB_Status", group1 = "Gender", value = "PD1",  label_name = "p.format", method = "t.test")
# compareMutPlot(luad_merge, group2="TMB_Status", group1 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")
# compareMutPlot(luad_merge, group2="TMB_Status", group1 = "Gender", value = "PDL2", label_name = "p.format", method = "t.test")



## test PDL1 expression prediction result, use Hellmann etc. and JCO Rizvi etc. datasets
pdl1_hellmann <- sampleInfo_Hellmann %>% filter(PDL1_Expression != "Unknown") %>% 
    mutate(PDL1_Expression = as.numeric(PDL1_Expression))
pdl1_hellmann %>% group_by(Gender) %>% summarise(N=n())
pdl1ROC_Hellmann <- plotROC(pdl1_hellmann, PDL1_Expression, Clinical_Benefit, Gender)
ggsave("pdl1ROC_Hellmann.pdf", plot=pdl1ROC_Hellmann, width=5, height=4)

pdl1_jco <- sampleInfo_JCO_Rizvi %>% filter(!is.na(PDL1_Score)) 
pdl1_jco %>% group_by(Gender) %>% summarise(N=n())
pdl1ROC_jco <- plotROC(pdl1_jco, PDL1_Score, Clinical_Benefit, Gender)
ggsave("pdl1ROC_JCO_Rizvi.pdf", plot=pdl1ROC_jco, width=5, height=4)


## expand to MELA datasets
MELA_Cell <- read_csv("MELA_cell.csv")
MELA_Science <- read_csv("MELA_science.csv")

MELA_Cell <- MELA_Cell %>% mutate(Gender=ifelse(Gender=="F", "Female", "Male"),
                                  Response=ifelse(Response=="R", "DCB", "NDB"))
MELA_Science <- MELA_Science %>% filter(group!="long-survival") %>% 
    mutate(gender=ifelse(gender=="female", "Female", "Male"),
           group=ifelse(group=="response", "DCB", "NDB"))
    
tmb_Mela_Cell <- compareMutPlot(MELA_Cell, group1 = "Gender", 
                                group2 = "Response", 
                                value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
tmb_ROC_Mela_Cell <- plotROC(MELA_Cell, mutation_load, Response, Gender)

tmb_Mela_Sci <- compareMutPlot(MELA_Science, group1 = "gender", group2 = "group",
                               value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
tmb_ROC_Mela_Sci <- plotROC(MELA_Science, mutation_load, group, gender)

tmb_mela <- cowplot::plot_grid(tmb_Mela_Cell, tmb_Mela_Sci, nrow=1, ncol=2, align = "v")
tmbRoc_mela <- cowplot::plot_grid(tmb_ROC_Mela_Cell, tmb_ROC_Mela_Sci, nrow=1, ncol=2, align = "v")

MELA_Cell %>% group_by(Gender) %>% summarise(N=n(), Median=median(mutation_load),
                                             Min=min(mutation_load), Max=max(mutation_load))
MELA_Science %>% group_by(gender) %>% summarise(N=n(), Median=median(mutation_load),
                                                Min=min(mutation_load), Max=max(mutation_load))
ggsave("tmb_MELA.pdf", plot = tmb_mela, width = 7, height = 3)
ggsave("tmbROC_MELA.pdf", plot = tmbRoc_mela, width = 10, height = 4)

MELA_Cell %>% filter(Gender == "Male") %>% summary()


compareBoxplot(MELA_Cell,
               x = "Gender", y = "mutation_load")
compareBoxplot(MELA_Science,
               x = "gender", y = "mutation_load")


# test smoke status
stmb_hellmann <- sampleInfo_Hellmann %>% mutate(Smoke=ifelse(Smoking_History == "never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
# sampleInfo_Hellmann %>% mutate(Smoke=ifelse(Smoking_History == "never", "Never", "Current/former")) %>% 
#     filter(Smoke == "Never") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
ggsave("stmbROC_Hellmann.pdf", plot=stmb_hellmann, width=5, height=4)


stmb_sci_Rizvi <- sampleInfo_Sci_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
# sampleInfo_Sci_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
#     filter(Smoke == "Never") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
ggsave("stmbROC_Sci_Rizvi.pdf", plot=stmb_sci_Rizvi, width=5, height=4)


stmb_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
# sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
#     filter(Smoke == "Never", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
ggsave("stmbROC_JCO_Rizvi_panel410.pdf", plot=stmb_JCO_Rizvi, width=5, height=4)

sampleInfo_Hellmann %>% mutate(Smoke=ifelse(Smoking_History == "never", "Never", "Current/former")) %>% 
     plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Smoke, positive = "DCB")
sampleInfo_Sci_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
   plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Smoke, positive = "DCB")
sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Smoke, positive = "DCB")

p <- rbind(sampleInfo_Forde %>% select(TMB_NonsynSNP, Gender, Clinical_Benefit), 
      sampleInfo_Hellmann %>% select(TMB_NonsynSNP, Gender, Clinical_Benefit),
      sampleInfo_Sci_Rizvi %>%  select(TMB_NonsynSNP, Gender, Clinical_Benefit)) %>% 
          plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit, group = Gender)

ggsave("ROC_combine_NSCLC_wes_.pdf", plot=p, width=5, height=4)

 # t <- plotROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
# 
# plotROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
# plotROC(sampleInfo_Hellmann, TMB_NonsynSNP, Clinical_Benefit, Gender)
# 
# ggsave("C:/Users/wangshx/Desktop/exam_ROC.pdf", plot=t, width=5, height=4)

# roc_Sci_Rizvi <- calcROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
# ggpubr::ggline(data = roc_Sci_Rizvi, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")
# 
# roc_JCO_Rizvi <- calcROC(sampleInfo_JCO_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
# ggline(data = roc_JCO_Rizvi, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")
# 
# roc_Hellmann <- calcROC(sampleInfo_Hellmann, TMB_NonsynSNP, Clinical_Benefit, Gender)
# ggline(data = roc_Hellmann, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")
# 
# roc_Forde <- calcROC(sampleInfo_Forde, TMB_NonsynSNP, Clinical_Benefit, Gender)
# ggline(data = roc_Forde, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")

# rm(list=ls(pattern = "^p|^roc"))
