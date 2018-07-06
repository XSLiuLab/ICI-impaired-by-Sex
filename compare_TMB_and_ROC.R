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
                                        names(auc)[3], " AUC =", round(auc[3], 3), "\n"),
                          size = 2)
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

# LUAD_TMB3 <- LUAD_TMB2 %>% mutate(TMB_Status = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP)[4],
#                                          "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP)[2],
#                                                       "Low", "Mid"))) %>% filter(TMB_Status != "Mid")

# use 5 as a cutoff

LUAD_TMB3 <- LUAD_TMB2 %>% mutate(TMB = TMB_NonsynSNP / 30 ) %>% 
    mutate(TMB_Status = ifelse(TMB > 4, "High", "Low"))



# LUAD_TMB3 <- LUAD_TMB2 %>% mutate(TMB_Status = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP, 0.75),
#                                                       "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP, 0.25),
#                                                                      "Low", "Mid"))) %>% filter(TMB_Status != "Mid")

LUAD_TMB3 %>% group_by(Gender) %>% 
    summarise(N=n(), N_DOWN=length(which(TMB_Status=="Low")), N_UP=length(which(TMB_Status=="High")))

rm(luad_maf); gc()
# load gene levels
load("C:/Users/wangshx/Desktop/data/summary_of_genes.RData")
geneList = c("AQP11", "ARSD", "BCORL2", "COL9A1", "CST5", "CYorf15A",
             "CYorf15B", "DDX3Y", "EIF1AX", "EIF1AY", "EIF2S3", "FAM150A",
             "FUNDC1", "GEMIN8", "HDHD1A", "KAL1", "KDM5D", "KDM6A",
             "NCRNA00183", "NCRNA00185", "NCRNA00230B", "NLGN4Y",
             "PART1", "PCDH11Y", "PNPLA4", "PRKY", "RPS4X", "RPS4Y1", "SOX11",
             "SRY", "STS", "TBL1Y", "TMSB4Y", "TRAPPC2", "TSIX", "TTTY14",
             "TTTY15", "USP9Y", "UTY", "XIST", "ZFX", "ZFY", "ZNF826", "ZRSR2")
selt_genes <- GeneSummary$Expr[GeneSymbol %in% geneList]
selt_genes <- selt_genes %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue)

# selt_genes <- GeneSummary$Expr[GeneSymbol %in% c("PDCD1", "CD274", "PDCD1LG2", "CTLA4")]
# selt_pros <- GeneSummary$Protein[GeneSymbol %in% c("PDCD1-M-E", "PD-L1-R-V")]
selt_genes <- selt_genes %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue) %>% 
    rename(PD1 = PDCD1, PDL1 = CD274, PDL2 = PDCD1LG2)

selt_pros <-  selt_pros %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue) %>% 
    rename(PD1_pro = `PDCD1-M-E`, PDL1_pro = `PD-L1-R-V`)

# 
# LUAD_TMB3 <- LUAD_TMB2
LUAD_TMB3$Tumor_Sample_Barcode <- paste0(LUAD_TMB3$Tumor_Sample_Barcode, "-01")

luad_merge <- dplyr::left_join(LUAD_TMB3, selt_genes, by="Tumor_Sample_Barcode")
luad_merge2 <- dplyr::left_join(luad_merge, selt_pros, by="Tumor_Sample_Barcode")

luad_merge2 %>% group_by(Gender, TMB_Status) %>% summarize(N_expr=n(), N_pro=length(which(!is.na(PDL1_pro))))

# genelist plot
df = luad_merge %>% tidyr::gather(key = Genes, expression, AQP11:ZRSR2)


ggplot(df, aes(Tumor_Sample_Barcode, Genes, group=Gender)) + 
    geom_tile(aes(fill = expression),colour = "white") + 
    scale_fill_gradient(low = "blue",high = "red") + scale_x_discrete(breaks=NULL)

annotation = df %>% select(Tumor_Sample_Barcode, TMB_Status, Gender) %>% unique.data.frame()
rownames(annotation) = annotation$Tumor_Sample_Barcode
annotation = annotation[, -1]
annotation = annotation %>% rownames_to_column(var="id") %>% 
    arrange(Gender, TMB_Status) %>% 
    column_to_rownames(var="id")
annotation = annotation %>% select(Gender, TMB_Status)

df = luad_merge[, -1:-7]
rownames(df) = luad_merge$Tumor_Sample_Barcode
df = t(df)
df = df[, rownames(annotation)]


pheatmap::pheatmap(df,cluster_cols=FALSE,cluster_rows=TRUE,legend=TRUE,
                   color=colorRampPalette(c("green","black","red"))(1000),
                   border_color=FALSE,fontsize=10,fontsize_row=12,fontsize_col=12,
                   annotation=annotation,annotation_legend=TRUE, show_colnames = FALSE)

compareSignature = function(geneList){
    selt_genes <- GeneSummary$Expr[GeneSymbol %in% geneList] %>% 
        gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue)
    
    # LUAD_TMB3$Tumor_Sample_Barcode <- paste0(LUAD_TMB3$Tumor_Sample_Barcode, "-01")
    geneNames = colnames(selt_genes)[-1]
    
    luad_merge <- dplyr::left_join(LUAD_TMB3, selt_genes, by="Tumor_Sample_Barcode")
    # df = luad_merge %>% tidyr::gather(key = Genes, expression, AQP11:ZRSR2)
    df = luad_merge %>% tidyr::gather_("Genes", "expression", geneNames)
    annotation = df %>% select(Tumor_Sample_Barcode, TMB_Status, Gender) %>% unique.data.frame()
    rownames(annotation) = annotation$Tumor_Sample_Barcode
    annotation = annotation[, -1]
    annotation = annotation %>% rownames_to_column(var="id") %>% 
        arrange(Gender, TMB_Status) %>% 
        column_to_rownames(var="id")
    annotation = annotation %>% select(Gender, TMB_Status)
    
    df = luad_merge[, -1:-7]
    rownames(df) = luad_merge$Tumor_Sample_Barcode
    df = t(df)
    df = df[, rownames(annotation)]
    
    
    pheatmap::pheatmap(df,cluster_cols=FALSE,cluster_rows=TRUE,legend=TRUE,
                       color=colorRampPalette(c("green","black","red"))(1000),
                       border_color=FALSE,fontsize=10,fontsize_row=12,fontsize_col=12,
                       annotation=annotation,annotation_legend=TRUE, show_colnames = FALSE)
}

geneList = read_csv("G:/biodata/datasets/Science_Davoli/Immune_Marker_GeneList.csv")
compareSignature(geneList %>% filter(Cell_Type=="CD4.mature") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="CD8.effector") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="NK.cells") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="B.cells") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="T.reg") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="Dendritic") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="T.reg") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="CD8.effector.NK.cells") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="Macrophages") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="Macrophages.M2") %>% select(Gene_Name) %>% unlist)
compareSignature(geneList %>% filter(Cell_Type=="Macrophages.M1") %>% select(Gene_Name) %>% unlist)

load("C:/Users/wangshx/Desktop/data/TCGA_LUAD_TIL.RData")
compareTIL = function(TIL_df){
    TIL_df = TIL$timer 
    TIL_df = TIL$quanTIseq
    TIL_df = left_join(LUAD_TMB3, TIL_df, by="Tumor_Sample_Barcode")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "B_cell",  label_name = "p.format", method = "wilcox.test")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "T_cell.CD8",  label_name = "p.format", method = "wilcox.test")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "T_cell.CD4",  label_name = "p.format", method = "wilcox.test")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "Neutrophil",  label_name = "p.format", method = "wilcox.test")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "Macrophage",  label_name = "p.format", method = "wilcox.test")
    compareMutPlot(TIL_df, group1="Gender", group2 = "TMB_Status", value = "DC",  label_name = "p.format", method = "wilcox.test")
}


compareBoxplot(luad_merge,
               x = "Gender", y = "CTLA4", method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PD1", method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PDL1",  method = "wilcox.test")
compareBoxplot(luad_merge,
               x = "Gender", y = "PDL2",  method = "wilcox.test")



compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PD1",  label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "PDL2", label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge, group1="TMB_Status", group2 = "Gender", value = "CTLA4", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge2, group1="TMB_Status", group2 = "Gender", value = "PDL1_pro", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge, group2 ="TMB_Status", group1 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")
compareMutPlot(luad_merge2, group2 ="TMB_Status", group1 = "Gender", value = "PDL1_pro", label_name = "p.format", method = "t.test")
#######
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
MELA_NEJM <- read_csv("MELA_nejm.csv")

MELA_Cell <- MELA_Cell %>% mutate(Gender=ifelse(Gender=="F", "Female", "Male"),
                                  Response=ifelse(Response=="R", "DCB", "NDB")) %>% 
    rename(Clinical_Benefit=Response)
MELA_Science <- MELA_Science %>% filter(group!="long-survival") %>% 
    mutate(gender=ifelse(gender=="female", "Female", "Male"),
           group=ifelse(group=="response", "DCB", "NDB")) %>% 
    rename(Gender=gender, Clinical_Benefit=group)

MELA_NEJM <- MELA_NEJM %>% rename(Gender=gender, mutation_load=Mutation_load) %>% 
    mutate(Gender=ifelse(Gender=="F", "Female", "Male"),
           Clinical_Benefit=ifelse(Clinical_Benefit=="LB", "DCB", "NDB"))

tmb_MELA_NEJM <- compareMutPlot(MELA_NEJM, group1 = "Gender", 
               group2 = "Clinical_Benefit", 
               value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
tmb_ROC_NEJM <- plotROC(MELA_NEJM, mutation_load, Clinical_Benefit, Gender, positive = "LB")

tmb_Mela_Cell <- compareMutPlot(MELA_Cell, group1 = "Gender", 
                                group2 = "Response", 
                                value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
tmb_ROC_Mela_Cell <- plotROC(MELA_Cell, mutation_load, Response, Gender)

tmb_Mela_Sci <- compareMutPlot(MELA_Science, group1 = "Gender", group2 = "group",
                               value = "mutation_load") + xlab("") +ylab("Tumor Mutation Burden")
tmb_ROC_Mela_Sci <- plotROC(MELA_Science, mutation_load, group, Gender)

tmb_mela <- cowplot::plot_grid(tmb_Mela_Cell, tmb_Mela_Sci, tmb_MELA_NEJM, nrow=1, ncol=3, align = "v")
tmbRoc_mela <- cowplot::plot_grid(tmb_ROC_Mela_Cell, tmb_ROC_Mela_Sci, tmb_ROC_NEJM, nrow=1, ncol=3, align = "v")

MELA_Cell %>% group_by(Gender) %>% summarise(N=n(), Median=median(mutation_load),
                                             Min=min(mutation_load), Max=max(mutation_load))
MELA_Science %>% group_by(gender) %>% summarise(N=n(), Median=median(mutation_load),
                                                Min=min(mutation_load), Max=max(mutation_load))
ggsave("tmb_MELA.pdf", plot = tmb_mela, width = 10, height = 3)
ggsave("tmbROC_MELA.pdf", plot = tmbRoc_mela, width = 14, height = 4)

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
t <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>%
    filter(Smoke == "Never", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = TMB_NonsynSNP, target = Clinical_Benefit,group = Gender, positive = "DCB")
ggsave("stmbROC_JCO_Rizvi_panel410.pdf", plot=stmb_JCO_Rizvi, width=5, height=4)
ggsave("stmbROC_JCO_Rizvi_panel410_NeverSmoker.pdf", plot=t, width=5, height=4)

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



all_f <- c(0.912, 0.715, 0.722, 0.67, 0.267, 0.817, 0.743)
all_m <- c(0.633, 0.629, 0.391, 0.565, 0.759, 0.604, 0.772)

t.test(all_f, all_m)
wilcox.test(all_f, all_m)

nsclc_f <- c(0.912, 0.715, 0.722, 0.67)
nsclc_m <- c(0.633, 0.629, 0.391, 0.565)
t.test(nsclc_f, nsclc_m, paired = TRUE) 
wilcox.test(nsclc_f, nsclc_m, paired=TRUE)

# remove cell paper
all_f2 <- c(0.912, 0.715, 0.722, 0.67, 0.817, 0.743)
all_m2 <- c(0.633, 0.629, 0.391, 0.565, 0.604, 0.772)
t.test(all_f2, all_m2, paired = TRUE)
wilcox.test(all_f2, all_m2, paired = TRUE)

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
