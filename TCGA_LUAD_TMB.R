p_load(tidyverse, ggpubr)

# compare samples using boxplot and add significant levels
compareBoxplot <- function(df, x=NULL, y=NULL, label_name=c("p.format", "p.signif"), 
                           method=c("wilcox.test", "t.test", "anova", "kruskai.test")){
    label_name <- match.arg(label_name)
    method <- match.arg(method)
    df <- as.data.frame(df)
    name_sort <- names(table(df[, x]))
    df[, x] <- factor(df[, x], levels = name_sort)
    p <- ggboxplot(df, x=x, y=y, ggtheme = theme_pubr(base_size = 8))
    p + stat_compare_means(label = label_name, label.x.npc = "center", method = method) + theme(plot.title = element_text(hjust = 0.5))
}

# compare mutation profile in two-level group variable
compareMutPlot <- function(dat, group1="Gender", group2="Clinical_Benefit", 
                           value="TMB_Total", label_name="p.format", method = "wilcox.test"){
    require(ggpubr)
    dat <- as.data.frame(dat)
    #my_comparisons <- combn(names(table(dat[, group2])), 2, simplify = FALSE)
    #my_comparisons  <- list(c("DCB", "NDB"))
    name_sort <- names(table(dat[, group1]))
    dat[, group1] <- factor(dat[, group1], levels = name_sort)
    p <- ggboxplot(dat, x = group1, y = value,
                   color = group2, palette = c("red", "blue"),
                   add = "jitter", shape = group2, font.label = list(size=6), 
                   add.params = list(size=2), 
                   ggtheme = theme_pubr()) + theme(plot.title = element_text(hjust = 0.5))
    
    p + stat_compare_means(aes_string(group=group2), label = label_name, 
                           method = method)        # Add pairwise comparisons p-value
    # stat_compare_means(label.x.npc = "center", 
    #                    label.y.npc = "top" )                   # Add global p-value
    
    # res <- summarySE(data = dat, measurevar = value, groupvars = c(group1, group2))
    # ggplot(res, aes_string(x=group1, y=value, fill=group2)) + 
    #      geom_bar(position = position_dodge() ,stat = "identity") +
    #      geom_errorbar(aes_string(ymin=paste0(value, "-se"), ymax=paste0(value, "+se")),
    #                    width=.2,
    #                    position = position_dodge(.9)) + 
    #      theme_bw()  + scale_fill_npg()
}


load("C:/Users/wangshx/Desktop/tcga.Rdata")

LUAD_TMB3 <- LUAD_TMB3 %>% mutate(TMB_Status2 = ifelse(TMB_NonsynSNP >= quantile(TMB_NonsynSNP)[4],
                                         "High", ifelse(TMB_NonsynSNP <= quantile(TMB_NonsynSNP)[2],
                                                      "Low", "Mid")))

selt_genes <- GeneSummary$Expr[GeneSymbol %in% c("PDCD1", "CD274", "PDCD1LG2", "CTLA4")]
selt_pros <- GeneSummary$Protein[GeneSymbol %in% c("PDCD1-M-E", "PD-L1-R-V")]
selt_genes <- selt_genes %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue) %>% 
    rename(PD1 = PDCD1, PDL1 = CD274, PDL2 = PDCD1LG2)

selt_pros <-  selt_pros %>% 
    gather(Tumor_Sample_Barcode, mvalue, starts_with("TCGA")) %>% spread(GeneSymbol, mvalue) %>% 
    rename(PD1_pro = `PDCD1-M-E`, PDL1_pro = `PD-L1-R-V`)

luad_merge <- dplyr::left_join(LUAD_TMB3, selt_genes, by="Tumor_Sample_Barcode")
luad_merge2 <- dplyr::left_join(luad_merge, selt_pros, by="Tumor_Sample_Barcode")

compareMutPlot(luad_merge, group2="TMB_Status", group1 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge2, group2="TMB_Status", group1 = "Gender", value = "PDL1_pro", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "PDL1", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge2 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "PDL1_pro", label_name = "p.format", method = "t.test")

#############
luad_merge %>% filter(TMB_Status2 != "Mid") %>% group_by(Gender, TMB_Status2) %>% 
    summarise(N=n(), Mean=mean(PDL1), SD=sd(PDL1))
luad_merge2 %>% filter(TMB_Status2 != "Mid") %>% group_by(Gender, TMB_Status2) %>% 
    summarise(N=sum(!is.na(PDL1_pro)), Mean=mean(PDL1_pro, na.rm=T), SD=sd(PDL1_pro, na.rm=T))


########### Immune Score
# load immune score genes
immuneScore_gene = read_csv(file="G:/biodata/GeneList/Immune_Score.csv", skip = 1)
all(immuneScore_gene$Gene_Name %in% GeneSummary$Expr$GeneSymbol)
immuneGene_expr = GeneSummary$Expr[GeneSymbol %in% immuneScore_gene$Gene_Name]

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
    if(any(x < 0, na.rm = TRUE)){
        return(NaN)
    }
    if(zero.propagate){
        if(any(x == 0, na.rm = TRUE)){
            return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
    } else {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
}

immuneGene_score = apply(immuneGene_expr[, -1], 2,  gm_mean)
immuneGene_score = tibble(Tumor_Sample_Barcode=names(immuneGene_score), imScore=immuneGene_score)

luad_merge3 = dplyr::left_join(luad_merge, immuneGene_score, by="Tumor_Sample_Barcode")

compareMutPlot(luad_merge3, group2="TMB_Status", group1 = "Gender", value = "imScore", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge3 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "imScore", label_name = "p.format", method = "t.test")

## Load TIL data
timer = read_tsv("G:/biodata/TCGA/TCGA_sixcell_population_TIMER.txt")
timer$Tumor_Sample_Barcode = substr(timer$sample, 1, 15)
luad_merge4 = timer %>% filter(Tumor_Sample_Barcode %in% luad_merge3$Tumor_Sample_Barcode) %>% 
    full_join(luad_merge3, by="Tumor_Sample_Barcode") %>% select(- sample)

# cutoff 5
compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "T_cell.CD8", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "T_cell.CD4", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "B_cell", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "Neutrophil", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "Macrophage", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4, group2="TMB_Status", group1 = "Gender", value = "DC", label_name = "p.format", method = "wilcox.test")

# 25% up and down
compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "T_cell.CD8", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "T_cell.CD4", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "B_cell", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "Neutrophil", label_name = "p.format", method = "wilcox.test")

compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "Macrophage", label_name = "p.format", method = "t.test")

compareMutPlot(luad_merge4 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "DC", label_name = "p.format", method = "wilcox.test")


## Load TCIA
### Load TCIA data
tcia = read_tsv("G:/biodata/TCGA/TCGA-all-cellTypeFractionsAll-20180702.tsv")
tcia = unique(tcia)
tcia$quanTIseq_lsei_TIL10 = as.numeric(tcia$quanTIseq_lsei_TIL10)
tcia$cibersort_LM22 = as.numeric(tcia$cibersort_LM22)
tcia$Tumor_Sample_Barcode = paste0(tcia$patientBarcode, "-01")
tcia_cibersort = tcia %>% select(-patientBarcode, -quanTIseq_lsei_TIL10) %>%
    filter(disease == "LUAD") %>% select(-disease) %>% filter(!is.na(cibersort_LM22)) %>% 
    spread(cell_type, cibersort_LM22) %>% rename(T_cell_CD8 = `CD8 T cells`)

luad_merge5 = tcia_cibersort %>% filter(Tumor_Sample_Barcode %in% luad_merge3$Tumor_Sample_Barcode) %>% 
    full_join(luad_merge3, by="Tumor_Sample_Barcode") 

compareMutPlot(luad_merge5, group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test") 
compareMutPlot(luad_merge5 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test")

tcia_quanTI = tcia %>% select(-patientBarcode, -cibersort_LM22) %>%
    filter(disease == "LUAD") %>% select(-disease) %>% filter(!is.na(quanTIseq_lsei_TIL10)) %>% 
    spread(cell_type, quanTIseq_lsei_TIL10) %>% rename(T_cell_CD8 = `CD8 T cells`)
luad_merge6 = tcia_quanTI %>% filter(Tumor_Sample_Barcode %in% luad_merge3$Tumor_Sample_Barcode) %>% 
    full_join(luad_merge3, by="Tumor_Sample_Barcode") 

compareMutPlot(luad_merge6, group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test") 
compareMutPlot(luad_merge5 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test")


### Load TIL cell score
TIL_score = read_csv("G:/biodata/TCGA/TCGA_TIL_cell_Score.csv", skip = 1)
TIL_score$Tumor_Sample_Barcode = substr(TIL_score$Tumor_Sample_Barcode, 1, 15)
TIL_score$Tumor_Sample_Barcode = gsub("\\.", '-', TIL_score$Tumor_Sample_Barcode)

luad_merge7 = TIL_score %>% filter(Tumor_Sample_Barcode %in% luad_merge3$Tumor_Sample_Barcode) %>% 
    full_join(luad_merge3, by="Tumor_Sample_Barcode") %>% rename(T_cell_CD8 = `CD8 T cells`)

compareMutPlot(luad_merge7, group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test") 

compareMutPlot(luad_merge7 %>% filter(TMB_Status2!="Mid"), group2="TMB_Status", group1 = "Gender", value = "T_cell_CD8", label_name = "p.format", method = "wilcox.test")


###### Meta analysis
p_load(altmeta)
data("aex")
