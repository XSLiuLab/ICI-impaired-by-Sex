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
CV.nonsyn$JCO_Rizvi <- calcCV(sampleInfo_JCO_Rizvi, target = TMB_NonsynSNP, Gender)
CV.nonsyn$Hellmann  <- calcCV(sampleInfo_Hellmann, target = TMB_NonsynSNP, Gender)
CV.nonsyn$Forde     <- calcCV(sampleInfo_Forde, target = TMB_NonsynSNP, Gender)

CV.nonsyn$LUAD      <- calcCV(LUAD_TMB2, target =TMB_NonsynSNP, Gender, TNM )

calcCV(sampleInfo_Sci_Rizvi, target = NeoCounts, Gender)
calcCV(sampleInfo_Hellmann, target = NeoCounts, Gender)
calcCV(sampleInfo_Forde, target = NeoCounts, Gender)

## ROC plot and analysis
calcROC <- function(.data, predict_var, target, group_var, positive="DCB"){
    # predic_var here must be a numeric value
    require(tidyverse)
    predict_var <- enquo(predict_var)
    target <- enquo(target)
    group_var <- enquo(group_var)
    
    groups <- .data %>% filter(!is.na(!! predict_var)) %>% select(!! group_var) %>% 
        unlist() %>% table() %>% names()
    
    total_res <- list()
    # process groups one by one
    j <- 1
    for (i in groups){
        df <- list()
        df <- .data %>% filter(!is.na(!! predict_var), !! group_var == i) %>%
            arrange(desc(!! predict_var)) %>% 
            mutate(isPositive = ifelse(!! target == positive, 1, 0))
        
        # select a threshold, calculate true positive and false positive value
        ths <- df %>% select(!! predict_var) %>% unlist
        
        mat <- base::sapply(ths, function(th){
            # true positive
            tp <- df %>% filter(!! predict_var >= th) %>% filter(isPositive == 1) %>% nrow
            # false positive
            fp <- df %>% filter(!! predict_var >= th) %>% filter(isPositive == 0) %>% nrow
            # true negative
            tn <- df %>% filter(!! predict_var < th) %>% filter(isPositive == 0) %>% nrow
            # false negative
            fn <- df %>% filter(!! predict_var < th) %>% filter(isPositive == 1) %>% nrow
            
            # true positive rate
            tpr <- tp / (tp + fn)
            # false positive rate
            fpr <- fp / (fp + tn)
            
            return(c(tp, fp, tn, fn, tpr, fpr))
            # combine
        })
        
        res <- t(mat)
        res <- data.frame(res)
        # fake a (0, 0) point
        res <- rbind(c(rep(NA, 4), 0, 0), res)
        colnames(res) <- c("tp", "fp", "tn", "fn", "tpr", "fpr")
        res$Group <- i
        total_res[[j]] <- res
        j <- j + 1
    }
    
    dat <- base::Reduce(rbind, total_res)
    return(dat)
}

roc_Sci_Rizvi <- calcROC(sampleInfo_Sci_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggline(data = roc_Sci_Rizvi, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")

roc_JCO_Rizvi <- calcROC(sampleInfo_JCO_Rizvi, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggline(data = roc_JCO_Rizvi, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")

roc_Hellmann <- calcROC(sampleInfo_Hellmann, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggline(data = roc_Hellmann, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")

roc_Forde <- calcROC(sampleInfo_Forde, TMB_NonsynSNP, Clinical_Benefit, Gender)
ggline(data = roc_Forde, x = "fpr", y = "tpr", linetype = "Group", shape = "Group")
