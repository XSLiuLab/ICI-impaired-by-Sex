#--------------------------------------------
# New analysis for new data from van ellen
#-------------------------------------------
library(tidyverse)
library(readxl)
library(data.table)

clinical = read_xlsx("/Volumes/data/biodata/datasets/NatGen2018_Van/van2018_clinical.xlsx")
mutation = fread("/Volumes/data/biodata/datasets/NatGen2018_Van/van2018_mutation.txt")

summary_mutation = mutation[Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation"),
                            .(Nonsyn=.N), by = .(pair_id)]
summary_mutation = as.tibble(summary_mutation)

clinical = full_join(clinical, summary_mutation)
clinical$TMB = clinical$Nonsyn / 30

write_csv(clinical, path = "../NatGen_2018_Van_Data.csv")

#--------- select new patients and analyzing

## clean data
clean_cli = clinical %>% 
    filter(cancer_type == "Lung", nchar(Tumor_Sample_Barcode) >= 20) 

clean_cli %>% 
    mutate(Response = case_when(
        va_response == "clinical benefit" ~ "R",
        va_response == "no clinical benefit" ~ "NR",
        TRUE ~ NA_character_
    ),
    Gender = case_when(
        sex == "FEMALE" ~ "Female",
        sex == "MALE" ~ "Male",
        TRUE ~ NA_character_
    )) -> clean_cli

clean_cli2 = clinical %>% 
    filter(cancer_type == "Lung") 

clean_cli2 %>% 
    mutate(Response = case_when(
        va_response == "clinical benefit" ~ "R",
        va_response == "no clinical benefit" ~ "NR",
        TRUE ~ NA_character_
    ),
    Gender = case_when(
        sex == "FEMALE" ~ "Female",
        sex == "MALE" ~ "Male",
        TRUE ~ NA_character_
    )) -> clean_cli2

#------------
clean_cli3 = clinical %>% 
    filter(cancer_type == "Melanoma") 

clean_cli3 %>% 
    mutate(Response = case_when(
        va_response == "clinical benefit" ~ "R",
        va_response == "no clinical benefit" ~ "NR",
        TRUE ~ NA_character_
    ),
    Gender = case_when(
        sex == "FEMALE" ~ "Female",
        sex == "MALE" ~ "Male",
        TRUE ~ NA_character_
    )) -> clean_cli3

#-----------------------------------
clean_cli4 = clinical %>% 
    filter(cancer_type == "Bladder") 

clean_cli4 %>% 
    mutate(Response = case_when(
        va_response == "clinical benefit" ~ "R",
        va_response == "no clinical benefit" ~ "NR",
        TRUE ~ NA_character_
    ),
    Gender = case_when(
        sex == "FEMALE" ~ "Female",
        sex == "MALE" ~ "Male",
        TRUE ~ NA_character_
    )) -> clean_cli4


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


my_palette <- c("black", "red", "blue")
add_modify <- theme(axis.line = element_line(colour = 'black', size = 1.0),
                    axis.ticks = element_line(colour = "black", size = 1.0))
#axis.text.y = element_blank(), 
#axis.text.x = element_blank()) 

plotROC(clean_cli, predict_col = TMB, target = Response, group = Gender, positive = "R") + 
    scale_color_manual(values = my_palette) + add_modify -> p1

ggsave("../CompareROC-NatGen-Lung.pdf", plot = p1, width = 4, height = 3)

plotROC(clean_cli2, predict_col = TMB, target = Response, group = Gender, positive = "R") + 
    scale_color_manual(values = my_palette) + add_modify

plotROC(clean_cli3, predict_col = TMB, target = Response, group = Gender, positive = "R") + 
    scale_color_manual(values = my_palette) + add_modify

plotROC(clean_cli4, predict_col = TMB, target = Response, group = Gender, positive = "R") + 
    scale_color_manual(values = my_palette) + add_modify

# HR
require(survival)
require(survminer)

surv_plot <- function(fit){
    p <- ggsurvplot(fit, 
                    font.tickslab=14,
                    # legend="right", 
                    legend.title = "",
                    legend.labs = c("Low", "High"),
                    palette = c("#0000FF", "#FF0000"),
                    risk.table = TRUE,
                    pval = TRUE,
                    pval.coord = c(1, 0.05),
                    surv.scale = "percent",
                    break.time.by = 5,
                    axes.offset = FALSE,
                    xlab="",
                    ylab="",
                    # ggtheme=theme_wsx,
                    # tables.theme = theme_wsx,
                    xlim = c(0, fit$time+5),
                    legend="none",
                    risk.table.height = 0.35,
                    fontsize=7,
                    pval.size = 6)
    ggpar(
        p,
        font.xtickslab = 20,
        font.ytickslab = 20 
    )
}

# use unified cutoff
clean_cli <- clean_cli %>% mutate(Cutoff = factor(ifelse(TMB > 4, "High", "Low"),
                                          levels = c("Low", "High")), 
                                  PFS_Months = pfs_days / 30) %>% 
    rename(PFS_Event = pfs_censor)
clean_cli$PFS_Event = as.integer(clean_cli$PFS_Event)

clean_cli3_u <- clean_cli3 %>% mutate(Cutoff = factor(ifelse(TMB > 8, "High", "Low"),
                                                  levels = c("Low", "High")), 
                                  PFS_Months = pfs_days / 30) %>% 
    rename(PFS_Event = pfs_censor)
clean_cli3_u$PFS_Event = as.integer(clean_cli3_u$PFS_Event)

fit3 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli3_u, Gender=="Male"))
surv_plot(fit3)
fit4 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli3_u, Gender=="Female"))
surv_plot(fit4)
fit5 =survfit(Surv(os_days, os_censor==0) ~ Cutoff,
              data = subset(clean_cli3_u, Gender=="Male"))
ggsurvplot(fit5, pval = TRUE)
fit6 =survfit(Surv(os_days, os_censor==0) ~ Cutoff,
              data = subset(clean_cli3_u, Gender=="Female"))
ggsurvplot(fit6, pval = TRUE)


fit1 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli, Gender=="Male"))
surv1 = surv_plot(fit1)
ggsurvplot(fit1, pval = TRUE)

fit2 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli, Gender=="Female"))
ggsurvplot(fit2, pval = TRUE) 
surv2 = surv_plot(fit2)

ggsave("../M-KMplot_NatGen_Male.pdf", plot = print(surv1, newpage = FALSE),
       width = 7, height = 5)
ggsave("../M-KMplot_NatGen_Female.pdf", plot = print(surv2, newpage = FALSE),
       width = 7, height = 5)

coxph(Surv(PFS_Months, PFS_Event==0) ~ Cutoff, data = subset(clean_cli, Gender=="Male"))


# separate Male and Female
clean_cli_sp = clean_cli %>% 
    mutate(Cutoff = factor(case_when(
        Gender == "Female" & TMB > median(.$TMB[.$Gender=="Female"], na.rm = TRUE) ~ "High",
        Gender == "Female" & TMB <= median(.$TMB[.$Gender=="Female"], na.rm = TRUE) ~ "Low",
        Gender == "Male" & TMB > median(.$TMB[.$Gender=="Male"], na.rm = TRUE) ~ "High",
        Gender == "Male" & TMB <= median(.$TMB[.$Gender=="Male"], na.rm = TRUE) ~ "Low"
    ),
        levels = c("Low", "High")), 
    PFS_Months = pfs_days / 30) %>% 
    rename(PFS_Event = pfs_censor)
clean_cli_sp$PFS_Event = as.integer(clean_cli_sp$PFS_Event)

fit1 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli_sp, Gender=="Male"))


surv1 = surv_plot(fit1)
ggsurvplot(fit1, pval = TRUE)

fit2 = survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
              data = subset(clean_cli_sp, Gender=="Female"))
ggsurvplot(fit2, pval = TRUE) 
surv2 = surv_plot(fit2)

ggsave("../M-KMplot_NatGen_Male.pdf", plot = print(surv1, newpage = FALSE),
       width = 7, height = 5)
ggsave("../M-KMplot_NatGen_Female.pdf", plot = print(surv2, newpage = FALSE),
       width = 7, height = 5)


# fit5 =survfit(Surv(PFS_Months, PFS_Event==0) ~ Cutoff,
#               data = clean_cli)
# ggsurvplot(fit5, pval = TRUE) 
# #-----------
# fit3 =survfit(Surv(os_days, os_censor==0) ~ Cutoff,
#               data = subset(clean_cli, Gender=="Male"))
# ggsurvplot(fit3, pval = TRUE)
# 
# fit4 =survfit(Surv(os_days, os_censor==0) ~ Cutoff,
#               data = subset(clean_cli, Gender=="Female"))
# ggsurvplot(fit4, pval = TRUE) 
