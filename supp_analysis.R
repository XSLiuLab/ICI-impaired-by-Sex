# Supp analysis
# Thu Feb 14 19:06:34 2019 ------------------------------
library(tidyverse)

# setwd("G:/biodata/immunotherapyDatasets/")
setwd("/Volumes/data/biodata/immunotherapyDatasets/")

# load NSCLC datasets
load("Rdata/sampleInfo_cache.RData")

################# Load data ###############################

# Unify the Tumor Mutation Burden to nonsynonymous mutation/MB
sampleInfo_Hellmann <-  sampleInfo_Hellmann %>% 
    mutate(sTMB = TMB_NonsynSNP / 30)
sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% 
    mutate(sTMB = TMB_NonsynSNP / 30)
# 0.98, 1.06 and 1.22 megabases in 
# the 341-, 410- and 468- gene panels, respectively
sampleInfo_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% 
    mutate(sTMB = ifelse(Gene_Panel == "IMPACT341", TMB_NonsynSNP / 0.98,
                         ifelse(Gene_Panel == "IMPACT410", 
                                TMB_NonsynSNP / 1.06, TMB_NonsynSNP / 1.22)))

sel_cols = c("Tumor_Sample_Barcode", "Age", "Gender", "Smoking_History",
             "PFS_Months", "PFS_Event", "Clinical_Benefit", "sTMB", "dataset")

nsclc =  bind_rows(
    select_at(mutate(sampleInfo_Hellmann, dataset="Hellmann 2018"), sel_cols),
    select_at(mutate(sampleInfo_JCO_Rizvi, dataset="Rizvi 2018"), sel_cols),
    select_at(mutate(sampleInfo_Sci_Rizvi, dataset="Rizvi 2015"), sel_cols)
) %>% mutate(
    Smoking_History = case_when(
        Smoking_History %in% c("never", "Never") ~ "Never",
        TRUE ~ "Current/former"
    )
)

table(nsclc$Gender)
table(nsclc$Smoking_History)
table(nsclc$PFS_Event)
table(nsclc$Clinical_Benefit)

# basic summary based on variable like sex
summary(filter(nsclc, Gender=="Female"))
summary(filter(nsclc, Gender=="Male"))

table(nsclc$Gender, nsclc$Clinical_Benefit)

# Cutoff analysis
library(survival)
library(survminer)

dyHR <- function(data, time="PFS_Months", event="PFS_Event"){
    data <- as.data.frame(data)
    male <- subset(data, Gender == "Male")
    female <- subset(data, Gender == "Female")
    ms_male <- sort(male[, "sTMB"])
    ms_female <- sort(female[, "sTMB"])
    
    fm <- as.formula(paste0("Surv(", time, ", ", event, ") ~ cutoff"))
    s_male <- c()
    for (i in ms_male){
        if (i > 0){
            male <- male %>% mutate(cutoff = factor(ifelse(sTMB > i, "High", "Low"), levels = c("Low", "High")))
            res_male <- summary(coxph(formula = fm, data=male))$coefficients[2]
            s_male <- rbind(s_male, c(i, res_male))
        }
    }
    
    s_female <- c()
    for (i in ms_female){
        if (i > 0){
            female <- female %>% mutate(cutoff = factor(ifelse(sTMB > i, "High", "Low"), levels = c("Low", "High")))
            res_female <- summary(coxph(formula = fm, data=female))$coefficients[2]
            s_female <- rbind(s_female, c(i, res_female))
        }
    }
    
    #res <- list(male=s_male, female=s_female)
    # res <- cbind(s_male, s_female)
    # colnames(res) <- c("x_male", "HR_male", "x_female", "HR_female")
    res <- data.frame(cutoff=c(s_male[,1], s_female[,1]), 
                      HR=c(s_male[,2], s_female[,2]),
                      Gender = c(rep("Male",nrow(s_male)),
                                 rep("Female",nrow(s_female))))
    return(res)
}

cs_nsclc <- dyHR(nsclc)

library(cowplot)
cs_nsclc2 = dplyr::filter(cs_nsclc, cutoff>=1 & cutoff<=20)
p = ggplot(cs_nsclc2, aes(x=cutoff, y=HR, color=Gender)) + 
    geom_line() + xlab("TMB Cutoff") + ylab("Hazard Ratio") +
    coord_cartesian(xlim = c(1,20)) +  
    scale_color_manual(values = c("red", "blue")) + 
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20))
p
save_plot("/Volumes/paper/backup/data/neoQ/HR_vs_TMBcutoff2.pdf", p, base_aspect_ratio = 1.4)
## use cutoff 17 to determine HR

# NSCLC
nsclc <- nsclc %>% mutate(Cutoff = factor(ifelse(sTMB > 16, "High", "Low"),
                                          levels = c("Low", "High") ))
nsclc <- nsclc %>% mutate(Cutoff = factor(ifelse(sTMB > 4, "High", "Low"),
                                          levels = c("Low", "High") ))

summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(nsclc, Gender=="Male")))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(nsclc, Gender=="Female")))

sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% mutate(Cutoff = factor(ifelse(sTMB > 16, 
                                                                               "High", "Low"),
                                                                        levels = c("Low", "High") ))
sampleInfo_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% mutate(Cutoff = factor(ifelse(sTMB > 16, 
                                                                               "High", "Low"),
                                                                        levels = c("Low", "High") ))
sampleInfo_Hellmann <- sampleInfo_Hellmann %>% mutate(Cutoff = factor(ifelse(sTMB > 16, 
                                                                             "High", "Low"),
                                                                      levels = c("Low", "High") ))

summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_Sci_Rizvi, Gender=="Male")))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_Sci_Rizvi, Gender=="Female")))

summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_JCO_Rizvi, Gender=="Male")))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_JCO_Rizvi, Gender=="Female")))

summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_Hellmann, Gender=="Male")))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(sampleInfo_Hellmann, Gender=="Female")))

tableData <- function(dataset){
    male <- dataset %>% 
        filter(Gender=="Male", Clinical_Benefit %in% c("DCB", "NDB")) %>% 
        select(Clinical_Benefit, Cutoff) %>% table()
    female <- dataset %>% 
        filter(Gender=="Female", Clinical_Benefit %in% c("DCB", "NDB")) %>% 
        select(Clinical_Benefit, Cutoff) %>% table()
    res <- list(male =male, female=female)
    return(res)
}

tb_sci_Rizvi <-  tableData(sampleInfo_Sci_Rizvi)
tb_jco_Rizvi <-  tableData(sampleInfo_JCO_Rizvi)
tb_Hellmann  <-  tableData(sampleInfo_Hellmann)

tb_sci_Rizvi$male %>% fisher.test()
tb_sci_Rizvi$female %>% fisher.test()

tb_jco_Rizvi$male %>% fisher.test()
tb_jco_Rizvi$female %>% fisher.test()

tb_Hellmann$male %>% fisher.test()
tb_Hellmann$female %>% fisher.test()


surv_plot <- function(fit){
    p <- ggsurvplot(fit, 
                    font.tickslab=12,
                    legend.title = "",
                    legend.labs = c("TMB-Low", "TMB-High"),
                    palette = c("#0000FF", "#FF0000"),
                    risk.table = TRUE,
                    pval = TRUE,
                    pval.coord = c(1, 0.05),
                    surv.scale = "percent",
                    break.time.by = 5,
                    axes.offset = TRUE,
                    xlab="Months",
                    ylab="Percent progression-free",
                    # ggtheme=theme_wsx,
                    # tables.theme = theme_wsx,
                    xlim = c(0, fit$time+5),
                    legend = c(0.8, 0.8),
                    risk.table.height = 0.35,
                    fontsize=7,
                    pval.size = 6)
    ggpar(
        p,
        font.xtickslab = 15,
        font.ytickslab = 15 
    )
}

### NSCLC
fit1 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=nsclc %>% filter(Gender=="Male"))
fit2 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=nsclc %>% filter(Gender=="Female"))


surv1_1 <- surv_plot(fit1)
surv1_2 <- surv_plot(fit2)

surv1_1
surv1_2


arrange_ggsurvplots(list(surv1_1, surv1_2))
res = arrange_ggsurvplots(list(surv1_1, surv1_2), print = FALSE)

ggsave("/Volumes/paper/backup/data/neoQ/KMplot_TMBcutoff16.pdf", res, width = 12, height = 6)

