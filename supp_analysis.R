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
ggplot(cs_nsclc, aes(x=cutoff, y=HR, color=Gender)) + 
    geom_point() + geom_line() + xlim(c(20,0)) + xlab("TMB Cutoff") 

## use cutoff 17 to determine HR

# NSCLC
nsclc <- nsclc %>% mutate(Cutoff = factor(ifelse(sTMB > 16, "High", "Low"),
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

### NSCLC
fit1 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_Sci_Rizvi %>% filter(Gender=="Male"))
fit2 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_Sci_Rizvi %>% filter(Gender=="Female"))
fit3 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_JCO_Rizvi %>% filter(Gender=="Male"))
fit4 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_JCO_Rizvi %>% filter(Gender=="Female"))
fit5 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_Hellmann %>% filter(Gender=="Male"))
fit6 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff, 
                data=sampleInfo_Hellmann %>% filter(Gender=="Female"))

surv1_1 <- surv_plot(fit1)
surv1_2 <- surv_plot(fit2)
surv1_3 <- surv_plot(fit3)
surv1_4 <- surv_plot(fit4)
surv1_5 <- surv_plot(fit5)
surv1_6 <- surv_plot(fit6)


p4_1 <- plotROC(nsclc, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-NSCLC-ALL.pdf", plot = p5_1, width = 4, height = 3)

