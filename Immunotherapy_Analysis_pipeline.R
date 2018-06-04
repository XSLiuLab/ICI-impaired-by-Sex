# tidy and re-construct analyses for immunotherapy on TMB and ROC etc..
# Shixiang Wang

library(tidyverse)

setwd("G:/biodata/immunotherapyDatasets/")

# load NSCLC datasets
load("Rdata/sampleInfo_cache.RData")

################# Load data ###############################
# load MELA datasets and transform them
MELA_NEJM <- read_csv("MELA/MELA_NEJM.csv")
MELA_Cell <- read_csv("MELA/MELA_Cell.csv")
MELA_Science <- read_csv("MELA/MELA_Science.csv")

MELA_NEJM <- MELA_NEJM %>% 
    mutate(Gender = ifelse(Gender == "F", "Female", "Male"), 
           OS = OS_Year * 12, 
           Event = OS_Status)

# MELA_Cell <- MELA_Cell %>% 
#     mutate(Gender = ifelse(Gender == "F", "Female", "Male"),
#            OS = OS_Day / 30, 
#            Event = ifelse(OS_Status=="Dead", 1, 0))

MELA_Science <- MELA_Science %>% 
    mutate(Gender = ifelse(Gender == "female", "Female", "Male"),
           OS = OS_Day / 30, 
           Event = OS_Status)

MELA_Science <- filter(MELA_Science, group != "long-survival")

############### Unify data ##########################
# Unify the Tumor Mutation Burden to nonsynonymous mutation/MB

sampleInfo_Forde    <-  sampleInfo_Forde %>% 
                        mutate(sTMB = TMB_NonsynSNP / 30, 
                        sNeo = NeoCounts / 30)
sampleInfo_Hellmann <-  sampleInfo_Hellmann %>% 
                        mutate(sTMB = TMB_NonsynSNP / 30, 
                        sNeo = NeoCounts / 30)
sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% 
                        mutate(sTMB = TMB_NonsynSNP / 30, 
                        sNeo = NeoCounts / 30)
# 0.98, 1.06 and 1.22 megabases in the 341-, 410- and 468- gene panels, respectively
sampleInfo_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% 
                        mutate(sTMB = ifelse(Gene_Panel == "IMPACT341", TMB_NonsynSNP / 0.98,
                                             ifelse(Gene_Panel == "IMPACT410", 
                                                    TMB_NonsynSNP / 1.06, TMB_NonsynSNP / 1.22)))

MELA_NEJM <- MELA_NEJM %>% mutate(sTMB = TMB / 30)
MELA_Science <- MELA_Science %>% mutate(sTMB = TMB / 30)

MELA_NEJM <- MELA_NEJM %>% mutate(Clinical_Benefit = ifelse(RESP == "LB", "DCB", "NDB"))
MELA_Science <- MELA_Science %>% mutate(Clinical_Benefit = ifelse(group == "response", "DCB", "NDB"))

# save(MELA_NEJM, MELA_Science, sampleInfo_Forde, sampleInfo_Hellmann, sampleInfo_JCO_Rizvi,
#      sampleInfo_Sci_Rizvi, file="Rdata/cache_Unify_sampleinfo.RData")

############### Compare sTMB (standardized Tumor Mutation Burden) ######################
# library(RColorBrewer)

p1_1 <- compareMutPlot(subset(sampleInfo_Sci_Rizvi, Clinical_Benefit %in% c("DCB", "NDB")), value = "sTMB") 
p1_2 <- compareMutPlot(sampleInfo_Hellmann, value="sTMB")
p1_3 <- compareMutPlot(sampleInfo_JCO_Rizvi, value="sTMB")
p1_4 <- compareMutPlot(sampleInfo_Forde, value="sTMB")
p1_5 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT341"), value = "sTMB") 
p1_6 <- compareMutPlot(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT410"), value = "sTMB") 

ggsave("Figures/CompareMut-Science-Rizvi.pdf", plot = p1_1, width = 4, height = 3)
ggsave("Figures/CompareMut-Hellmann.pdf", plot = p1_2, width = 4, height = 3)
ggsave("Figures/CompareMut-JCO-Rizvi.pdf", plot = p1_3, width = 4, height = 3)
ggsave("Figures/CompareMut-Forde.pdf", plot = p1_4, width = 4, height = 3)
ggsave("Figures/CompareMut-JCO-Rizvi-341.pdf", plot = p1_5, width = 4, height = 3)
ggsave("Figures/CompareMut-JCO-Rizvi-410.pdf", plot = p1_6, width = 4, height = 3)

p1_7 <- compareMutPlot(MELA_NEJM, value="sTMB")
p1_8 <- compareMutPlot(MELA_Science, value="sTMB")
ggsave("Figures/CompareMut-MELA-NEJM.pdf", plot = p1_7, width = 4, height = 3)
ggsave("Figures/CompareMut-MELA-Science.pdf", plot = p1_8, width = 4, height = 3)

############ Compare ROC ################################
my_palette <- c("black", "red", "blue")
add_modify <- theme(axis.line = element_line(colour = 'black', size = 1.0),
                    axis.ticks = element_line(colour = "black", size = 1.0))
                    #axis.text.y = element_blank(), 
                    #axis.text.x = element_blank()) 

p2_1 <- plotROC(sampleInfo_Sci_Rizvi, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p2_2 <- plotROC(sampleInfo_Hellmann, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p2_3 <- plotROC(sampleInfo_JCO_Rizvi, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p2_4 <- plotROC(sampleInfo_Forde, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p2_5 <- plotROC(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT341"), sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p2_6 <- plotROC(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT410"), sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-Science-Rizvi.pdf", plot = p2_1, width = 4, height = 3)
ggsave("Figures/CompareROC-Hellmann.pdf", plot = p2_2, width = 4, height = 3)
ggsave("Figures/CompareROC-JCO-Rizvi.pdf", plot = p2_3, width = 4, height = 3)
ggsave("Figures/CompareROC-Forde.pdf", plot = p2_4, width = 4, height = 3)
ggsave("Figures/CompareROC-JCO-Rizvi-341.pdf", plot = p2_5, width = 4, height = 3)
ggsave("Figures/CompareROC-JCO-Rizvi-410.pdf", plot = p2_6, width = 4, height = 3)

######  PDL1 predict ROC
pdl1_hellmann <- sampleInfo_Hellmann %>% filter(PDL1_Expression != "Unknown") %>% 
    mutate(PDL1_Expression = as.numeric(PDL1_Expression))
pdl1_hellmann %>% group_by(Gender) %>% summarise(N=n())
pdl1ROC_Hellmann <- plotROC(pdl1_hellmann, PDL1_Expression, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/pdl1ROC_Hellmann.pdf", plot=pdl1ROC_Hellmann, width=4, height=3)

pdl1_jco <- sampleInfo_JCO_Rizvi %>% filter(!is.na(PDL1_Score)) 
pdl1_jco %>% group_by(Gender) %>% summarise(N=n())
pdl1ROC_jco <- plotROC(pdl1_jco, PDL1_Score, Clinical_Benefit, Gender)  + 
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/pdl1ROC_JCO_Rizvi.pdf", plot=pdl1ROC_jco, width=4, height=3)

###### TMB predict in smoker or non-smoker
stmb_hellmann <- sampleInfo_Hellmann %>% mutate(Smoke=ifelse(Smoking_History == "never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = sTMB, target = Clinical_Benefit,group = Gender, positive = "DCB") +
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/stmbROC_Hellmann.pdf", plot=stmb_hellmann, width=4, height=3)


stmb_sci_Rizvi <- sampleInfo_Sci_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former") %>% plotROC(predict_col = sTMB, target = Clinical_Benefit,group = Gender, positive = "DCB") +
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/stmbROC_Sci_Rizvi.pdf", plot=stmb_sci_Rizvi, width=4, height=3)


stmb_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>% 
    filter(Smoke == "Current/former", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = sTMB, target = Clinical_Benefit,group = Gender, positive = "DCB") +
    scale_color_manual(values = my_palette) + add_modify
t <- sampleInfo_JCO_Rizvi %>% mutate(Smoke=ifelse(Smoking_History == "Never", "Never", "Current/former")) %>%
    filter(Smoke == "Never", Gene_Panel=="IMPACT410") %>% plotROC(predict_col = sTMB, target = Clinical_Benefit,group = Gender, positive = "DCB") +
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/stmbROC_JCO_Rizvi_panel410.pdf", plot=stmb_JCO_Rizvi, width=4, height=3)
ggsave("Figures/stmbROC_JCO_Rizvi_panel410_NeverSmoker.pdf", plot=t, width=4, height=3)


############# Cutoff analysis
library(survival)
library(survminer)
colm <- c("sTMB", "Clinical_Benefit", "Gender")
nsclc <- rbind(sampleInfo_Hellmann[, c(colm, "PFS_Months", "PFS_Event")],
               sampleInfo_JCO_Rizvi[, c(colm, "PFS_Months", "PFS_Event")], 
               sampleInfo_Sci_Rizvi[, c(colm, "PFS_Months", "PFS_Event")])

quantile(nsclc[nsclc$Gender=="Female", "sTMB"], probs = seq(0, 1, 0.1))
quantile(nsclc[nsclc$Gender=="Male", "sTMB"], probs = seq(0, 1, 0.1))
quantile(nsclc[, "sTMB"], probs = seq(0, 1, 0.1))

mela <- rbind(MELA_NEJM[, c(colm, "OS", "Event")],
              MELA_Science[, c(colm, "OS", "Event")])
quantile(mela[mela$Gender=="Female", "sTMB"] %>% unlist, probs = seq(0, 1, 0.1))
quantile(mela[mela$Gender=="Male", "sTMB"] %>% unlist, probs = seq(0, 1, 0.1))
quantile(mela[, "sTMB"] %>% unlist, probs = seq(0, 1, 0.1))


dyHR <- function(data, time="PFS_Months", event="PFS_Event"){
    data <- as.data.frame(data)
    male <- subset(data, Gender=="Male")
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


ggplot(cs_nsclc, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth() 

cs_nsclc_sciRizvi <- dyHR(sampleInfo_Sci_Rizvi)
ggplot(cs_nsclc_sciRizvi, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth()

cs_nsclc_JCORizvi <- dyHR(sampleInfo_JCO_Rizvi)
ggplot(cs_nsclc_JCORizvi, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth()

cs_nsclc_Hellmann <- dyHR(sampleInfo_Hellmann)
ggplot(cs_nsclc_JCO, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth()

save(cs_nsclc, cs_nsclc_Hellmann, cs_nsclc_JCORizvi, cs_nsclc_sciRizvi, file="Rdata/cache_HR_nsclc.RData")

cs_mela <- dyHR(mela, time = "OS", event = "Event")
ggplot(cs_mela, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth() 

cs_mela_nejm <- dyHR(MELA_NEJM, time = "OS", event = "Event")
ggplot(cs_mela_nejm, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth() 

cs_mela_science <- dyHR(MELA_Science, time = "OS", event = "Event")
ggplot(cs_mela_science, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth() 


## use cutoff to determine HR

# NSCLC
nsclc <- nsclc %>% mutate(Cutoff = factor(ifelse(sTMB > 4, "High", "Low"),
                                          levels = c("Low", "High") ))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(nsclc, Gender=="Male")))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = subset(nsclc, Gender=="Female")))

sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% mutate(Cutoff = factor(ifelse(sTMB > 4, 
                                                                "High", "Low"),
                                          levels = c("Low", "High") ))
sampleInfo_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% mutate(Cutoff = factor(ifelse(sTMB > 4, 
                                                                               "High", "Low"),
                                                                        levels = c("Low", "High") ))
sampleInfo_Hellmann <- sampleInfo_Hellmann %>% mutate(Cutoff = factor(ifelse(sTMB > 4, 
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

# MELA
mela<- mela %>% mutate(Cutoff = factor(ifelse(sTMB > 4, "High", "Low"),
                                          levels = c("Low", "High") ))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(mela, Gender=="Male")))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(mela, Gender=="Female")))

MELA_NEJM <- MELA_NEJM %>% mutate(Cutoff = factor(ifelse(sTMB > 4, "High", "Low"),
                                       levels = c("Low", "High") ))
MELA_Science <- MELA_Science %>% mutate(Cutoff = factor(ifelse(sTMB > 4, "High", "Low"),
                                                  levels = c("Low", "High") ))

summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(MELA_NEJM, Gender=="Male")))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(MELA_NEJM, Gender=="Female")))

summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(MELA_Science, Gender=="Male")))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = subset(MELA_Science, Gender=="Female")))

tb_mela_nejm <- tableData(MELA_NEJM)
tb_mela_science <- tableData(MELA_Science)

tb_mela_nejm$male %>% fisher.test()
tb_mela_nejm$female %>% fisher.test()

tb_mela_science$male %>% fisher.test()
tb_mela_science$female %>% fisher.test()


#### plot K-M curve
# fm_nsclc <- as.formula("Surv(PFS_Months, PFS_Event) ~ Cutoff") 
# fm_mela <- as.formula("Surv(OS, Event) ~ Cutoff") 

# theme_wsx <- theme_survminer(base_size = 8, font.x = 8, font.y = 8,
#                              font.main = 10, font.submain = 8,
#                              font.tickslab = 8)

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

# par(mar=rep(0.1, 4))
ggsave("Figures/M-KMplot_sci_Rizvi_Male.pdf", plot = print(surv1_1, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_sci_Rizvi_Female.pdf", plot = print(surv1_2, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_JCO_Rizvi_Male.pdf", plot = print(surv1_3, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_JCO_Rizvi_Female.pdf", plot = print(surv1_4, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_Hellmann_Male.pdf", plot = print(surv1_5, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_Hellmann_Female.pdf", plot = print(surv1_6, newpage = FALSE),
       width = 7, height = 5)


#### MELA
fit2_1 <- survfit(Surv(OS, Event) ~ Cutoff, 
                data=MELA_NEJM %>% filter(Gender=="Male"))
fit2_2 <- survfit(Surv(OS, Event) ~ Cutoff, 
                data=MELA_NEJM %>% filter(Gender=="Female"))
fit2_3 <- survfit(Surv(OS, Event) ~ Cutoff, 
                  data=MELA_Science %>% filter(Gender=="Male"))
fit2_4 <- survfit(Surv(OS, Event) ~ Cutoff, 
                  data=MELA_Science %>% filter(Gender=="Female"))


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
                    break.time.by = 12,
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

surv2_1 <- surv_plot(fit2_1)
surv2_2 <- surv_plot(fit2_2)
surv2_3 <- surv_plot(fit2_3)
surv2_4 <- surv_plot(fit2_4)

ggsave("Figures/M-KMplot_MELA_NEJM_Male.pdf", plot = print(surv2_1,  newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_MELA_NEJM_Female.pdf", plot = print(surv2_2, newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_MELA_Science_Male.pdf", plot = print(surv2_3,  newpage = FALSE),
       width = 7, height = 5)
ggsave("Figures/M-KMplot_MELA_Science_Female.pdf", plot = print(surv2_4,  newpage = FALSE),
       width = 7, height = 5)


## add ROC analysis for MELA
p3_1 <- plotROC(MELA_NEJM, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p3_2 <- plotROC(MELA_Science, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

p4_1 <- plotROC(nsclc, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
p5_1 <- plotROC(mela, sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-MELA-NEJM.pdf", plot = p3_1, width = 4, height = 3)
ggsave("Figures/CompareROC-MELA-Science.pdf", plot = p3_2, width = 4, height = 3)

ggsave("Figures/CompareROC-MELA-ALL.pdf", plot = p4_1, width = 4, height = 3)
ggsave("Figures/CompareROC-NSCLC-ALL.pdf", plot = p5_1, width = 4, height = 3)

 # plotROC(sampleInfo_Sci_Rizvi, sTMB, Clinical_Benefit, Gender)
# plotROC(sampleInfo_JCO_Rizvi, sTMB, Clinical_Benefit, Gender)
# plotROC(sampleInfo_Hellmann, sTMB, Clinical_Benefit, Gender) + 
#     scale_color_manual(values = c("black", "red", "blue")) + theme(axis.line = element_line(colour = 'black', size = 1.0),
#                                                                    axis.ticks = element_line(colour = "black", size = 1.0),
#                                                                    axis.text.y = element_blank(), 
#                                                                    axis.text.x = element_blank())


# compareMutPlot(MELA_Science, group1 = "Gender", group2 = "Clinical_Benefit",
#                value = "sTMB") + xlab("") +ylab("Tumor Mutation Burden") 
# 
# compareMutPlot(MELA_NEJM, group1 = "Gender", group2 = "Clinical_Benefit",
#                value = "sTMB") + xlab("") +ylab("Tumor Mutation Burden") 
