# tidy and re-construct analyses for immunotherapy on TMB and ROC etc..
# Shixiang Wang

library(tidyverse)

# setwd("G:/biodata/immunotherapyDatasets/")
setwd("/Volumes/data/biodata/immunotherapyDatasets/")

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

setwd("/home/wsx/Downloads/neoQ")
############### Compare sTMB (standardized Tumor Mutation Burden) ######################
# library(RColorBrewer)

################# New Figure

library(cowplot)

benefit_info = 
    dplyr::bind_rows(
        dplyr::select(dplyr::mutate(sampleInfo_Sci_Rizvi, Study="Rizvi 2015"), Gender, Clinical_Benefit, sTMB, Study) %>% dplyr::filter(Clinical_Benefit %in% c("DCB", "NDB")),
        dplyr::select(dplyr::mutate(sampleInfo_Hellmann, Study="Hellmann 2018"), Gender, Clinical_Benefit, sTMB, Study),
        dplyr::select(dplyr::mutate(sampleInfo_JCO_Rizvi, Study="Rizvi 2018"), Gender, Clinical_Benefit, sTMB, Study),
        dplyr::filter(dplyr::select(dplyr::mutate(sampleInfo_JCO_Rizvi, Study="Rizvi 2018 -341"), Gender, Clinical_Benefit, sTMB, Gene_Panel, Study), Gene_Panel=="IMPACT341") %>% dplyr::select(-Gene_Panel),
        dplyr::filter(dplyr::select(dplyr::mutate(sampleInfo_JCO_Rizvi, Study="Rizvi 2018 -410"), Gender, Clinical_Benefit, sTMB, Gene_Panel, Study), Gene_Panel=="IMPACT410") %>% dplyr::select(-Gene_Panel)
    )

benefit_info$Study = factor(benefit_info$Study, 
                            levels = c("Rizvi 2015", "Hellmann 2018", "Rizvi 2018",
                                       "Rizvi 2018 -341", "Rizvi 2018 -410"))
dodge <- position_dodge(width = 0.9)
p = ggplot(benefit_info, aes(x=Gender, y=sTMB, fill=Clinical_Benefit)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_brewer(palette="RdBu") + 
    stat_compare_means(aes_string(group="Clinical_Benefit"), label = "p.format", 
                       method = "wilcox.test") + ylab("TMB (mutations per MB)") + 
    facet_wrap(~Study, scales = "free_y") + theme(legend.position = c(0.7, 0.2), axis.title.x = element_blank()) + 
    labs(fill = "Clinical Benefit")
save_plot(
    "boxplot_vs.pdf", plot = p, base_height = 6, base_aspect_ratio = 1.6
)

#########################



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
ggplot(cs_nsclc_Hellmann, aes(x=cutoff, y=HR, color=Gender)) + geom_point() + geom_smooth()

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
nsclc$dataset = c(rep("Hellmann 2018", nrow(sampleInfo_Hellmann)),
                  rep("Rizvi 2018", nrow(sampleInfo_JCO_Rizvi)),
                  rep("Rizvi 2015", nrow(sampleInfo_Sci_Rizvi)))

## NEW figures for ROC comparison ################################
library(pROC)

roc_df = dplyr::select(nsclc, c(Gender, Smoking_History, Clinical_Benefit, sTMB, dataset))
roc_df$response = ifelse(roc_df$Clinical_Benefit=="DCB", 1, 0)

do_roc = function(df, study) {
    library(pROC)
    df_f = subset(df, Gender == "Female")
    df_m = subset(df, Gender == "Male")
    roc_f = roc(df_f$response, df_f$sTMB, direction = "<")
    roc_m = roc(df_m$response, df_m$sTMB, direction = "<")
    roc_a = roc(df$response, df$sTMB, direction = "<")
    roc_l = list(Female=roc_f, Male=roc_m, ALL=roc_a)
    
    message("AUC for female: ", auc(roc_f))
    message("AUC for male: ", auc(roc_m))
    message("AUC for all: ", auc(roc_a))
    
    ggroc(roc_l, legacy.axes = TRUE) -> p1
    
    return(data.frame(specificity = p1$data$specificity,
                      sensitivity = p1$data$sensitivity,
                      group = p1$data$name,
                      study = study))
}

roc_list = list()
roc_list[["Rizvi_2015"]] = do_roc(subset(roc_df, dataset=="Rizvi 2015"), study = "Rizvi 2015")
roc_list[["Hellmann_2018"]] = do_roc(subset(roc_df, dataset=="Hellmann 2018"), study = "Hellmann 2018")
roc_list[["Rizvi_2018"]] = do_roc(subset(roc_df, dataset=="Rizvi 2018"), study = "Rizvi 2018")
roc_list[["NSCLC_combined"]] = do_roc(roc_df, study = "NSCLC combined")
# This data need to be corrected,
# ALL  0.629
# Male 0.563
# Female 0.686

sampleInfo_JCO_Rizvi$response = ifelse(sampleInfo_JCO_Rizvi$Clinical_Benefit=="DCB", 1, 0)
roc_list[["Rizvi_2018_341"]] = do_roc(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT341"), study = "Rizvi 2018 -341")
roc_list[["Rizvi_2018_410"]] = do_roc(subset(sampleInfo_JCO_Rizvi, Gene_Panel=="IMPACT410"), study = "Rizvi 2018 -410")


# ROC ALL
roc_all = dplyr::rbind_list(roc_list)
roc_all$study = factor(roc_all$study, 
                       levels = c("Rizvi 2015", "Hellmann 2018", "Rizvi 2018",
                                "Rizvi 2018 -341", "Rizvi 2018 -410", "NSCLC combined"))
roc_all$group = factor(roc_all$group, 
                       levels = c("Male", "ALL", "Female"))
ggplot(roc_all, aes(x=1-specificity, y=sensitivity, linetype=group, color=group)) +
    geom_line() + scale_color_manual(values = c("blue", "black", "red")) +  
    facet_wrap(~study) + labs(linetype="", color="") -> p2

save_plot(
    "ROC_ALL.pdf", plot = p2, base_height = 6, base_aspect_ratio = 1.6
)

# ROC Smoking
# sampleInfo_Sci_Rizvi$response = ifelse(sampleInfo_Sci_Rizvi$Clinical_Benefit=="DCB", 1, 0)
# sampleInfo_Hellmann$response = ifelse(sampleInfo_Hellmann$Clinical_Benefit=="DCB", 1, 0)

roc_list2 = list()
roc_list2[["Rizvi_2015"]] = do_roc(subset(roc_df, 
                                          dataset=="Rizvi 2015" & 
                                              Smoking_History=="Current/former"),
                                   study = "Rizvi 2015 smoker")
roc_list2[["Hellmann_2018"]] = do_roc(subset(roc_df, 
                                          dataset=="Hellmann 2018" & 
                                              Smoking_History=="Current/former"),
                                   study = "Hellmann 2018 smoker")
roc_list2[["Rizvi_2018"]] = do_roc(sampleInfo_JCO_Rizvi %>% 
                                       dplyr::mutate(Smoke=ifelse(Smoking_History == "Never", 
                                                           "Never", "Current/former")) %>% 
                                       dplyr::filter(Smoke == "Current/former", 
                                                     Gene_Panel=="IMPACT410"),
                                   study = "Rizvi 2018 smoker")
roc_list2[["Rizvi_2018_nonsmoker"]] = do_roc(sampleInfo_JCO_Rizvi %>% 
                                       dplyr::mutate(Smoke=ifelse(Smoking_History == "Never", 
                                                                  "Never", "Current/former")) %>% 
                                       dplyr::filter(Smoke == "Never", 
                                                     Gene_Panel=="IMPACT410"),
                                   study = "Rizvi 2018 non-smoker")


roc_smoke = dplyr::rbind_list(roc_list2)
roc_smoke$study = factor(roc_smoke$study, 
                       levels = c("Rizvi 2015 smoker", "Hellmann 2018 smoker", 
                                  "Rizvi 2018 smoker", "Rizvi 2018 non-smoker"))
roc_smoke$group = factor(roc_smoke$group, 
                       levels = c("Male", "ALL", "Female"))
ggplot(roc_smoke %>% dplyr::filter(study!="Rizvi 2018 non-smoker"),
       aes(x=1-specificity, y=sensitivity, linetype=group, color=group)) +
    geom_line() + scale_color_manual(values = c("blue", "black", "red")) +  
    facet_wrap(~study, ncol = 3) + labs(linetype="", color="") +
    theme(legend.position = "none") -> p3

save_plot(
    "ROC_smoker.pdf", plot = p3, base_height = 3, base_aspect_ratio = 2.5
)

ggplot(roc_smoke %>% dplyr::filter(study=="Rizvi 2018 non-smoker"),
       aes(x=1-specificity, y=sensitivity, linetype=group, color=group)) +
    geom_line() + scale_color_manual(values = c("blue", "black", "red")) +  
    facet_wrap(~study) + labs(linetype="", color="") +
    theme(legend.position = "none") -> p4

save_plot(
    "ROC_nonsmoker.pdf", plot = p4, base_height = 3, base_aspect_ratio = 1.1
)

###############################################

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

tb_sci_Rizvi <-  tableData(sampleInfo_Sci_Rizvi) # 3 with NR removed
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

ggsave("Figures/CompareROC-MELA-ALL.pdf", plot = p5_1, width = 4, height = 3)
ggsave("Figures/CompareROC-NSCLC-ALL.pdf", plot = p4_1, width = 4, height = 3)

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


### Oncoplot for Maf data
impGenes =  c("TP53", "KRAS", "STK11", "EGFR", "ALK", "BRAF", "RET", "ROS1", "PETEN", "POLE")

oncoplot(Science_Rizvi_maf, fontSize = 12)
Science_Rizvi_maf_F = subsetMaf(Science_Rizvi_maf, tsb=subset(Science_Rizvi_maf@clinical.data, 
                                                              Gender=="F")$Tumor_Sample_Barcode,
                                mafObj = TRUE)
oncoplot(Science_Rizvi_maf_F, fontSize = 12)
Science_Rizvi_maf_M = subsetMaf(Science_Rizvi_maf, tsb=subset(Science_Rizvi_maf@clinical.data, 
                                                              Gender=="M")$Tumor_Sample_Barcode,
                                mafObj = TRUE)
oncoplot(Science_Rizvi_maf_M, fontSize = 12, top = 50)

coOncoplot(m1 = Science_Rizvi_maf_M, m2 = Science_Rizvi_maf_F, m1Name = "Rizvi et al (2015) - Male",
           m2Name = "Rizvi et al (2015) - Female", genes = impGenes)

JCO_Rizvi_maf_F = subsetMaf(JCO_Rizvi_maf, tsb=subset(JCO_Rizvi_maf@clinical.data, 
                                                              Gender=="Female")$Tumor_Sample_Barcode,
                                mafObj = TRUE)
JCO_Rizvi_maf_M = subsetMaf(JCO_Rizvi_maf, tsb=subset(JCO_Rizvi_maf@clinical.data, 
                                                      Gender=="Male")$Tumor_Sample_Barcode,
                            mafObj = TRUE)
coOncoplot(m1 = JCO_Rizvi_maf_M, m2 = JCO_Rizvi_maf_F, m1Name = "Rizvi et al (2018) - Male",
           m2Name = "Rizvi et al (2018) - Female", genes = impGenes)

Hellmann_maf_F = subsetMaf(Hellmann_maf, tsb=subset(Hellmann_maf@clinical.data, 
                                                     Gender=="female")$Tumor_Sample_Barcode,
                           mafObj = TRUE)
Hellmann_maf_M = subsetMaf(Hellmann_maf, tsb=subset(Hellmann_maf@clinical.data, 
                                                    Gender=="male")$Tumor_Sample_Barcode,
                           mafObj = TRUE)
coOncoplot(m1 = Hellmann_maf_M, m2 = Hellmann_maf_F, m1Name = "Hellmann et al (2018) - Male",
           m2Name = "Hellmann et al (2018) - Female", genes = impGenes)


### Compare ROC considering Gene Mutation
my_palette <- c("black", "red", "blue")
add_modify <- theme(axis.line = element_line(colour = 'black', size = 1.0),
                    axis.ticks = element_line(colour = "black", size = 1.0))
#axis.text.y = element_blank(), 
#axis.text.x = element_blank()) 

sample_jcoRizvi = as.character(JCO_Rizvi_maf@clinical.data$Tumor_Sample_Barcode)

####### Jco Rizvi 
## TP53
sample_jcoRizvi_TP53_plus = JCO_Rizvi_maf@data %>% filter(Hugo_Symbol == "TP53") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_jcoRizvi_TP53_minus = setdiff(sample_jcoRizvi, sample_jcoRizvi_TP53_plus)

m_1 = plotROC(sampleInfo_JCO_Rizvi %>% filter(Tumor_Sample_Barcode %in% sample_jcoRizvi_TP53_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

m_2 = plotROC(sampleInfo_JCO_Rizvi %>% filter(Tumor_Sample_Barcode %in% sample_jcoRizvi_TP53_plus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-Mutation-JCO-Rizvi-TP53minus.pdf", plot = m_1, width = 4, height = 3)
ggsave("Figures/CompareROC-Mutation-JCO-Rizvi-TP53plus.pdf", plot = m_2, width = 4, height = 3)
## KRAS
sample_jcoRizvi_KRAS_plus = JCO_Rizvi_maf@data %>% filter(Hugo_Symbol == "KRAS") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_jcoRizvi_KRAS_minus = setdiff(sample_jcoRizvi, sample_jcoRizvi_KRAS_plus)

m_3 = plotROC(sampleInfo_JCO_Rizvi %>% filter(Tumor_Sample_Barcode %in% sample_jcoRizvi_KRAS_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

m_4 = plotROC(sampleInfo_JCO_Rizvi %>% filter(Tumor_Sample_Barcode %in% sample_jcoRizvi_KRAS_plus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-Mutation-JCO-Rizvi-KRAS_minus.pdf", plot = m_3, width = 4, height = 3)
ggsave("Figures/CompareROC-Mutation-JCO-Rizvi-KRASplus.pdf", plot = m_4, width = 4, height = 3)
## EGFR/ALK
sample_jcoRizvi_EGFR_plus = JCO_Rizvi_maf@data %>% filter(Hugo_Symbol == "EGFR") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_jcoRizvi_ALK_plus = JCO_Rizvi_maf@data %>% filter(Hugo_Symbol == "ALK") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()

sample_jcoRizvi_EGFR_ALK_minus = setdiff(sample_jcoRizvi, union(sample_jcoRizvi_EGFR_plus, sample_jcoRizvi_ALK_plus))

m_5 = plotROC(sampleInfo_JCO_Rizvi %>% filter(Tumor_Sample_Barcode %in% sample_jcoRizvi_EGFR_ALK_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/CompareROC-Mutation-JCO-Rizvi-EGFR_ALK_minus.pdf", plot = m_5, width = 4, height = 3)


sample_Hellmann = as.character(Hellmann_maf@clinical.data$Tumor_Sample_Barcode)
#######  Hellmann
## TP53
sample_Hellmann_TP53_plus = Hellmann_maf@data %>% filter(Hugo_Symbol == "TP53") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_Hellmann_TP53_minus = setdiff(sample_Hellmann, sample_Hellmann_TP53_plus)
# sample_hellmann_TP53_minus = sample_hellmann[! sample_hellmann %in% sample_hellmann_TP53_plus]


m_6 = plotROC(sampleInfo_Hellmann %>% filter(Tumor_Sample_Barcode %in% sample_Hellmann_TP53_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

m_7 = plotROC(sampleInfo_Hellmann %>% filter(Tumor_Sample_Barcode %in% sample_Hellmann_TP53_plus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-Mutation-Hellmann-TP53minus.pdf", plot = m_6, width = 4, height = 3)
ggsave("Figures/CompareROC-Mutation-Hellmann-TP53plus.pdf", plot = m_7, width = 4, height = 3)
## KRAS
sample_Hellmann_KRAS_plus = Hellmann_maf@data %>% filter(Hugo_Symbol == "KRAS") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_Hellmann_KRAS_minus = setdiff(sample_Hellmann, sample_Hellmann_KRAS_plus)

m_8 = plotROC(sampleInfo_Hellmann %>% filter(Tumor_Sample_Barcode %in% sample_Hellmann_KRAS_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

m_9 = plotROC(sampleInfo_Hellmann %>% filter(Tumor_Sample_Barcode %in% sample_Hellmann_KRAS_plus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify

ggsave("Figures/CompareROC-Mutation-Hellmann-KRAS_minus.pdf", plot = m_8, width = 4, height = 3)
ggsave("Figures/CompareROC-Mutation-Hellmann-KRASplus.pdf", plot = m_9, width = 4, height = 3)
## EGFR/ALK
sample_Hellmann_EGFR_plus = Hellmann_maf@data %>% filter(Hugo_Symbol == "EGFR") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()
sample_Hellmann_ALK_plus = Hellmann_maf@data %>% filter(Hugo_Symbol == "ALK") %>% select(Tumor_Sample_Barcode) %>% unlist %>% as.character() %>% unique()

sample_Hellmann_EGFR_ALK_minus = setdiff(sample_Hellmann, union(sample_Hellmann_EGFR_plus, sample_Hellmann_ALK_plus))

m_10 = plotROC(sampleInfo_Hellmann %>% filter(Tumor_Sample_Barcode %in% sample_Hellmann_EGFR_ALK_minus),
        sTMB, Clinical_Benefit, Gender) + 
    scale_color_manual(values = my_palette) + add_modify
ggsave("Figures/CompareROC-Mutation-Hellmann-EGFR_ALK_minus.pdf", plot = m_10, width = 4, height = 3)
