# Survival and HR analysis for immunotherapy datasets

# create function used for visualize Hazard Ratio (forest plot)
plotHR <- function(df, time="PFS_Months", status="PFS_Event", expr,
                   main = "Hazard ratio", fontsize = 0.7,
                   refLabel = "reference", noDigits = 3){
    df <- as.data.frame(df)
    #expr <- substitute(expr)
    #print(expr)
    expr <- paste0("Surv(", time, ",", status,")", " ~ ", expr)
    expr <- formula(expr)
    if(is.null(expr)){
        stop("Please specify an expression used in coxph() function.")
    }
    
    res <- list()
    res$model <- coxph(expr, data=df)
    res$ggforest <- survminer::ggforest(model = res$model, data = df, main = main,
                                        fontsize = fontsize, refLabel = refLabel,
                                        noDigits = noDigits)
    res$diagno <- survminer::ggcoxdiagnostics(res$model, type="schoenfeld", os.scale="time")
    print(res$model)
    res$ggforest
    
    return(res)
}


setwd("G:/biodata/immunotherapyDatasets/")
load("Rdata/sampleInfo_cache.RData")

require(survival)
require(survminer)

res.cut <- surv_cutpoint(sampleInfo_Sci_Rizvi,
                         time="PFS_Months", 
                         event = "PFS_Event",
                         variables = "TMB_NonsynSNP")
summary(res.cut)
plot(res.cut, "TMB_NonsynSNP", palette="jco")

res.cat <- surv_categorize(res.cut)
head(res.cat)
fit <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_NonsynSNP,
               data = res.cat)
ggsurvplot(fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE)


buildCModel <- function(.data, time, status, gender, measure){
    require(tidyverse)
    time      <- enquo(time)
    status    <- enquo(status)
    gender    <- enquo(gender)
    measure   <- enquo(measure)
    
    # transform quos to strings
    timeN    <- quo_name(time)
    statusN  <- quo_name(status)
    genderN  <- quo_name(gender)
    measureN <- quo_name(measure)
    
    df_male  <- .data %>% filter(!!gender == "Male")
    df_male  <- df_male %>% 
        mutate(TMB_Status = ifelse(!! measure >= median(!! measure), "High", "Low"),
               Cutoff = ifelse(!! measure >= as.numeric(surv_cutpoint(df_male,
                                                                      time=timeN, 
                                                                      event = statusN,
                               variables = measureN)$cutpoint[1]),
                               "High", "Low")) %>% 
        mutate(Cutoff = factor(Cutoff, levels=c("Low", "High")), 
               TMB_Status = factor(TMB_Status, levels=c("Low", "High")))
    
    df_female  <- .data %>% filter(!!gender == "Female") 
    df_female  <- df_female %>% 
        mutate(TMB_Status = ifelse(!! measure >= median(!! measure), "High", "Low"),
               Cutoff = ifelse(!! measure >= as.numeric(surv_cutpoint(df_female,
                                                                      time=timeN, 
                                                                      event = statusN,
                                variables = measureN)$cutpoint[1]),
                               "High", "Low")) %>% 
        mutate(Cutoff = factor(Cutoff, levels=c("Low", "High")), 
               TMB_Status = factor(TMB_Status, levels=c("Low", "High")))
    
    # male <- list()
    # fit_median <- formula(paste0("Surv(", timeN, ",", statusN,")", " ~ ", "TMB_Status"))
    # fit_cutoff <- formula(paste0("Surv(", timeN, ",", statusN,")", " ~ ", "Cutoff"))
    # male$surv_model1 <- survfit(fit_median, data=df_male)
    # male$surv_model2 <- survfit(fit_cutoff, data=df_male)
    # male$coxph_model1 <- coxph(fit_median, data=df_male)
    # male$coxph_model2 <- coxph(fit_cutoff, data=df_male)
    # male$data <- df_male
    # 
    # female <- list()
    # female$surv_model1 <- survfit(fit_median, data=df_female)
    # female$surv_model2 <- survfit(fit_cutoff, data=df_female)
    # female$coxph_model1 <- coxph(fit_median, data=df_female)
    # female$coxph_model2 <- coxph(fit_cutoff, data=df_female)
    # female$data <- df_female
    # 
    # res <- list(male, female)
    
    res <- list()
    res$male <- df_male
    res$female <- df_female
    return(res)
}

df_sci_Rizvi <- buildCModel(sampleInfo_Sci_Rizvi, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
View(df_sci_Rizvi$male)

fit1 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff,
               data = df_sci_Rizvi$male)
fit2 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff,
                data = df_sci_Rizvi$female)
ggsurvplot(fit1, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit2, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
        data = df_sci_Rizvi$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
      data = df_sci_Rizvi$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_sci_Rizvi$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_sci_Rizvi$female))


df_JCO_Rizvi <- buildCModel(sampleInfo_JCO_Rizvi, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit3 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff,
                data = df_JCO_Rizvi$male)
fit4 <- survfit(Surv(PFS_Months, PFS_Event) ~ Cutoff,
                data = df_JCO_Rizvi$female)
ggsurvplot(fit3, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit4, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
      data = df_JCO_Rizvi$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
      data = df_JCO_Rizvi$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi$female))


df_Hellmann <- buildCModel(sampleInfo_Hellmann, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)

fit5 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = df_Hellmann$male)
fit6 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = df_Hellmann$female)
ggsurvplot(fit5, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit6, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
      data = df_Hellmann$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
      data = df_Hellmann$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_Hellmann$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_Hellmann$female))

# JCO datasets have different gene panels
df_JCO_Rizvi_341 <- buildCModel(sampleInfo_JCO_Rizvi %>% filter(Gene_Panel=="IMPACT341"),
                            PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit7_1 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = df_JCO_Rizvi_341$male)
fit7_2 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = df_JCO_Rizvi_341$female)
ggsurvplot(fit7_1, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit7_2, risk.table = TRUE, conf.int = TRUE, pval = TRUE)



summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_341$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_341$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_341$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_341$female))

df_JCO_Rizvi_410 <- buildCModel(sampleInfo_JCO_Rizvi %>% filter(Gene_Panel=="IMPACT410"),
                                PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit7_3 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                  data = df_JCO_Rizvi_410$male)
fit7_4 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                  data = df_JCO_Rizvi_410$female)
ggsurvplot(fit7_3, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit7_4, risk.table = TRUE, conf.int = TRUE, pval = TRUE)



summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_410$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_410$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_410$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_410$female))

df_JCO_Rizvi_468 <- buildCModel(sampleInfo_JCO_Rizvi %>% filter(Gene_Panel=="IMPACT468"),
                                PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit7_5 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                  data = df_JCO_Rizvi_468$male)
fit7_6 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                  data = df_JCO_Rizvi_468$female)
ggsurvplot(fit7_5, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit7_6, risk.table = TRUE, conf.int = TRUE, pval = TRUE)



summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_410$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
              data = df_JCO_Rizvi_410$female))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_410$male))
summary(coxph(Surv(PFS_Months, PFS_Event) ~ Cutoff,
              data = df_JCO_Rizvi_410$female))



################## For MELA dataset
MELA_NEJM <- read_csv("MELA/MELA_NEJM.csv")
MELA_Cell <- read_csv("MELA/MELA_Cell.csv")
MELA_Science <- read_csv("MELA/MELA_Science.csv")

MELA_NEJM <- MELA_NEJM %>% 
    mutate(Gender = ifelse(Gender == "F", "Female", "Male"), 
           OS = OS_Year * 12, 
           Event = OS_Status)

MELA_Cell <- MELA_Cell %>% 
    mutate(Gender = ifelse(Gender == "F", "Female", "Male"),
           OS = OS_Day / 30, 
           Event = ifelse(OS_Status=="Dead", 1, 0))

MELA_Science <- MELA_Science %>% 
    mutate(Gender = ifelse(Gender == "female", "Female", "Male"),
           OS = OS_Day / 30, 
           Event = OS_Status)

MELA_Science <- filter(MELA_Science, group != "long-survival")

df_mela_NEJM <- buildCModel(MELA_NEJM, OS, Event, Gender, TMB)

fit8_1 <- survfit(Surv(OS, Event) ~ Cutoff,
                  data = df_mela_NEJM$male)
fit8_2 <- survfit(Surv(OS, Event) ~ Cutoff,
                  data = df_mela_NEJM$female)
ggsurvplot(fit8_1, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit8_2, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

summary(coxph(Surv(OS, Event) ~ TMB_Status,
              data = df_mela_NEJM$male))
summary(coxph(Surv(OS, Event) ~ TMB_Status,
              data = df_mela_NEJM$female))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = df_mela_NEJM$male))
summary(coxph(Surv(OS, Event) ~ Cutoff,
              data = df_mela_NEJM$female))



df_mela_Cell <- buildCModel(MELA_Cell, OS, Event, Gender, TMB)

fit8_3 <- survfit(Surv(OS, Event) ~ TMB_Status,
                  data = df_mela_Cell$male)
fit8_4 <- survfit(Surv(OS, Event) ~ TMB_Status,
                  data = df_mela_Cell$female)
ggsurvplot(fit8_3, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit8_4, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

df_mela_Science <- buildCModel(MELA_Science, OS, Event, Gender, TMB)

fit8_5 <- survfit(Surv(OS, Event) ~ Cutoff,
                  data = df_mela_Science$male)
fit8_6 <- survfit(Surv(OS, Event) ~ Cutoff,
                  data = df_mela_Science$female)
ggsurvplot(fit8_5, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
ggsurvplot(fit8_6, risk.table = TRUE, conf.int = TRUE, pval = TRUE)

summary(coxph(Surv(OS, Event) ~ TMB_Status,
              data = df_mela_Science$male))
summary(coxph(Surv(OS, Event) ~ TMB_Status,
              data = df_mela_Science$female))


############# build another model - divide samples into 4 parts
buildCModel2 <- function(.data, time, status, gender, measure){
    require(tidyverse)
    time      <- enquo(time)
    status    <- enquo(status)
    gender    <- enquo(gender)
    measure   <- enquo(measure)
    
    # transform quos to strings
    timeN    <- quo_name(time)
    statusN  <- quo_name(status)
    genderN  <- quo_name(gender)
    measureN <- quo_name(measure)
    
    qGroup <- function(x){
        qn <- quantile(x)
        ifelse(x < qn[2], "[0, 25%)",
               ifelse(x < qn[3], "[25%, 50%)", 
                      ifelse(x < qn[4], "[50%, 75%)", "[75%, 1]" 
                             )))
    }
    
    df_male  <- .data %>% filter(!!gender == "Male") %>% as.data.frame()
    df_male$TMB_Status <- qGroup(df_male[, measureN])
    df_male$TMB_Status <- factor(df_male$TMB_Status, levels=c("[0, 25%)","[25%, 50%)",
                                                              "[50%, 75%)", "[75%, 1]" ))

    df_female  <- .data %>% filter(!!gender == "Female") %>% as.data.frame()
    df_female$TMB_Status <- qGroup(df_female[, measureN])
    df_female$TMB_Status <- factor(df_female$TMB_Status, levels=c("[0, 25%)","[25%, 50%)",
                                                              "[50%, 75%)", "[75%, 1]" ))    

    res <- list()
    res$male <- df_male
    res$female <- df_female
    return(res)
}


buildCModel3 <- function(.data, time, status, gender, measure){
    require(tidyverse)
    time      <- enquo(time)
    status    <- enquo(status)
    gender    <- enquo(gender)
    measure   <- enquo(measure)
    
    # transform quos to strings
    timeN    <- quo_name(time)
    statusN  <- quo_name(status)
    genderN  <- quo_name(gender)
    measureN <- quo_name(measure)
    
    qGroup <- function(x){
        qn <- quantile(x, c(0.2,0.4,0.6,0.8,1.0))
        ifelse(x < qn[1], "[0, 20%)",
               ifelse(x < qn[2], "[20%, 40%)", 
                      ifelse(x < qn[3], "[40%, 60%)", 
                             ifelse(x < qn[4], "[60%, 80%)", "[80%, 1]" )
                      )))
    }
    
    df_male  <- .data %>% filter(!!gender == "Male") %>% as.data.frame()
    df_male$TMB_Status <- qGroup(df_male[, measureN])
    df_male$TMB_Status <- factor(df_male$TMB_Status, levels=c("[0, 20%)","[20%, 40%)",
                                                              "[40%, 60%)", "[60%, 80%)", "[80%, 1]" ))
    
    df_female  <- .data %>% filter(!!gender == "Female") %>% as.data.frame()
    df_female$TMB_Status <- qGroup(df_female[, measureN])
    df_female$TMB_Status <- factor(df_female$TMB_Status, levels=c("[0, 20%)","[20%, 40%)",
                                                                  "[40%, 60%)", "[60%, 80%)", "[80%, 1]" ))    
    
    res <- list()
    res$male <- df_male
    res$female <- df_female
    return(res)
}

mela_nejm <-  buildCModel3(MELA_NEJM, OS, Event, Gender, TMB)
fit10_1 <- survfit(Surv(OS, Event) ~ TMB_Status,
                  data = mela_nejm$male)
fit10_2 <- survfit(Surv(OS, Event) ~ TMB_Status,
                  data = mela_nejm$female)
ggsurvplot(fit10_1, risk.table = TRUE, pval = TRUE)
ggsurvplot(fit10_2, risk.table = TRUE, pval = TRUE)

mela_science <- buildCModel2(MELA_Science, OS, Event, Gender, TMB)
fit10_3 <- survfit(Surv(OS, Event) ~ TMB_Status,
                   data = mela_science$male)
fit10_4 <- survfit(Surv(OS, Event) ~ TMB_Status,
                   data = mela_science$female)
ggsurvplot(fit10_3, risk.table = TRUE, pval = TRUE)
ggsurvplot(fit10_4, risk.table = TRUE, pval = TRUE)

mela_cell <- buildCModel2(MELA_Cell, OS, Event, Gender, TMB)
fit10_11 <- survfit(Surv(OS, Event) ~ TMB_Status,
                   data = mela_cell$male)
fit10_12 <- survfit(Surv(OS, Event) ~ TMB_Status,
                   data = mela_cell$female)
ggsurvplot(fit10_11, risk.table = TRUE, pval = TRUE)
ggsurvplot(fit10_12, risk.table = TRUE, pval = TRUE)

nsclc_hellmann <-  buildCModel2(sampleInfo_Hellmann, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit10_5 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = nsclc_hellmann$male)
fit10_6 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                data = nsclc_hellmann$female)
ggsurvplot(fit10_5, risk.table = TRUE,  pval = TRUE)
ggsurvplot(fit10_6, risk.table = TRUE, pval = TRUE)

nsclc_sci_Rizvi <-  buildCModel2(sampleInfo_Sci_Rizvi, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit10_7 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                   data = nsclc_sci_Rizvi$male)
fit10_8 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                   data = nsclc_sci_Rizvi$female)
ggsurvplot(fit10_7, risk.table = TRUE,  pval = TRUE)
ggsurvplot(fit10_8, risk.table = TRUE, pval = TRUE)

nsclc_JCO_Rizvi <-  buildCModel2(sampleInfo_JCO_Rizvi, PFS_Months, PFS_Event, Gender, TMB_NonsynSNP)
fit10_9 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                   data = nsclc_JCO_Rizvi$male)
fit10_10 <- survfit(Surv(PFS_Months, PFS_Event) ~ TMB_Status,
                   data = nsclc_JCO_Rizvi$female)
ggsurvplot(fit10_9, risk.table = TRUE,  pval = TRUE)
ggsurvplot(fit10_10, risk.table = TRUE, pval = TRUE)
