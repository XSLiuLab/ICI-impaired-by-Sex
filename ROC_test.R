# test ROC and calculate I square

library(pROC)

setwd("G:/biodata/immunotherapyDatasets/")
load("Rdata/sampleInfo_cache.RData")

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

# add numerical response column
#MELA_Cell <- MELA_Cell %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))
MELA_Science <- MELA_Science %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))
MELA_NEJM <- MELA_NEJM %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))

sampleInfo_Forde <- sampleInfo_Forde %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))
sampleInfo_Hellmann <- sampleInfo_Hellmann %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))
sampleInfo_Sci_Rizvi <- sampleInfo_Sci_Rizvi %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))
sampleInfo_JCO_Rizvi <- sampleInfo_JCO_Rizvi %>% mutate(response=ifelse(Clinical_Benefit=="DCB", 1, 0))

sampleInfo_JCO_Rizvi341 <- sampleInfo_JCO_Rizvi %>% filter(Gene_Panel=="IMPACT341")
sampleInfo_JCO_Rizvi410 <- sampleInfo_JCO_Rizvi %>% filter(Gene_Panel=="IMPACT410")


test_roc <- function(response, predictor, data){
    
    form <- formula(paste0(response, ' ~ ', predictor))
    roc1 <- roc(form, data=data, subset=(Gender == "Female"), plot=FALSE, ci=TRUE)
    roc2 <- roc(form, data=data, subset=(Gender == "Male"), plot=FALSE, ci=TRUE)
    
    CI1 <- roc1$ci
    CI2 <- roc2$ci
    var1 <- var(roc1)
    var2 <- var(roc2)
    variance <- var1 + var2
    
    tr <- roc.test(roc1, roc2)
    statistic <- tr$statistic
    pvalue <- tr$p.value
    estimate <- tr$estimate
    
    return(list(statistic=statistic, variance=variance, p.value=pvalue,
                AUC=estimate, CI1=CI1, CI2=CI2, Var1=var1, Var2=var2))
}

r_Sci_Rizvi <- test_roc("response", "TMB_NonsynSNP", sampleInfo_Sci_Rizvi)
r_JCO_Rizvi <- test_roc("response", "TMB_NonsynSNP", sampleInfo_JCO_Rizvi)
r_Hellmann  <- test_roc("response", "TMB_NonsynSNP", sampleInfo_Hellmann)
# r_Forde     <- test_roc("response", "TMB_NonsynSNP", sampleInfo_Forde)

r_JCO_Rizvi341 <- test_roc("response", "TMB_NonsynSNP", sampleInfo_JCO_Rizvi341)
r_JCO_Rizvi410 <- test_roc("response", "TMB_NonsynSNP", sampleInfo_JCO_Rizvi410)

r_MELA_Science <- test_roc("response", "mutation_load", MELA_Science)
r_MELA_NEJM    <- test_roc("response", "mutation_load", MELA_NEJM)

st <- c(r_Sci_Rizvi$statistic, r_JCO_Rizvi341$statistic, 
        r_JCO_Rizvi410$statistic, r_Hellmann$statistic, 
        r_MELA_NEJM$statistic, r_MELA_Science$statistic)
va <- c(r_Sci_Rizvi$variance, r_JCO_Rizvi341$variance, 
        r_JCO_Rizvi410$variance, r_Hellmann$variance, 
        r_MELA_NEJM$variance, r_MELA_Science$variance)

roc_df <- data.frame(statistic=st, variance=va)

library(altmeta)
roc_meta <- metahet(roc_df$statistic, s2=roc_df$variance)
roc_meta2 <- metahet(roc_df[1:4,]$statistic, s2=roc_df[1:4, ]$variance)

# data(aSAH)
# t1 = roc(aSAH$outcome, aSAH$s100b,
#     levels=c("Good", "Poor"))
# t = roc(aSAH$outcome, aSAH$s100b,
#     percent=TRUE, plot=TRUE, ci=TRUE)
# # Using subset (only with formula)
# roc(outcome ~ s100b, data=aSAH, subset=(gender == "Male"))
# roc(outcome ~ s100b, data=aSAH, subset=(gender == "Female"))
