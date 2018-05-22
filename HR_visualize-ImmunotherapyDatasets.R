# Analyze clinical data and visualize by HR forest plots
library(survival)
library(survminer)

# lung$age <- ifelse(lung$age > 70, ">70", "<=70")
# fit <- coxph(Surv(time, status) ~ sex+ph.ecog+age, data=lung)
# ggforest(fit) 


setwd("G:/biodata/immunotherapyDatasets/")
load(file ="Rdata/cache_data_code.RData")

summaryNeos <- function(neo_path, nq_path){
    require(tidyverse)
    # read neoantigen result
    neo <- read_tsv(neo_path)
    # read neoantigen quality result
    nq  <- read_tsv(nq_path)
    colnames(neo) <- make.names(colnames(neo))
    colnames(nq)  <- make.names(colnames(nq))
    
    # exclude Excluded == 1
    nq_filter <- filter(nq, Excluded != 1)
    nq_sm <- nq_filter %>% 
        group_by(Sample) %>% 
        summarise(SingleNCounts = n(), # SingleNCounts means SNV neoantigen
                  NQuality = max(NeoantigenRecognitionPotential)) %>%
        select(Sample, SingleNCounts, NQuality)
    
    print("Summary Neoantigen Quality results:\n")
    print(nq_sm)
    
    neo_sm <- neo %>% 
        group_by(Sample.Name) %>%
        summarize(NeoCounts=n(),
                  NeoC_Strong=length(which(Best.MT.Score < 50))) %>% 
        select(Sample.Name, NeoCounts, NeoC_Strong)
    print("Summary Neoantigen (pvacseq) results:\n")
    print(neo_sm)
    
    dplyr::full_join(nq_sm, neo_sm, by=c("Sample"="Sample.Name"))
        
}

# summary neoantigen results
Neos_Forde    <- summaryNeos(neo_path = "Neo_result/merged_neoantigens_Forde.tsv", 
                    nq_path = "Neo_result/neoantigen_Quality_Forde.txt")
Neos_Hellmann <- summaryNeos(neo_path = "Neo_result/merged_neoantigens_Hellmann.tsv", 
                    nq_path = "Neo_result/neoantigen_Quality_Hellmann.txt")
Neos_Sci_Rizvi<- summaryNeos(neo_path = "Neo_result/merged_neoantigens_Sci_Rizvi.tsv", 
                    nq_path = "Neo_result/neoantigen_Quality_Sci_Rizvi.txt")


# add neoantigen results to sampleInfo
sampleInfo_Forde <- dplyr::full_join(sampleInfo_Forde, Neos_Forde, by=c("Tumor_Sample_Barcode"="Sample"))
sampleInfo_Hellmann <- dplyr::full_join(sampleInfo_Hellmann, Neos_Hellmann, 
                                        by=c("Tumor_Sample_Barcode"="Sample"))
sampleInfo_Sci_Rizvi <- dplyr::full_join(sampleInfo_Sci_Rizvi, Neos_Sci_Rizvi,
                                         by=c("Tumor_Sample_Barcode"="Sample"))

#save(sampleInfo_Forde, sampleInfo_Hellmann, sampleInfo_JCO_Rizvi, sampleInfo_Sci_Rizvi, file = "Rdata/sampleInfo_cache.RData")
# remove Neoantigen load data column from original paper
# we use Neoantigen load computed by consistent method (NetMHC 4.0) here
sampleInfo_Hellmann  <- select(sampleInfo_Hellmann, -Nonsyn, -Neoantigen)
sampleInfo_Sci_Rizvi <- select(sampleInfo_Sci_Rizvi, -Neoantigen, -Nonsyn)


fit <- coxph(Surv(time, status) ~ sex+ph.ecog+age, data=lung)
# Now, create function used for visualize Hazard Ratio (forest plot)
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

TidyInfo_Hellmann <- sampleInfo_Hellmann %>% mutate(Performance_status = as.factor(ifelse(ECOG_performance_status==0, "ECOG 0", "ECOG 1")),
                                Gender = as.factor(Gender), 
                                Age = Age, 
                               Smoking_Status = as.factor(ifelse(Smoking_History=="never", "Never", "Current/former")),
                               Histology = as.factor(Histology), 
                               Best_overall_response = as.factor(ifelse(Treatment_Best_Response%in%c("CR","PR"), "CR/PR",
                                                                ifelse(Treatment_Best_Response%in%c("SD","PD"), "SD/PD", "NE"))), 
                               Clinical_Benefit = as.factor(Clinical_Benefit), 
                               PDL1_Expression_Status = as.factor(ifelse(PDL1_Expression%in%c("1", "Unknown"),
                                                            "Unknow/<1%", ">1%")))
                               #TMB_Status = ifelse(TMB_NonsynVariants>median(TMB_NonsynVariants, na.rm = TRUE), "High", "Low"),
                               #Neoantigen_Load_Status = ifelse(NeoCounts>median(NeoCounts, na.rm=TRUE), "High", "Low"),
                               #Neoantigen_Quality_Status = ifelse(NQuality>median(NQuality, na.rm=TRUE), "High", "Low"))

HR_Hellmann_1 <- plotHR(TidyInfo_Hellmann,
       expr="Gender+Histology+Age+Smoking_Status+Performance_status+PDL1_Expression_Status+Best_overall_response+Clinical_Benefit+TMB_Total+TMB_NonsynSNP+TMB_NonsynVariants+SingleNCounts+NQuality+NeoCounts+NeoC_Strong")
HR_Hellmann_2 <- plotHR(TidyInfo_Hellmann, expr="TMB_Total+TMB_NonsynSNP+TMB_NonsynVariants+SingleNCounts+NQuality+NeoCounts+NeoC_Strong")
ggcoxzph(cox.zph(HR_Hellmann_2$model))

save(sampleInfo_Sci_Rizvi, file = "test.RData")

HR_Sci_Rizvi_1 <- plotHR(sampleInfo_Sci_Rizvi,
                         expr="Gender+Histology+Age+Smoking_History+Smoking_Years+Smoke_Signature+Best_Overall_Response+Clinical_Benefit")


# compare mutation and neoantigens
p1_1 <- compareBoxplot(sampleInfo_Forde, x="Gender", y="TMB_Total") 
p1_1 <- ggpar(p1_1 , title="Forde etc. dataset", ylab="Total Mutation", font.main = 14)
p1_2 <- compareBoxplot(sampleInfo_Hellmann, x="Gender", y="TMB_Total")
p1_2 <- ggpar(p1_2 , title="Hellmann etc. dataset", ylab="Total Mutation", font.main = 14)
p1_3 <- compareBoxplot(sampleInfo_Sci_Rizvi, x="Gender", y="TMB_Total")
p1_3 <- ggpar(p1_3 , title="Science Rizvi etc. dataset", ylab="Total Mutation", font.main = 14)
p1_4 <- compareBoxplot(sampleInfo_JCO_Rizvi, x="Gender", y="TMB_Total")
p1_4 <- ggpar(p1_4 , title="JCO Rizvi etc. dataset", ylab="Total Mutation", font.main = 14)
library(cowplot)
p <- plot_grid(p1_1, p1_2, p1_3, p1_4, nrow=2, ncol=2, align = "v")
ggsave("Compare-Total-between-F-and-M.pdf", plot = p, width = 6, height = 5)


CorrMut_Hellmann <-  ggstatsplot::ggcorrmat(
        data = sampleInfo_Hellmann,
        corr.method = "spearman",                 # correlation method
        sig.level = 0.05,                        # threshold of significance
        cor.vars = c(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants, NQuality, NeoCounts, NeoC_Strong),         # a range of variables can be selected  
        cor.vars.names = c("Total Mutation", "Nonsynonymous Mutation", "Nonsynonymous Mutation + INDEL",
                           "Neoantigen Quality", "Neoantigen Load (ic50<500)", "Neoantigen Load (ic50<50)"),
        title = "",
        #subtitle = "Iris dataset by Anderson",
        caption = expression(
            paste(
                italic("Note"),
                ": X denotes correlation non-significant at ",
                italic("p "),
                "< 0.005; adjusted alpha"
            )
        )
    )
CorrMut_Forde <-  ggstatsplot::ggcorrmat(
    data = sampleInfo_Forde,
    corr.method = "spearman",                 # correlation method
    sig.level = 0.05,                        # threshold of significance
    cor.vars = c(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants, NQuality, NeoCounts, NeoC_Strong),         # a range of variables can be selected  
    cor.vars.names = c("Total Mutation", "Nonsynonymous Mutation", "Nonsynonymous Mutation + INDEL",
                       "Neoantigen Quality", "Neoantigen Load (ic50<500)", "Neoantigen Load (ic50<50)"),
    title = "",
    #subtitle = "Iris dataset by Anderson",
    caption = expression(
        paste(
            italic("Note"),
            ": X denotes correlation non-significant at ",
            italic("p "),
            "< 0.005; adjusted alpha"
        )
    )
)

CorrMut_Sci_Rizvi <-  ggstatsplot::ggcorrmat(
    data = sampleInfo_Sci_Rizvi,
    corr.method = "spearman",                 # correlation method
    sig.level = 0.05,                        # threshold of significance
    cor.vars = c(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants, NQuality, NeoCounts, NeoC_Strong),         # a range of variables can be selected  
    cor.vars.names = c("Total Mutation", "Nonsynonymous Mutation", "Nonsynonymous Mutation + INDEL",
                       "Neoantigen Quality", "Neoantigen Load (ic50<500)", "Neoantigen Load (ic50<50)"),
    title = "",
    #subtitle = "Iris dataset by Anderson",
    caption = expression(
        paste(
            italic("Note"),
            ": X denotes correlation non-significant at ",
            italic("p "),
            "< 0.005; adjusted alpha"
        )
    )
)

CorrMut_JCO_Rizvi <-  ggstatsplot::ggcorrmat(
    data = sampleInfo_Hellmann,
    corr.method = "spearman",                 # correlation method
    sig.level = 0.05,                        # threshold of significance
    cor.vars = c(TMB_Total, TMB_NonsynSNP, TMB_NonsynVariants),         # a range of variables can be selected  
    cor.vars.names = c("Total Mutation", "Nonsynonymous Mutation", "Nonsynonymous Mutation + INDEL"),
    title = "",
    #subtitle = "Iris dataset by Anderson",
    caption = expression(
        paste(
            italic("Note"),
            ": X denotes correlation non-significant at ",
            italic("p "),
            "< 0.005; adjusted alpha"
        )
    )
)


CorrMut_Forde <- CorrMut_Forde + theme_void(base_size=6)
ggsave(filename = "CorrMut_Forde.pdf", plot = CorrMut_Forde, width = 8, height = 7)
CorrMut_Hellmann <- CorrMut_Hellmann + theme_void(base_size=6)
ggsave(filename = "CorrMut_Hellmann.pdf", plot = CorrMut_Hellmann, width = 8, height = 7)
CorrMut_Sci_Rizvi <- CorrMut_Sci_Rizvi + theme_void(base_size=6)
ggsave(filename = "CorrMut_Sci_Rizvi.pdf", plot = CorrMut_Sci_Rizvi, width = 8, height = 7)
CorrMut_JCO_Rizvi <- CorrMut_JCO_Rizvi + theme_void(base_size=6)
ggsave(filename = "CorrMut_JCO_Rizvi.pdf", plot = CorrMut_JCO_Rizvi, width = 8, height = 7)
