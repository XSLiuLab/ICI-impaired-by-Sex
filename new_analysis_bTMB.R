# bTMB as supple
# <https://www.nature.com/articles/s41591-018-0134-3#Sec19>

p_load(tidyverse, readxl)

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
                                    color = groupN)) + geom_roc(show.legend = TRUE)
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


# bTMB.df = read_xlsx("G:/biodata/datasets/bTMB_Gandara.etc.xlsx", sheet = 1)
cli.poplar = read_xlsx("G:/biodata/datasets/bTMB_Gandara.etc.xlsx", sheet = 2)
cli.oak = read_xlsx("G:/biodata/datasets/bTMB_Gandara.etc.xlsx", sheet = 3)

# bTMB = bTMB.df %>% 
#     filter(effect %in% c("missense", "nonsense"))

cli.poplar2 = cli.poplar %>% 
    filter(TRT01P != "Docetaxel" & BEP == "Y")
cli.oak2 = cli.oak %>% 
    filter(TRT01P != "Docetaxel" & BEP == "Y")

# cli.oak2 = cli.oak2 %>% filter(!BCOR %in% c(".", "NE")) 
# cli.poplar2 = cli.poplar2 %>% filter(!BCOR %in% c(".", "NE")) 

cli.poplar2$trial = "poplar"
cli.oak2$trial = "oak"

cli = bind_rows(select(cli.poplar2, trial, PtID, btmb, BAGE, SEX, HIST, TOBHX, KRASMUT:OS.CNSR, BCOR),
                select(cli.oak2, trial, PtID, btmb, BAGE, SEX, HIST, TOBHX, KRASMUT:OS.CNSR, BCOR))

table(cli$BCOR)


# validate results of paper 
cli.poplar3 = cli.poplar %>% 
    filter(BEP == "Y")
cli.oak3 = cli.oak %>% 
    filter(BEP == "Y")

cli.poplar3$trial = "poplar"
cli.oak3$trial = "oak"

cli_paper = bind_rows(cli.poplar3 %>% select(trial, PtID, btmb, BAGE, SEX, HIST, TOBHX, KRASMUT:OS.CNSR, BCOR, TRT01P),
                      cli.oak3 %>% select(trial, PtID, btmb, BAGE, SEX, HIST, TOBHX, KRASMUT:OS.CNSR, BCOR, TRT01P))


cli_paper = cli_paper %>% 
    mutate(MutationStatus = ifelse(as.numeric(btmb) > 16, "High", "Low"))

fit = surv_fit(Surv(PFS, event = PFS.CNSR == 0) ~ TRT01P, 
               data = cli_paper %>% filter(as.numeric(btmb) > 16))
ggsurvplot(fit, pval = TRUE)
fit = surv_fit(Surv(OS, event = OS.CNSR == 0) ~ TRT01P, 
               data = cli_paper %>% filter(as.numeric(btmb) > 16))
ggsurvplot(fit, pval = TRUE)


fit = surv_fit(Surv(PFS, event = PFS.CNSR == 0) ~ MutationStatus,
               data = cli_paper)

ggsurvplot(fit, pval = TRUE)
cli_paper = cli_paper %>% 
    mutate(MutationStatus = ifelse(as.numeric(btmb) > 4*1.1, "High", "Low"))

# cli_paper = cli_paper %>% 
#     mutate(MutationStatus = ifelse(as.numeric(btmb) > 0, "High", "Low"))
fit = surv_fit(Surv(PFS, event = PFS.CNSR == 0) ~ MutationStatus,
               data = cli_paper %>% filter(TRT01P!="Docetaxel"))

ggsurvplot(fit, pval = TRUE)

# forest plot
cli_cox = cli_paper %>% 
    mutate(nbTMB = as.numeric(btmb) / 1.1,
           Age = BAGE,
           Sex = factor(SEX, levels = c("M", "F")),
           HIST = factor(HIST, levels = c("SQUAMOUS", "NON-SQUAMOUS")),
           Smoking = factor(TOBHX, levels = c("NEVER", "PREVIOUS", "CURRENT")),
           Treat = factor(TRT01P, levels = c("Docetaxel", "MPDL3280A")),
           Mutation = factor(MutationStatus, levels = c("Low", "High")))

fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + Sex + HIST + Smoking + Treat + nbTMB, data = cli_cox)
ggforest(fit)    
fit = coxph(Surv(PFS, PFS.CNSR==0) ~ Age + Sex + HIST + Smoking + Treat + Mutation, data = cli_cox)
ggforest(fit)

fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + Sex + HIST + Smoking + nbTMB,
            data = cli_cox %>% filter(TRT01P!="Docetaxel"))
ggforest(fit)    

fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + Sex + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P!="Docetaxel"))
ggforest(fit)  


fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + Sex + HIST + Smoking + nbTMB,
            data = cli_cox %>% filter(TRT01P=="Docetaxel"))
ggforest(fit)    

fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + Sex + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel"))
ggforest(fit)  

# male vs female
fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "M"))
ggforest(fit)  

fit = coxph(Surv(OS, OS.CNSR==0) ~ Age + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "F"))
ggforest(fit)  

fit = coxph(Surv(PFS, PFS.CNSR==0) ~ Age + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "M"))
ggforest(fit)  

fit = coxph(Surv(PFS, PFS.CNSR==0) ~ Age + HIST + Smoking + Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "F"))
ggforest(fit)  

# 
fit = coxph(Surv(OS, OS.CNSR==0) ~ Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "M"))
ggforest(fit)  

fit = coxph(Surv(OS, OS.CNSR==0) ~ Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "F"))
ggforest(fit)  

#
fit = coxph(Surv(PFS, PFS.CNSR==0) ~ Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "M"))
ggforest(fit)  

fit = coxph(Surv(PFS, PFS.CNSR==0) ~ Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "F"))
ggforest(fit) 


fit = survfit(Surv(PFS, PFS.CNSR==0) ~ Mutation,
            data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "M"))
ggsurvplot(fit, pval = TRUE)  

fit = survfit(Surv(PFS, PFS.CNSR==0) ~ Mutation,
              data = cli_cox %>% filter(TRT01P=="Docetaxel", Sex == "F"))
ggsurvplot(fit, pval = TRUE) 

########### DCB definition
# way1 
df = cli %>% 
    mutate(nbTMB = as.numeric(btmb)/1.1, 
           Clinical_Benefit = case_when(
               BCOR == "CR" ~ "DCB",
               (BCOR %in% c("PR", "SD")) & PFS > 6 ~ "DCB",
               (BCOR == "PD") & PFS <= 6 ~ "NDB",
               (BCOR %in% c("PR", "SD")) & PFS <= 6 ~"NDB",
               TRUE ~ NA_character_
               # TRUE ~ "NDB"
           ),
           Clinical_Benefit = factor(Clinical_Benefit))
# way2
df = cli %>% 
    mutate(nbTMB = as.numeric(btmb)/1.1, 
           Clinical_Benefit = case_when(
               BCOR %in% c("CR", "PR") ~ "DCB",
               TRUE ~ "NDB"
           )) 


# df_oak %>% select(OS, OS.CNSR, PFS, PFS.CNSR, Clinical_Benefit, BCOR, bTMB, SEX) %>% arrange(desc(PFS)) %>% View()
library(ggpubr)

compareMutPlot(df, group1 = "SEX", value = "nbTMB")

# DCB and NDB
df %>% filter(!is.na(Clinical_Benefit)) %>% 
    ggboxplot(x="SEX", y="nbTMB",  color="Clinical_Benefit",
              xlab = "Clinical Benefit",  palette = c("red", "blue"),
              ylab = "blood TMB") + 
    rotate_x_text(angle = 45) + 
    stat_compare_means(aes_string(group="Clinical_Benefit"), label = "p.format", 
                       method = "wilcox.test")
    # geom_hline(yintercept = mean(df_oak$bTMB, na.rm=TRUE), linetype=2) + 


df %>% 
    ggboxplot(x="BCOR", y="nbTMB",  color="BCOR", add="jitter", xlab = "Response Status", 
              ylab = "blood TMB", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df$nbTMB, na.rm=TRUE), linetype=2) + 
    stat_compare_means(method = "anova", label.y=70, label.x = 1)




df_2 = df %>% mutate(Cutoff = factor(ifelse(nbTMB > 4, "High", "Low"),
                                     levels = c("Low", "High") )) %>% 
    mutate(OS.event = ifelse(OS.CNSR == 0, 1, 0),
        PFS.event = ifelse(PFS.CNSR == 0, 1, 0),
        Gender = case_when(
        SEX == "F" ~ "Female",
        SEX == "M" ~ "Male",
        TRUE ~ NA_character_)) %>% mutate(
            Gender = factor(Gender, levels = c("Male", "Female"))
        )

plotROC(df_2, nbTMB, Clinical_Benefit, Gender, all = FALSE)

library(survival)
summary(coxph(Surv(OS, OS.event) ~ Cutoff,
              data = subset(df_2, Gender=="Male")))
summary(coxph(Surv(OS, OS.event) ~ Cutoff,
              data = subset(df_2, Gender=="Female")))

summary(coxph(Surv(PFS, PFS.event) ~ Cutoff,
              data = subset(df_2, Gender=="Male")))
summary(coxph(Surv(PFS, PFS.event) ~ Cutoff,
              data = subset(df_2, Gender=="Female")))

surv_plot <- function(fit, ...){
    p <- ggsurvplot(fit, 
                    font.tickslab=14,
                    #legend="right", 
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
                    pval.size = 6, ...)
    ggpar(
        p,
        font.xtickslab = 20,
        font.ytickslab = 20 
    )
}

### NSCLC
# library(survminer)
# # fit = surv_fit(Surv(OS, OS.event) ~ Cutoff, data = as.data.frame(df_2), group.by = "Gender")
# # ggsurvplot_list(fit, data = df_2, pval = TRUE)
# 
# 
# fit1 <- survfit(Surv(OS, OS.event) ~ Cutoff,
#                       data = subset(df_2, Gender=="Male"))
# fit2 <- survfit(Surv(OS, OS.event) ~ Cutoff,
#                 data = subset(df_2, Gender=="Female"))
# surv_plot(fit1, data = subset(df_2, Gender=="Male"))
# surv_plot(fit2, data = subset(df_2, Gender=="Female"))
# 
# 
# 
# select_cf = function(x, method = "OS"){
#     df_2 = df %>%
#         mutate(Cutoff = factor(ifelse(nbTMB > x, "High", "Low"),
#                                          levels = c("Low", "High") )) %>% 
#         mutate(OS.event = ifelse(OS.CNSR == 0, 1, 0),
#                PFS.event = ifelse(PFS.CNSR == 0, 1, 0),
#                Gender = case_when(
#                    SEX == "F" ~ "Female",
#                    SEX == "M" ~ "Male",
#                    TRUE ~ NA_character_)) %>% mutate(
#                        Gender = factor(Gender, levels = c("Male", "Female"))
#                    )
#     df_male = df_2 %>% filter(Gender == "Male")
#     df_female = df_2 %>% filter(Gender == "Female")
#     
#     if(method == "OS"){
#         fit1 <- survfit(Surv(OS, OS.event) ~ Cutoff,
#                         data = df_male)
#         fit2  <- survfit(Surv(OS, OS.event) ~ Cutoff,
#                          data = df_female)
#     }else{
#         fit1 <- survfit(Surv(PFS, PFS.event) ~ Cutoff,
#                         data = df_male)
#         fit2  <- survfit(Surv(PFS, PFS.event) ~ Cutoff,
#                          data = df_female)
#     }
# 
#     p1 = list()
#     p1[[1]] = ggsurvplot(fit1, data = df_male, pval = TRUE)
#     p1[[2]] = ggsurvplot(fit2, data = df_female, pval = TRUE)
#     # p1[[1]] = surv_plot(fit1, data = df_male)
#     # p1[[2]] = surv_plot(fit2, data = df_female)
#     arrange_ggsurvplots(p1, print = TRUE,  ncol = 2, nrow = 1)
#     
# }
# 
# select_cf(4, "PFS")
# select_cf(15, "PFS")
# 
# 
# select_cf(4)
# select_cf(15)
# 
# surv_cutpoint(df_2 %>% filter(Gender=="Female"), "PFS", "PFS.event", variables = "nbTMB")
# surv_cutpoint(df_2 %>% filter(Gender=="Male"), "PFS", "PFS.event", variables = "nbTMB")
