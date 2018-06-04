# Functions

#' @function Summary functions for error bar plot.
#' @references http://www.jianshu.com/p/003138ac593b

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    #library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    library(plyr)
    
    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                           .fun = function(xx, col, na.rm) {
                               c(subjMean = mean(xx[,col], na.rm=na.rm))
                           },
                           measurevar,
                           na.rm
    )
    
    # Put the subject means with original data
    data <- merge(data, data.subjMean)
    
    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
        mean(data[,measurevar], na.rm=na.rm)
    
    # Remove this subject mean column
    data$subjMean <- NULL
    
    return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
    
    # Ensure that the betweenvars and withinvars are factors
    factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                         FUN=is.factor, FUN.VALUE=logical(1))
    
    if (!all(factorvars)) {
        nonfactorvars <- names(factorvars)[!factorvars]
        message("Automatically converting the following non-factors to factors: ",
                paste(nonfactorvars, collapse = ", "))
        data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
    }
    
    # Get the means from the un-normed data
    datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                       na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
    
    # Drop all the unused columns (these will be calculated with normed data)
    datac$sd <- NULL
    datac$se <- NULL
    datac$ci <- NULL
    
    # Norm each subject's data
    ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
    
    # This is the name of the new column
    measurevar_n <- paste(measurevar, "_norm", sep="")
    
    # Collapse the normed data - now we can treat between and within vars the same
    ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                        na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
    
    # Apply correction from Morey (2008) to the standard error and confidence interval
    #  Get the product of the number of conditions of within-S variables
    nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                    FUN.VALUE=numeric(1)))
    correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
    
    # Apply the correction factor
    ndatac$sd <- ndatac$sd * correctionFactor
    ndatac$se <- ndatac$se * correctionFactor
    ndatac$ci <- ndatac$ci * correctionFactor
    
    # Combine the un-normed means with the normed results
    merge(datac, ndatac)
}




##################################
# Functions used for prepare data#
##################################
##> Mutation Signature Analysis

autoMutSig <- function(maf, ref_genome, prefix="chr", add=TRUE, ignoreChr="chrM", useSyn=TRUE, n=6, nTry=6, plotBestFitRes=FALSE, minMut=5, useCNV=FALSE){
    tnm <- maftools::trinucleotideMatrix(maf=maf, ref_genome = ref_genome, ignoreChr = ignoreChr, prefix = prefix, add = add, useSyn = useSyn)
    signature <- maftools::extractSignatures(tnm, nTry = nTry, plotBestFitRes = plotBestFitRes)
    se <- maftools::signatureEnrichment(maf=maf, sig_res = signature)
    mut_signatures <- c(tnm, signature, se)
    return(mut_signatures)
}

##> Infer Heterogeneity for Tumor MAF
autoTumorHeter <- function(maf, tsb="ALL", top=5, vafCol=NULL, segFile=NULL, ignChr=NULL, minVaf=0, maxVaf=1, useSyn=FALSE){
    if(tsb=="ALL"){
        tsb=maftools::getSampleSummary(maf)$Tumor_Sample_Barcode
    }
    het <- maftools::inferHeterogeneity(maf=maf, tsb=tsb, top=top, segFile=segFile, ignChr=ignChr, minVaf=minVaf, maxVaf=maxVaf, useSyn=useSyn)
    math <- het$clusterData %>% select(Tumor_Sample_Barcode, MATH) %>% distinct()
    res <- list(heter_data=het, MATH=math)
    return(res)
}


# function summary data from Neoantigen Quality file
summaryNQres <- function(input=NULL){
    options(stringsAsFactors=FALSE)
    
    if(is.null(input)){
        stop("Error: please give your result file of neoantigen quality as input!")
    }
    
    neoQuality <- read.delim(file=input, header=TRUE)
    NQ_sm <- neoQuality %>% 
        group_by(Sample) %>% 
        summarize(MissenNCounts=n(),
                  NQuality=max(NeoantigenRecognitionPotential[Excluded != 1]),    
                  NQualityNotExculd=max(NeoantigenRecognitionPotential), 
                  NQMean=mean(NeoantigenRecognitionPotential[Excluded != 1]), 
                  NQMeanNotExlcud=mean(NeoantigenRecognitionPotential), 
                  NQSum=sum(NeoantigenRecognitionPotential[Excluded != 1]), 
                  NQSumNotExlcud=sum(NeoantigenRecognitionPotential))
    
    res <- list()
    res$rawData   <- neoQuality
    res$NQsummary <- NQ_sm
    
    return(res)
}

# function summary data from merged neoantigen file
summaryMergedNeos <- function(input=NULL){
    if(is.null(input)){
        stop("Error: please give your result file of merged neoantigens as input!")    
    }
    neoantigens <- read.delim(file=input, header=TRUE)
    Neo_sm <- neoantigens %>% 
        group_by(Sample.Name, Best.MT.Score.Method) %>%
        summarize(NeoCounts=n())
    
    res <- list()
    res$rawData <- neoantigens
    res$Neosummary <- Neo_sm
    
    return(res)
}

# get Tumor Mutation Burden of samples 
getSampleTMB <- function(maf){
    maf.silent <- maf@maf.silent
    sample.silent <- maf.silent[,.N, .(Tumor_Sample_Barcode)]
    sample.nonsilent <- getSampleSummary(maf)
    res <- dplyr::full_join(sample.silent, sample.nonsilent, by="Tumor_Sample_Barcode")
    res %>% mutate(TMB_Total=ifelse(!is.na(N), N+total, total), 
                   TMB_NonsynSNP=Missense_Mutation+Nonsense_Mutation,
                   TMB_NonsynVariants=total) %>% select(TMB_Total:TMB_NonsynVariants, Tumor_Sample_Barcode)
}

##> retrieve HLA information to compute neoantigen load/quality
generateHLAfile <- function(df, path, tsb="Tumor_Sample_Barcode", HLA="HLA"){
    df <- as.data.frame(df)
    Cols <- c(tsb, HLA)
    Allcols <- colnames(df)
    if(all(Cols %in% Allcols)){
        df <- df[, Cols]
        write_tsv(df, path=path)
    }else{
        stop("Please check your colnames.")
    }
}

###############################################
# Functions used for analyze or visulize data #
###############################################

# compare samples using boxplot and add significant levels
compareBoxplot <- function(df, x=NULL, y=NULL, label_name=c("p.format", "p.signif"), 
                           method=c("wilcox.test", "t.test", "anova", "kruskai.test")){
    label_name <- match.arg(label_name)
    method <- match.arg(method)
    df <- as.data.frame(df)
    name_sort <- names(table(df[, x]))
    df[, x] <- factor(df[, x], levels = name_sort)
    p <- ggboxplot(df, x=x, y=y, ggtheme = theme_pubr(base_size = 8))
    p + stat_compare_means(label = label_name, label.x.npc = "center", method = method) + theme(plot.title = element_text(hjust = 0.5))
}

# compare mutation profile in two-level group variable
compareMutPlot <- function(dat, group1="Gender", group2="Clinical_Benefit", 
                           value="TMB_Total", label_name="p.format", method = "wilcox.test"){
    require(ggpubr)
    dat <- as.data.frame(dat)
    #my_comparisons <- combn(names(table(dat[, group2])), 2, simplify = FALSE)
    #my_comparisons  <- list(c("DCB", "NDB"))
    name_sort <- names(table(dat[, group1]))
    dat[, group1] <- factor(dat[, group1], levels = name_sort)
    p <- ggboxplot(dat, x = group1, y = value,
                   color = group2, palette = c("red", "blue"),
                   add = "jitter", shape = group2, font.label = list(size=6), 
                   add.params = list(size=2), 
                   ggtheme = theme_pubr()) + theme(plot.title = element_text(hjust = 0.5))
    
    p + stat_compare_means(aes_string(group=group2), label = label_name, 
                           method = method)        # Add pairwise comparisons p-value
    # stat_compare_means(label.x.npc = "center", 
    #                    label.y.npc = "top" )                   # Add global p-value
    
    # res <- summarySE(data = dat, measurevar = value, groupvars = c(group1, group2))
    # ggplot(res, aes_string(x=group1, y=value, fill=group2)) + 
    #      geom_bar(position = position_dodge() ,stat = "identity") +
    #      geom_errorbar(aes_string(ymin=paste0(value, "-se"), ymax=paste0(value, "+se")),
    #                    width=.2,
    #                    position = position_dodge(.9)) + 
    #      theme_bw()  + scale_fill_npg()
}

# summary data by group variables, include min, max and median
groupSummary <- function(df, summarise_var=NULL, ...){
    summarise_var  <- enquo(summarise_var)
    if(summarise_var != quo(NULL)){
        group_var <- quos(...)
        #print(summarise_var)
        df %>% 
            group_by(!!! group_var) %>% 
            dplyr::summarize(n = n(), 
                             min = min(!!summarise_var), 
                             max = max(!!summarise_var), 
                             median = median(!!summarise_var))
    }else{
        stop("summarise_var can not be null!")
    }
}

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
