# This code file used to clean all four immunotherapy datasets which collected from four difference papers.
# Author: ShixiangWang
# mail: <wangshx@shanghaitech.edu.cn>


setwd("G:/biodata/immunotherapyDatasets/")
library(tidyverse)
library(data.table)
#library(NMF)
library(maftools)

# ref_38genome = "G:/biodata/reference/hg38.fa"
# ref_19genome = "G:/biodata/reference/hg19.fa"
# load functions
source("C:/Users/wangshx/Desktop/data/neoQ/Allfunctions.R")

##################################
# clean the datasets one by one ##
##################################
# NOTE: To better clean the data, some raw data files been processed by hand.


#<<<<<<<<<<<<<<<<<<<<<< #1 begin
# datasets from <Rizvi, N. A., et al. "Cancer immunology. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. " Science 348.6230(2015):124-128.>
# Files download from cBioPortal and Scicence 

Science_Rizvi_clin <- read_csv("Science_Rizvi_clinical_info2.csv") # the clinical info file from cBioportal has something wrong with the smoking signature result, so use data from paper
Science_Rizvi_mut  <- fread("Science_Rizvi_mutations.txt")
Science_Rizvi_neos <- read_csv("Science_Rizvi_neoantigens.csv")

Science_Rizvi_mut$Tumor_Sample_Barcode[Science_Rizvi_mut$Tumor_Sample_Barcode == "R7495_2"] <- "R7495"
Science_Rizvi_clin_clean <- Science_Rizvi_clin %>% rename(Tumor_Sample_Barcode=Sample_ID)

Science_Rizvi_HLA <- Science_Rizvi_neos %>% select(Sample, HLA) %>% distinct() %>% filter(Sample!="SI9946_2")
Science_Rizvi_HLA$Sample[Science_Rizvi_HLA$Sample=="R7495_2"] <- "R7495"

Science_Rizvi_sampleInfo <- full_join(Science_Rizvi_clin_clean ,
                                      Science_Rizvi_HLA, 
                                      by=c("Tumor_Sample_Barcode"="Sample")) %>% 
    mutate(HLA = gsub(pattern = "([A-C])([0-9]{2})([0-9]{2})", 
                      replacement = "HLA-\\1*\\2:\\3", HLA)) %>% 
    mutate(HLA = gsub(pattern = ",,", replacement = ",", HLA), Tumor_Stage = "IV") %>% 
    setDT()

Science_Rizvi_maf <- maftools::read.maf(maf=Science_Rizvi_mut, clinicalData = Science_Rizvi_sampleInfo, isTCGA = FALSE)

#dir.create(path="Rdata")
save(Science_Rizvi_maf, file="Rdata/Science_Rizvi_maf.RData")
save(Science_Rizvi_sampleInfo, file="Rdata/Science_Rizvi_sampleInfo.RData")
# Science_Rizvi_MutSig <- autoMutSig(maf = Science_Rizvi_maf, ref_genome = ref_19genome)
# Science_Rizvi_Heter  <- autoTumorHeter(maf = Science_Rizvi_maf )

rm(Science_Rizvi_clin, Science_Rizvi_clin_clean, Science_Rizvi_Heter, Science_Rizvi_sampleInfo, Science_Rizvi_HLA, Science_Rizvi_mut, Science_Rizvi_neos, Science_Rizvi_MutSig, Science_Rizvi_maf);gc()
#<<<<<<<<<<<<<<<<<<<<<< #1 end


#<<<<<<<<<<<<<<<<<<<<<< #2 begin
# datasets from <Rizvi et al. "Molecular Determinants of Response to Anti–Programmed Cell Death (PD)-1 and Anti–Programmed Death-Ligand (PD-L)-Ligand 1 Blockade in Patients With Non–Small-Cell Lung Cancer Profiled With Targeted Next-Generation Sequencing. " JCO.>
# Files download from cBioPortal

JCO_Rizvi_clin <- read_csv("JCO_Rizvi_clinical_info.csv")
JCO_Rizvi_mut  <- fread("JCO_Rizvi_mutations.txt")

JCO_Rizvi_clin_clean <- JCO_Rizvi_clin %>% rename(Tumor_Sample_Barcode=SAMPLE_ID, 
                                                  Treatment_Type=TREATMENT_TYPE,
                                                  Age=AGE, Gender=SEX, Smoking_History=SMOKER,
                                                  Lines_of_TX=LINES_OF_TX,
                                                  PFS_Months=PFS_MONTHS, PFS_Event=PFS_STATUS, 
                                                  Clinical_Benefit=DURABLE_CLINICAL_BENEFIT, 
                                                  Patient_ID=PATIENT_ID, Gene_Panel=GENE_PANEL,
                                                  Mutation_Rate = MUTATION_RATE,
                                                  PDL1_Score=PDL1_SCORE,
                                                  Histology=CANCER_TYPE_DETAILED) %>% select(c(-13,-17))
JCO_Rizvi_clin_clean$Histology[JCO_Rizvi_clin_clean$Histology == "Large Cell Neuroendocrine Carcinoma"] <- "LCNC"
JCO_Rizvi_clin_clean$Histology[JCO_Rizvi_clin_clean$Histology == "Lung Adenocarcinoma"] <- "LUAD"
JCO_Rizvi_clin_clean$Histology[JCO_Rizvi_clin_clean$Histology == "Lung Squamous Cell Carcinoma"] <- "LUSC"
JCO_Rizvi_clin_clean$Histology[JCO_Rizvi_clin_clean$Histology == "Non-Small Cell Lung Cancer"] <- "NSCLC"
JCO_Rizvi_clin_clean <- JCO_Rizvi_clin_clean %>% mutate(
    Clinical_Benefit = ifelse(Clinical_Benefit=="YES", "DCB", "NDB"),
    PFS_Event = ifelse(PFS_Event == "Progressed", 1, 0)
)


JCO_Rizvi_sampleInfo <- JCO_Rizvi_clin_clean %>% mutate(Tumor_Stage = "IIIB-IV") %>%  setDT()
JCO_Rizvi_maf <- maftools::read.maf(maf=JCO_Rizvi_mut, clinicalData = JCO_Rizvi_sampleInfo, isTCGA = FALSE)

#dir.create(path="Rdata")
save(JCO_Rizvi_maf, file="Rdata/JCO_Rizvi_maf.RData")
save(JCO_Rizvi_sampleInfo, file="Rdata/JCO_Rizvi_sampleInfo.RData")
#autoMutSig(maf = JCO_Rizvi_maf, ref_genome = ref_19genome)
#autoTumorHeter(maf = JCO_Rizvi_maf )

rm(JCO_Rizvi_clin, JCO_Rizvi_clin_clean, JCO_Rizvi_maf, JCO_Rizvi_mut, JCO_Rizvi_sampleInfo);gc()
#<<<<<<<<<<<<<<<<<<<<<< #2 end

#<<<<<<<<<<<<<<<<<<<<<< #3 begin
# datasets from <Hellmann et al. "Genomic Features of Response to Combination Immunotherapy in Patients with Advanced NonSmall-Cell Lung Cancer. " Cancer Cell.>
# Files download from paper supplement files

Hellmann_clin <- read_csv("hellmann_clinical_info.csv")
Hellmann_mut  <- read_csv("hellmann_mutation.csv")
Hellmann_neos <- read_csv("hellmann_neoantigen.csv")

Hellmann_clin_clean <- Hellmann_clin %>% rename(Tumor_Sample_Barcode=UniqueSubjectIdentifier,
                                                Patient_ID=PatientID, Gender=Sex, 
                                                Smoking_History=Smoking_Status,
                                                PDL1_Expression=`PD-L1_expression_percent`,
                                                PFS_Months = PFS_time, PFS_Event=PFS_status,
                                                Treatment_Best_Response=BestOverallResponse,
                                                Clinical_Benefit=ClinicalBenefit,
                                                Nonsyn=Nonsynonymous_tumor_mutation,
                                                Neoantigen=Predicted_neoantigen_burden) %>% 
    mutate(HLA = paste(paste0("HLA-", HLA_A1), 
                       paste0("HLA-",HLA_A2), 
                       paste0("HLA-",HLA_B1), 
                       paste0("HLA-",HLA_B2), 
                       paste0("HLA-", HLA_C1), 
                       paste0("HLA-", HLA_C2), 
                       sep=","))

Hellmann_sampleInfo <- Hellmann_clin_clean %>% arrange(Patient_ID) %>% mutate(Tumor_Stage="IV") %>% setDT() 

Hellmann_mut$PatientID <- as.factor(Hellmann_mut$PatientID)
Hellmann_mut$Tumor_Sample_Barcode <- Hellmann_sampleInfo$Tumor_Sample_Barcode[Hellmann_mut$PatientID]
Hellmann_mut$NCBI_Build <- "GRCh37"
Hellmann_mut$Chromosome <- Hellmann_mut$chr
Hellmann_mut$Hugo_Symbol <- sub("([A-Za-z0-9]+);.*$","\\1", Hellmann_mut$gene_name)
Hellmann_mut$Strand <- "+"
Hellmann_mut$Start_Position <- Hellmann_mut$start
Hellmann_mut$Gene <- sub("([A-Za-z0-9]+);.*$","\\1", Hellmann_mut$gene_id)
Hellmann_mut$t_depth <- Hellmann_mut$tumor_depth
Hellmann_mut$t_alt_count <- Hellmann_mut$tumor_alt_depth
Hellmann_mut$t_vaf <- Hellmann_mut$tumor_vaf
Hellmann_mut$Reference_Allele <- Hellmann_mut$ref
Hellmann_mut$Tumor_Seq_Allele2 <- Hellmann_mut$alt
Hellmann_maf_min <- Hellmann_mut %>% select(Tumor_Sample_Barcode:Tumor_Seq_Allele2) %>% 
    mutate(Reference_Allele=ifelse(is.na(Reference_Allele),'-',Reference_Allele), 
           Tumor_Seq_Allele2=ifelse(is.na(Tumor_Seq_Allele2), "-", Tumor_Seq_Allele2))
write_tsv(Hellmann_maf_min, path = "Rdata/Hellmann_maf_min.maf")

Hellmann_maf_min2 <- Hellmann_maf_min %>% select(-Gene) %>% 
    mutate(End_Position = ifelse( Reference_Allele == '-', 
                                  Start_Position + 1,
                                  ifelse(Tumor_Seq_Allele2 == '-', 
                                         Start_Position + nchar(Reference_Allele) - 1, 
                                         Start_Position)))
write_tsv(Hellmann_maf_min2, path = "Rdata/Hellmann_maf_min2.maf")
# we need re-annotate this file with maf2maf.pl <https://github.com/mskcc/vcf2maf> and then process with maftools
# commands: ~/ProjectsManager/biotools/ensembl-vep$ ../vcf2maf/maf2maf.pl --input-maf ~/Desktop/Hellmann_maf_min.maf --output-maf ~/Desktop/Hellmann_reannaoted.maf --ref-fasta ~/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa  --vep-path=./ --tmp-dir $HOME/tmp

Hellmann_inputMAF <- fread(input = "Rdata/Hellmann.maf")
Hellmann_maf <- maftools::read.maf(maf=Hellmann_inputMAF, clinicalData = Hellmann_sampleInfo, isTCGA = FALSE)

save(Hellmann_maf, file="Rdata/Hellmann_maf.RData")
save(Hellmann_sampleInfo, file="Rdata/Hellmann_sampleInfo.RData")
#<<<<<<<<<<<<<<<<<<<<<< #3 end

#<<<<<<<<<<<<<<<<<<<<<< #4 begin
# datasets from <Forde et al. "Molecular Determinants of Response to Anti–Programmed Cell Death (PD)-1 and Anti–Programmed Death-Ligand (PD-L)-Ligand 1 Blockade in Patients With Non–Small-Cell Lung Cancer Profiled With Targeted Next-Generation Sequencing. " NEJM.>
# Files download from paper supplement files

Forde_clin <- read_csv("Forde_clinical_info.csv")
Forde_mut  <- read_csv("Forde_mutation.csv")
Forde_neos <- read_csv("Forde_neoantigens.csv")

Forde_clin_clean <- Forde_clin %>% dplyr::rename(Tumor_Sample_Barcode=SampleID, 
                                          Patient_ID=PatientID, Smoking_History=`Smoking History`,
                                          Smoking_Pack_Years = PackYears, 
                                          Histology=HistopathologicDiagnosis,
                                          Tumor_Stage=TNM)
Forde_neos <- Forde_neos[1:19584, -24]
Forde_neos_sum <- Forde_neos %>% dplyr::group_by(SampleID) %>% 
                        dplyr::summarise(Neoantigen_Burden=n(), 
                        HLA=paste(unique(HLA_Allele), collapse = ",")) %>% 
                        dplyr::mutate(SampleID=substr(SampleID, 1, 7)) %>% 
                        dplyr::rename(Tumor_Sample_Barcode=SampleID)

Forde_sampleInfo <- dplyr::full_join(Forde_clin_clean, Forde_neos_sum, by="Tumor_Sample_Barcode")%>% setDT()
Forde_sampleInfo <- Forde_sampleInfo %>% 
    dplyr::rename(Neoantigen = Neoantigen_Burden, Smoking_Years = Smoking_Pack_Years)


pattern1 <- "^chr[0-9XYM]+_[0-9]+-[0-9]+_([ACGT]{0,})_([ACGT]{0,})$"
Forde_maf_min <- Forde_mut %>% mutate(SampleID=substr(SampleID, 1, 7),
                                      NCBI_Build="GRCh37",
                                      Hugo_Symbol=Gene,
                                      Strand = "+",
                                      t_vaf=as.numeric(sub("([0-9]+)%", "\\1",DistinctMutantReads))/100,
                                      Chromosome=sub("^chr([0-9XYM]+)_[0-9]+.*$", "\\1", NucleotidePosition),
                                      Start_Position=sub("^chr[0-9XYM]+_([0-9]+)-[0-9]+.*$", "\\1", NucleotidePosition),
                                      End_Position=sub("^chr[0-9XYM]+_[0-9]+-([0-9]+).*$", "\\1", NucleotidePosition),
                                      Reference_Allele=ifelse(nchar(sub(pattern1, "\\1", NucleotidePosition))>0, sub(pattern1, "\\1", NucleotidePosition), "-"),
                                      Tumor_Seq_Allele2=ifelse(nchar(sub(pattern1, "\\2", NucleotidePosition))>0, sub(pattern1, "\\2", NucleotidePosition), "-")
                                      ) %>% rename(Tumor_Sample_Barcode=SampleID) %>% 
                                    select(Tumor_Sample_Barcode, NCBI_Build:Tumor_Seq_Allele2)

write_tsv(Forde_maf_min, path = "Rdata/Forde_maf_min.maf")

Forde_inputMaf <- fread(input="Rdata/Forde.maf")
Forde_maf <- maftools::read.maf(maf=Forde_inputMaf, clinicalData = Forde_sampleInfo, isTCGA = FALSE)

save(Forde_maf, file="Rdata/Forde_maf.RData")
save(Forde_sampleInfo, file="Rdata/Forde_sampleInfo.RData")
# 
#<<<<<<<<<<<<<<<<<<<<<< #4 end

rm(list = ls()); gc()
