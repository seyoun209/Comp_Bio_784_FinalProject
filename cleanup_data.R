library(limma)
library(data.table)
library(ggplot2)
library(bumphunter)



rsem <- fread("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/01.data/data_mrna_seq_v2_rsem.txt.gz") |> as.data.frame()
sample <- fread("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/01.data/data_clinical_sample.txt.gz") |> as.data.frame()
patient <- fread("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/01.data/data_clinical_patient.txt.gz") |> as.data.frame()
methyl <- fread("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/01.data/TCGA-BLCA_DNAmethylation450k_Methyl.txt_v2.gz") |> as.data.frame()
methly_beta <- methyl[,-c(1:4)]
methyl_nodup <- as.data.frame(sapply(unique(names(methly_beta)),function(col) rowMeans(methly_beta[names(methly_beta) == col])))
                              
colnames(rsem) <- c(colnames(rsem)[1:2],gsub('.{3}$',"",colnames(rsem)[-c(1:2)]))
tumor_status <- c("TUMOR FREE","WITH TUMOR")
patient_subset <- patient[,c(2,21)][is.element(patient[,21],tumor_status),]

sub_methyl <- cbind(methyl[,1:4],methyl_nodup[,(is.element(colnames(methyl_nodup), patient_subset[,1]))])
sub_rsem <- cbind(rsem[,1:2],rsem[,(is.element(colnames(rsem), patient_subset[,1]))])

subset_nm <- colnames(sub_methyl[,!is.element(colnames(sub_methyl),colnames(sub_rsem))])[-c(1:4)]
#three samples are not overlapping with methylation and gene expression. Therefore, to match up, we subtract "TCGA-XF-A9SG" "TCGA-ZF-A9RG" "TCGA-XF-AAMF".

sub_methyl_2 <-sub_methyl[,!is.element(colnames(sub_methyl),subset_nm)]











