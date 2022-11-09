library(limma)
library(dplyr)
library(readr)


# Prepare data ------------------------------------------------------------

# # saved these data from cleanup_data.R
# saveRDS(pro.met.re, 'promoter_regions.rds')
# saveRDS(sub_methyl_3, 'sub_methyl.rds')

# load data saved from end of cleanup_data.R
pro.met.re = readRDS('promoter_regions.rds')
sub_methyl_3 = readRDS('sub_methyl.rds')

# subset methylation data according to pro.met.re
# and remove old, unnecessary columns
working.data = left_join(pro.met.re, sub_methyl_3, by=c('methyl'='Hybridization-REF'))
working.data = working.data %>% dplyr::select(-Gene_Symbol, -Chromosome, -Genomic_Coordinate)

# only need working.data
rm(list=setdiff(ls(), "working.data"))

# stash row information
meta.data = working.data[,1:4]

# average methylation by gene
working.data = working.data %>%
  dplyr::select(-methyl, -chr, -pos) %>%
  group_by(gene) %>%
  summarize_all(mean)

# load patient data
clinical.patient = read_delim('data_clinical_patient.txt', skip=4)
clinical.sample = read_delim('data_clinical_sample.txt', skip=4)

# combine the data and subset to those in working.data
patient = inner_join(clinical.patient, clinical.sample, by='PATIENT_ID')
patient = patient %>%
  filter(PATIENT_ID %in% colnames(working.data)) %>%
  group_by(PATIENT_ID) %>%  # one patient has two samples, remove second one
  dplyr::slice(n=1) %>%
  ungroup()
rm(list=c('clinical.patient', 'clinical.sample'))

# set column names lower case
colnames(patient) = tolower(colnames(patient))

# put patient data in same order as methylation data
patient = patient[match(colnames(working.data[,-1]), patient$patient_id), ]


# Limma -------------------------------------------------------------------


# patient covariate matrix
design.matrix = model.matrix(~ tumor_status + sex + age,
                             data = patient)
colnames(design.matrix) = c('(Intercept)', 'tumor', 'male', 'age')

# fit models
limma.fit = lmFit(working.data[,-1], design=design.matrix)
limma.fit.eb = eBayes(limma.fit)

# look at results
top.table = topTable(limma.fit.eb, coef='tumor', sort.by='p', p.value=0.2, number=10000)

# extract gene names and meta information
genes.select = working.data$gene[as.numeric(row.names(top.table))]
genes.final = meta.data %>% filter(gene %in% genes.select)

# save it
write_csv(genes.final, 'diff_meth_genes.csv')
