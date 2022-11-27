library(limma)
library(dplyr)
library(readr)
library(ggplot2)

theme_set(theme_bw())


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
genes.select.up = working.data$gene[as.numeric(row.names(top.table %>% filter(logFC>0)))]
genes.select.down = working.data$gene[as.numeric(row.names(top.table %>% filter(logFC<0)))]
genes.final = meta.data %>% filter(gene %in% genes.select)

top.table$gene = genes.select

# save it
write_csv(genes.final, 'diff_meth_genes.csv')


# looking at residual plots
row_num = 5
fitted.values = limma.fit.eb$design %*% limma.fit.eb$coefficients[row_num,]
plot(fitted.values, working.data[row_num, 2:ncol(working.data)])


# Plots -------------------------------------------------------------------

# volcano plot
volcanoplot(limma.fit.eb, coef='tumor', names=working.data$gene,
            highlight = 135)

# do manual volcano plot

# get all results
all.res = topTable(limma.fit.eb, coef='tumor', sort.by='p', number=100000)
all.res$highlight = all.res$adj.P.Val < 0.2

# plot
gg_volcano = ggplot(all.res %>% mutate(P.Value = -log10(P.Value)), 
                    aes(x=logFC, y=P.Value, color=highlight)) +
  geom_point(size = 0.7) +
  scale_color_manual(values = c('black', 'blue')) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  scale_y_continuous(limits = c(0, 6.2)) +
  labs(x = 'beta coefficient',
       y = '-log10(p-value)') +
  theme(legend.position = 'none')

# pca
working.data.select = working.data %>% filter(gene %in% genes.select)
pc.0 = prcomp(t(working.data.select[,2:ncol(working.data.select)]))
dat = data.frame(pc.0$x[,1:2], tumor = patient$tumor_status)
ggplot(dat, aes(x=PC1, y=PC2, col=tumor)) +
  geom_point() +
  labs(title='pca on 135 identified genes')


# mRNA --------------------------------------------------------------------

# read in data
mrna = read_delim('data_mrna_seq_v2_rsem.txt')

# subset to match working.data
mrna = mrna %>%
  rename(gene = Hugo_Symbol) %>%
  rename_all(function(x) gsub('-01$', '', x)) %>%
  select_at(colnames(working.data))

# do log(1+x) transform
mrna[,2:ncol(mrna)] = log(1+mrna[,2:ncol(mrna)])

# do pca
pc.1 = prcomp(t(mrna[,2:ncol(mrna)]))
dat = data.frame(pc.1$x[,1:2], tumor = patient$tumor_status)
ggplot(dat, aes(x=PC1, y=PC2, col=tumor)) +
  geom_point() +
  labs(title='pca on transformed gene expression, all genes')

# subset to match 135 identified genes
# only 120 in the dataset
mrna = mrna %>%
  filter(gene %in% genes.select)

# pca
pc.2 = prcomp(t(mrna[,2:ncol(mrna)]))
dat = data.frame(pc.2$x[,1:2], tumor = patient$tumor_status)
ggplot(dat, aes(x=PC1, y=PC2, col=tumor)) +
  geom_point() +
  labs(title='pca on transformed gene expression, 135 identified genes')


# correlation b/t gene expression and diff methylation
mrna.avg = cbind(mrna$gene, rowMeans(mrna[,2:ncol(mrna)])) %>% as.data.frame()
top.table.2 = top.table %>%
  filter(gene %in% mrna$gene)
top.table.2 = top.table.2[match(mrna.avg$V1, top.table.2$gene), ]
cor(top.table.2$logFC, as.numeric(mrna.avg$V2))

dat = data.frame(diff_meth = top.table.2$logFC, 
                 gene_exp = as.numeric(mrna.avg$V2))
ggplot(data=dat, aes(x=diff_meth, y=gene_exp)) +
  geom_point() +
  labs(x = 'Differential methylation (beta)',
       y = 'Average gene expression, log(1+x)')
