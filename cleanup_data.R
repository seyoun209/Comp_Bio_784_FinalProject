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

# This step is removing probes with row containing missing values: 102,702 probes are 

library(tidyr)
library(dplyr)
sub_methyl_3 <- sub_methyl_2 %>% drop_na()


##################################
##promoter region genes###########
##################################

library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(tidyverse)

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)


sub_methyl_3[,"Chromosome"] <- sub("^","chr",sub_methyl_3[,"Chromosome"])
methyl.ran <- GRanges(Rle(sub_methyl_3[,"Chromosome"]),IRanges(sub_methyl_3[,"Genomic_Coordinate"],sub_methyl_3[,"Genomic_Coordinate"]))


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter.re <- promoters(txdb,upstream=1500,downstream=200)
#promoter.methyl <- as.matrix(findOverlaps(methyl.ran,promoter.re))

fi.cns <- c("TXID","TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
txTable <- gsub(" ","",as.matrix(biomaRt::select(txdb,keys=names(promoter.re),columns=fi.cns,keytype="TXNAME")))
txTable <- txTable[!is.na(txTable[,2]),]

map.id <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters="entrezgene_id", values=txTable[,"GENEID"],mart=ensembl)
mer.txTable <- unique(as.matrix(merge(txTable,map.id,by.x="GENEID",by.y="entrezgene_id")))
rownames(mer.txTable) <- mer.txTable[,"TXNAME"]


inter.pro.methyl <- as.matrix(findOverlaps(methyl.ran,promoter.re))
inter.pro.methyl <- inter.pro.methyl[is.element(names(promoter.re)[inter.pro.methyl[,2]],rownames(mer.txTable)),]
pro.met.re <- cbind(mer.txTable[names(promoter.re)[inter.pro.methyl[,2]],"hgnc_symbol"],sub_methyl_3[inter.pro.methyl[,1],c(1,3,4)])
colnames(pro.met.re) <- c("gene","methyl","chr","pos")
pro.met.re <- unique(pro.met.re[which(pro.met.re[,"gene"] != ""),])

#head(pro.met.re)


# for (i in 1:length(pro.met.re[,1])){
#   ea.met.id <- pro.met.re[i,1]
#   ea.exp.id <- pro.met.re[i,2]
#   ea.methyl <- sub_methyl_3[sub_methyl_3[,1] == ea.exp.id,]
#   ea.exp <- ea.exp.id

##################################
##bumphunters###########
##################################



pos <- list(pos1=seq(1,1000,35),pos2=seq(2001,3000,35),pos3=seq(1,1000,50))
chr <- rep(paste0("chr",c(1,1,2)), times = sapply(pos,length))
pos <- unlist(pos, use.names=FALSE)
cl <- clusterMaker(chr, pos, maxGap = 300)

