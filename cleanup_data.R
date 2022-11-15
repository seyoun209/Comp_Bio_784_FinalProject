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
library("ChAMP")

#ensembl=useMart("ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#save(ensembl,file="ensembl")
load("ensembl")

sub_methyl_3[,"Chromosome"] <- sub("^","chr",sub_methyl_3[,"Chromosome"])
methyl.ran <- GRanges(Rle(sub_methyl_3[,"Chromosome"]),IRanges(sub_methyl_3[,"Genomic_Coordinate"],sub_methyl_3[,"Genomic_Coordinate"]))


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter.re <- promoters(txdb,upstream=1500,downstream=500)
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

##################################
##bumphunters###########
##################################

library(minfi)
library(MEAL)
library(MultiDataSet)

promoter_met_mat <- sub_methyl_3[(is.element(sub_methyl_3[,"Hybridization-REF"],pro.met.re[,"methyl"])),]
rownames(promoter_met_mat) <- promoter_met_mat[,1]
mat <- promoter_met_mat[,-c(1:4)] %>% as.matrix()
grset <- makeGenomicRatioSetFromMatrix(mat)
save(grset,file="grset")
patient_coldata <- patient[is.element(patient[,21],tumor_status),]
patient_coldata <- patient_coldata[(!is.element(patient_coldata[,2],subset_nm)),]
colnames(patient_coldata) <- gsub("\\s", "", colnames(patient_coldata))
patient_coldata[,21] <- gsub("\\s", "", patient_coldata[,21])
patient_coldata$NoninvasiveBladderCancerTreatmentType %<>% gsub("\\(|)", "", .) %>% gsub("\\[|]","",.) %>% gsub("\\s", "",.)
colData(grset) <- cbind(colData(grset),patient_coldata[ order(match(patient_coldata[,2],rownames(colData(grset)))), ])

status <- pData(grset)$PersonNeoplasmStatus
treatment <- pData(grset)$NoninvasiveBladderCancerTreatmentType


deseignmatrix <- model.matrix(~status+treatment)
deseignmatrix2 <- model.matrix(~status)
dmrs <- bumphunter(grset, design =deseignmatrix, type="Beta",cutoff=0.1)
dmrs2 <- bumphunter(grset, design =deseignmatrix2, type="Beta",cutoff=0.1)
#save(grset,file="grset")
res <- runPipeline(set = grset, variable_names = "PersonNeoplasmStatus")
#save(res,file="res")

targetRange <- GRanges("chr22:45608465- 45608516")

myDMR <- champ.DMR(beta=assays(grset)$Beta,pheno=grset$PersonNeoplasmStatus,method="Bumphunter")
#head(myDMR$DMRcateDMR)
DMR.GUI(DMR=myDMR,beta=assays(grset)$Beta,pheno=grset$PersonNeoplasmStatus,arraytype="450K")
knitr::include_graphics("Figure/DMR.GUI.jpg")
CpG.GUI(CpG=rownames(assays(grset)$Beta),arraytype="450K")
QC.GUI(beta=assays(grset)$Beta,arraytype="450K")
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)


#####
geneset <- c("PRKCZ","CDH5","CTNNB1","APC2","TLE2")
tumor_free <- assay(grset)[,which(grset$PersonNeoplasmStatus == "TUMORFREE" )]
withtumor <- assay(grset)[,which(grset$PersonNeoplasmStatus == "WITHTUMOR" )]
ppi_gene_mat <- promoter_met_mat[is.element(promoter_met_mat[,2],geneset),]

for (i in 1:length(geneset)){
  print(i)
  onegene <- ppi_gene_mat[which(ppi_gene_mat[,2] == geneset[i]),]
  cpgid <- rownames(onegene)
  tfree_subset <- tumor_free[is.element(rownames(tumor_free),cpgid),]
  wtiht_subset <- withtumor[is.element(rownames(withtumor),cpgid),]
  rowmean_control <- cbind(t(rbind(rowMeans(tfree_subset))),rep("tumor_free"))
  rowmean_case <- cbind(t(rbind(rowMeans(wtiht_subset))),rep("with_tumor"))
  combine.mat <- rbind(rowmean_control,rowmean_case) %>% as.data.frame()
  colnames(combine.mat) <- c("mean_of_beta","group")
  p <- ggplot(combine.mat, aes(x=combine.mat$group, y=as.double(combine.mat$mean_of_beta),fill=group)) + geom_boxplot()+scale_fill_brewer(palette="Blues")+theme_classic()
  png(paste("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/03.plots","/",geneset[i],".png",sep=""))
  print(p)
  dev.off()
}

for (i in 1:length(geneset)){
  print(i)
  rsem_subset <- rsem[is.element(rsem[,1],geneset[i]),]
  rsem_test <- t(rsem_subset[-c(1:2)])
  onegene <- ppi_gene_mat[which(ppi_gene_mat[,2] == geneset[i]),]
  cpgid <- rownames(onegene)
  tfree_subset <- tumor_free[is.element(rownames(tumor_free),cpgid),]
  wtiht_subset <- withtumor[is.element(rownames(withtumor),cpgid),]
  colmean_control <- cbind(t(rbind(colMeans(tfree_subset))),rep("tumor_free"))
  colmean_case <- cbind(t(rbind(colMeans(wtiht_subset))),rep("with_tumor"))
  combine.mat <- rbind(colmean_control,colmean_case) %>% as.data.frame()
  combinedall <- merge(combine.mat,rsem_test,by='row.names',all.x=FALSE) %>% as.data.frame()
  colnames(combinedall) <- c("sample","beta","group","expression")
  p <- ggplot(combinedall, aes(x=as.double(combinedall$beta), y=as.double(combinedall$expression),color=combinedall$group)) + geom_point()+geom_smooth(method=lm, se=FALSE)+scale_color_brewer(palette="Paired")+theme_minimal()
  png(paste("/pine/scr/s/e/seyoun/project_784/Comp_Bio_784_FinalProject/03.plots","/",geneset[i],"corr",".png",sep=""))
  print(p)
  dev.off()
}

