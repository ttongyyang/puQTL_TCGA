rm(list = ls())
GTEx_TCGA <- read.delim("../puQTL/5.analysis/normal/GTEx-TCGA.txt",header = F)
all_share <- c()
all_TCGA_specific_gene <- c()
all_GTEx_specific_gene <- c()
promoter_anno <- read.delim("../puQTL/GTEX/tmp/proActiv/data/promoterCoordinates.gencode.38.bed",header = F)
promoter_anno <- promoter_anno[,c(3,4,5,6)]
promoter_anno$V4 <- gsub(promoter_anno$V4,pattern = "prmtr\\.",replacement = "")
promoter_anno$V5 <- gsub(promoter_anno$V5,pattern = "\\..*",replacement = "")
hg19_promoter <- read.delim("../puQTL/GTEX/08.analysis/promoterCoordinates.gencode.v19.csv",header = T)
hg19_promoter <- hg19_promoter[,c(2,6)]
hg19_promoter$promoterId<- gsub(hg19_promoter$promoterId,pattern = "prmtr\\.",replacement = "")
colnames(hg19_promoter)[2] <- "promoter"
cancer_specific <- c()
share <- c()
normal_specific <- c()
all_result <- data.frame()
for (i in 1:nrow(GTEx_TCGA)) {
  GTEx <- read.delim(paste0("../puQTL/GTEX_new_sample/4.MatrixEQTL/1.significant/relative/",GTEx_TCGA[i,1]),header = F) 

  TCGA <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/TCGA-",GTEx_TCGA[i,2],".cis"),header = T)
  TCGA <- TCGA[,c(1,3,6,7,8,10)]
  colnames(promoter_anno)[2] <- "promoter"
  promoter_anno$promoter <- as.numeric(promoter_anno$promoter)
  TCGA <- dplyr::inner_join(TCGA,promoter_anno)
  hg19_promoter$promoter <- as.numeric(hg19_promoter$promoter)
  TCGA <- dplyr::inner_join(TCGA,hg19_promoter)
  TCGA$promoter <- paste0(TCGA$V3,"-",TCGA$promoter_gene)
  GTEx <- GTEx[,c(1,5,6)]
  GTEx$promoter <- paste0(str_split_fixed(GTEx$V5, "-",4)[,4],"-",str_split_fixed(GTEx$V5, "-",4)[,2])
  GTEx$QTL <- paste0(GTEx$V1,":",GTEx$promoter)
  TCGA$QTL <- paste0(TCGA$SNP_id,":",TCGA$promoter)
  share_QTL <- TCGA[TCGA$QTL%in%GTEx$QTL,]
  share <- c(share,share_QTL$QTL)
  cancer_specific <- c(cancer_specific,TCGA$QTL)
  normal_specific <- c(normal_specific,GTEx$QTL)
  part_result <- data.frame(share_QTL=nrow(share_QTL),TCGA_specific=nrow(TCGA)-nrow(share_QTL),GTEX_specific=nrow(GTEx)-nrow(share_QTL),cancer = paste0(GTEx_TCGA[i,2],":",GTEx_TCGA[i,1]))
  all_result <- rbind(all_result,part_result)
}
all_result_share <- all_result[,c(1,4)]
colnames(all_result_share) <- c("type","cancer")
all_result_share$tmp <- "all_result_share"
TCGA_specific <- all_result[,c(2,4)]
colnames(TCGA_specific) <- c("type","cancer")
TCGA_specific$tmp <- "TCGA_specific"
GTEx_specific <- all_result[,c(3,4)]
colnames(GTEx_specific) <- c("type","cancer")
GTEx_specific$tmp <- "GTEx_specific"
all_plot <- rbind(all_result_share,TCGA_specific)
length(unique(share))
6066
length(unique(cancer_specific))
131257-6066
length(unique(normal_specific))
480711-6066

all_result

ggplot(all_plot, aes(type,cancer, fill=tmp))+geom_bar(stat='identity',position='fill') +
  labs(x = 'tissue',y = 'Fraction') +theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20)) + coord_flip()



ggplot(all_result,aes(share_QTL,cancer))+
  geom_bar(stat="identity")


ggplot(all_result,aes(share_QTL,cancer))+
  geom_col(width = 0.7)+
  geom_text(aes(x=share_QTL+6,label=share_QTL))+theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20)) 

pdf("../puQTL/5.analysis/tissue_specific/share_QTL.pdf",width=5 ,height = 5)
p <- ggplot(all_result,aes(share_QTL,cancer))+
  geom_col(width = 0.7)+
  geom_text(aes(x=share_QTL+6,label=share_QTL))+theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20)) 
print(p)
dev.off()


pdf("../puQTL/5.analysis/tissue_specific/TCGA_specific.pdf",width=5 ,height = 5)
p <- ggplot(all_result,aes(TCGA_specific,cancer))+
  geom_col(width = 0.7)+
  geom_text(aes(x=TCGA_specific+6,label=TCGA_specific))+theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20)) 
print(p)
dev.off()



pdf("../puQTL/5.analysis/tissue_specific/GTEX_specific.pdf",width=5 ,height = 5)
p <- ggplot(all_result,aes(GTEX_specific,cancer))+
  geom_col(width = 0.7)+
  geom_text(aes(x=GTEX_specific+6,label=GTEX_specific))+theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        # panel.border = element_line(colour = "black"),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20)) 
print(p)
dev.off()








rm(list = ls())
GTEx_TCGA <- read.delim("../puQTL/5.analysis/normal/GTEx-TCGA.txt",header = F)
all_share <- c()
all_TCGA_specific_gene <- c()
all_GTEx_specific_gene <- c()
promoter_anno <- read.delim("../puQTL/GTEX/tmp/proActiv/data/promoterCoordinates.gencode.38.bed",header = F)
promoter_anno <- promoter_anno[,c(3,4,5,6)]
promoter_anno$V4 <- gsub(promoter_anno$V4,pattern = "prmtr\\.",replacement = "")
promoter_anno$V5 <- gsub(promoter_anno$V5,pattern = "\\..*",replacement = "")
hg19_promoter <- read.delim("../puQTL/GTEX/08.analysis/promoterCoordinates.gencode.v19.csv",header = T)
hg19_promoter <- hg19_promoter[,c(2,6)]
hg19_promoter$promoterId<- gsub(hg19_promoter$promoterId,pattern = "prmtr\\.",replacement = "")
colnames(hg19_promoter)[2] <- "promoter"
cancer_specific <- c()
share <- c()
normal_specific <- c()
all_result <- data.frame()

setwd("../puQTL/web/data/survival-table")
survival_dir <- list.files()
survival_dir <- gsub(survival_dir,pattern = "TCGA-",replacement = "")
survival_dir <- gsub(survival_dir,pattern = "\\.survival.Rdata",replacement = "")
GTEx_TCGA <- GTEx_TCGA[GTEx_TCGA$V2 %in% survival_dir,]
for (i in 1:nrow(GTEx_TCGA)) {
  GTEx <- read.delim(paste0("../puQTL/GTEX_new_sample/4.MatrixEQTL/1.significant/relative/",GTEx_TCGA[i,1]),header = F) 
  # GTEx_vcf <- read.delim(paste0("../puQTL/GTEX/06.puQTL/vcf/",GTEx_TCGA[i,1],".vcf"),header = F)
  # GTEx_vcf <- GTEx_vcf[,c(2,3)]
  # GTEx_vcf <- unique(GTEx_vcf)
  TCGA <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/TCGA-",GTEx_TCGA[i,2],".cis"),header = T)
  TCGA <- TCGA[,c(1,3,6,7,8,10)]
  colnames(promoter_anno)[2] <- "promoter"
  promoter_anno$promoter <- as.numeric(promoter_anno$promoter)
  TCGA <- dplyr::inner_join(TCGA,promoter_anno)
  hg19_promoter$promoter <- as.numeric(hg19_promoter$promoter)
  TCGA <- dplyr::inner_join(TCGA,hg19_promoter)
  TCGA$promoter <- paste0(TCGA$V3,"-",TCGA$promoter_gene)
  GTEx <- GTEx[,c(1,5,6)]
  GTEx$promoter <- paste0(str_split_fixed(GTEx$V5, "-",4)[,4],"-",str_split_fixed(GTEx$V5, "-",4)[,2])
  GTEx$QTL <- paste0(GTEx$V1,":",GTEx$promoter)
  TCGA$QTL <- paste0(TCGA$SNP_id,":",TCGA$promoter)
  share_QTL <- TCGA[TCGA$QTL%in%GTEx$QTL,]
  specific_QTL <- TCGA[!(TCGA$QTL%in%GTEx$QTL),]
  load(paste0("../puQTL/web/data/survival-table/TCGA-",GTEx_TCGA[i,2],".survival.Rdata"))
  specific_QTL_survival <- specific_QTL[specific_QTL$SNP_id%in%survival$SNPID,]
  print(nrow(specific_QTL_survival)/nrow(specific_QTL))
  share_QTL_survival <- share_QTL[share_QTL$SNP_id%in%survival$SNPID,]
  print(nrow(share_QTL_survival)/nrow(share_QTL))
}






rm(list = ls())
GTEx_TCGA <- read.delim("../puQTL/5.analysis/normal/GTEx-TCGA.txt",header = F)
all_share <- c()
all_TCGA_specific_gene <- c()
all_GTEx_specific_gene <- c()
promoter_anno <- read.delim("../puQTL/GTEX/tmp/proActiv/data/promoterCoordinates.gencode.38.bed",header = F)
promoter_anno <- promoter_anno[,c(3,4,5,6,1)]
promoter_anno$V4 <- gsub(promoter_anno$V4,pattern = "prmtr\\.",replacement = "")
promoter_anno$V5 <- gsub(promoter_anno$V5,pattern = "\\..*",replacement = "")
hg19_promoter <- read.delim("../puQTL/GTEX/08.analysis/promoterCoordinates.gencode.v19.csv",header = T)
hg19_promoter <- hg19_promoter[,c(2,6)]
hg19_promoter$promoterId<- gsub(hg19_promoter$promoterId,pattern = "prmtr\\.",replacement = "")
colnames(hg19_promoter)[2] <- "promoter"
cancer_specific <- c()
share <- c()
normal_specific <- c()
all_result <- data.frame()
for (i in 1:nrow(GTEx_TCGA)) {
  GTEx <- read.delim(paste0("../puQTL/GTEX_new_sample/4.MatrixEQTL/1.significant/relative/",GTEx_TCGA[i,1]),header = F) 
  TCGA <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/TCGA-",GTEx_TCGA[i,2],".cis"),header = T)
  TCGA <- TCGA[,c(1,3,6,7,8,10,4,5)]
  colnames(promoter_anno)[2] <- "promoter"
  promoter_anno$promoter <- as.numeric(promoter_anno$promoter)
  TCGA <- dplyr::inner_join(TCGA,promoter_anno)
  hg19_promoter$promoter <- as.numeric(hg19_promoter$promoter)
  TCGA <- dplyr::inner_join(TCGA,hg19_promoter)
  TCGA$promoter <- paste0(TCGA$V3,"-",TCGA$promoter_gene)
  GTEx <- GTEx[,c(1,5,6)]
  GTEx$promoter <- paste0(str_split_fixed(GTEx$V5, "-",4)[,4],"-",str_split_fixed(GTEx$V5, "-",4)[,2])
  GTEx$QTL <- paste0(GTEx$V1,":",GTEx$promoter)
  TCGA$QTL <- paste0(TCGA$SNP_id,":",TCGA$promoter)
  share_QTL <- TCGA[TCGA$QTL%in%GTEx$QTL,]
  specific_QTL <- TCGA[!(TCGA$QTL%in%GTEx$QTL),]
  specific_QTL_bed <- specific_QTL[,c(12,2,1,7,8)]
  specific_QTL_bed$Alt <- gsub(specific_QTL_bed$Alt,pattern = "\\,.*",replacement = "")
  specific_QTL_bed$V1 <- gsub(specific_QTL_bed$V1,pattern = "_.*",replacement = "")
  share_QTL_bed <- share_QTL[,c(12,2,1,7,8)]
  share_QTL_bed <- unique(share_QTL_bed)
  share_QTL_bed <- unique(share_QTL_bed)
  share_QTL_bed$Alt <- gsub(share_QTL_bed$Alt,pattern = "\\,.*",replacement = "")
  share_QTL_bed$V1 <- gsub(share_QTL_bed$V1,pattern = "_.*",replacement = "")
  write.table(specific_QTL_bed,paste0("../puQTL/5.analysis/normal/deepSEA/0.vcf/cancer/",GTEx_TCGA[i,2],".vcf"),col.names = F,row.names = F,sep = "\t",quote = F)
  write.table(share_QTL_bed,paste0("../puQTL/5.analysis/normal/deepSEA/0.vcf/share/",GTEx_TCGA[i,2],".vcf"),col.names = F,row.names = F,sep = "\t",quote = F)
}






#!/bin/bash
for i in `ls ../puQTL/5.analysis/normal/deepSEA/0.vcf/cancer`
do 
k=${i%%.*}
./apps/anaconda2/bin/python ../puQTL/DeepSEA/DeepSEA-v0.94c/rundeepsea.py ../puQTL/5.analysis/normal/deepSEA/0.vcf/cancer/${k}.vcf  ../puQTL/5.analysis/normal/deepSEA/1.rawresult/cancer/${k}
done    
