rm(list = ls())
setwd("../puQTL/TF/result")
cancer_list <- list.files(pattern = "cis$")
all_enrich <- read.delim("../puQTL/5.analysis/enrichment/enrichment.txt")
all_enrich <- all_enrich[,c(1,2)]
for (i in 1:length(cancer_list)) {
  TF_result <- read.delim(paste0("../puQTL/TF/result/",cancer_list[i]),header = T)
  TF_result <- TF_result[,c(6,11,13)]
  TF_all <- unique(TF_result$geneSymbol)
  all_result_TF <- read.delim(paste0("../puQTL/5.analysis/TF/0.rawdata/",gsub(cancer_list[i],pattern = "cis",replacement = "bed"))
                              ,header = F)  
  all_result_TF <- all_result_TF[!(all_result_TF$V1 %in%TF_result$SNP_id), ]
  tissue <- gsub(cancer_list[i],pattern = "\\.cis",replacement = "")
  
  all_puQTL_TF <- data.frame()
  for (k in 1:length(TF_all)) {
    a <- TF_all[k]
    tmp_puQTL <- TF_result[TF_result$geneSymbol == a ,]
    no_puQTL <- all_result_TF[all_result_TF$V2 == a ,]
    tmp_puQTL_nobreak <- all_enrich[i,2] - nrow(tmp_puQTL)
    no_puQTL_nobreak <- all_enrich[i,1] - nrow(no_puQTL)
    part_result <- fisher.test(matrix(c(nrow(tmp_puQTL),nrow(no_puQTL),tmp_puQTL_nobreak,no_puQTL_nobreak),nrow=2))
    if (part_result$p.value < 0.05) {
      puQTL_TF_part <- data.frame(TF=a,OR =part_result$estimate,p = part_result$p.value,low=part_result$conf.int[1],upper=part_result$conf.int[2])
      write.table(puQTL_TF_part,paste0("../puQTL/5.analysis/TF/1.rawresult/fin/",cancer_list[i]),append = T,col.names = T,row.names = F,sep = "\t",quote = F)  
    }
    puQTL_TF_part <- data.frame(TF=a,OR =part_result$estimate,p = part_result$p.value,low=part_result$conf.int[1],upper=part_result$conf.int[2])
    all_puQTL_TF <- rbind(all_puQTL_TF,puQTL_TF_part)
  }
  all_puQTL_TF <- all_puQTL_TF[order(all_puQTL_TF$p),]  
  write.table(all_puQTL_TF,paste0("../puQTL/5.analysis/TF/1.rawresult/",cancer_list[i]),col.names = T,row.names = F,sep = "\t",quote = F)
}






all_TF <- read.delim("../puQTL/5.analysis/TF/1.rawresult/all_TF.txt",header = F)
colnames(all_TF) <- "TF"
setwd("../puQTL/5.analysis/TF/1.rawresult")
cancer_list <- list.files(pattern = "cis$")

for (i in 1:length(cancer_list)) {
  cancer_data <- read.delim(paste0("../puQTL/5.analysis/TF/1.rawresult/",cancer_list[i]),header = T)
  cancer_data <- cancer_data[,c(1,3)]
  cancer_data$p <- -log10(cancer_data$p)
  all_TF <- dplyr::inner_join(all_TF,cancer_data)
  colnames(all_TF)[i+1] <- gsub(gsub(cancer_list[i],pattern = "\\.cis",replacement = ""),pattern = "TCGA-",replacement = "")
}



rownames(all_TF)<- all_TF$TF
all_TF <- all_TF[,-1]
path <- paste0("../puQTL/5.analysis/enrichment/plot/TF.pdf")
pdf(path,width=6 ,height = 6.5)
pheatmap(as.matrix(all_TF),cluster_cols = F,
         cluster_rows = T,
         color = c(colorRampPalette(colors = c("#1781b5","#f6dcce"))(98/5),colorRampPalette(colors = c("#f6dcce","#d11a2d"))(98/0.2))
)
dev.off()





bk <- c(seq(0,5,by=1),seq(6,96,by=1)) 

plot <- pheatmap::pheatmap( 
  all_TF, 
  scale = 'row', 
  color = c(colorRampPalette(colors = c("blue","white"))(98/5),colorRampPalette(colors = c("white","red"))(98/92)), 
  legend_breaks=seq(0,98,1),breaks=bk, 
  cluster_cols = F,cluster_rows = F, 
  show_colnames = F,fontsize = 14, 
  #gaps_col = c(123,186,323,493), ###res=5 
  gaps_col = c(105,171,327,492), 
  border_color = 'white' 
)




all_TF <- read.delim("../puQTL/5.analysis/TF/1.rawresult/all_TF.txt",header = F)
colnames(all_TF) <- "TF"
setwd("../puQTL/5.analysis/TF/1.rawresult")
cancer_list <- list.files(pattern = "cis$")

for (i in 1:length(cancer_list)) {
  cancer_data <- read.delim(paste0("../puQTL/5.analysis/TF/1.rawresult/",cancer_list[i]),header = T)
  cancer_data <- cancer_data[,c(1,3)]
  cancer_data$p <- -log10(cancer_data$p)
  all_TF <- dplyr::inner_join(all_TF,cancer_data)
  colnames(all_TF)[i+1] <- gsub(gsub(cancer_list[i],pattern = "\\.cis",replacement = ""),pattern = "TCGA-",replacement = "")
}




plot_TF <- all_TF$TF

all_TF <- read.delim("../puQTL/5.analysis/TF/1.rawresult/all_TF.txt",header = F)
colnames(all_TF) <- "TF"
setwd("../puQTL/5.analysis/TF/1.rawresult")
cancer_list <- list.files(pattern = "cis$")

for (i in 1:length(cancer_list)) {
  cancer_data <- read.delim(paste0("../puQTL/5.analysis/TF/1.rawresult/",cancer_list[i]),header = T)
  cancer_data <- cancer_data[,c(1,3)]
  cancer_data$p <- -log10(cancer_data$p)
  all_TF <- dplyr::left_join(all_TF,cancer_data)
  colnames(all_TF)[i+1] <- gsub(gsub(cancer_list[i],pattern = "\\.cis",replacement = ""),pattern = "TCGA-",replacement = "")
}



all_TF_tmp <- all_TF[!(all_TF$TF%in%plot_TF),]
rownames(all_TF_tmp)<- all_TF_tmp$TF
all_TF_tmp <- all_TF_tmp[,-1]

all_TF_tmp[is.na(all_TF_tmp)] <- 0
expSamples = apply(all_TF_tmp, 1, function(x) { length(x[x < -log10(0.05/362)]) } )
all_TF_tmp <- all_TF_tmp[expSamples<33,]
# path <- paste0("../puQTL/5.analysis/enrichment/plot/TF.pdf")
# pdf(path,width=6 ,height = 6.5)
p <- pheatmap(as.matrix(all_TF_tmp),cluster_cols = F,
              cluster_rows = T,
              color = c(colorRampPalette(colors = c("#1781b5","#f6dcce"))(98/5),colorRampPalette(colors = c("#f6dcce","#d11a2d"))(98/0.2))
)
# dev.off()
all_TF_tmp <-  all_TF_tmp[p$tree_row$order,]

TF_plot_order <- rownames(all_TF_tmp)
plot_TF <- all_TF$TF
all_TF <- read.delim("../puQTL/5.analysis/TF/1.rawresult/all_TF.txt",header = F)
colnames(all_TF) <- "TF"
setwd("../puQTL/5.analysis/TF/1.rawresult")
cancer_list <- list.files(pattern = "cis$")

for (i in 1:length(cancer_list)) {
  cancer_data <- read.delim(paste0("../puQTL/5.analysis/TF/1.rawresult/",cancer_list[i]),header = T)
  cancer_data <- cancer_data[,c(1,3)]
  cancer_data$p <- -log10(cancer_data$p)
  all_TF <- dplyr::left_join(all_TF,cancer_data)
  colnames(all_TF)[i+1] <- gsub(gsub(cancer_list[i],pattern = "\\.cis",replacement = ""),pattern = "TCGA-",replacement = "")
}



all_TF_tmp <- all_TF
rownames(all_TF_tmp)<- all_TF_tmp$TF
all_TF_tmp <- all_TF_tmp[,-1]
all_TF_tmp <- all_TF_tmp[TF_plot_order,]
p <- pheatmap(as.matrix(all_TF_tmp),cluster_cols = F,
              cluster_rows = F,
              color = c(colorRampPalette(colors = c("#1781b5","#f6dcce"))(98/5),colorRampPalette(colors = c("#f6dcce","#d11a2d"))(98/0.2))
)


path <- paste0("../puQTL/5.analysis/enrichment/plot/all_rest_TF.pdf")
pdf(path,width=10 ,height = 20)
p <- pheatmap(as.matrix(all_TF_tmp),cluster_cols = F,
              cluster_rows = F,
              color = c(colorRampPalette(colors = c("#1781b5","#f6dcce"))(98/5),colorRampPalette(colors = c("#f6dcce","#d11a2d"))(98/0.2))
)
dev.off()
