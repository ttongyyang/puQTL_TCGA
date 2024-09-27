rm(list = ls())
setwd("../puQTL/5.analysis/coloc/0.rawdata/")
gwas_list <- list.files(pattern = "ma")
setwd("../puQTL/5.analysis/coloc/1.puQTL")
cancer_list <- list.files(pattern = "")
cancer_sample <- read.delim("../puQTL/coloc/sample.txt",header = F)
gwas_anno <- read.delim("../puQTL/5.analysis/Mendelian_randomization/filter_gene")
gwas_anno <- gwas_anno[grepl("finn",gwas_anno$ID),]
gwas_anno$cancer <- gsub(gwas_anno$cancer,pattern = " ",replacement = "")
gwas_anno$ID <- gsub(gwas_anno$ID,pattern = " ",replacement = "")
gwas_anno$ID <- str_split_fixed(gwas_anno$ID, "-",3)[,3]

gwas_sample <- read.delim("../puQTL/5.analysis/coloc/sample.txt",header = F)

for (k in 1:length(cancer_list)) {
  QTL_data <- read.delim(paste0("../puQTL/5.analysis/coloc/1.puQTL/",cancer_list[k]),header = F)
  colnames(QTL_data) <- c("rs_id","promoter_ID","gene","beta","pval","MAF")
  QTL_data$se=sqrt(((QTL_data$beta)^2)/qchisq(QTL_data$pval,1,lower.tail=F))
  QTL_data$varbeta <- (QTL_data$se)^2
  gwas_cancer <- gwas_anno[gwas_anno$cancer == gsub(cancer_list[k],pattern = "TCGA-",replacement = ""),]
  if (nrow(gwas_cancer)!=0) {
    for (i in 1:nrow(gwas_cancer)) {
      if (length(gwas_list[grepl(gwas_cancer[i,2],gwas_list)]) == 1) {
        GWAS_data <- read.delim(paste0("../puQTL/5.analysis/Mendelian_randomization/SMR/gwas/finngen_analysis/",gwas_list[grepl(gwas_cancer[i,2],gwas_list)]),header = F)
        if (nrow(GWAS_data) >0 ) {
          GWAS_data <- GWAS_data[,c(1,5,6,7,8)]
          colnames(GWAS_data) <- c("rs_id","GWAS_beta","GWAS_se","GWAS_P","N") 
          fin_data <- dplyr::inner_join(GWAS_data,QTL_data)
          if (nrow(fin_data) > 0) {
            fin_data$varbeta <- (as.numeric(fin_data$GWAS_se))^2
            cancer_size <- cancer_sample[gsub(cancer_sample$V1,pattern = " ",replacement = "") ==gsub(cancer_list[k],pattern = "TCGA-",replacement = ""),]$V2
            coloc_list <- split(fin_data,f=fin_data$gene)
            for (m in 1:length(coloc_list)) {
              input <- coloc_list[[m]]
              input <- input[!duplicated(input$rs_id),]
              gwas_number <- gwas_sample[gwas_sample$V1 == gsub(gwas_list[grepl(gwas_cancer[i,2],gwas_list)],pattern = "\\.ma",replacement = ""),]$V2
              input$pval <- as.numeric(input$pval)
              input$GWAS_P <- as.numeric(input$GWAS_P)
              input$GWAS_beta <- as.numeric(input$GWAS_beta)
              input$varbeta <- as.numeric(input$varbeta)
              result <- coloc.abf(dataset1=list(snp = input$rs_id,pvalues=input$pval, type="quant", N=cancer_size),
                                  dataset2=list(snp = input$rs_id,pvalues=input$GWAS_P,beta=input$GWAS_beta,varbeta=input$varbeta, type="cc", N=gwas_number),
                                  MAF=input$MAF)
              summary <- data.frame(result$summary)
              summary <- data.frame(t(summary))
              summary$traits <- gwas_list[grepl(gwas_cancer[i,2],gwas_list)]
              summary$cancer <- cancer_list[k]
              summary$promoter <- paste0(input$promoter_ID[1],":",input$gene[1])
              write.table(summary,paste0("../puQTL/5.analysis/coloc/2.result/0.rawdata/new/sum_",cancer_list[k]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
              result <- result$results
              part_result <- result
              if (nrow(part_result) > 0) {
                part_result <- part_result[,c(1,2,14)]
                part_result$traits <- gwas_list[grepl(gwas_cancer[i,2],gwas_list)]
                part_result$cancer <- cancer_list[k]
                part_result$first_p4 <- summary$PP.H4.abf
                part_result$promoter <- paste0(input$promoter_ID[1],":",input$gene[1])
                write.table(part_result,paste0("../puQTL/5.analysis/coloc/2.result/0.rawdata/new/add_",cancer_list[k]),sep = "\t",quote = F,col.names = F,row.names = F,append = T)
              }
            }
          }
        }
      }
      
    }
  }
  
}