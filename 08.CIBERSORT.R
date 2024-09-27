rm(list = ls())
CIBERSORT <- read.delim("../puQTL/5.analysis/immune/CIBERSORT/format.txt",header = T)
setwd("../puQTL/QTL_result/0.rawdata/significant/cis/correct")
cancer_list <- list.files(pattern = "cis")
for (i in 1:length(cancer_list)) {
  puQTL <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/",cancer_list[i]),header = T)
  cancer <- gsub(cancer_list[i],pattern = "TCGA-",replacement = "")
  cancer <- gsub(cancer,pattern = "\\.cis",replacement = "")
  promoter_usage <- read.delim(paste0("../puQTL/1.combine_data/expression/TCGA_",cancer,".exp"))
  promoter_usage$gene <- gsub(promoter_usage$gene,pattern = "\\.ENSG.*",replacement = "")
  promoter_usage$gene <- gsub(promoter_usage$gene,pattern = "prmtr\\.",replacement = "")
  promoter_usage <- promoter_usage[promoter_usage$gene %in% puQTL$promoter,]
  colnames(promoter_usage) <-  gsub(colnames(promoter_usage),pattern = "\\.",replacement = "-")
  test_data <- t(promoter_usage)
  all_data_result <- data.frame()
  for (k in 1:ncol(test_data)) {
    part_data <- as.data.frame(test_data[,c(k)])
    part_data$SampleID <- rownames(part_data)
    colnames(part_data) <- c("promoter_usage","SampleID")
    colnames(CIBERSORT)[1] <- "SampleID"
    part_data <- dplyr::inner_join(part_data,CIBERSORT)
    for (f in 4:ncol(part_data)) {
      analysis <- part_data[,c(1,f)]
      analysis <- analysis[!is.na(analysis$promoter_usage),]
      if (sd(as.numeric(analysis[,2])) != 0) {
        part_result <- cor.test(as.numeric(analysis$promoter_usage),as.numeric(analysis[,2]))
        fin_part_result <- data.frame(pomoter=colnames(test_data)[k],
                                      celltype = colnames(part_data)[f],
                                      p_value = part_result$p.value,
                                      cor=part_result$estimate)
        all_data_result <- rbind( all_data_result,fin_part_result)
      }
      
      
    }
  }
  write.table(all_data_result,paste0("../puQTL/5.analysis/immune/CIBERSORT/0.result/",cancer,".txt"),col.names = T,row.names = F,sep = "\t",quote = F)
}




rm(list = ls())
setwd("../puQTL/5.analysis/immune/CIBERSORT/0.result")
cancer_list <- list.files(pattern = "txt")
all_result <- data.frame()
for (i in 1:length(cancer_list)) {
  all_data <- read.delim(paste0("../puQTL/5.analysis/immune/CIBERSORT/0.result/",cancer_list[i]))
  all_data <- all_data[all_data$p_value <0.05,]
  cancer <- gsub(cancer_list[i],pattern = "\\.txt",replacement = "")
  puQTL <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/TCGA-",cancer,".cis"),header = T)
  all_list <- split(all_data,f=all_data$celltype)
  for (k in 1:length(all_list)) {
    cell_type <- gsub(names(all_list[k]),pattern = "\\.",replacement = "-")
    cell_type_number <- nrow(all_list[[k]])
    cell_type_fraction <- length(unique(all_list[[k]]$pomoter))/length(unique(puQTL$promoter))
    part_result <- data.frame(cancer=cancer,cell_type = cell_type,number =cell_type_number,fraction=cell_type_fraction)
    all_result <- rbind(all_result,part_result)
  }
}










rm(list = ls())
setwd("../puQTL/5.analysis/immune/CIBERSORT/0.result")
cancer_list <- list.files(pattern = "txt")
all_result <- data.frame()
for (i in 1:length(cancer_list)) {
  all_data <- read.delim(paste0("../puQTL/5.analysis/immune/CIBERSORT/0.result/",cancer_list[i]))
  all_data <- all_data[all_data$p_value <0.05,]
  cancer <- gsub(cancer_list[i],pattern = "\\.txt",replacement = "")
  puQTL <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/TCGA-",cancer,".cis"),header = T)
  cell_type_number <- length(unique(all_data$pomoter))
  cell_type_fraction <- length(unique(all_data$pomoter))/length(unique(puQTL$promoter))
  part_result <- data.frame(cancer=cancer,number =cell_type_number,fraction=cell_type_fraction)
  all_result <- rbind(all_result,part_result)
}

path <- paste0("../puQTL/5.analysis/enrichment/plot/drug.pdf")
pdf(path,width=10 ,height = 6.5)
print(p)
dev.off()

p <- ggplot(all_result)+
  geom_col(aes(x= cancer,
               y= number),
           color='black',width=0.6,fill='#d5a478')+
  labs(x=NULL,y='Number of immune-related puQTL-gene',
       title = 'a')+
  theme_test(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(color = 'black',face = 'bold'),
        plot.margin = margin(1,0.5,0.5,2.5,'cm'),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold',
                                  size=13,hjust = 0.5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1100),
                     sec.axis = sec_axis(~./1100,
                                         name = 'Proportions of immune-related puQTL-gene',
                                         breaks = seq(0,1,0.1)))+
  geom_line(aes(x= cancer,
                y=fraction*1100,
                group=1),
            linetype=3,cex=1)+
  geom_point(aes(x= cancer,
                 y=fraction*1100),
             color='#589c47',size=3.5)+
  geom_text(aes(x= cancer,
                y=number,
                label=number),
            vjust=-0.5,size=3.5,fontface='bold')+
  geom_rect(aes(xmin=8,xmax=11,ymin=1030,ymax=1060),
            fill='#d5a478',color='black')+
  annotate('text',x=13,y=1045,label='Gene Count',
           fontface='bold',size=4.5)+
  annotate('segment',x=16,xend = 20,y=1045,yend = 1045,
           linetype=3,cex=1)+
  annotate(geom='text',x=18,y=1045,label='•',
           size=12,color='#589c47')+
  annotate('text',x=21,y=1050,label='Proportions',
           fontface='bold',size=4.5)+
  ggplot(data)+
  geom_col(aes(x= reorder(GO_terms,-`-Log10(P value)`),
               y=Gene_count),
           color='black',width=0.6,fill='#d5a478')+
  labs(x=NULL,y='Gene count',
       title = 'Enriched GO Terms in Biological Process (expanded gene families)')+
  theme_test(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(color = 'black',face = 'bold'),
        plot.margin = margin(1,0.5,0.5,2.5,'cm'),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold',
                                  size=13,hjust = 0.5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,6000),
                     sec.axis = sec_axis(~./42,
                                         name = '-Log10(P value)',
                                         breaks = seq(0,140,20)))+
  geom_line(aes(x= reorder(GO_terms,-`-Log10(P value)`),
                y=`-Log10(P value)`*42,
                group=1),
            linetype=3,cex=1)+
  geom_point(aes(x= reorder(GO_terms,-`-Log10(P value)`),
                 y=`-Log10(P value)`*42),
             color='#589c47',size=3.5)+
  geom_text(aes(x= reorder(GO_terms,-`-Log10(P value)`),
                y=Gene_count,
                label=Gene_count),
            vjust=-0.5,size=3.5,fontface='bold')+
  geom_rect(aes(xmin=4,xmax=5.4,ymin=5550,ymax=5790),
            fill='#d5a478',color='black')+
  annotate('text',x=7.5,y=5680,label='Gene Count',
           fontface='bold',size=4.5)+
  annotate('segment',x=10.9,xend = 12.6,y=5670,yend = 5670,
           linetype=3,cex=1)+
  annotate(geom='text',x=11.7,y=5700,label='•',
           size=12,color='#589c47')+
  annotate('text',x=15.4,y=5680,label='-Log10(P value)',
           fontface='bold',size=4.5)+
  geom_line(data = df,aes(a,b),cex=0.5)


