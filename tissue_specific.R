cat *cis | cut -f1,6,7 | grep -v promoter | sort | uniq -c | awk -F " " '{print $1}' | sort | uniq -c > ../puQTL/5.analysis/tissue_specific/summary
all_summary <- read.delim("../puQTL/5.analysis/tissue_specific/summary",header = F)
colnames(all_summary) <- c("number","overlap")
all_summary <- all_summary[order(all_summary$overlap),]
all_summary$overlap <- paste0(all_summary$overlap)
mouse <- all_summary

mouse$overlap  <- factor(mouse$overlap, levels=mouse$overlap)
dt = mouse[order(mouse$number, decreasing = TRUE),]
dt$fraction<-dt$number/sum(dt$number)
dt$ymax<-cumsum(dt$fraction)
dt$ymin<-c(0,head(dt$ymax,n=-1))
dt$labelPosition<-(dt$ymax + dt$ymin)/2
dt$label<- dt$V2
color1 <-  c("#BEBADA","#CCEBC5","#FB8072","#377EB8","#4DAF4A","#984EA3",
             "#FF7F00","#FFFF33","#A65628","#F781BF","#999999",
             "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854",
             "#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3",                                                                                                   
             "#80B1D3","#FDB462","#B3DE69","#FCCDE5",                                                                                                 
             "#D9D9D9","#BC80BD","#FFED6F","#1B9E77",                                                                                                 
             "#D95F02","#7570B3","#E7298A","#E41A1C")
pdf("../puQTL/5.analysis/tissue_specific/summary.pdf",width=12 ,height = 10)
p <- ggplot(dt,aes(ymax=round(ymax,4),ymin=round(ymin,4),
                   xmax=4,xmin=3))+
  geom_rect(aes(fill=overlap))+
  #geom_text(aes(x = 3.8, y = ((ymin+ymax)/2),label =label),size=5)+
  coord_polar(theta = "y")+
  scale_fill_manual(values=color1)+
  xlim(2,4)+
  theme_void()+theme(legend.title=element_blank())+ 
  theme(legend.text=element_text(size=15),legend.key.size =unit(30, "pt"))+
  theme(plot.title = element_text(hjust = 0.5,size = 40))
print(p)
dev.off()






rm(list = ls())
all_summary <- read.delim("../puQTL/5.analysis/tissue_specific/puQTL_summary.txt",header = F)
all_summary$tmp <- paste0(all_summary$V2,"_",all_summary$V3,"_",all_summary$V4)
all_summary <- all_summary[,c(1,5)]
colnames(all_summary) <- c("number","tmp")
setwd("../puQTL/QTL_result/0.rawdata/significant/cis/correct/")
cancer_list <- list.files(pattern = "cis")
all <- data.frame()
for (i in 1:33) {
  data <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/",cancer_list[i]),header = T)
  data <- data[,c(1,6,7,8,10)]
  data$tmp <- paste0(data$SNP_id,"_",data$promoter,"_",data$promoter_gene)
  data <- dplyr::inner_join(data,all_summary)
  data <- data[,c(4,5,7)]
  all <- rbind(all,data)
}


all$number <- factor(all$number,levels = c(1:33))
all$beta <- abs(all$beta)
pdf(paste0("../puQTL/5.analysis/tissue_specific/beta_boxplot.pdf"),width=15,height = 8)
p = ggplot(all, aes(x=number, y=beta)) + 
  geom_boxplot(fill=color1)+
  scale_fill_manual(values=color1)+labs(x = 'The number of tissues shared by the same puQTL',y = 'Absolute effect size of puQTL') +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20))
print(p)
dev.off()


pdf(paste0("../puQTL/5.analysis/tissue_specific/beta_boxplot.pdf"),width=15,height = 8)
p = ggplot(all, aes(x=number, y=beta)) + 
  geom_boxplot(fill=color1)+
  scale_fill_manual(values=color1)+labs(x = 'The number of tissues shared by the same puQTL',y = 'Absolute effect size of puQTL') +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20))
print(p)
dev.off()


pdf(paste0("../puQTL/5.analysis/tissue_specific/P_boxplot.pdf"),width=15,height = 8)
p = ggplot(all, aes(x=number, y=-log10(p.value))) + 
  geom_boxplot(fill=color1)+
  scale_fill_manual(values=color1)+labs(x = 'The number of tissues shared by the same puQTL',y = 'P-value of puQTL') +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank())+ theme(plot.title = element_text(hjust = 0.5,size = 20))
print(p)
dev.off()



rm(list = ls())
all_summary <- read.delim("../puQTL/5.analysis/tissue_specific/puQTL_summary.txt",header = F)
all_summary$tmp <- paste0(all_summary$V2,"_",all_summary$V3,"_",all_summary$V4)
all_summary <- all_summary[,c(1,5)]
colnames(all_summary) <- c("number","tmp")
setwd("../puQTL/QTL_result/0.rawdata/significant/cis/correct/")
cancer_list <- list.files(pattern = "cis")
all <- data.frame()
for (i in 1:33) {
  data <- read.delim(paste0("../puQTL/QTL_result/0.rawdata/significant/cis/correct/",cancer_list[i]),header = T)
  data <- data[,c(1,6,7,8,10)]
  data$tmp <- paste0(data$SNP_id,"_",data$promoter,"_",data$promoter_gene)
  data <- dplyr::inner_join(data,all_summary)
  data <- data[,c(4,5,7)]
  all <- rbind(all,data)
}

all$type <- ifelse(all$number>11,"share","uniq")
all$type <- ifelse(all$number<11&all$number>1,"Inter",all$type)
all$type <- factor(all$type,levels = c("share","Inter","uniq"))


all$beta[abs(all$beta) > 2.5] <- 2.5
p = ggplot(all, aes(x=type, y=abs(beta),color = type)) + 
  geom_boxplot()+labs(x = '',y = 'Absolute effect size of puQTL') +
  scale_color_manual(values=c("#4778a0", "#d76a6d", "#8bc086"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme_classic()
pdf(paste0("../puQTL/5.analysis/tissue_specific/figure_2B.pdf"),width=10,height = 10)

print(p)
dev.off()
all$tmp <- -log10(all$p.value)
all$tmp[all$tmp > 30] <- 30


p <- ggplot(all, aes(x = tmp, fill = type)) +
  scale_fill_manual(values=c("#4778a0", "#d76a6d", "#8bc086"))+
  geom_density(alpha = 0.5)+theme_bw()+theme_classic()

pdf(paste0("../puQTL/5.analysis/tissue_specific/figure_2C.pdf"),width=10,height = 10)

print(p)
dev.off()


rm(list = ls())
source("../puQTL/GTEX/08.analysis/MASH/test/gtexresults-master/code/normfuncs.R")
thresh <- 0.05
out      <- readRDS("../puQTL/MASH/test/gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxbeta <- out$strong.b
maxz    <- out$strong.z
standard.error <- out$strong.s
out      <-readRDS("../puQTL/MASH/test/gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
pm.mash.beta   <- pm.mash*standard.error
lfsr           <- out$lfsr
lfsr[lfsr < 0] <- 0

sigmat <- (lfsr <= thresh)
nsig   <- rowSums(sigmat)
pdf("../puQTL/5.analysis/tissue_specific/mag_his.pdf",width=8,height = 5)
hist((het.func(het.norm(effectsize=pm.mash.beta[nsig>0,]),threshold=0.5)),
     main="",xlab="",breaks=0.5:33.5,col="grey",freq=FALSE,ylim=c(0,0.5),
     xaxt="n")
axis(1,at = seq(1, 33, by=1),labels = c(1:33))
mtext("All Tissues")
dev.off()



sign.func <- function (normeffectsize)
  apply(normeffectsize,1,function(x)(sum(x>0)))
sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
pdf("../puQTL/5.analysis/tissue_specific/sign_his.pdf",width=8,height = 5)
hist(sign.func(het.norm(effectsize=pm.mash.beta[nsig>0,])),main="",xlab="",
     breaks=0.5:33.5,col="grey",freq=FALSE,xaxt="n",ylim=c(0,0.4))
axis(1, at=seq(1, 33, by=1), labels=c(1:33))
mtext("Number of tissues shared by sign")
dev.off()





rm(list = ls())
source("../puQTL/GTEX/08.analysis/MASH/test/gtexresults-master/code/normfuncs.R")
thresh <- 0.05
out      <- readRDS("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/1.MASH/2.analysis/all_tissue/gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxbeta <- out$strong.b
maxz    <- out$strong.z
standard.error <- out$strong.s
out      <-readRDS("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/1.MASH/2.analysis/test/gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
pm.mash.beta   <- pm.mash*standard.error
aa <- colnames(pm.mash.beta)[!colnames(pm.mash.beta)%in%c("Brain-Amygdala","Brain-Spinal-cord--cervical-c-1","Brain-Substantia-nigra","Kidney-Cortex","Minor-Salivary-Gland")]
pm.mash.beta <- pm.mash.beta[,aa]
lfsr           <- out$lfsr
lfsr <- lfsr[,aa]
lfsr[lfsr < 0] <- 0

sigmat <- (lfsr <= thresh)
nsig   <- rowSums(sigmat)
pdf("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/mag_his.pdf",width=8,height = 5)
hist((het.func(het.norm(effectsize=pm.mash.beta[nsig>0,]),threshold=0.5)),
     main="",xlab="",breaks=0.5:44.5,col="grey",freq=FALSE,ylim=c(0,0.8),
     xaxt="n")
axis(1,at = seq(1, 44, by=1),labels = c(1:44))
mtext("All Tissues")
dev.off()



sign.func <- function (normeffectsize)
  apply(normeffectsize,1,function(x)(sum(x>0)))
sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
pdf("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/sign_his.pdf",width=8,height = 5)
hist(sign.func(het.norm(effectsize=pm.mash.beta[nsig>0,])),main="",xlab="",
     breaks=0.5:44.5,col="grey",freq=FALSE,xaxt="n",ylim=c(0,0.2))
axis(1, at=seq(1, 44, by=1), labels=c(1:44))
mtext("Number of tissues shared by sign")
dev.off()
















rm(list = ls())
out      <- readRDS("../puQTL/MASH/test/gtexresults-master/fastqtl_to_mash_output/tissue.mash.rds")
maxb     <- out$strong.b
maxz     <- out$strong.z
out      <-readRDS("../puQTL/MASH/test/gtexresults-master/mashr_flashr_workflow_output/tissue.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
lfsr.all       <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta   <- pm.mash*standard.error
thresh       <- 0.05
pm.mash.beta <- pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash    <- lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size <- matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size) <- rownames(shared.fold.size) <- colnames(maxz)
for (i in 1:ncol(lfsr.mash))
  for (j in 1:ncol(lfsr.mash)) {
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j] = mean(quotient > 0.5 & quotient < 2)
  }

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(64)

lat=shared.fold.size[cancer_order,cancer_order]
p <- corrplot(lat,
              type = "upper",
              is.corr = FALSE,
              order = "hclust",
              col.lim = c(0,1), 
              col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                           "#E0F3F8","#91BFDB","#4575B4")))(64),
              method = "shade",
              tl.col = "black",
              diag = TRUE)


Squamous_cancer.order <- c("TCGA-HNSC.cis","TCGA-LUSC.cis","TCGA-UCEC.cis","TCGA-BLCA.cis")


Squamous_lat=shared.fold.size[Squamous_cancer.order,Squamous_cancer.order]
corrplot(Squamous_lat,
         type = "upper",
         is.corr = FALSE,
         order = "hclust",
         col.lim = c(0,1), 
         col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                      "#E0F3F8","#91BFDB","#4575B4")))(64),
         method = "shade",
         tl.col = "black",
         diag = TRUE)

Neuronal_cancer.order <- c("TCGA-GBM.cis","TCGA-LGG.cis","TCGA-PCPG.cis")
Neuronal_cancer.order <- c("TCGA-UVM.cis","TCGA-SKCM.cis")
Neuronal_cancer.order <- c("TCGA-KICH.cis","TCGA-KIRC.cis","TCGA-KIRP.cis")
Neuronal_cancer.order <- c("TCGA-LUAD.cis","TCGA-LUSC.cis")
Neuronal_lat=shared.fold.size[Neuronal_cancer.order,Neuronal_cancer.order]
corrplot(Neuronal_lat,
         type = "upper",
         is.corr = FALSE,
         order = "hclust",
         col.lim = c(0,1), 
         col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                      "#E0F3F8","#91BFDB","#4575B4")))(64),
         method = "shade",
         tl.col = "black",
         diag = TRUE)


#sign
rm(list = ls())
out      <- readRDS("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/1.MASH/2.analysis/relative/brain/gtexresults-master/fastqtl_to_mash_output/brain.mash.rds")
maxb     <- out$strong.b
maxz     <- out$strong.z
out      <-readRDS("../puQTL/GTEX_new_sample/5.analysis/tissue_specific/1.MASH/2.analysis/relative/brain/gtexresults-master/mashr_flashr_workflow_output/brain.mash.EZ.posterior.rds")
pm.mash        <- out$PosteriorMean
lfsr.all <- out$lfsr
standard.error <- maxb/maxz
pm.mash.beta <- pm.mash*standard.error

thresh=0.05
pm.mash.beta=pm.mash.beta[rowSums(lfsr.all<0.05)>0,]
lfsr.mash=lfsr.all[rowSums(lfsr.all<0.05)>0,]
shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(maxz)
for(i in 1:ncol(lfsr.mash)){
  for(j in 1:ncol(lfsr.mash)){
    sig.row=which(lfsr.mash[,i]<thresh)
    sig.col=which(lfsr.mash[,j]<thresh)
    a=(union(sig.row,sig.col))
    quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
    shared.fold.size[i,j]=mean(quotient > 0)
  }
}




clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(64)

lat=shared.fold.size
p <- corrplot(lat,
              type = "upper",
              title = "aa",
              is.corr = FALSE,
              order = "hclust",
              col.lim = c(0.8,1), 
              col = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                           "#E0F3F8","#91BFDB","#4575B4")))(64),
              method = "shade",
              tl.col = "black",
              diag = TRUE)





































../puQTL/QTL_result/0.rawdata/significant/cis/correct$ grep ENSG00000142208 *cis > ../puQTL/QTL_result/0.rawdata/ENSG00000142208_tissue
rm(list = ls())
gene_puQTL <- read.delim("../puQTL/QTL_result/0.rawdata/ENSG00000142208_tissue",header = F)
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "\\.cis.*",replacement = "")
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "TCGA-",replacement = "")

gene_puQTL <- arrange(gene_puQTL,desc(-log(gene_puQTL$V5)))
gene_puQTL$V1  <- factor(gene_puQTL$V1, levels=gene_puQTL$V1)
p <- ggplot(gene_puQTL, aes(x=V1, y=-log(V5))) +
  geom_segment( aes(x=V1, xend=V1, y=0, yend=-log(V5))) +
  geom_point( size=2.5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5)+
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12))+
  theme(legend.position='none')+  ##去除legend+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
pdf("../puQTL/5.analysis/tissue_specific/ENSG00000142208_P.pdf",width=7,height = 4)
print(p)
dev.off()



../puQTL/QTL_result/0.rawdata/significant/cis/correct$ grep ENSG00000108654 *cis > ../puQTL/QTL_result/0.rawdata/ENSG00000108654
rm(list = ls())
gene_puQTL <- read.delim("../puQTL/QTL_result/0.rawdata/ENSG00000108654",header = F)
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "\\.cis.*",replacement = "")
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "TCGA-",replacement = "")

gene_puQTL <- arrange(gene_puQTL,desc(-log(gene_puQTL$V5)))
gene_puQTL$V1  <- factor(gene_puQTL$V1, levels=gene_puQTL$V1)
p <- ggplot(gene_puQTL, aes(x=V1, y=-log(V5))) +
  geom_segment( aes(x=V1, xend=V1, y=0, yend=-log(V5))) +
  geom_point( size=2.5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5)+
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12))+
  theme(legend.position='none')+  ##去除legend+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
pdf("../puQTL/5.analysis/tissue_specific/ENSG00000108654_P.pdf",width=4,height = 4)
print(p)
dev.off()


../puQTL/QTL_result/0.rawdata/significant/cis/correct$ grep ENSG00000119888 *cis > ../puQTL/QTL_result/0.rawdata/ENSG00000119888
rm(list = ls())
gene_puQTL <- read.delim("../puQTL/QTL_result/0.rawdata/ENSG00000141646.txt",header = F)
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "\\.cis.*",replacement = "")
gene_puQTL$V1 <- gsub(gene_puQTL$V1,pattern = "TCGA-",replacement = "")

gene_puQTL <- arrange(gene_puQTL,desc(-log(gene_puQTL$V5)))
gene_puQTL$V1  <- factor(gene_puQTL$V1, levels=gene_puQTL$V1)
p <- ggplot(gene_puQTL, aes(x=V1, y=-log(V5))) +
  geom_segment( aes(x=V1, xend=V1, y=0, yend=-log(V5))) +
  geom_point( size=2.5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1.5)+
  theme(axis.title.y = element_text(size = 16,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12))+
  theme(legend.position='none')+  ##去除legend+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
pdf("../puQTL/5.analysis/tissue_specific/ENSG00000141646_P.pdf",width=5,height = 4)
print(p)
dev.off()









snp <- read.delim("../puQTL/3.figure/boxplot/pairs/TCGA-COAD.snpmatrix")
expression <-  read.delim("../puQTL/3.figure/boxplot/pairs/TCGA-COAD.exp") 
rs <- "rs920784:48318120:C:T"  
erna <- "prmtr.34503.ENSG00000141646.9"
# main 
genotype = snp[snp$snpid == rs,]
genotype =genotype[,-1]
exp = expression[expression$gene == erna,]
exp = exp[,-1]
exp = exp[,colnames(genotype)]
genotype[which(genotype==0)] = "AA"
genotype[which(genotype==1)] = "Aa"
genotype[which(genotype==2)] = "aa"
genotype = as.vector(as.matrix(genotype))

dat = data.frame(genotype=genotype,expression=as.numeric(exp))
dat = dat[which(!is.na(genotype)),]
dat$genotype = factor(dat$genotype,levels=c("AA","Aa","aa"))

library(EnvStats)
p <- ggplot(dat,aes(x=genotype,y=expression,fill=genotype)) +
  geom_boxplot() +
  stat_n_text() +
  labs(x="Genotype", y="Normalized promoter usage") +
  theme_bw() +theme(legend.position="top")+theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(colour="black", size = rel(1.6)),
        axis.text.y = element_text(colour="black", size = rel(1.6)),
        axis.title.x = element_text(colour = "black", size = rel(1.2)),
        axis.title.y = element_text(colour = "black", size = rel(1.2)),
        strip.text.x = element_text(colour = "black", size = rel(1.6)),
        legend.text = element_text(colour = "black",size = rel(1.0)),
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=1, linetype=1),
        plot.title = element_text(hjust=0.5))
pdf(file= paste0("../puQTL/5.analysis/tissue_specific/boxplot_ENSG00000141646.pdf"),width = 6,height = 5)
print(p)
dev.off()


