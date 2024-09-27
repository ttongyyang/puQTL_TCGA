#choose promoter
mkdir expression
for i in `cat cancer.txt`;
do
echo $i
Rscript ../puQTL/pu_act/3.sample_name/bin/filter_promoter.R ../puQTL/pu_act/3.sample_name/expression/TCGA-${i}.exp ../puQTL/pu_act/3.sample_name/promoter_filter/TCGA_${i}.exp_tmp
done

#get promoter list
for i in `ls ../puQTL/pu_act/3.sample_name/promoter_filter`
do
cut -f1 ../puQTL/pu_act/3.sample_name/promoter_filter/$i > ../puQTL/1.combine_data/location/promoter_list/$i
done

#impute normalize expression
mkdir normalize_impute
for i in `cat samples`;
do
echo $i
Rscript ../puQTL/pu_act/3.sample_name/bin/normalize_impute.R ../puQTL/pu_act/3.sample_name/promoter_filter/TCGA_${i}.exp_tmp ../puQTL/pu_act/3.sample_name/normalize_impute/TCGA_${i}.exp 
done

#omit 80% sample
#2)impute promoter
add_cancer <- c("KIRC","THCA","UCEC")
for (i in 1:length(add_cancer)) {
  exp = read.table(paste0("../puQTL/pu_act/3.sample_name/promoter_filter/TCGA_",add_cancer[i],".exp_tmp"),header=T,row.names=1,sep="\t",check.names=F)
  expSamples = apply(exp, 2, function(x) { length(x[is.na(x)==T]) } )
  exp = exp[,which(expSamples < nrow(exp)*0.8)]
  samples = colnames(exp)
  library(impute)
  .libPaths("/usr/lib64/R/library")
  library(preprocessCore)
  z_exp <- scale(exp, center = TRUE, scale = TRUE)
  library(preprocessCore)
  expNorm = normalize.quantiles(as.matrix(z_exp))
  rownames(expNorm) = rownames(z_exp)
  colnames(expNorm) = colnames(z_exp)
  expr <- impute.knn(as.matrix(expNorm))
  expr <- expr$data
  write.table(expr,file=paste0("../puQTL/pu_act/3.sample_name/normalize_impute/TCGA_",add_cancer[i],".exp"),row.names=T,col.names=T,sep="\t",quote=F)
}

#normalize expression
mkdir normalize_impute
for i in `cat samples`;
do
echo $i
Rscript ../puQTL/pu_act/3.sample_name/bin/normalize.R ../puQTL/pu_act/3.sample_name/promoter_filter/TCGA_${i}.exp_tmp ../puQTL/pu_act/3.sample_name/normalize_expression/TCGA_${i}.exp 
done

#omit 80% sample
#1)normalize promoter
rm(list = ls())
library(peer)
cancer_list <- c("KIRC","THCA","UCEC")
for (i in 1:3) {
  exp = read.table(paste0("../puQTL/pu_act/3.sample_name/promoter_filter/TCGA_",cancer_list[i],".exp_tmp"), row.names=1, header=T, check.names=F)
  library(preprocessCore)
  expSamples = apply(exp, 2, function(x) { length(x[is.na(x)==T]) } )
  exp = exp[,which(expSamples < nrow(exp)*0.8)]
  z_exp <- scale(exp, center = TRUE, scale = TRUE)
  library(preprocessCore)
  expNorm = normalize.quantiles(as.matrix(z_exp))
  rownames(expNorm) = rownames(z_exp)
  colnames(expNorm) = colnames(z_exp)
  fout <- paste0("../puQTL/pu_act/3.sample_name/normalize_expression/TCGA_",cancer_list[i],".exp")
  write.table(expNorm,file=fout,row.names=T,col.names=T,sep="\t",quote=F)
}














# ensure same samples order, get order from genotype
mkdir sampleIDs
for i in `cat samples`;
do
echo $i
head -1 ../puQTL/pu_act/3.sample_name/genotype/matrix/TCGA-${i}.snpmatrix | sed 's/\t/\n/g' | sed '1d' > sampleIDs/TCGA-${i}.samples
done

# #80% NA KIRC THCA UCEC
# for i in KIRC THCA UCEC;
# do
# echo $i
# head -1 ../puQTL/pu_act/3.sample_name/normalize_expression/TCGA_${i}.exp | sed 's/\t/\n/g'  > sampleIDs/TCGA-${i}.samples
# done

mkdir genotype
mkdir snp
cp ../puQTL/pu_act/3.sample_name/genotype/matrix/*  ../puQTL/1.combine_data/snp




# make eigen
mkdir eigen
for i in `cat samples`;
do
echo $i
cut -f1,3- ../puQTL/pu_act/3.sample_name/eigen/TCGA-${i}.eigenvec > eigen/TCGA-${i}.tmp
done

for i in `cat samples`;
do
echo $i
Rscript bin/format.eigen.R eigen/TCGA-${i}.tmp sampleIDs/TCGA-${i}.samples eigen/TCGA-${i}.eigen
done
rm -rf eigen/TCGA-*.tmp


# make covariates
mkdir covariates
for i in `cat samples`;
do
echo $i
Rscript bin/format.covariates.R ../puQTL/pu_act/3.sample_name/covariates/TCGA-${i}.covariates sampleIDs/TCGA-${i}.samples covariates/TCGA-${i}.tmp
done

for i in `cat samples`;
do
echo $i
sed 's/FEMALE/0/g' covariates/TCGA-${i}.tmp | sed 's/MALE/1/g' | sed 's/\[Not Available\]/NA/g' > covariates/TCGA-${i}.covariates
rm -rf covariates/TCGA-${i}.tmp
done

#make  expression
mkdir expression
for i in `cat samples`;
do
echo $i
Rscript bin/format.exp.R ../puQTL/pu_act/3.sample_name/normalize_expression/TCGA_${i}.exp sampleIDs/TCGA-${i}.samples expression/TCGA_${i}.exp
done

# make peer
mkdir peer_nor
for i in `cat samples`;
do
echo $i
Rscript bin/format.peer.R ../puQTL/pu_act/peer/factors/TCGA_${i}.exp sampleIDs/TCGA-${i}.samples peer_nor/TCGA-${i}.peer
done

# mkdir location
mkdir location
mkdir location/snp
mkdir location/promoter
# 1) snp locations
for i in `cat samples`;
do
echo $i
grep -v "#" ../puQTL/pu_act/3.sample_name/genotype/vcf/TCGA-${i}.vcf | awk '{print $3"\tchr"$1"\t"$2}' | sed "1isnpid\tchr\tpos" > location/snp/TCGA-${i}.snp.loc
done

# 2) promoter locations
for i in `ls ../puQTL/1.combine_data/expression`;
do
cut -f1 ../puQTL/1.combine_data/expression/$i > ../puQTL/1.combine_data/location/promoter_list/$i;
done

rm(list = ls())
setwd("../puQTL/1.combine_data/location/promoter_list/")
cancer_list <- list.files(pattern = "TCGA")
result <- data.frame()
loc <- read.delim("./share/TCGA/Pro/00.rawdata/promoterCoordinates.gencode.v19.csv",header = T)
loc <- loc[,c(1:3,5,6)]
for (i in 1:length(cancer_list)) {
  raw_data <- read.delim(paste0("../puQTL/1.combine_data/location/promoter_list/",cancer_list[i]),header = T)
  promoter_loc <- as.data.frame(raw_data)
  colnames(promoter_loc) <- "promoter"
  promoter_loc$pro <- gsub(promoter_loc$promoter,pattern="\\.ENSG.*",replacement="")
  colnames(promoter_loc) <- c("promoter","promoterId")
  result <- dplyr::inner_join(promoter_loc,loc,by="promoterId")
  result <- result[,c(1,3,5,4)]
  colnames(result) <- c("geneid","chr","left","right")
  write.table(result,paste0("../puQTL/1.combine_data/location/promoter/",cancer_list[i]),sep = "\t",quote=F,col.names = T,row.names = F)
}


# concate covariates
mkdir combine
for i in `cat samples`;
do
echo $i
cp covariates/TCGA-${i}.covariates combine/TCGA-${i}.covariates.tmp
sed '1d' eigen/TCGA-${i}.eigen > combine/TCGA-${i}.eigen.tmp
sed '1d' peer_nor/TCGA-${i}.peer > combine/TCGA-${i}.peer.tmp
cat combine/TCGA-${i}.covariates.tmp combine/TCGA-${i}.eigen.tmp combine/TCGA-${i}.peer.tmp > combine/TCGA-${i}.covariates
done

# as some types of cancer limited to single gender, MALE or FEMALE, such as CESC, OV, UCS, TGCT, PRAD, UCEC
for i in CESC OV UCS TGCT PRAD UCEC;
do
echo $i
mv combine/TCGA-${i}.covariates combine/TCGA-${i}.tmp
grep -v gender combine/TCGA-${i}.tmp > combine/TCGA-${i}.covariates
rm -rf combine/TCGA-${i}.tmp
done



