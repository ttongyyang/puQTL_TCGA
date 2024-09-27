! /bin/bash
#get data format data
cp ./share/TCGA/Pro/01.cleandata/step2/* ../puQTL/pu_act/raw_data

for i in `ls ../puQTL/pu_act/raw_data`;
do
head -1 ../puQTL/pu_act/raw_data/$i | sed 's/\t/\n/g' | awk -F "-" '{print $1"."$2"."$3}' |awk '{printf "%s ",$0}' | sed 's/ /\t/g' | awk -F "\t" '{printf $0"\n"}' > colname_${i};
sed 1d ../puQTL/pu_act/raw_data/$i > no_col_${i}; 
cat colname_${i} no_col_${i} > ../puQTL/pu_act/2.re_name/$i; 
rm colname_${i}; rm no_col_${i}; 
done

#get and add sample state and quality(A,B,C) and label(TCGA.OR.A5KO) 
for i in `ls ../puQTL/pu_act/raw_data `
do 
head -1 ../puQTL/pu_act/raw_data/$i | sed 's/\t/\n/g' | awk -F "-" '{print $1"."$2"."$3}' |awk '{printf "%s ",$0}' | sed 's/ /\t/g' | awk -F "\t" '{printf $0"\n"}' |sed 's/Promoter../gene/g' | sed 's/\t$//g' > colname_${i}
head -1 ../puQTL/pu_act/raw_data/$i  | sed 's/\t/\n/g' | awk -F "-" '{print $4}' | awk -F "" '{print $1$2}'  | sed '1,1s/^/tumor/' | awk '{printf "%s ",$0}' | sed 's/ /\t/g' | awk -F "\t" '{printf $0"\n"}' | sed 's/\t$//g'  > tumor_${i}
head -1 ../puQTL/pu_act/raw_data/$i  | sed 's/\t/\n/g' | awk -F "-" '{print $4}' | awk -F "" '{print $3}'  | sed '1,1s/^/quality/' | awk '{printf "%s ",$0}' | sed 's/ /\t/g' | awk -F "\t" '{printf $0"\n"}'| sed 's/\t$//g' > quality_${i}
sed 1d ../puQTL/pu_act/raw_data/$i > no_col_${i}
head -1 ../puQTL/pu_act/raw_data/$i  > col_${i}
cat col_${i} tumor_${i} quality_${i} colname_${i} no_col_${i} > ../puQTL/pu_act/2.re_name/$i
rm tumor_${i}
rm quality_${i}
rm no_col_${i}
rm col_${i}
rm colname_${i}
done


#choose tumor delete health sample
rm(list = ls())
setwd("../puQTL/pu_act/2.re_name/")
cancer_list <- list.files(pattern = "TCGA")
result <- data.frame()
for (i in 1:length(cancer_list)) {
  re_name <- read.delim(paste0("../puQTL/pu_act/2.re_name/",cancer_list[i]))
  rownames(re_name) <- re_name[,1]
  result <- re_name[,re_name[1,]<10] #choose tumor
  result <- result[3,]
  write.table(result,paste0("../puQTL/pu_act/1.pu_sample_name/",cancer_list[i]),sep = "\t",col.names = F,row.names = F,quote = F)
}

for i in `ls ../puQTL/pu_act/1.pu_sample_name | grep TCGA`; do sed -i 's/\t/\n/g' $i; done
for i in `ls ../puQTL/pu_act/1.pu_sample_name | grep TCGA`; do sort ../puQTL/pu_act/1.pu_sample_name/$i | uniq > ../puQTL/pu_act/1.pu_sample_name/uniq/$i; done


# Same samples for everything(overlap)
# 1) genotype
mkdir sampleIDs
for i in `cat samples`;
do
echo $i
cut -f1 ./project/QTL/01.imputation/filter/step2/TCGA-${i}.ped | awk -F "-" '{print $1"-"$2"-"$3}' > sampleIDs/TCGA-${i}.snp.sampleIDs
sort sampleIDs/TCGA-${i}.snp.sampleIDs | uniq > sampleIDs/TCGA-${i}.snp.samples
done


# 3) age and gender 
for i in `cat samples`;
do
echo $i
cut -f1 ./QTL/04.otherCovarites/covariates/TCGA-${i}.covariates | sed "1d" | sort | uniq > sampleIDs/TCGA-${i}.covariates.samples
done

# 3) promoter activity 
for i in `ls ../puQTL/pu_act/2.re_name`;
do 
sed -n '4p' ../puQTL/pu_act/2.re_name/$i > ../puQTL/pu_act/1.pu_sample_name/$i;
sed -i 's/gene\t//g' ../puQTL/pu_act/1.pu_sample_name/$i; 
done

for i in `ls ../puQTL/pu_act/1.pu_sample_name`;
do 
sed -i 's/\t/\n/g' $i; 
done

for i in `ls ../puQTL/pu_act/1.pu_sample_name | grep TCGA`; 
do 
sort ../puQTL/pu_act/1.pu_sample_name/$i | uniq > ../puQTL/pu_act/1.pu_sample_name/uniq/$i;
done

sed -i 's/\./-/g' ../puQTL/pu_act/1.pu_sample_name/uniq/*
  
  cp ../puQTL/pu_act/1.pu_sample_name/uniq/* ../puQTL/pu_act/3.sample_name/sampleIDs/
  
  
  
  # 4) overlap
  mkdir overlap
for i in `cat samples`;
do
echo $i
cat sampleIDs/TCGA-${i}.snp.samples sampleIDs/TCGA-${i} sampleIDs/TCGA-${i}.covariates.samples | sort | uniq -c | awk '$1==3{print $2}' > overlap/TCGA-${i}.samples
done


# construct expression matrix
rm(list = ls())
setwd("../puQTL/pu_act/2.re_name/")
cancer_list <- list.files(pattern = "TCGA")
result <- data.frame()
for (i in 1:length(cancer_list)) {
  re_name <- read.delim(paste0("../puQTL/pu_act/2.re_name/",cancer_list[i]))
  rownames(re_name) <- re_name[,1]
  re_name <- re_name[,-1]
  #tumor sample 
  result <- re_name[,re_name[1,]<10]
  #tumor sample high quality
  result <- result[,result[2,]=="A"]
  #tumor sample B quality
  cand_B <- as.data.frame(re_name[,re_name[2,]!="A"]) 
  colnames(cand_B) <- colnames(re_name)[which(re_name[2,]!="A")]
  
  add <- as.data.frame(cand_B[,!(cand_B[3,] %in% result[3,])])
  colnames(add) <- colnames(cand_B)[which(!(cand_B[3,] %in% result[3,]))]
  if (nrow(add)!=0) {
    add_B <- as.data.frame(add[,add[2,]=="B"])
    colnames(add_B) <- colnames(add)[which(add[2,]=="B")]
    if (nrow(add_B)!=0) {
      result <- cbind(result,add_B) 
    }
    cand_C <- as.data.frame(add[,add[2,]!="B"])
    if (nrow(cand_C)!=0) {
      colnames(cand_C) <- colnames(add)[which(add[2,]!="B")]
      add_C <- as.data.frame(cand_C[,!(cand_C[3,] %in% add_B[3,])])
      colnames(add_C) <- colnames(cand_C)[which(!(cand_C[3,] %in% add[3,]))]
      result <- cbind(result,add_C)
    }
  }
  # print(length(unique(colnames(result))))
  # 
  if (length(unique(as.character(result[3,])))!=length(as.character(result[3,]))) {
    dup_list <- as.character(result[3,])[duplicated(as.character(result[3,]))]
    new_result <- result[,!(as.character(result[3,])%in% dup_list)]
    dep_result <- result[,as.character(result[3,])%in% dup_list]
    for (k in 1:length(dup_list)) {
      part_result <- dep_result[,dep_result[3,] %in% dup_list[k]]
      new_result <- cbind(new_result,part_result[,which.min(part_result[1,])])
    }
    result <- new_result
  } 
  colnames(result) <- result[3,]
  result <- result[-c(1,2,3),]
  print(ncol(result))
  print(length(unique(colnames(result))))
  write.table(result,paste0("../puQTL/pu_act/4.pu_act/",cancer_list[i]),sep = "\t",col.names = T,row.names = T,quote = F)
}

# rm(list = ls())
# setwd("../puQTL/pu_act/2.re_name/")
# cancer_list <- list.files(pattern = "TCGA")
# result <- data.frame()
# for (i in 1:length(cancer_list)) {
#   re_name <- read.delim(paste0("../puQTL/pu_act/2.re_name/",cancer_list[i]))
#   rownames(re_name) <- re_name[,1]
#   re_name <- re_name[,-1]
#   #tumor sample 
#   result <- re_name[,re_name[1,]<10]
#   #tumor sample high quality
#   result <- result[,result[2,]=="A"]
#   #tumor sample B quality
#   cand_B <- as.data.frame(re_name[,re_name[2,]!="A"]) 
#   colnames(cand_B) <- colnames(re_name)[which(re_name[2,]!="A")]
#   
#   add <- as.data.frame(cand_B[,!(cand_B[3,] %in% result[3,])])
#   colnames(add) <- colnames(cand_B)[which(!(cand_B[3,] %in% result[3,]))]
#   if (nrow(add)!=0) {
#   add_B <- as.data.frame(add[,add[2,]=="B"])
#   colnames(add_B) <- colnames(add)[which(add[2,]=="B")]
#   if (nrow(add_B)!=0) {
#     result <- cbind(result,add_B) 
#   }
#   cand_C <- as.data.frame(add[,add[2,]!="B"])
#   if (nrow(cand_C)!=0) {
#     colnames(cand_C) <- colnames(add)[which(add[2,]!="B")]
#     add_C <- as.data.frame(cand_C[,!(cand_C[3,] %in% add_B[3,])])
#     colnames(add_C) <- colnames(cand_C)[which(!(cand_C[3,] %in% add[3,]))]
#     result <- cbind(result,add_C)
#   }
#   }
#   # print(length(unique(colnames(result))))
#   # 
#   if (length(unique(as.character(result[3,])))!=length(as.character(result[3,]))) {
#     dup_list <- as.character(result[3,])[duplicated(as.character(result[3,]))]
#     new_result <- result[,!(as.character(result[3,])%in% dup_list)]
#     dep_result <- result[,as.character(result[3,])%in% dup_list]
#     for (k in 1:length(dup_list)) {
#       part_result <- dep_result[,dep_result[3,] %in% dup_list[k]]
#       new_result <- cbind(new_result,part_result[,which.min(part_result[1,])])
#     }
#     result <- new_result
#   } 
#   colnames(result) <- result[3,]
#   result <- result[-c(1,2,3),]
#   print(ncol(result))
#   print(length(unique(colnames(result))))
#   write.table(result,paste0("../puQTL/pu_act/4.pu_act/",cancer_list[i]),sep = "\t",col.names = T,row.names = T,quote = F)
# }
for i in `ls ../puQTL/pu_act/4.pu_act`; do sed -i '1,1s/^/gene\t/g' $i; done
for i in `ls`; do sed -i '1s/\./-/g' $i; done

mkdir expression
for i in `cat samples`;
do
echo $i
Rscript ../puQTL/pu_act/3.sample_name/bin/expression.format.R ../puQTL/pu_act/4.pu_act/TCGA-${i} overlap/TCGA-${i}.samples expression/TCGA-${i}.exp
done

#run factors








#because omit 80% sample overlap sample change
for i in `cat samples`;
do
echo $i
head -1 ../puQTL/pu_act/3.sample_name/normalize_expression/TCGA_${i}.exp | sed 's/\t/\n/g' > ../puQTL/pu_act/3.sample_name/new_overlap/TCGA-${i}.samples
done

# construct covariates
#../puQTL/pu_act/3.sample_name/
mkdir covariates
for i in `cat samples`;
do
echo $i
Rscript bin/covariates.format.R ./project/QTL/04.otherCovarites/covariates/TCGA-${i}.covariates ../puQTL/pu_act/3.sample_name/new_overlap/TCGA-${i}.samples covariates/TCGA-${i}.covariates
done


# construct snp
mkdir genotype
mkdir genotype/raw
for i in `cat samples`;
do
echo $i
grep -F -w -f new_overlap/TCGA-${i}.samples ./project/QTL/01.imputation/filter/step2/TCGA-${i}.ped | awk -F " " '{split($1,a,"-");print a[1]"-"a[2]"-"a[3]" "a[1]"-"a[2]"-"a[3]" "$0}' | cut -d " " -f1,2,5- > genotype/raw/TCGA-${i}.ped
cp ./project/QTL/01.imputation/filter/step2/TCGA-${i}.map genotype/raw
done

#2_12
# clean snp samples
mkdir genotype/ped
for i in `cat samples`;
do
echo $i
python bin/snp.clean.py genotype/raw/TCGA-${i}.ped genotype/ped/TCGA-${i}.ped
done


# convert snp in plink format to vcf format using plink v1.9
mkdir genotype/vcf
for i in `cat samples`;
do
echo $i
/home/tongyang/apps/plink/plink --file genotype/ped/TCGA-${i} --recode vcf --out genotype/vcf/TCGA-${i}
done


# make matrix for snp
mkdir genotype/matrix
for i in `cat samples`;
do
echo $i
grep -v "##" genotype/vcf/TCGA-${i}.vcf > genotype/vcf/TCGA-${i}.tmp
python bin/parse.py -v genotype/vcf/TCGA-${i}.tmp -o genotype/matrix/TCGA-${i}.snpmatrix
done  


mkdir eigen
mkdir eigen/raw
for i in `cat samples`;
do
echo $i
awk '{split($1,a,"-");print a[1]"-"a[2]"-"a[3]"\t"a[1]"-"a[2]"-"a[3]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ./project/QTL/03.eigen/pca/TCGA-${i}.eigenvec | grep -F -w -f new_overlap/TCGA-${i}.samples > eigen/raw/TCGA-${i}.eigenvec
done


for i in `cat samples`;
do
echo $i
python bin/eigen.clean.py eigen/raw/TCGA-${i}.eigenvec eigen/TCGA-${i}.eigenvec
done
