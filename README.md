# puQTL_TCGA
Abnormal promoter usage has been shown to be a key factor in cancer development and progression. In this research, we systematically analyzed large-scale cancer omics data to explore the impact of genetic variants on promoter usage, referred to as promoter usage quantitative trait loci (puQTLs).


## 0. Requirements
Software: vcftools version 0.1.16, plink version 1.9, PEER version 1.0, MatrixEQTL version 2.3, SnpEff version 5.0, bedtools version 2.29.0, mashr version 0.2.49,  motifBreakR version 2.4.0,  Coloc version 5.2.2, R version 4.0.3,  LDlinkR version 1.3.0.

## 1. Data preprocess for puQTL analysis
The shell scripts 01.preprocess.sh and 02.QTL_preprocess.sh were employed to preprocess individual genotype data and promoter usage for puQTL analysis. Key variables, including principal components of the genotype, PEER factors, gender, and age, were all accounted for in preparation for identifying puQTLs.

## 2. puQTL mapping for each cancer
The R script 03.MatrixeQTL.R was used to identify puQTLs for each tissue type. The MatrixeQTL R package was employed to test associations between genotype and promoter usage. 

## 3. Cancer specificity analysis for puQTLs
The 04.tissue_specific.R utilized the mash method to elucidate the heterogeneity of puQTL effect sizes across different cancer types.

## 7. puQTL analysis in tumorigenesis 
The R script 05.normal_cancer.R was used to compare the puQTL in normal and cancer sample. 

## 4. Enrichment analysis for puQTLs
The script 06.puQTL_SNPeff.sh  06.TF.R was used to evaluate enrichment of puQTLs in experimentally annotated epigenomic regulatory features and transcription factors binding region.

## 5. colocalization analysis for puQTL in cancer
The R script 07.coloc.R was used for colocalization analysis of puQTLs with GWAS associated cancer diseases.

## 6. causal analysis for puQTL in cancer
The shell script 08.SMR.sh was used to investigate the potential causal puQTLs and promoters in cancer.

## 7. The role of puQTLs in immunotherapy
The R script 09.CIBERSORT.R was used to investigate the potential effect of puQTLs on immune cell.
