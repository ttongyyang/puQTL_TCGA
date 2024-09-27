#!/bin/bash
while read line
do
	for i in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 MT X Y
	do
		./packages/smr-1.3.1-linux-x86_64/smr-1.3.1 --bfile ../ref_seq/EUR_bfile/1000G.EUR.QC.${i} --gwas-summary ../puQTL/5.analysis/Mendelian_randomization/SMR/gwas/fin/${line}.ma --beqtl-summary ../puQTL/5.analysis/SMR/1.puQTL/TCGA-LUAD --out ../puQTL/5.analysis/SMR/2.result/${i}/${i}.txt --thread-num 10
	done
done < trait.txt
