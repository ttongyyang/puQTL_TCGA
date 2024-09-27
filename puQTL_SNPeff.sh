#/bin/bash 
for i in `ls /data/tongyang/puQTL/QTL_result/0.rawdata/significant/cis/correct/vcf`
do 
	        k=${i%.*}
		java -Xmx8g -jar ./packages/biosoft/snpeff/snpEff/snpEff.jar -c ./packages/biosoft/snpeff/snpEff/snpEff.config GRCh37.75 ../puQTL/QTL_result/0.rawdata/significant/cis/correct/vcf/${k}.vcf -verbose -stats ../puQTL/5.analysis/snpEff/0.rawdata/${k}.html -csvStats ../puQTL/5.analysis/snpEff/0.rawdata/${k}.csv > ../puQTL/5.analysis/snpEff/0.rawdata/${i}.vcf
	done

