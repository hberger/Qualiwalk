BAM_FOLDER=/data_genome1/SharedSoftware/inhouse_development/QualiWalker/test
INPUT_FILES="1168_A_merged_dedup_realign_STXBP6_region.bam 1168_C_merged_dedup_realign_STXBP6_region.bam"
REF_GENOME=/data_genome1/References/Human/Sequences/Genome/human_g1k_v37/human_g1k_v37.fasta.bgz
KNOWN_SNP=/data_genome1/References/Human/Variation/dbSNP/dbSNP_v142/All_20150217.vcf.gz
#chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

export REF_GENOME
export BAM_FOLDER
export KNOWN_SNP

chroms=14:25281306-25519095

for f in $INPUT_FILES; do 
		FF=${BAM_FOLDER}/${f}
		sample_name=$(basename $f)
		echo ${BAM_FOLDER}/${FF}

		#function quali_walk { c=$1;  python /data_genome1/SharedSoftware/inhouse_development/QualiWalker/QualiWalker.py -f $REF_GENOME -r $c -s $KNOWN_SNP -c ${sample_name}_${c} -o ${sample_name}_${c}_qual_windows.txt $FF; }
		function quali_walk { c=$1;  python /data_genome1/SharedSoftware/inhouse_development/QualiWalker/QualiWalker.py -f $REF_GENOME -r $c -s $KNOWN_SNP -o ${sample_name}_${c}_qual_windows.txt $FF; }
		
		export -f quali_walk
		export FF
		export sample_name

		echo $chroms | sed -e 's/ /\n/g' | parallel -t -j 6 quali_walk

		#head -n 1 ${sample_name}_1_qual_windows.txt > ${sample_name}_all_qual_windows.txt
		#for i in ${sample_name}_{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}_qual_windows.txt; do tail -n +2 $i  >> ${sample_name}_all_qual_windows.txt; done

		#bcftools concat -O z ${sample_name}_{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}.vcf > ${sample_name}_qw_varcalls_all.vcf.bgz
		#tabix -p vcf ${sample_name}_qw_varcalls_all.vcf.bgz

		#rm ${sample_name}_{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}.vcf*
		#rm ${sample_name}_{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}_qual_windows.txt
done


