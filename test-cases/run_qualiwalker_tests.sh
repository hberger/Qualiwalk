INPUT_FOLDER=./input_data
REF_GENOME=/data_genome1/References/Human/Sequences/Genome/human_g1k_v37/human_g1k_v37.fasta.bgz
KNOWN_SNP=/data_genome1/References/Human/Variation/dbSNP/dbSNP_v142/All_20150217.vcf.gz
R_binary=/data_genome1/SharedSoftware/R/3.4/bin/Rscript
R_compscript=/data_genome1/SharedSoftware/inhouse_development/R/compare_tables.R
QW_BINARY=/data_genome1/SharedSoftware/inhouse_development/QualiWalker/QualiWalker.py

export REF_GENOME
export KNOWN_SNP
export R_binary
export R_compscript
export QW_BINARY

function quali_walk_test { 
	options=$3; 
	input=$2; 
	test_id=$1; 
	test ! -d output_data/${test_id} && mkdir output_data/${test_id}
	python $QW_BINARY -f $REF_GENOME -s $KNOWN_SNP -o output_data/${test_id}/${test_id}_qual_windows.txt $options $input &> output_data/${test_id}/${test_id}.log; 
	$R_binary $R_compscript expected_results/${test_id} output_data/${test_id} &> test_status/${test_id}.status
}

function quali_walk_test_w_vcf { 
	options=$3; 
	input=$2; 
	test_id=$1; 
	test ! -d output_data/${test_id} && mkdir output_data/${test_id}
	python $QW_BINARY -f $REF_GENOME -s $KNOWN_SNP -c output_data/${test_id}/${test_id} -o output_data/${test_id}/${test_id}_qual_windows.txt $options $input &> output_data/${test_id}/${test_id}.log; 
	$R_binary $R_compscript expected_results/${test_id} output_data/${test_id} &> test_status/${test_id}.status
}

export -f quali_walk_test
export -f quali_walk_test_w_vcf

# TEST 1 - normal region, strict 300bp 
region=1:1009400-1009699
#if [ ! -d ouput_data/test_1 ]; then mkdir output_data/test_1; fi
#python $QW_BINARY -f $REF_GENOME -r $region -s $KNOWN_SNP -o output_data/test_1/NA128787A_normal_region_only.txt ${INPUT_FOLDER}/NA12878A_normal_region.bam &> output_data/test_1/NA12787A_normal_region.log
#$R_binary $R_compscript expected_results/test_1 output_data/test_1 &> test_status/test_1.status
quali_walk_test test_1 ${INPUT_FOLDER}/NA12878A_normal_region.bam "-r $region"
		
# TEST 2 - low coverage region, strict 300bp 
region=1:1000300-1000599
#if [ ! -d ouput_data/test_2 ]; then mkdir output_data/test_2; fi
#python $QW_BINARY -f $REF_GENOME -r $region -s $KNOWN_SNP -o output_data/test_2/NA128787A_low_cov_region_only.txt ${INPUT_FOLDER}/NA12878A_low_cov_region.bam &> output_data/test_2/NA12787A_low_cov_region.log
#$R_binary $R_compscript expected_results/test_2 output_data/test_2 &> test_status/test_2.status
quali_walk_test test_2 ${INPUT_FOLDER}/NA12878A_low_cov_region.bam "-r $region"

# TEST 3 - normal region, strict 300bp; with VCF
region=1:1009400-1009699
#if [ ! -d ouput_data/test_3 ]; then mkdir output_data/test_3; fi
#python $QW_BINARY -f $REF_GENOME -r $region -c output_data/test_3/test_3 -s $KNOWN_SNP -o output_data/test_3/NA128787A_normal_region_only.txt ${INPUT_FOLDER}/NA12878A_normal_region.bam &> output_data/test_3/NA12787A_normal_region.log
#$R_binary $R_compscript expected_results/test_3 output_data/test_3 &> test_status/test_3.status
quali_walk_test_w_vcf test_3 ${INPUT_FOLDER}/NA12878A_normal_region.bam "-r $region"

