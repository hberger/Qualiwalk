INPUT_FOLDER=./input_data
REF_GENOME=/data_genome1/References/Human/Sequences/Genome/human_g1k_v37/human_g1k_v37.fasta.bgz
KNOWN_SNP=/data_genome1/References/Human/Variation/dbSNP/dbSNP_v142/All_20150217.vcf.gz
R_binary=/data_genome1/SharedSoftware/R/3.4/bin/Rscript
R_compscript=/data_genome1/SharedSoftware/inhouse_development/R/compare_tables.R
PW_BINARY=/data_genome1/SharedSoftware/inhouse_development/QualiWalker/PositionWalker.py

export REF_GENOME
export KNOWN_SNP
export R_binary
export R_compscript
export PW_BINARY

function position_walk_test { 
	options=$3; 
	input=$2; 
	test_id=$1; 
	test ! -d output_data/${test_id} && mkdir output_data/${test_id}
	python $PW_BINARY -f $REF_GENOME -s $KNOWN_SNP -o output_data/${test_id}/${test_id}_qual_windows.txt $options $input &> output_data/${test_id}/${test_id}.log; 
	$R_binary $R_compscript expected_results/${test_id} output_data/${test_id} &> test_status/${test_id}.status
}

function position_walk_test_w_vcf { 
	options=$3; 
	input=$2; 
	test_id=$1; 
	test ! -d output_data/${test_id} && mkdir output_data/${test_id}
	python $PW_BINARY -f $REF_GENOME -s $KNOWN_SNP -c output_data/${test_id}/${test_id} -o output_data/${test_id}/${test_id}_qual_windows.txt $options $input &> output_data/${test_id}/${test_id}.log; 
	$R_binary $R_compscript expected_results/${test_id} output_data/${test_id} &> test_status/${test_id}.status
}

export -f position_walk_test
export -f position_walk_test_w_vcf

		
# TEST PW001 - normal region, strict 300bp, PositionWalker
region=1:1009400-1009699
TEST=test_PW001
position_walk_test $TEST ${INPUT_FOLDER}/NA12878A_normal_region.bam "-r $region"


