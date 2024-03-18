#! usr/env/bin bash

## Create environment
conda env create -f ./environments/teama.yml
conda activate teama

input_dir="$1"
output_dir="$2"
mkdir $output_dir/trim
## Trimming using fastp
cat ../listx.txt | xargs -n 1 -P 24 -I {} sh -c " fastp -i $output_dir/filteredreads/{}_filtered_R1.fastq.gz -I $output_dir/filteredreads/{}_filtered_R2.fastq.gz -o $output_dir/trim/{}_R1.fq.gz -O $output_dir/trim/{}_R2.fq.gz -3 -M 30" 

conda deactivate