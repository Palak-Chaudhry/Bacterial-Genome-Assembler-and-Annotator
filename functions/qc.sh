#! usr/env/bin bash

## Create environment
conda env create -f ./environments/teama.yml
conda activate teama
input_dir="$2"
output_dir="$3"
if [ "$1" == "pre" ];
then
    ## QC raw reads
    mkdir $input_dir/pre_qc
    find $input_dir | grep ".fastq.gz" | xargs -n 1 -P 12 fastqc -o $input_dir/pre_qc
    multiqc $input_dir/pre_qc/ -o $input_dir/pre_qc/
else
    ## QC trmmed reads
    mkdir $input_dir/post_qc
    find "${input_dir}/trim/" | grep ".fastq.gz" | xargs -n 1 -P 12 fastqc -o $input_dir/post_qc
    multiqc $input_dir/post_qc/ -o $input_dir/post_qc/
fi

conda deactivate
