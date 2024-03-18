#! usr/env/bin bash

## Create environment
conda env create -f ./environments/teama.yml
conda activate teama
input_dir="$1"
output_dir="$2"
mkdir $output_dir/spades
##Assemble using spades
cat ../listx.txt | xargs -n 1 -P 46 -I {} sh -c "spades --phred-offset 33 -o $output_dir/spades/{}/ -1 $output_dir/trim/{}_R1.fq.gz -2 $output_dir/trim/{}_R2.fq.gz"
for i in $(cat ../list.txt); do mv $output_dir/spades/${i}/contigs.fasta $output_dir/spades/${i}_contigs.fasta ; done

conda deactivate