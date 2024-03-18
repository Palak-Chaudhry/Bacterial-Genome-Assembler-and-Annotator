#! usr/env/bin bash

## Create environment
conda env create -f ./environments/filtercontig.yml
conda activate filtercontig
input_dir="$1"
output_dir="$2"
mkdir $output_dir/filtered_spades
## Filter contigs using filter.contig.py
cat listx.txt | xargs -n 1 -P 46 -I {} sh -c " python filter.contigs.py -i $output_dir/spades/{}_contigs.fasta -o $output_dir/filtered_spades/{}.fna > $output_dir/filtered_spades/{}.stdout "
cat listx.txt | xargs -n 1 -P 46 -I {} sh -c " python filter.contigs.py -i $output_dir/spades/{}_contigs.fasta -o $output_dir/filtered_spades/{}_1000.fna > $output_dir/filtered_spades/{}_1000.stdout "
cat listx.txt | xargs -n 1 -P 46 -I {} sh -c " python filter.contigs.py -i $output_dir/spades/{}_contigs.fasta -o $output_dir/filtered_spades/{}_300.fna > $output_dir/filtered_spades/{}_300.stdout "

conda deactivate