#!/bin/bash

##Extract 16S from files
##Extract 16SRNA 
# ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "cat $output_dir/barrnap/{}.gff | grep "Name=16S_rRNA;product=16S ribosomal RNA" > $output_dir/barrnap/{}_16S.gff"
# ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "bedtools getfasta -fi $input_dir/{}.fna -bed $output_dir/barrnap/{}_16S.gff -fo $output_dir/barrnap/{}_16S.fna"
# ##Strain_Elizabethkingia_anophelis_R26

# # create reference database for blasting from ensemble ref seqs
# #this one is from the .dna.toplevel.fa
# wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-58/bacteria//fasta/bacteria_99_collection/elizabethkingia_anophelis_r26_gca_002023665/dna/Elizabethkingia_anophelis_r26_gca_002023665.ASM202366v2_.dna.toplevel.fa.gz
# gunzip Elizabethkingia_anophelis_r26_gca_002023665.ASM202366v2_.dna.toplevel.fa.gz
# makeblastdb -in Elizabethkingia_anophelis_r26_gca_002023665.ASM202366v2_.dna.toplevel.fa -dbtype nucl -out  ./tools/Elizabethkingia_Anophelis_R26_refdb_dna

# benchmarking
# bash ./misc/evaluation.sh $input_dir/ ./tools/Elizabethkingia_Anophelis_R26_refdb_dna prodigal
# Usage: ./misc/evaluation.sh <input_files_directory> <input reference db directory>

# Step 1: Define variables
input_files_dir=$1
reference_db=$2
tool=$3

# Ensure output directory for BLAST results exists
output_dir="blast_results"
mkdir -p $output_dir

# Loop over all .faa files in the input directory
for predicted_genes in "$input_files_dir"/*.fa; do
    # Extract the base filename without the directory path and extension
    base_name=$(basename "$predicted_genes" .fa)
    
    # Define output file name based on input file
    blast_output="$output_dir/${base_name}_${tool}.txt"
    
    # Run tBLASTN
    # tblastn \
    # -query "$predicted_genes" \
    # -db "$reference_db" \
    # -outfmt '6 std qlen' \
    # -task tblastn-fast \
    # | awk '$3>80 && $4>0.9*$13' > "$blast_output"
    blastn \
    -query "$predicted_genes" \
    -db "$reference_db" \
    -outfmt '6 std qlen slen' \
    -task blastn \
    | awk '$3>80 && $4>0.9*$13' > "$blast_output"

    # Count the number of unique hits
    true_positives=$(wc -l < "$blast_output")
    
    # Count the number of predicted genes
    total_predictions=$(grep "^>" "$predicted_genes" | wc -l)

    # Calculate false positives (assuming all predictions should have a hit)
    false_positives=$(($total_predictions - $true_positives))

    #Calculate sensitivity
    real_genes=3790
    sensitivity=$(echo "scale=4; $true_positives / $real_genes" | bc)
    
    
    echo "File: $predicted_genes"
    echo "Total Predictions: $total_predictions"
    echo "Unique Hits (True Positives): $true_positives"
    echo "False Positives: $false_positives"
    echo "Sensitivity: $sensitivity"
    echo "-------------------------------------"

    echo "$true_positives" >> $output_dir/${tool}_tp.txt
    echo "$total_predictions" >> $output_dir/${tool}_total.txt
    echo "$false_positives" >> $output_dir/${tool}_fp.txt
    echo "$sensitivity" >> $output_dir/${tool}_s.txt
done