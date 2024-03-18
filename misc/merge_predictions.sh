#! usr/bin/env bash

##Merge Predictions
mkdir gene_prediction/merged_results
#Find the genes that overlap between the prodigal and glimmer predictions (with a minimum fraction of overlap of 0.95)
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/prodigal/{}.gff -b $output_dir/glimmer{}.gff -r -f 0.95 -wa > $output_dir/merged/{}_prodigal_glimmer.gff"
echo "prodigal and glimmer predictions completed"
#Find the genes that overlap between the prodigal and gms2 predictions (with a minimum fraction of overlap of 0.95)
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/prodigal/{}.gff -b $output_dir/gms2/{}.gff -r -f 0.95 -wa > $output_dir/merged/{}_prodigal_gms2.gff"
echo "prodigal and gms2 predictions completed"
#Find the genes that overlap between the gms2 and glimmer predictions (with a minimum fraction of overlap of 0.95)
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/gms2/{}.gff -b $output_dir/glimmer{}.gff -r -f 0.95 -wa > $output_dir/merged/{}_gms2_glimmer.gff"
echo "gms2 and glimmer predictions completed"
#Find the genes that are predicted by all three tools with a minimum fraction overlap of 0.95 (progigal, genemark, and glimmer)
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/merged/{}_prodigal_glimmer.gff -b $output_dir/gms2/{}.gff -r -f 0.95 -wa > $output_dir/merged/{}_prodigal_gms2_glimmer.gff"
echo "gms2 and prodigal and glimmer predictions completed"
#Find the genes predicted by prodigal and genemark that are not in glimmer
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/merged/{}_prodigal_gms2.gff -b $output_dir/merged/{}_prodigal_gms2_glimmer.gff -v > $output_dir/merged/{}_prodigal_gms2_only.gff"
echo "no glimmer predictions completed"
#Find the genes predicted by glimmer and genemark that are not in prodigal
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/merged/{}_gms2_glimmer.gff -b $output_dir/merged/{}_prodigal_gms2_glimmer.gff -v > $output_dir/merged/{}_glimmer_gms2_only.gff"
echo "no prodigal predictions completed"
#Find the genes predicted by prodigal and glimmer that are not in genemark
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools intersect -a $output_dir/merged/{}_prodigal_glimmer.gff -b $output_dir/merged/{}_prodigal_gms2_glimmer.gff -v > $output_dir/merged/{}_prodigal_glimmer_only.gff"
echo "no gms2 predictions completed"
#Merge all 4 DNA prediction files (prodigal_genemark_only, prodigal_glimmer_only, genemark_glimmer_only, and prodigal_genemark_glimmer)
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat $output_dir/merged/{}_prodigal_gms2_glimmer.gff $output_dir/merged/{}_prodigal_glimmer_only.gff $output_dir/merged/{}_prodigal_gms2_only.gff $output_dir/merged/{}_glimmer_gms2_only.gff | sort -n -k 4 > $output_dir/merged/{}_cds.gff"
echo "merge completed"
#Get fasta file for the coding regions
ls $input_dir| cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools getfasta -s -fi ./raw_reads_p/{}.fna -bed $output_dir/merged/{}_cds.gff -fo $output_dir/merged/{}_cds.fna"
echo "fasta files completed"
mv $output_dir/merged/*_cds.fna ./gene_prediction/
