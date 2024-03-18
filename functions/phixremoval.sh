## Create environment
conda env create -f ./environments/phixremoval.yml
conda activate phixremoval
input_dir="$1"
output_dir="$2"
mkdir $output_dir/{phix,filteredreads}
#PhiX removal
bowtie2-build GCF_000819615.1_ViralProj14015_genomic.fna PhixIndex

cat ../list.txt | xargs -n 1 -P 24 -I {} sh -c "bowtie2 --local -t -x PhixIndex -1 $input_dir/{}_R1.fastq.gz -2 $input_dir/{}_R2.fastq.gz --un ./ > $output_dir/phix/{}_filtered.sam"
cat ../listx.txt | xargs -n 1 -P 24 -I {} samtools fastq $output_dir/phix/{}_filtered.sam -1 $output_dir/filteredreads/{}_filtered_R1.fastq.gz -2 $output_dir/filteredreads/{}_filtered_R2.fastq.gz

conda deactivate