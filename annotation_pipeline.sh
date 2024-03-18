#! usr/bin/env bash

# ##Create Folders
# mkdir gene_prediction/{prodigal,glimmer,mga,barrnap,aragorn} gene_annotation/{pilercr,signalp,tmhmm,eggnog,vfdb,card} tools

# ##Create environment
# conda install -c bioconda -c conda-forge prodigal metagene_annotator barrnap aragorn bedtools python3 pigz piler-cr eggnog-mapper rgi -y

parse_args() {

    # Function to parse arguments
    # Specifying usage message
    # -c : Compare tools (set -t to all for this mode)
    usage="Usage: bash annotation_pipeline.sh -i <input directory> -o <output directory> -[OPTIONS]
              Bacterial Gene Prediction for Illumina short-reads. The options available are:
                        -i : Input Directory for genome prediction (directory containing assembled fna)[required]
                        -o : Output directory for all results [required]
                        -t : Tool to choose from [prodigal/metageneannotator/glimmer] (defaults to Prodigal)
                        -v : Flag to turn on verbose mode
                        -h : Print usage instructions"

  #Specifying default Arguments
  t=prodigal
  v=0
  c=0

  #Getopts block, will take in the arguments as inputs and assign them to variables
  while getopts "i:o:t:cvh" option; do
          case $option in
                  i) input_dir=$OPTARG;;
                  o) output_dir=$OPTARG;;
                  t) tool=$OPTARG;;
                #   c) c=1;;
                  v) v=1;;
                  h) echo "$usage"
                        exit 0;;
                 \?) echo "Invalid option."
                    "$usage"
                             exit 1;;
          esac
  done


  #Check for presence of required arguments
  if [ ! "$input_dir" ] || [ ! "$output_dir" ]
  then
    echo "ERROR: Required arguments not provided!"
    echo "$usage"
    exit 1
  fi

  if [ ! -d "$input_dir" ]
  then 
    echo "ERROR: Not a valid input reads directory"
    echo "$usage"
    exit 1
  fi

#   if [ "$c" == 1 and "$t" != "all"]
#   then
#     echo "ERROR: set -t to all"
#     echo "$usage"
#     exit 1
#   fi
}

##GENE PREDICTION
run_prodigal(){
    if  [ $v == 1 ]
    then
	   echo "Running Prodigal"
    fi

    mkdir $output_dir/prodigal
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "prodigal -m -c -f gff -i $i/{}.fna -o ./$output_dir/prodigal/{}.gff "
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "prodigal -m -c -f gff -i $i/{}.fna -a ./$output_dir/prodigal_results/{}.faa "
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools getfasta -fi $i/{}.fna -bed ./$output_dir/prodigal/{}.gff -name -fo ./$output_dir/prodigal/{}.fna"

    if  [ $v == 1 ]
    then
	   echo "Prodigal completed. Results will be available in $output_dir/prodigal (.gff, .fna, .faa)"
    fi
}

run_glimmer(){
    if  [ $v == 1 ]
    then
	   echo "Downloading Glimmer"
    fi
    wget https://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz
    tar xzf glimmer302b.tar.gz
    rm glimmer302b.tar.gz
    cd glimmer3.02/src
    make
    cd ..
    mkdir temp $output_dir/glimmer
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "./glimmer3.02/bin/build-icm ./temp/{}.icm < $input_dir/{}.fna"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "./glimmer3.02/bin/glimmer3 $input_dir/{}.fna ./temp/{}.icm ./temp/{}"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bash ./functions/glimmer_to_gff.sh ./$output_dir/glimmer/{}.predict ./$output_dir/glimmer/{}.gff"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools getfasta -fi $input_dir/{}.fna -bed ./$output_dir/glimmer/{}.gff -name -fo ./$output_dir/glimmer/{}.fna"
    rm -r glimmer3.02 
    rm ./$output_dir/glimmer/*.detail
    rm ./$output_dir/glimmer/*.predict
    if  [ $v == 1 ]
    then
	   echo "Glimmer completed. Results will be available in $output_dir/glimmer (.gff, .fna)"
    fi
}

run_mga(){
    if  [ $v == 1 ]
    then
	   echo "starting metagene annotator"
    fi
    mkdir $output_dir/mga
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "mga $input_dir/{}.fna -s > ./$output_dir/mga/{}.out"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bash ./functions/mga_to_gff3.sh ./$output_dir/mga/{}.out ./$output_dir/mga/{}.gff"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "bedtools getfasta -fi $input_dir/{}.fna -bed ./$output_dir/mga/{}.gff -name -fo ./$output_dir/mga/{}.fna"
    rm $output_dir/mga/*.out
    if  [ $v == 1 ]
    then
	   echo "Metagene Annotator completed. Results will be available in $output_dir/mga (.gff, .fna)"
    fi
}

# run_all(){
#     run_prodigal $input_dir $output_dir
#     run_glimmer $input_dir $output_dir
#     run_mga $input_dir $output_dir
# }

run_barrnap(){
    if  [ $v == 1 ]
    then
	   echo "starting barrnap"
    fi
    mkdir $output_dir/barrnap
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "barrnap --outseq $output_dir/barrnap/{}.fna $input_dir/{}.fna > ./$output_dir/barrnap/{}.gff"
    if  [ $v == 1 ]
    then
	   echo "Barrnap completed. Results will be available in $output_dir/barrnap (.gff, .fna)"
    fi
}

run_aragorn(){
    if  [ $v == 1 ]
    then
	   echo "starting aragorn"
    fi
    mkdir $output_dir/aragorn
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "aragorn -v -fons -o ./$output_dir/aragorn/{}.fna $input_dir/{}.fna"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "bash ./functions/aragorn_to_gff.sh ./$output_dir/aragorn/{}.out.txt ./$output_dir/aragorn/{}.gff"
    rm ./$output_dir/aragorn/*.out.txt
    if  [ $v == 1 ]
    then
	   echo "Aragorn completed. Results will be available in $output_dir/aragorn (.gff, .fna)"
    fi

}

run_merge_files(){
    mkdir $output_dir/merged
    if  [ $v == 1 ]
    then
	   echo "Merging coding and non coding gene files"
    fi
    if [ "$tool" == "prodigal" ]
    then
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/prodigal/{}.fna ./$output_dir/barrnap/{}.fna ./$output_dir/aragorn/{}.fna > ./$output_dir/merged/{}.fna"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/prodigal/{}.gff ./$output_dir/barrnap/{}.gff ./$output_dir/aragorn/{}.gff > ./$output_dir/merged/{}.gff"
    fi
    if [ "$tool" == "glimmer" ]
    then
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/glimmer/{}.fna ./$output_dir/barrnap/{}.fna ./$output_dir/aragorn/{}.fna > ./$output_dir/merged/{}.fna"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/glimmer/{}.gff ./$output_dir/barrnap/{}.gff ./$output_dir/aragorn/{}.gff > ./$output_dir/merged/{}.gff"
    fi
    if [ "$tool" == "metageneannotator" ]
    then
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/mga/{}.fna ./$output_dir/barrnap/{}.fna ./$output_dir/aragorn/{}.fna > ./$output_dir/merged/{}.fna"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat ./$output_dir/mga/{}.gff ./$output_dir/barrnap/{}.gff ./$output_dir/aragorn/{}.gff > ./$output_dir/merged/{}.gff"
    fi
    if  [ $v == 1 ]
    then
	   echo "Merging completed. Results will be available in $output_dir/merged (.gff, .fna)"
    fi
}


##GENE ANNOTATION
run_eggnog(){

    if  [ $v == 1 ]
    then
	   echo "Eggnog started. This will take time as downloading database"
    fi
    mkdir temp $2/eggnog
    download-eggnog-data.py -dbname Bacteria --data_dir
    create_dbs.py --dbname bacteria --taxa Bacteria --data-dir
    if  [ $v == 1 ]
    then
	   echo "Database downloaded. Running annotation"
    fi
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c \
    "emapper.py --cpu 64 --data-dir --dmnd-db -i $input_dir/{} -o $2/eggnog/{}/" 
    rm bacteria*
    if  [ $v == 1 ]
    then
	   echo "annotation completed. Results will be available in $output_dir/eggnog"
    fi
}

##Convert file
run_pilercr(){
    if  [ $v == 1 ]
    then
	   echo "Pilercr started"
    fi
    mkdir $2/pilercr
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c \
    "pilercr -in ./$output_dir/merged/{}.fna -out ./$output_dir/pilercr/{}.log"
    if  [ $v == 1 ]
    then
	   echo "pilercr completed. Results will be available in $output_dir/pilercr"
    fi
}

run_signalp(){
    if  [ $v == 1 ]
    then
	   echo "Signalp started"
    fi
    for i in $(ls $input_dir | cut -d "." -f1); 
    do signalp6 -ff ./$output_dir/prodigal/${i}.faa -od ./$output_dir/signalp/${i}/ -fmt none;
    done
    for i in $(ls $input_dir | cut -d "." -f1); 
    do mv ./$output_dir/signalp/${i}/output.gff3 ./$output_dir/signalp/${i}/${i}.gff;
    done
    if  [ $v == 1 ]
    then
	   echo "Signalp completed. Results will be available in $output_dir/signalp"
    fi
}

##Convert output
run_card(){
    if  [ $v == 1 ]
    then
	   echo "Card rgi downloading database"
    fi
    mkdir ./functions/card
    cd ./functions/card
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
    mkdir -p wildcard
    tar -xjf wildcard_data.tar.bz2 -C wildcard
    gunzip wildcard/*.gz
    rgi card_annotation -i ./card.json > card_annotation.log 2>&1
    rgi wildcard_annotation -i wildcard --card_json ./card.json -v 4.0.0 > wildcard_annotation.log 2>&1
    #Loading all databases 
    rgi load   --card_json ./card.json   --debug --local   --card_annotation card_database_v3.2.9.fasta   \
    --card_annotation_all_models card_database_v3.2.9_all.fasta   --wildcard_annotation wildcard_database_v4.0.0.fasta   \
    --wildcard_annotation_all_models wildcard_database_v4.0.0_all.fasta   --wildcard_index ./wildcard/index-for-model-sequences.txt   \
    --wildcard_version 4.0.0  --amr_kmers ./wildcard/all_amr_61mers.txt   --kmer_database ./wildcard/61_kmer_db.json   --kmer_size 61  
    #comparing all files against the databases
    cd ~
    if  [ $v == 1 ]
    then
	   echo "Card rgi database downloaded. Running annotation now"
    fi
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c \
    "rgi main --input_sequence ./$output_dir/merged/${i}.fna --output_file ./$output_dir/card/${i}.tsv --local --clean"
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c \
    "bash ./functions/card_to_gff.sh $output_dir/card/{}.tsv ./$output_dir/card/{}.gff "
    if  [ $v == 1 ]
    then
	   echo "card rgi completed. Results will be available in $output_dir/card"
    fi
}

##Convert file
run_vfdb(){
    if  [ $v == 1 ]
    then
	   echo "Virulence factor database downloading"
    fi
    mkdir -p functions/vfdb $output_dir/vfdb
    cd functions/vfdb
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    gunzip VFDB_setA_pro.fas.gz
    rm VFDB_setA_pro.fas.gz
    makeblastdb -in VFDB_setA_pro.fas -dbtype prot -out VFDB_setA_pro
    cd ~
    if  [ $v == 1 ]
    then
	   echo "VFDB annotation started"
    fi
    ls $input_dir | cut -d "." -f1 | xargs -I {} sh -c \
    "blastp -query ./$output_dir/prodgial/${i}.faa -db ./functions/vfdb/VFDB_setA_pro -out \
    ./$output_dir/vfdb/${i}_results.tsv -evalue 0.001 -outfmt 6" 
    # run awk command which add virulence factor info to the tsv files (matches based on VF number)
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c \
    "bash add_vir_info.sh ./$output_dir/vfdb/${i}_results.tsv ./$output_dir/vfdb/{}.add.tsv"
    # remove extra first line (header)
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c \
    "awk 'NR > 1' ./$output_dir/vfdb/{}.add.tsv > ./$output_dir/vfdb/{}.add_nohead.tsv"
    # run file to convert to gff
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c \
    "bash blastpOut_to_gff.sh ./$output_dir/vfdb/{}.add_nohead.tsv ./$output_dir/vfdb/{}.gff"
    rm  $output_dir/card/*.tsv
    if  [ $v == 1 ]
    then
	   echo "vfdb annotation completed. Results will be available in $output_dir/vfdb"
    fi
}

run_tmhmm(){

    pip3 install pybiolib
    
    mkdir $output_dir/tmhmm
    if  [ $v == 1 ]
    then
	   echo "tmhmm started"
    fi
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c \
    "biolib run DTU/DeepTMHMM --fasta $output_dir/merged/{}.tsv"
    mv *.fna $output_dir/tmhmm/
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c \
    "bash ./functions/tmhmm_to_gff.sh $output_dir/tmhmm/{}.tsv ./$output_dir/tmhmm/{}.gff "
    if  [ $v == 1 ]
    then
	   echo "tmhmm annotation completed. Results will be available in $output_dir/tmhmm"
    fi

}

run_merge_annotations(){
    if  [ $v == 1 ]
    then
	   echo "Merging all the different annotation files"
    fi
    mkdir $output_dir/results
    ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 12 -I {} sh -c "cat $output_dir/eggnog/{}.gff $output_dir/signalp/{}/{}.gff $output_dir/card/{}.gff $output_dir/vfdb/{}.gff $output_dir/tmhmm/{}.gff $output_dir/pilercr/{}.gff > $output_dir/results/{}.gff"
    if  [ $v == 1 ]
    then
	   echo "Merging completes. Results will be available in $output_dir/results"
    fi
}
# run_trial(){
#     if  [ $v == 1 ]
#     then
# 	   echo "starting trial"
#     fi
#     ls $input_dir | cut -d "." -f1 | xargs -n 1 -P 46 -I {} sh -c "cat $input_dir/{}.fna | head -1"

# }
main(){
    parse_args "$@"
    
    if  [ $v == 1 ]
    then
	   echo "Gene Prediction and Annotation starts"
    fi

    mkdir -p "${output_dir}"
    # run_trial

    ##Gene Prediction
    if [ "$tool" == "prodigal" ]
    then
        run_prodigal
    elif [ "$tool" == "glimmer" ]
    then
        run_glimmer
    elif [ "$tool" == "metageneannotator" ]
    then
        run_mga
    elif [ "$tool" == "all" ]
    then
        run_all
    fi
    run_barrnap
    run_aragorn

    run_merge_files

    ##Gene Annotation
    run_eggnog
    run_pilercr
    run_signalp
    run_card
    run_vfdb
    run_tmhmm
    run_merge_annotations
}

main "$@"