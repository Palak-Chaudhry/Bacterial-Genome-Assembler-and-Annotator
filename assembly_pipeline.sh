#! usr/env/bin bash

## To create list of raw read file names

parse_args() {

    # Function to parse arguments
    # Specifying usage message
    # -c : Compare tools (set -t to all for this mode)
    usage="Usage: bash assembly_pipeline.sh -i <input directory> -o <output directory> -[OPTIONS]
              Bacterial gene assembly for Illumina short-reads. The options available are:
                        -i : Input Directory for genome assembly (directory containing raw reads in fastq.gz format)[required]
                        -o : Output directory for all results [required]
                        -v : Flag to turn on verbose mode
                        -h : Print usage instructions"

  #Specifying default Arguments
  v=0

  #Getopts block, will take in the arguments as inputs and assign them to variables
  while getopts "i:o:t:vh" option; do
          case $option in
                  i) input_dir=$OPTARG;;
                  o) output_dir=$OPTARG;;
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
}

main(){
    parse_args "$@"
    python3 getlist.py

    if  [ $v == 1 ]
    then
	   echo "Starting quality check"
    fi
    bash ./functions/qc.sh pre $input_dir $output_dir

    if  [ $v == 1 ]
    then
	   echo "Quality check complete results will be found in outputdir/preqc, starting phix removal"
    fi
    
    bash ./functions/phixremoval.sh $input_dir $output_dir

    if  [ $v == 1 ]
    then
	   echo "Phix successfully removed. Trimming starting"
    fi
    
    bash ./functions/trim.sh $input_dir $output_dir
    
    if  [ $v == 1 ]
    then
	   echo "Trimming complete, checking quality"
    fi
    
    bash ./functions/qc.sh post $input_dir $output_dir
    
    if  [ $v == 1 ]
    then
	   echo "Quality check complete, results can be found at outputdir/postqc. Starting SPAdes assembly. Might take a while"
    fi
    
    bash ./functions/spades.sh $input_dir $output_dir
    
    if  [ $v == 1 ]
    then
	   echo "SPAdes assembly complete, filtering high confidence contigs"
    fi
    
    bash ./functions/filterconfig.sh $input_dir $output_dir
    
    if  [ $v == 1 ]
    then
	   echo "Filtering complete. Final results available at outputdir/filtered_spades. You can choose your desired filtered reads"
    fi
}

main "$@"