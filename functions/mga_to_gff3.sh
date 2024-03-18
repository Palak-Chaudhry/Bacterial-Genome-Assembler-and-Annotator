#! usr/bin/env bash

input_file=$1
output_file=$2

# Initialize variables to store contig and GC/RBS values
current_contig=""
current_gc_rbs=""

# Read the input file line by line
while IFS= read -r line; do
    if [[ $line =~ ^# ]]; then
        # Extract contig information from comment lines
        if [[ $line =~ NODE_ ]]; then
            current_contig=$(echo "$line" | awk '{print $2}')
            echo "$line"
        elif [[ $line =~ "gc =" ]]; then
            current_gc_rbs=$(echo "$line")
        fi
    else
        # Extract gene information from non-comment lines
        fields=($line)
        gene_name=${fields[0]}
        start_pos=${fields[1]}
        end_pos=${fields[2]}
        score=${fields[6]}
        strand=${fields[3]}
        phase=${fields[4]}
        
        # Print the formatted output
        echo -e "$current_contig\tMGA\tCDS\t$start_pos\t$end_pos\t$score\t$strand\t$phase\tID=$gene_name;$current_gc_rbs;"
    fi
done < $input_file > $output_file