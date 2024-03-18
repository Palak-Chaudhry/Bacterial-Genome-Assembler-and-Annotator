#! usr/bin/env bash

input_file=$1
output_file=$2
query=""
while IFS= read -r line; do
    if [[ $line =~ ^\> ]]; then
        query=$( echo "$line" | awk '{print $1}')
    elif [[ $line =~ found ]]; then
        continue
    else
        fields=($line)
        geneid=${fields[0]}
        gene_name=${fields[1]}
        pos=${fields[2]}
        start=$(echo "$pos" | grep -oP '\b\d+\b' | sed -n '1p')
        stop=$(echo "$pos" | grep -oP '\b\d+\b' | sed -n '2p')
        gc=${fields[3]}
        seq=${fields[4]}
        if [[ "$start" -lt "$stop" ]]; then
            echo -e "$query\t$gene_name\ttRNA\t$start\t$stop\t.\t-\t.\tID=$geneid;gc=$gc;seq=$seq;"
        else
            echo -e "$query\t$gene_name\ttRNA\t$stop\t$start\t.\t+\t.\tID=$geneid;gc=$gc;seq=$seq;"
        fi
    fi
done < $input_file > $output_file