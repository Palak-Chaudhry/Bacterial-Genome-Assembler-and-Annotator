#! usr/bin/env bash

awk 'BEGIN {OFS="\t"} /^>/ {contig=$1; sub(">", "", contig); next} {
    strang = ($4 < 0) ? "-" : "+"; gsub(/[+-]/,"", $4);
    if ($3 - $2 < 0 ) print contig, "GLIMMER", $1, $3, $2, $5, strang, $4, "ID="$1"; NOTE:GLIMMER ORF prediction;";
    else print contig, "GLIMMER", $1, $2, $3, $5, strang, $4, "ID="$1"; NOTE:GLIMMER ORF prediction;";
}' $1 > $2