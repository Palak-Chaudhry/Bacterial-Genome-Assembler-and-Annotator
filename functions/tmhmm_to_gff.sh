#!/bin/bash

awk 'BEGIN{OFS="\t"} {
    split($2, len, "="); 
    split($3, expAA, "="); 
    split($4, first60, "="); 
    split($5, predHel, "="); 
    split($6, topology, "="); 
    print $1, "custom_source", "gene", 1, len[2], ".", ".", ".", "ID=" $1 ";Length=" len[2] ";ExpAA=" expAA[2] ";First60=" first60[2] ";PredHel=" predHel[2] ";Topology=" topology[2] 
    }'
$1 > $2
