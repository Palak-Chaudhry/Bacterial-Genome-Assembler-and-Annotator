#! usr/bin/env bash

awk -F"\t" 'BEGIN {OFS="\t"}
 !/^#/ && NF {print $1, "CARD", "gene", $2, $3, ".", "+", ".", "Name=" $7 ";Resistance=" $4 ";Accession=" $5}' $1 > $2
