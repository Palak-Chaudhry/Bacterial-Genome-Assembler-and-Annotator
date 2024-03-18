#! usr/bin/env bash

awk 'BEGIN {OFS="\t"} /^#/ {next} 
     {if($3 == "CDS") print $1, $2, $3, $4, $5, $6, $7, $8, $9$10$11$12$13;}' $1 > $2
