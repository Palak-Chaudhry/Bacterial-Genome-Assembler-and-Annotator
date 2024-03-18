#! usr/bin/env bash

# gff: 1(seqname) 2(source) 3(feature type (gene,variation)) 4(start) 5(end) 6(score) 7(strand) 8(frame-indicated that the first base of feature is first base of codon) 9(attr-additional info)
# blastp: 1(query seqid) 2(subject seqid) 3(% identical matches) 4(align length) 5(mismatches) 6(gap openings) 7(query start) 8(qend) 9(subject align start) 10(subject align end) 11(e-val) 12(bitscore)
awk 'BEGIN {FS="\t"; OFS="\t"} {
    print $2, "BlastP", "alignment", $7, $8, $11, ".", ".", "ID=" $1 ";Percent_identity=" $3 ";Alignment_length=" $4 ";Mismatches=" $5 ";Bit score=" $12 ";Virulence_info=" $13
}' $1 > $2

