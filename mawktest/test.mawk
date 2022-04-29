#!/bin/sh


# >>>filter10xPB<<<
# v0.9
# takes a fastq(gz) file preprocessed with <extractPBAdapters
# and filters for sequences containing 25N[TSO]> / <[TSO]25N

# USAGE: 
# gunzip < fastq.gz | extractPBadapters [options] | filter10xPB




mawk '
BEGIN {
    OFS="\t";
}
$3~ /[ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]\[TSO\]>/ {
    # remove exessive >
    gsub("\]>>+", "]>", $3);
    gsub("\]>>+", "]>", $4);
    
    print($1,$2,$3,$4);
    }
$3~ /<\[TSO\][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG][ACTG]/ {
    # remove exessive >
    gsub("<+<\[", "<[", $3);
    gsub("<+<\[", "<[", $4);
    
    print($1,$2,$3,$4);
    }
'