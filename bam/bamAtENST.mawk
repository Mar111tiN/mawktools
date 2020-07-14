#!/bin/sh

# filter out the bam reads covering a given transcript provided as ENST
# USAGE:
# $ bamAtENST bam_file TransID APPRIS_list

ENST=$2;
APPRIS=$3;

samtools view $1 $(cat $APPRIS | awk '$6 == "'$ENST'" {printf("%s:%s-%s\n",$1,$2,$3)}')