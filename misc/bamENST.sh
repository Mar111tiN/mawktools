#!/bin/sh
# $ bamCoverage bamfile ENST APPRISlist


# get the coords of the ENST (and the appropriate Transcript List)

ENST=$2;
APPRIS=$3;
coords=$(cat $APPRIS | awk '$6 == "'$ENST'" {printf("%s:%s-%s\n",substr($1,4),$2,$3)}');
if [[ -n "$coords" ]]; then
  echo "Reading" $coords
  samtools view $1 $coords;
else
  echo $ENST "not found."
fi