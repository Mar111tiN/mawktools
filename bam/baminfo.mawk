#!/bin/sh

# adds two columns to a bam csv
# NCount: the number of Ns in a read sequence
# Qcount: the accumulated Quality Count for the read

mawk '
NR==1 { # HEADER column
  printf("%s\t%s\t%s\t%s\n", $0,"NCount","Qsum","Qmean");
  # make lookup table for phred qualities
  for (i=33; i++<126;) {
    char=sprintf("%c", i-1);
    asc[char]=i - 34;
  }
  # get the coords for the SEQ and QUAL field
  n = split("SEQ,QUAL",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  } 
}
##### PROCESS LINES #################
NR>1 {
  seq=$COORD["SEQ"];
  seqL=length(seq);
  qual=$COORD["QUAL"];
  countN=gsub("N", "", seq);
  # calculate the sum of base qualities
  qualcount=0;
  for (i=0; i++ < seqL;) {
    q=asc[substr(qual,i,1)];
    qualcount+=q;
  }
  qualmean = int((qualcount * 10 / seqL) + 0.5) / 10;
  printf("%s\t%s\t%s\t%s\n", $0, countN, qualcount,qualmean);
}'