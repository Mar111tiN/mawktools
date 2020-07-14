#!/bin/sh

# 3 CORES

### cleans samtools mpileup output 
mawk '
NR == 1 {
    samples = (NF - 3) / 3;
}
{
  # cycle through the reads
  for (i = 0; i++ < samples;) {
    col = (i*3) + 2;
    # remove position traces from all $col fields
    gsub(/\^.|\$/,"",$col);

    # remove the indel traces
    while (match($col,/[+-][0-9]+/)) {
      pre = substr($col,1,RSTART-2);
      indel = substr($col,RSTART,1); # + or -
      base = substr($col,RSTART-1,1); # base before the deletion
      l = substr($col,RSTART+1,RLENGTH-1);
      post = substr($col,RSTART+RLENGTH+l);
      if (indel == "-") {
        if (match(base,/[ACGT.]/)) {
          base = "D";
        } else {
          base = "d";
        }
      } else {
        if (match(base,/[ACGT.]/)) {
          base = "I";
        } else {
          base = "i";
        }            
      }
      $col = pre base post;
    } 
  }
  printf("%s",$1);
  for (i=1; i++<3;) {
    printf("\t%s",$i);
  }
  for (s=0;s++<samples;) {
    depth=(s*3)+1;
    read=(s*3)+2;
    printf("\t%s\t%s",$depth,$read);
  }
  printf("\n");
}'  | mawk ' # cat - <(echo "STOP\n") $SNP |   ###### if the SNP markers are necessary, they have to be read in here 
## $5 ~ /[.,]*[ACTG]+[.,]*/ | mawk    # if you only want to output mutant positions (complete the quotes!!)
NR == 1 { 
  samples = (NF - 3) / 2;
  ###### QUERY ##############
  # get the letters to look for
  len = split("Aa,Gg,Cc,Tt,Dd,Ii",Letters,",")
  ###### HEADER ################
  for (l=0;l++<len;) {
    LETPAT[l] = "[" Letters[l] "]"
  }
  printf("%s\t%s\t%s", "Chr","Start","Ref");
  ####### HAVE TO ADJUST ###########
  # loop through the letters
  for (s=0;s++<samples;) {
      printf("\t%s\t%s", "Depth" s, "Alt" s);
  }
  printf("\n");
}
{  ######### LINES #############
  printf("%s\t%s\t%s",$1,$2,$3);
  
  # loop through the reads
  for (s=0;s++<samples;) {
    col=3+(2*s);
    # loop through the letters
    for (l = 0;l++< len;) {
        COUNT[s "-" l] = gsub(LETPAT[l], "", $col);
    }
  }
  ######### OUTPUT #############
  # loop through the letters
  for (s=0;s++<samples;) {
    maxcount = 0;
    base = "";
    depthcol=2+(2*s);
    printf("\t%s\t",$depthcol); 
    # first line extra for pretty
    for (l=0;l++<len;) {
      if (COUNT[s "-" l] > maxcount) {
        maxcount = COUNT[s "-" l];
        base = substr(Letters[l],1,1);
      } 
    } 
    printf("%s%s", base,maxcount);
  }
  printf("\n");
}'