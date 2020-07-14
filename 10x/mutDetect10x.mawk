#!/bin/sh

# used in the MutDetect2 pipeline to sum up all mutation calls on CB-UB level

# USAGE:
# mutDetect10x < bamPileup_output "chr1:123124-123124:A>T"
# $ samtools view bam_file [region] | bam2csv | bamPileup region | bamPileup_output "chr1:123124-123124:A>T"

mawk '
NR == 1 {

  # get the Ref and Alt values from the CLI params
  param="'$1'";
  paramCount = split(param,PARAMS,">");
  Ref=PARAMS[1];
  Alt=PARAMS[2];
  # get the coords for the fields
  # for flexible data structures
  n = split("SEQ,QUAL,CIGAR,START",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  }

  # get the coords for the fields
  # for flexible data structures
  n = split("CB,UB,QueryBase,QueryQual",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  }

  printf "Summing mutation calls for Ref %s and Alt %s\n", Ref,Alt  > "/dev/stderr";
  # print Header
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","CB","UB","UBcount","Ref","RefQ","Alt","AltQ","Call")
}

# look at all hits
(NR > 1) && ($COORD["QueryBase"] != "-1") && ($COORD["CB"] != ".") {

  CB = $COORD["CB"];
  UB = $COORD["UB"];
  QueryBase = $COORD["QueryBase"];
  QueryQual = $COORD["QueryQual"];

  # build up the arrays
  CBA[CB] ++;
  
  CBUBA[CB "|" UB] ++;
  if (QueryBase == Alt) {
    CBUBAlt[CB "|" UB] ++
    CBUBAltQ[CB "|" UB] = CBUBAltQ[CB "|" UB] QueryQual;
    
  } else {
    CBUBRef[CB "|" UB] ++
    CBUBRefQ[CB "|" UB] = CBUBRefQ[CB "|" UB] QueryQual;
  } 
}

END {
  for (cbub in CBUBA) {
    split(cbub,CBUB,"|")
    CB = CBUB[1];
    UB = CBUB[2];
    UBcount = CBUBA[cbub];
    Ref = (cbub in CBUBRef) ? CBUBRef[cbub] : 0;
    RefQ = (cbub in CBUBRefQ) ? CBUBRefQ[cbub] : ".";
    Alt = (cbub in CBUBAlt) ? CBUBAlt[cbub] : 0;
    AltQ = (cbub in CBUBAltQ) ? CBUBAltQ[cbub] : ".";
    Call = (Ref > Alt) ? 1 : 2;
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",CB,UB,UBcount,Ref,RefQ,Alt,AltQ,Call)
  }
}
 '