#!/bin/sh

# used in the MutDetect2 pipeline to collapse all mutation calls from CB-UB level to CB level

# USAGE:
# collapseUB < mutDetect10x_output
# $ samtools view bam_file [region] | bam2csv | bamPileup region | bamPileup_output "chr1:123124-123124:A>T" | collapseUB


mawk '
NR == 1 {
  param="'$1'";
  # make lookup table for phred qualities
  for (i=33; i++<126;) {
    char=sprintf("%c", i-1);
    asc[char]=i - 34;
  }
    # get the coords for the fields
  # for flexible data structures
  n = split("CB,UB,Ref,RefQ,Alt,AltQ,Call",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  }
  # print Header
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n","CB","CBcount","Call","Ref","RefQ","Alt","AltQ");
}

# look at all hits
NR > 1 {
  q=asc[substr(qual,i,1)]-34;
  CB = $COORD["CB"];
  UB = $COORD["UB"];
  Ref = $COORD["Ref"];;
  RefQs = $COORD["RefQ"];;
  Alt = $COORD["Alt"];;
  AltQs = $COORD["AltQ"];;
  Call = $COORD["Call"];;

# get average quality for RefQ and AltQ
############# RefQ ######################
  RefQsum = 0;
  if (RefQs == ".") {
    RefQ = 0;
    RefP = 1;
  } else {
    for (i=0; i++< Ref;) {
      Q = substr(RefQs,i,1);
      RefQsum += asc[Q];
    }
    RefQ = RefQsum / Ref;
    # divide the P-value of Error with number of supporting reads
    # that should be rough estimate for read support
    RefP = exp((-RefQ/10)*log(10)) / Ref; 
  }
  ########## AltQ ##################
  AltQsum = 0;
  if (AltQs == ".") {
    AltQ = 0;
    AltP = 1;
  } else {
    for (i=0; i++< Alt;) {
      Q = substr(AltQs,i,1)
      AltQSum += asc[Q];
    }
    AltQ = AltQSum / Alt;
    AltP = exp((-AltQ/10)*log(10)) / Alt;  
  }
  # case Ref Call
  if (Call == 1) {
    if (AltP == 1) {
      CallP = 60;
    } else {
      CallP = int(-10 * log(RefP * (1- AltP)) / log(10) + 0.5);
    }
    
  # case Alt Call
  } else {
    if (RefP = 1) {
      CallP = 60;
    } else {
      CallP = int(-10 * log(AltP * (1 - RefP)) / log(10) + 0.5);
    }
  }
  
  # CBcount (how many UB per CB)
  CBA[CB] ++;
  if (Call == 1) {
    CBRef[CB] ++;
    CBRefQ[CB] = CBRefQ[CB] ":" CallP;
  # Alt Call
  } else {
    CBAlt[CB] ++;
    CBAltQ[CB] = CBAltQ[CB] ":" CallP;
  }
}

END {
  for (CB in CBA) {
    CBcount = CBA[CB];
    Ref = (CB in CBRef) ? CBRef[CB] : 0;
    RefQ = (CB in CBRefQ) ? substr(CBRefQ[CB],2) : ".";
    Alt = (CB in CBAlt) ? CBAlt[CB] : 0;
    AltQ = (CB in CBAltQ) ? CBAltQ[CB] : ".";
  
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",CB,CBcount,Call,Ref,RefQ,Alt,AltQ);
  }
}
 '