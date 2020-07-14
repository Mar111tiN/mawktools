#!/bin/sh

# used in the Coverage10x2 pipeline
# accumulates all coverage info both on total molecules as well as on CB-UB level
# used on individual transcripts providing the length of the transcript

# USAGE:
# covPileup < 10xbam2cov-output translength
# $ samtools view bam_file region | 10xbam2cov TransID + | covPileup10x 10422

mawk '

NR == 1 {
  TL="'$1'";
    # get the coords for the fields
  # for flexible data structures
  n = split("CB,UB,TransID,Start,End,Strand",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  }
  # init pointer
  pos = 0;
  cov = 0;
  UBcov = 0;
  # print Header
  printf("%s\t%s\t%s\t%s\n","5Pos","3Pos","CovMol","CovUB");
}

# look at all hits
NR > 1 {
  CB = $COORD["CB"];
  UB = $COORD["UB"];
  start = $COORD["Start"];
  end = $COORD["End"];
  strand = $COORD["Strand"];

  # with every read, total coverage increases by one between start and end
  # those positions are stored in COVUP and COVDOWN arrays
  COVup[start] ++;
  COVdown[end + 1] ++;

  # UBcollapsed coverage, overlapping reads from one CB-UB are translated into expanded end coords

  # UBUP and UBDOWN are the UB-collapsed coords of coverage increase and decrease
  if (UBend[CB "|" UB]  < start ) { # read has new CB-UB or a new overlap island of the same CB-UB
    if (UBend[CB "|" UB] > 0) { # if overlap (previous CB-UB has an actual end coord)..
      UBdown[UBend[CB "|" UB] + 1] ++;  # ..flush the old end coords into UBdown
    }
    UBup[start] ++; # increment UBup with the new CB-UB coords

  } 
  UBend[CB "|" UB] = end; # if overlap or not, the end coord for that CB-UB will be updated

  ########### OUTPUT ################
  if (pos < start) {
    for (p=pos; p++< start - 1;) { # walk over last coords
      cov = cov + COVup[p] - COVdown[p];
      delete COVup[p];
      delete COVdown[p];
      # accumulate the UBends for that position:
      for (CBUB in UBend) {
        if (UBend[CBUB] == p) {
          UBdown[p + 1] ++;
          delete UBend[CBUB];
        }
      } 
      UBcov = UBcov + UBup[p] - UBdown[p];
      delete UBup[p];
      delete UBdown[p];
      printf("%s\t%s\t%s\t%s\n", p, TL-p+1,cov,UBcov);
    }
    pos=start-1;
    stop=end;
  }
}
END { # spill out all the remaining data arrays
  for (p=pos; p++ < stop;){
    cov = cov + COVup[p] - COVdown[p];
    UBcov = UBcov + UBup[p] - UBdown[p];
    printf("%s\t%s\t%s\t%s\n", p, TL-p+1,cov,UBcov);
  }
}
'