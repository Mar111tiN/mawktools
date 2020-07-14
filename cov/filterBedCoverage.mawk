#!/bin/sh

# takes bedfile and chromosome as parameter
# bedfile is expected to be reduced to one chromosome <(cat bedFile | awk '/chrXX/')
# !!! does only work on single chromosome!!!
bedFile=$1;

cat ${bedFile} - | mawk '
BEGIN {
  ## INIT ######
  readBed=1;
  bedCount=0;
  bedPointer=1;
  bedStep=1;
}

readBed {
  if ($3 !~ "Coverage") { # reading bedFile line with right chromosome
    # store the bed regions as blocks in BEDSTART and BEDEND
    bedCount++;
    BEDSTART[bedCount] = $2;
    BEDEND[bedCount] = $3;
  } else { # reached header of coverage data
    ########## HEADER #######################
    # write HEADER and switch to coverage reads

    # get the coords for the fields
    # for flexible data structures
    n = split("Chr,Pos,Coverage",FIELDS,",");
    for (col in FIELDS) {
      for (i = 0; i++ < NF;) {
        if ($i == FIELDS[col]) {
          COORD[FIELDS[col]]=i;
        }
      }
    }
    # print Header
    printf("%s\t%s\t%s\t%s\n","#","Chr","Pos","Coverage"); 
    # switch to coverage mode
    readBed = 0;
    readCov = 1;
  }
  next;
}

readCov {  # switching to coverage
  # get data
  chrom=$COORD["Chr"];
  pos=$COORD["Pos"];
  cov=$COORD["Coverage"];

  while (pos >= BEDSTART[bedPointer] ) { # if pos below current bedRegion, drop line silently
    #print(bedPointer, BEDSTART[bedPointer], BEDEND[bedPointer], pos);
    if (bedPointer > bedCount) { # if all bedRegions have been read, stop output
      exit 0;
    }
    if (pos > BEDEND[bedPointer]) {
      bedPointer++;
    } else { # step to next bedRegion
      printf("%s\t%s\t%s\t%s\n", bedStep++, chrom, pos, cov);
      break; 
    }
  }
}'