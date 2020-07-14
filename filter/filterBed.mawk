#!/bin/sh

# filters any position-based file to positions included in a bed file
# takes bedfile and as parameter
# works only on stdin in a pipe
bedFile=$1;
useCounter=${2-0};

echo "STOPBED" > stop.file

cat ${bedFile} stop.file - | mawk '
BEGIN {
  ## INIT ######
  readBed=1;
  bedCount=0;
  bedPointer=1;
  bedStep=1;
  chromCount = 1;
}

readBed {
  if ($0 !~ "STOPBED") { # reading bedFile line with right chromosome
    # store the bed regions as blocks in BEDSTART and BEDEND
    if ( $1 !~ /chr/ ) {
      next;
    }
    bedCount++;
    BEDSTART[bedCount] = $2;
    BEDEND[bedCount] = $3;
    # next chrom
    if ($1 != currentChrom) {
      CHROMEND[currentChrom] = bedCount - 1 # fix the last bedPointer for each chrom
      currentChrom = $1;
      CHROMSTART[currentChrom] = bedCount; # for every newly encountered chromosome, create a marker for the bedPointer
      CHROMNAME[bedCount] = currentChrom; # for that position, remember the chrom name
      CHROMCOUNT[chromCount++] = currentChrom;
    }
    
  } else { # reached end of bedfile
    # switch to HEADER mode
    readBed = 0;
    writeHeader = 1;
  }
  next;
}

########## HEADER #######################
writeHeader { # check for existence of header and print out
  CHROMEND[currentChrom] = bedCount
  writeHeader = 0;
  readData = 1;
  chromCount = 1;
  currentChrom = CHROMNAME[1]; # reset the current chrom to the first one

  # for (cc in CHROMSTART) {
  #   print(cc, CHROMSTART[cc]);
  # }
  if ($0 ~ "Chr\t") { 
    print($0);
    next;
  } else { # move on to readData on same line
    print("<filterBedAll> No header detected") > "/dev/stderr";
  }
}

############# DATA #########################
readData {  # switching to data
  # get data
  chrom=$1;
  pos=$2;
  
  # check if chrom is contained in bedfile
  if (chrom in CHROMSTART == 0) next;
  # move to the right chromosome in the read file
  while (chrom != currentChrom) {
    currentChrom = CHROMCOUNT[++chromeCount];
    bedPointer = CHROMSTART[currentChrom];
  }
  
  # cycle through the bedRegions
  while (pos >= BEDSTART[bedPointer] ) { # if pos downstream of current bedRegion, drop line silently
    #print(bedPointer, BEDSTART[bedPointer], BEDEND[bedPointer], pos);
    if (pos > BEDEND[bedPointer]) { # step to next bedRegion
      bedPointer++;
      # if all bedRegions have been read, stop output
      if (bedPointer > bedCount) exit 0;
      # if bedpointer is in next chrom, skip
      if (bedPointer > CHROMEND[currentChrom]) next;
    } else { # write to file
      if (useCounter) { 
        printf("%s\t%s\n", bedStep++, $0);
      } else {
        print($0);
      }  
      break; 
    }
  }
}'

# cleanup
rm stop.file;