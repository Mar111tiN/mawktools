#!/bin/sh

# file with genomic pos in Pos | filterBed <bedfile> [chrom ..which chrom to use] [=0|1..useExonicCoords]
# marks any position-based file as belonging to positions included in a bed file
# takes bedfile and as parameter
# works only on stdin in a pipe

bedFile=$1;
# create random filename
suffix=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 10);
stopfile=stop.${suffix}.txt;
echo "STOPBED" > $stopfile

cat ${bedFile} $stopfile - | mawk '
BEGIN {
    ## INIT ######
    # get filter chrom as arg3 for matching in filterBed file
    # "7" --> chr7
    # default ""  --> chr
    filterChrom="'${2-""}'";
    if (filterChrom !~ "chr") {
        filterChrom = "chr" filterChrom;
    }
    if (filterChrom != "chr"){
        printf("filterBed>> Filtering chromosome %s\n", filterChrom) > "/dev/stderr";
    }
    readBed=1;
    bedCount=0; 
    bedPointer=1;
    bedStep=1;
    chromCount = 1;
}

readBed {
    if ($0 !~ "STOPBED") { # reading bedFile line with right chromosome
        # store the bed regions as blocks in BEDSTART and BEDEND
        if ( $1 !~ filterChrom ) {
            next;
        }
        bedCount++;
        # next chrom in case filterChrom is "chr"
        if ($1 != currentChrom) {
            CHROMEND[currentChrom] = bedCount - 1 # fix the last bedPointer for each chrom
            currentChrom = $1;
            CHROMSTART[currentChrom] = bedCount; # for every newly encountered chromosome, create a marker for the bedPointer
            CHROMNAME[bedCount] = currentChrom; # for that position, remember the chrom name
            CHROMCOUNT[chromCount++] = currentChrom;
        } 
        BEDSTART[bedCount] = $2;
        BEDEND[bedCount] = $3;
        # #####
        # print(bedCount,BEDSTART[bedCount], BEDEND[bedCount], exonicCoord, BEDSUM[bedCount]);
        # #####

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
        printf($0);
            printf("\tonTarget\n");
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
    if (chrom in CHROMSTART == 0) {
        # if not --> offTarget
        print($0, "0");
        next;
    }
    # move to the right chromosome in the read file
    # should best be already done in the samtools file
    while (chrom != currentChrom) {
        currentChrom = CHROMCOUNT[++chromeCount];
        bedPointer = CHROMSTART[currentChrom];
        # reset the chromosome flush if new chromosome has been reached
        chromflush = 0;
    }
    if (chromflush) {
        print($0, "0");
        next;        
    }

    # cycle through the bedRegions as long as they are in current chrom
    while (bedPointer <= CHROMEND[currentChrom]) {
        # get everything upstream of first bedPointer as offTarget
        if (pos < BEDSTART[bedPointer]) {
            print($0, "0");
            next;
        }

        if (pos <= BEDEND[bedPointer]) {
            # within current BEDcoords is onTarget
            print($0, "1");
            next;
        }
        # if pos downstream of current bedRegion, go to next bedRegion
        if (pos > BEDEND[bedPointer]) { # step to next bedRegion
            bedPointer++;
            # if all bedRegions have been read, stop output
            if (bedPointer > bedCount) {
                # stop analysing data
                print($0, "0");
                readData = 0;
                # print out remaining data
                flush = 1;
                next;
            }
        }
    }
    # if bedpointer is in next chrom, flush rest of the data
    print($0, "0");
    chromflush = 1;
}

flush { # after last bedPointer has been reached, flush out the remaining data
    print($0, "0");
}'
# cleanup
rm $stopfile;