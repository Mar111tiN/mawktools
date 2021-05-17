#!/bin/sh

# converts a bedfile into an index file containing all positions


####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--use-exon-pos)
        useExonicCoords=1
        shift
        ;;
        # snp_file output
        -c|--chrom)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            filterChrom=$2
            shift 2
        else
            echo "<bed2index> Error: chromosome argument is missing\n[-c|--chrom (default=all)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<bed2index> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done

# I have no positional args
# # set positional arguments in their proper place
eval set -- "$PARAMS"

useExonicCoords=${useExonicCoords-0};
filterChrom=${filterChrom-""};
bedFile=$1;

# echo "bedFile:" $bedFile;
# echo "useExonicCoords:" $useExonicCoords;
# echo "filterChrom:" $filterChrom;

cat $bedFile | mawk '
BEGIN {
    ## INIT ######
    # useExonicCoords
    useExonicCoords='$useExonicCoords';
    if (useExonicCoords == 1) {
        printf("<bed2index> Printing out exonic coordinates\n") > "/dev/stderr";
    }
    # get filter chrom as arg3 for matching in bed2index file
    # "7" --> chr7
    # default ""  --> chr

    printf("<bed2index> Reading bed file %s\n", "'$bedFile'") > "/dev/stderr";
    filterChrom="'$filterChrom'";
    # handle the filterChrom
    # convert integer chrom to chr-chrom
    if (filterChrom !~ "chr") {
        filterChrom = "chr" filterChrom;
    }
    
    if (filterChrom != "chr"){ # no filterChrom given --> filterChrom == "chr"
        # adjust the chrom-pattern for exact chrom match
        filterChrom = "^" filterChrom "\t";
        printf("<bed2index> Filtering chromosome %s\n", filterChrom) > "/dev/stderr";
    }
    bedCount=0; 
    bedPointer=1;
    bedStep=1;
    chromCount = 1;
}

$0 ~ /^chr/ {
    # print(filterChrom);
    if  ($0 !~ "Chr\t") { # reading bedFile line with right chromosome
        # store the bed regions as blocks in BEDSTART and BEDEND
        if ( $0 !~ filterChrom ) {
            next;
        }
        bedCount++;
        # next chrom in case filterChrom is "chr"
        if ($1 != currentChrom) {
            exonicCoord=0; # the current accumulated exon base count 
            CHROMEND[currentChrom] = bedCount - 1 # fix the last bedPointer for each chrom
            currentChrom = $1;
            CHROMSTART[currentChrom] = bedCount; # for every newly encountered chromosome, create a marker for the bedPointer
            CHROMNAME[bedCount] = currentChrom; # for that position, remember the chrom name
            CHROMCOUNT[chromCount++] = currentChrom;
        }
        BEDSTART[bedCount] = $2;
        BEDEND[bedCount] = $3;
        # exonicCoord will be increased by the length of the bed fragment
        # BEDSUM marks the exonic coords of the bed fragments end (defined by bedCount)
        exonicCoord = exonicCoord + $3 - $2 +1;
        # #####
        # print(bedCount,BEDSTART[bedCount], BEDEND[bedCount], exonicCoord, BEDSUM[bedCount]);
        # #####

        BEDSUM[bedCount] = exonicCoord;
        next;
    } else { # reached end of bedfile
            # switch to HEADER mode
            readBed = 0;
            writeHeader = 1;
    }
}

############# OUTPUT #########################
END {  
    CHROMEND[currentChrom] = bedCount
    chromCount = 1;
    currentChrom = CHROMNAME[1];
    # write HEADER
    printf("Chr\tPos")
    if (useExonicCoords == 1) {
        printf("\tExonPos");
    }
    printf("\n");


    for (i=0;i++<length(BEDSTART);) {
        if (i in CHROMNAME) {
            chrom=CHROMNAME[i]
        }
        start=BEDSTART[i];
        end=BEDEND[i];
        exonshift=BEDSUM[i]-end+start;
        for (p=start-1;p++<end;) {
            printf("%s\t%s",chrom,p);
            if (useExonicCoords == 1) {
                exonicCoord=p - start + exonshift;
                printf("\t%s", exonicCoord);
            }
            printf("\n");
        } 
    }
}'
