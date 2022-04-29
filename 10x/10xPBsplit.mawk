#!/bin/sh


# >>>10xPBsplit<<<
# v0.9
# takes a fastq(gz) file preprocessed with <10xPBfilter>
# and splits sequences of facing 10x structures { N[TSO]>N<[TSO]N } into separate fastq

# USAGE: 
# gunzip < fastq.gz |  10xPBextract [options] | 10xPBfilter -i | 10xPBsplit
# [     -m | --minSize             <INT=40>                     minimum size for middle sequence               ]
# [      O | --Overlap             <INT=74>                     sequence overlap between left and right split  ]

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # minSize
        -m|--minSize)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minSize=$2
            shift 2
        else
            echo "<10xPBsplit> Error: Param minSize is missing\n[-m|--minSize]" >&2
            exit 1
        fi
        ;;
        # Overlap
        -O|--Overlap)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            Overlap=$2
            shift 2
        else
            echo "<10xPBsplit> Error: Param Overlap is missing\n[-O|--Overlap]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<10xPBsplit> Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
done




mawk '
BEGIN {
    OFS="\t";
    minSize='${minSize-40}';
    Overlap='${Overlap-37}';
    # restrict maximum overlap
    if (Overlap > 100) {
        Overlap =100;
    }
    perc="%";
    printf("<10xPBsplit> Splitting overlapping sequences [minSize=%s|Overlap=%s%s]\n", minSize, Overlap, perc) >> "/dev/stderr";
}

$3~ /N\[TSO\]>N<\[TSO\]N/ {
    info=$2;
    shortInfo=$3;
    seq=$4;
    qual=$5;

    # get the center sequence and its size
    # check for minSize

    if (match(seq, /TSO\]>[ACTG]+<\[TSO/)) {
        # extract the center length
        # center=substr(seq,RSTART+5,RLENGTH-10);
        centerLength=RLENGTH-10;
        if (centerLength>minSize-1) {

            # split the sequences and qualites
            overlapLength = int(centerLength/2*(1+Overlap/100));
            seqL=substr(seq,1,RSTART+4+overlapLength);
            qualL=substr(qual,1,RSTART+4+overlapLength);

            seqR=substr(seq,RSTART+RLENGTH-5-overlapLength);
            qualR=substr(qual,RSTART+RLENGTH-5-overlapLength);

            # split the info fields
            match(info, /TSO\]>[0-9]+N<\[TSO/);
            infoL=substr(info,1,RSTART-1) "TSO]>" overlapLength "N"
            infoR= overlapLength "N<[TSO" substr(info,RSTART+RLENGTH) 

            # split the info fields
            match(shortInfo, /TSO\]>N<\[TSO/);
            shortInfoL=substr(shortInfo,1,RSTART-1) "TSO]>N"
            shortInfoR= "N<[TSO" substr(shortInfo,RSTART+RLENGTH) 
        }
        ####### OUTPUT #########
        ## printL
        print($1 "L", infoL, shortInfoL, seqL, qualL);
        print($1 "R", infoR, shortInfoR, seqR, qualR);
    }
    next;
}
{
    # for all other structures, just pass through the data
    print($0);
}
'