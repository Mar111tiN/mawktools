#!/bin/sh

# >>>10xPB2fastq<<<
# v0.9
# takes a fastq(gz) file preprocessed with the 10xPB toolchain components:
#     <10xPBextract>   |      transforms into one-line data and extracts the relevant oligo adapter sequences
# --> <10xPBinfo>      |      adds a long and a short info field discribing the sequence structure with regard to known adapters
# and filters for sequences containing 25N[TSO]> / <[TSO]25N

# USAGE: 
# gunzip < fastq.gz | 10xPBextract [options] | 10xPBinfo | 10xPBfilter
# [     -i | --fullInfo            <Flag=False>                 set if full structure should be output    ]
# [     -C | --CBlength            <INT=16>                     set length of cellular barcode            ]
# [     -U | --UMIlength           <INT=10>                     set length of UMI                         ]
# [     -d | --del                 <INT=1>                      how many deletions at most in CB+UMI      ]

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # mismatch
        -i|--fullInfo)
        fullInfo=1
        shift
        ;;
        # CBlength
        -C|--CBlength)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            CBlength=$2
            shift 2
        else
            echo "<10xPBfilter> Error: Param CBlength is missing\n[-C|--CBlength]" >&2
            exit 1
        fi
        ;;
        # UMIlength
        -U|--UMIlength)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            UMIlength=$2
            shift 2
        else
            echo "<10xPBfilter> Error: Param UMIlength is missing\n[-U|--UMIlength]" >&2
            exit 1
        fi
        ;;
        # deletions
        -d|--del)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            Del=$2
            shift 2
        else
            echo "<10xPBfilter> Error: Param del is missing\n[-d|--del]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<10xPBfilter> Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
done






mawk '
############# BEGIN #########################
BEGIN {
    OFS="\t";

    # params
    fullInfo='${fullInfo=0}';
    CBlength='${CBlength-16}';
    UMIlength='${UMIlength-10}';
    Del='${Del-1}';
    LMAX=CBlength + UMIlength - Del;
    printf("<l0xPBfilter> Filtering for 10x signature [CBlength=%s|UMIlength=%s|fullInfo=%s]\n", CBlength, UMIlength, fullInfo) >> "/dev/stderr";
}

############# READ LINES #########################
{
    has10xSignature=0;
}

############# READ MATCHING LINES #########################
$3~ /N\[TSO\]>/ {
    info=$2;
    # check for N=>25
    lmax=0
    while (match(info, /[1-9][0-9]+N\[TSO\]>/)) {
        l = substr(info, RSTART,RLENGTH-7);
        info = substr(info, RSTART+RLENGTH) 
        if (l > lmax) {
            lmax=l;
        }
    }
    if (lmax>=LMAX) {
        # remove excessive >
        gsub("\]>>+", "]>", $4);
        gsub("\]>>+", "]>", $5);
        has10xSignature=1;
    }
}

############# READ MATCHING LINES #########################
$3~ /<\[TSO\]N/ {
        info=$2;
    # check for N=>25
    lmax=0
    while (match(info, /<\[TSO\][1-9][0-9]+N/)) {
        l = substr(info, RSTART+6,RLENGTH-7);
        info = substr(info, RSTART+RLENGTH) 
        if (l > lmax) {
            lmax=l;
        }
    }
    if (lmax>=LMAX) {
        # remove excessive <
        gsub("<+<\[", "<[", $4);
        gsub("<+<\[", "<[", $5);
        has10xSignature=1; 
    }
}

############# READ LINES @CONDITION #########################
has10xSignature {
    if (fullInfo==1) {
            print($1,$2,$3,$4,$5);
        } else {
            print($1,$3,$4,$5);
        }
}
'