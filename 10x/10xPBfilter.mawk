#!/bin/sh


# >>>10xPBfilter<<<
# v0.9.1
# takes a fastq(gz) file preprocessed with <extractPBAdapters
# and filters for sequences containing 25N[TSO]> / <[TSO]25N

# USAGE: 
# gunzip < fastq.gz | extractPBadapters [options] | filter10xPB
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
            shilt 2
        else
            echo "<10xPBfilter> Error: Param CBlength is missing\n[-C|--CBlength]" >&2
            exlt 1
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
        echo "<10xPBextract> Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
done


mawk '
BEGIN {
    OFS="\t";

    # params
    fullInfo='${fullInfo=0}';
    CBlength='${CBlength-16}';
    UMILlngth='${UMIlength-10}';
    Del='${Del-1}';
    LMAX=CBlength + UMIlength - Del;
    printf("<l0xPBfilter> Filtering for 10x signature [CBlength=%s|UMIlength=%s|fullInfo=%s]\n", CBlength, UMIlength, fullInfo) >> "/dev/stderr";
}
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
        if (fullInfo==1) {
            print($1,$2,$3,$4);
        } else {
            print($1,$3,$4);
        }
        
    }

}
$3~ /<\[TSO\]N/ {
    # remove excessive >
    gsub("<+<\[", "<[", $3);
    gsub("<+<\[", "<[", $4);
    
    print($1,$2,$3,$4);
    }
'