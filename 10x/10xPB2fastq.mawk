#!/bin/sh


# >>>10xPBfilter<<<
# v0.9
# takes a fastq(gz) file preprocessed with the 10xPB toolchain components:
#     <10xPB2fastq>    |      transforms into one-line data and extracts the relevant oligo adapter sequences
# --> <10xPBinfo>      |      adds a long and a short info field discribing the sequence structure with regard to known adapters
# --> <10xPBfilter>    |      filters for sequences containing 25N[TSO]> / <[TSO]25N
# --> <10xPBsplit>     |      splits sequences with facing 10xSignature components into two separate lines
# --> <10xPBclean>     |      isolates 10x data { N[TSO]>N or N<[TSO]N } into clean data lines
# and filters for sequences starting with proper CB-UMI-TSO stretch and splits this into separate files

# USAGE: 
# gunzip < fastq.gz |  10xPBextract [options] | 10xPBinfo | 10xPBfilter -i | 10xPBsplit | 10xPBclean | 10xPB2fastq
# [     -o | --output_prefix       <Path to file>               output prefix for _R1.fastq and _R2.fastq file    ]
# [     -m | --max_mismatch        <Int=10>                     exit after N non-conform structures               ]

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # oligoFile
        -o|--outPrefix)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            outPrefix=$2
            shift 2
        else
            echo "<10xPB2fastq> Error: Path to outPrefix is missing\n[-o|--out_prefix]" >&2
            exit 1
        fi
        ;;
        # number of hits
        -h|--Nhits)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            Nhits=$2
            shift 2
        else
            echo "<10xPB2fastq> Error: Value for param max_mismatch is missing\n[-m|--max_mismatch]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<10xPB2fastq> Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
done


mawk '
############# BEGIN #########################
BEGIN {
    OFS="\t";

    # params
    maxMismatch='${maxMismatch-10}';
    outPrefix="'${outPrefix-""}'";
    if (outPrefix == "") {
        printf("<10xPB2fastq> No output prefix provided. Using \"PBextracted_R1|2.fastq\"!\n") >> "/dev/stderr";
        outPrefix = "PBextracted";
    }
    # set the file
    read1File=outPrefix "_R1.fastq";
    read2File=outPrefix "_R2.fastq";
    printf("<10xPB2fastq> Splitting 10xPB seq into %s (CB+UMI) and _2.fastq\n", read1File) >> "/dev/stderr";
    misMatchCount=0;
}
############# READ MATCHING LINES #########################
$2 ~ /N\[TSO\]>/ {
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
############# EXIT @NON-MATCHING LINES #########################
$2 !~ /N\[TSO\]>/ {
    misMatchCount++;
    MMREAD[misMatchCount]=$1;
    MMINFO[misMatchCount]=$2;
    if (misMatchCount>maxMismatch) {
        print("<10xPB2fastq> Detected read with un-cleaned data (please process your data with 10xPBclean tool!)\n") >> "/dev/stderr";
        for (i=0;i++<misMatchCount;) {
            printf("Read%s\t[%s]\n", $1,$2)
        }
    }
}
############# END  #########################
END { # OUTPUT NON-MATCHING LINES
    printf("<10xPB2fastq> Detected very few (%s) read with un-cleaned data (please check reads!)\n", misMatchCount) >> "/dev/stderr";
    for (i=0;i++<misMatchCount;) {
        printf("@%s\t[%s]\n", MMREAD[i],MMINFO[i]);
    }
}
'