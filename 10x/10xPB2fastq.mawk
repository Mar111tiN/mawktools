#!/bin/sh


# >>>10xPBfilter<<<
# v0.9
# v0.9.1 -- fix bug causing non-similar fastq output
# --> v.0.9.2 --> make CB+UMI length flexible (now hardcoded as 26)
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
# [     -s | --min_seq_length      <Int=25>                     minimum length of sequence after TSO              ]
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
        # maxMismatch
        -s|--min_seq_length)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minSeqLength=$2
            shift 2
        else
            echo "<10xPB2fastq> Error: Value for param min_seq_length is missing\n[-s|--min_seq_length]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<10xPB2fastq> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        # maxMismatch
        -m|--max_mismatch)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            maxMismatch=$2
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
    minSeqLength='${minSeqLength-25}';
    maxMismatch='${maxMismatch-10}';
    outPrefix="'${outPrefix-""}'";
    if (outPrefix == "") {
        printf("<10xPB2fastq> No output prefix provided. Using \"PBextracted_R1|2.fastq\"!\n") >> "/dev/stderr";
        outPrefix = "PBextracted";
    }
    # set the file
    read1File=outPrefix "_R1.fastq";
    read2File=outPrefix "_R2.fastq";
    printf("<10xPB2fastq> Splitting 10xPB seq into %s (CB+UMI) and _R2.fastq\n", read1File) >> "/dev/stderr";
    misMatchCount=0;
}
############# READ MATCHING LINES #########################
$2 ~ /N\[TSO\]>/ {
    info=$2;
    seq=$3;
    qual=$4;
    # extract the lengths

    match(info, /[0-9]+N\[TSO\]>/);
    xLength = substr(info,1,RLENGTH-7);
    match(info, /\[TSO\]>[0-9]+N/);
    seqLength = substr(info, RSTART+6,RLENGTH-7);

    # check for minSeqLength
    if (int(seqLength)<int(minSeqLength)) {
        next;
    }

    # add one random base for xLength=25
    if (xLength == 25) {
        seq= "G" seq;
        qual="X" qual;
        xLength = 26;
    }
    
    # extract the xSeq
    xSeq = substr(seq, xLength-25,26);
    xQual = substr(qual, xLength-25,26);
    if (xLength>25) {
        printf("@%s\n%s\n+\n%s\n", $1,xSeq, xQual) >> read1File;
        # extract the dataSeq
        dataStart = xLength + 7;
        dataSeq = substr(seq, dataStart);
        dataQual = substr(qual, dataStart);
        printf("@%s\n%s\n+\n%s\n", $1,dataSeq, dataQual) >> read2File;
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
            printf("Read%s\t[%s]\n", $1,$2) >> "/dev/stderr";
        }
    }
}
############# END  #########################
END { # OUTPUT NON-MATCHING LINES
    if (misMatchCount > 0) {
        printf("<10xPB2fastq> Detected very few (%s) read with un-cleaned data (please check reads!)\n", misMatchCount) >> "/dev/stderr";
        for (i=0;i++<misMatchCount;) {
            printf("@%s\t[%s]\n", MMREAD[i],MMINFO[i]) >> "/dev/stderr";
        }
    }

}'