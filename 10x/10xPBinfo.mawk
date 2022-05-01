#!/bin/sh

# >>>10xPBinfo<<<
# v0.9.1 fixed bug with single bases between tags
# takes a fastq(gz) file preprocessed with the 10xPB toolchain component:
#     <10xPBextract>   |      transforms into one-line data and extracts the relevant oligo adapter sequences
# and adds a long and a short info field discribing the sequence structure with regard to known adapters

# USAGE: gunzip < fastq.gz | 10xPBextract [options] | 10xPBinfo
mawk '
############# BEGIN #########################
BEGIN {  ### GET/WRITE HEADER
    # set field separators
    FS="\t";
    OFS="\t";
}
############ READ LINES #####################
{   
    seq=$2;
    qual=$3;
    ### retrieve the short string
    info=seq;
    gsub(">+",">",info);
    gsub("<+","<", info);
    # turn base stretches to N
    # starting bases
    if (match(info, "^[ACTG]+")) {
        
        info = RLENGTH "N" substr(info, RSTART+RLENGTH);
    }
    # middle bases L
    while (match(info, ">[ACTG]+")) {
        info = substr(info,1,RSTART) RLENGTH-1 "N" substr(info, RSTART+RLENGTH);
    }
    # middle bases R
    while (match(info, "[ACTG]+<")) {
        info = substr(info,1,RSTART-1) RLENGTH-1 "N" substr(info, RSTART+RLENGTH-1);
    }
    # end bases
    if (match(info, "[ACTG]+$")) {
        info = substr(info,1,RSTART-1) RLENGTH "N";
    }
    # remaining single bases like [Tag]>C[Tag]


    # make short info
    shortInfo = info;
    while (match(shortInfo, "[0-9]+N")) {
        gsub("[0-9]+N","N", shortInfo);
    }
    print($1,info,shortInfo,seq,qual);
}
'