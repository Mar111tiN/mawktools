#!/bin/sh

# >>>10xPBinfo<<<
# v0.9
# takes a fastq(gz) file preprocessed with <extractPBAdapters
# and retrieves an info about the general structure

# USAGE: 
# gunzip < fastq.gz | extractPBadapters [options] | filter10xPB             unique starting string for read name    ]
mawk '

############# BEGIN #########################
BEGIN {  ### GET/WRITE HEADER
    # set field separators
    FS="\t";
    OFS="\t";
}

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
    # middle bases
    while (match(info, "[ACTG][ACTG]+")) {
        # print("======", info);
        info = substr(info,1,RSTART-1) RLENGTH "N" substr(info, RSTART+RLENGTH);
    }
    # end bases
    if (match(info, "[ACTG]+$")) {
        info = substr(info,1,RSTART-1) RLENGTH "N";
    }
    # make short info
    shortInfo = info;
    while (match(shortInfo, "[0-9]+N")) {
        gsub("[0-9]+N","N", shortInfo);
    }
    print($1,info,shortInfo,seq,qual);
}
'