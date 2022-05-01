#!/bin/sh


# >>>10xPBclean<<<
# v0.9
# v0.9.1 -- fix counter bug (always incrementing)
# takes a fastq(gz) file preprocessed with the 10xPB toolchain components:
#     <10xPBextract>   |      transforms into one-line data and extracts the relevant oligo adapter sequences
# --> <10xPBinfo>      |      adds a long and a short info field discribing the sequence structure with regard to known adapters
# --> <10xPBfilter>    |      filters for sequences containing 25N[TSO]> / <[TSO]25N
# --> <10xPBsplit>     |      splits sequences with facing 10xSignature components into two separate lines
# and isolates 10x data { N[TSO]>N or N<[TSO]N } into clean data lines
# !!! is very computation-intensive -->v0.9.1: outsource the revcomp functionality

# USAGE: 
# gunzip < fastq.gz |  10xPBextract [options] | 10xPBinfo | 10xPBfilter -i | 10xPBsplit | 10xPBclean



mawk '
################ BEGIN #########################
BEGIN {
    OFS="\t";
    printf("<10xPBclean> Isolating 10x sequences of structure  N[TSO]>N | N<[TSO]N \n") >> "/dev/stderr";
    # set reverser as global
    REV["A"]="T";
    REV["T"]="A";
    REV["G"]="C";
    REV["C"]="G";
    REV["|"]="|";
}

############## FUNCTIONS ######################
function rev(qual) {
    gsub(/<\[TSO\]/, ">]OST[", qual);
    revqual="";
    len=split(qual,QUAL,"");
    for (l=0;l++<len;) {
        revqual=revqual QUAL[len-l+1];
    }
    return revqual
}

function revcomp(seq) {
    # mask the TSO tag
    gsub(/<\[TSO\]/, "|", seq);
    revseq="";
    len=split(seq,SEQ,"");
    for (l=0;l++<len;) {
        revseq=revseq REV[SEQ[len-l+1]];
    }
    # unmask TSO tag
    gsub(/\|/, "[TSO]>", revseq);
    return revseq
}

function revertInfo(info) {
        match(info, /[0-9]+N<\[TSO/);
        N2=substr(info,1,RLENGTH-5);
        match(info, /<\[TSO\][0-9]+N/);
        N1=substr(info,RSTART+6,RLENGTH-6);
        info=N1 "[TSO]>" N2;
        return info
}
{
    count=1;
}
############# READ MATCHING LINES #########################
$3~ /N\[TSO\]\>/ {
    info=$2;
    shortInfo=$3;
    seq=$4;
    qual=$5;

    # use fullInfo field for extracting the sequences
    while (match(info, /[0-9]+N\[TSO\]>[0-9]+N/)) {
        # get starting coords in seq from N-expansion
        # adding to rstart by all <INT>N
        seqStart=RSTART;
        seqLength=RLENGTH;
        # get the upstream info for seq-coord update
        preInfo=substr(info, 1,RSTART-1);
        # get the output info field (the match)
        outInfo=substr(info,RSTART,RLENGTH);
        # update the info field (remove left-side+match)
        info=substr(info,RSTART+RLENGTH);

        # >>>DEBUG
        # printf("%s\t", preInfo);
        # printf("%s\t%s\t", $1, seq);
        # >>>DEBUG

        while (match(preInfo, /[0-9]+N/)) {
            # increase seqStart by N
            seqStart = seqStart + substr(preInfo, RSTART,RLENGTH - 1) - RLENGTH;
            preInfo = substr(preInfo, RSTART+RLENGTH);
        }

        match(outInfo, /^[0-9]+N/)
        # increase seqLength by N from outInfo
        seqLength = seqLength + substr(outInfo, RSTART,RLENGTH-1) - RLENGTH;
        match(outInfo, /[0-9]+N$/)
        # increase seqLength by N from outInfo
        seqLength = seqLength + substr(outInfo, RSTART,RLENGTH-1) - RLENGTH;
        # print(seqStart, seqLength);

        # update seq and qual
        outSeq=substr(seq,seqStart, seqLength);
        outQual=substr(qual,seqStart, seqLength);
        seq=substr(seq,seqStart+seqLength);
        qual=substr(qual,seqStart+seqLength);

        ##### OUTPUT ###################
        print($1 count++, outInfo, outSeq, outQual);
    }
}

############# READ MATCHING LINES #########################
$3~ /<\[TSO\]N/ {
    info=$2;
    shortInfo=$3;
    seq=$4;
    qual=$5;

    # use fullInfo field for extracting the sequences
    while (match(info, /[0-9]+N<\[TSO\][0-9]+N/)) {
        # get starting coords in seq from N-expansion
        # adding to rstart by all <INT>N
        seqStart=RSTART;
        seqLength=RLENGTH;
        # get the upstream info for seq-coord update
        preInfo=substr(info, 1,RSTART-1);
        # get the output info field (the match)
        outInfo=substr(info,RSTART,RLENGTH);
        # update the info field (remove left-side+match)
        info=substr(info,RSTART+RLENGTH);

        # >>>DEBUG
        # printf("%s\t", preInfo);
        # printf("%s\t%s\t", $1, seq);
        # >>>DEBUG
        
        while (match(preInfo, /[0-9]+N/)) {
            # increase seqStart by N
            seqStart = seqStart + substr(preInfo, RSTART,RLENGTH - 1) - RLENGTH;
            preInfo = substr(preInfo, RSTART+RLENGTH);
        }

        match(outInfo, /^[0-9]+N/)
        # increase seqLength by N from outInfo
        seqLength = seqLength + substr(outInfo, RSTART,RLENGTH-1) - RLENGTH;
        match(outInfo, /[0-9]+N$/)
        # increase seqLength by N from outInfo
        seqLength = seqLength + substr(outInfo, RSTART,RLENGTH-1) - RLENGTH;
        # print(seqStart, seqLength);

        # update seq and qual
        outSeq=substr(seq,seqStart, seqLength);
        outQual=substr(qual,seqStart, seqLength);
        seq=substr(seq,seqStart+seqLength);
        qual=substr(qual,seqStart+seqLength);
#         ### --> maybe outsource this component to the next process
        reverse complement
        outInfo=revertInfo(outInfo);
        outSeq=revcomp(outSeq);
        outQual=rev(outQual);
        ##### OUTPUT ###################
        print($1 count++, outInfo, outSeq, outQual);
    }
}'