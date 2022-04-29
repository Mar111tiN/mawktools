#!/bin/sh


# >>>10xPBclean<<<
# v0.9
# takes a fastq(gz) file preprocessed with <10xPBsplit>
# and isolates 10x data { N[TSO]>N or N<[TSO]N } into clean data lines
# !!! is very computation-intensive --> outsource the revcomp functionality

# USAGE: 
# gunzip < fastq.gz |  10xPBextract [options] | 10xPBfilter -i | 10xPBsplit | 10xPBclean








mawk '
################ INIT #########################
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

############## functions ######################
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


($3~ /N\[TSO\]\>/) || ($3~ /<\[TSO\]N/) {
    info=$2;
    shortInfo=$3;
    seq=$4;
    qual=$5;

    # init Counter
    count=1;

    # isolate L-data
    while (match(seq, /[ATGC]+\[TSO\]>[ATGC]+/)) {
        outSeq=substr(seq,RSTART,RLENGTH);
        outQual=substr(qual,RSTART,RLENGTH);
        seq=substr(seq, 1,RSTART-1) substr(seq,RSTART+RLENGTH);
        qual=substr(qual, 1,RSTART-1) substr(qual,RSTART+RLENGTH);

        (match(info, /[0-9]+N\[TSO\]>[0-9]+N/));
        outInfo=substr(info,RSTART,RLENGTH);
        info=substr(info, 1,RSTART-1) substr(info,RSTART+RLENGTH);
        print($1 count++, outInfo, outSeq, outQual);
    }
    
    # isolate R-data
    while (match(seq, /[ATGC]+<\[TSO\][ATGC]+/)) {
        outSeq=substr(seq,RSTART,RLENGTH);
        outQual=substr(qual,RSTART,RLENGTH);
        seq=substr(seq, 1,RSTART-1) substr(seq,RSTART+RLENGTH);
        qual=substr(qual, 1,RSTART-1) substr(qual,RSTART+RLENGTH);

        (match(info, /[0-9]+N<\[TSO\][0-9]+N/));
        outInfo=substr(info,RSTART,RLENGTH);
        info=substr(info, 1,RSTART-1) substr(info,RSTART+RLENGTH);

        ### --> maybe outsource this component to the next process
        # reverse complement
        outInfo=revertInfo(outInfo);
        outSeq=revcomp(outSeq);
        outQual=rev(outQual);
        print($1 "<" count++, outInfo, outSeq, outQual);
    }
}
'