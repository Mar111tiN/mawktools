#!/bin/sh

## function that outputs the reverse complement of a sequence

mawk '

function revcomp(seq) {
    len=split(seq,SEQ,"");
    for (l=0;l++<len;) {
        revseq=revseq REV[SEQ[len-l+1]];
    }
    return revseq
}


BEGIN {
    # set reverser as global
    REV["A"]="T";
    REV["T"]="A";
    REV["G"]="C";
    REV["C"]="G";

    seq="'$1'";
    print(revcomp(seq));
}
'