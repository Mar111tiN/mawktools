#!/bin/sh

## function that looks for pattern in sequence and replaces it with a tag

mawk '
function search_replace(pattern,tag) {
    sl=length(seq);
    while (match(seq, pattern)) { 
        l=substr(seq,1,RSTART-1);
        r=substr(seq,RSTART+RLENGTH);
        seq=l tag r;
        hits++;
    }
    return hits
}

BEGIN {
    seq="'$1'";
    pattern="'$2'";
    tag="'$3'";
    DATA["seq"]=seq;
    DATA["hits"]=0;
    hits = search_replace(pattern, tag);
    print(seq, hits);
}'