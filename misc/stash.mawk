#!/bin/sh


mawk '
NR == 1 {
    LEN=10;
    LINE=50;

    # make the split pattern from the desired length
    dot=".";
    for (i=0;i++<LEN;) {
        pattern = pattern dot;
    }
}
/[AGTCagtc]/{
    gsub(pattern,"&\n", $0);
    A=gsub(/[AaTt]/,"",$0);
    T=gsub(/[GgCc]/,"C",$0);
    print;
    next;
}
{
    print("50N");
}' | mawk '
NR ==1 {
    cols="Start,GC/AT";
    split(cols,COLS,",");
    for (i = 0; i++ < length(COLS)-1;) {
        printf("%s\t", COLS[i]);
    }
    printf("%s\n",COLS[length(COLS)]);
    basePos=1;
    next;
}
/50N/ {
    basePos+=50;
    next;
}
{
    GC = gsub("C", "", $0);
    N = gsub("N", "", $0);
    ratio = (GC + (N/2) / 10;
    printf("%s\t%s\n",basePos, ratio);
    basePos+=10;
}'