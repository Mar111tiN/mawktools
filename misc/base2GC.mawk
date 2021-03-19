#!/bin/sh

mawk '
/[AGTCagtc]/{
    gsub(/........../,"&\n", $0);
    gsub(/[AaTt]/,"",$0);
    gsub(/[GgCc]/,"C",$0);
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
$0 ~ /C/ {
    GC=0;
    GC = gsub("C", "", $0);
    ratio = GC / 10;
    printf("%s\t%s\n",basePos, ratio);
    basePos+=10;
}'