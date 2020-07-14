#!/bin/sh
# works like matrix2EBinput but treats inserts and deletions separately

outputBase=$2;
splitFactor=$1;

# gunzip < $1 | 
mawk '
BEGIN {
    file_name = "'$2'";
    print file_name;
    splits = '$splitFactor'*4;
    for (x = 0; x++ < splits;) {
        toFile[x%splits]= int((x-1) / 4);
    }  
}
{
    file_number=toFile[NR%splits];
    print >> file_name"."file_number;
}

END {
    for (x = -1; ++x < splits;) {
        printf("NR=%s:file%s.\n", x, toFile[x]);
        }
}'
