#!/bin/sh

# short version of cleanpileup
### cleans output from samtools mpileup where only Chr, Start, and the read data has been left
# Chr   Start   Read1 Read2 Read3
# this can be done with pon2cols tool or deliberate cutting
# removes from all reads the traces of base location and indel lengths


mawk '
{
    read = $0;
    # remove position traces from all read fields
    gsub(/\^[^\t]|\$/,"",read);
    # 
    while (match(read,/[+-][0-9]+/)) {
        pre = substr(read,1,RSTART-2);
        indel = substr(read,RSTART,1); # + or -
        base = substr(read,RSTART-1,1); # base before the deletion
        l = substr(read,RSTART+1,RLENGTH-1);
        indel_bases = substr(read,RSTART+RLENGTH,l);
        post = substr(read,RSTART+RLENGTH+l);
        if (indel == "-") {
            if (match(indel_bases,/[ACGT]/)) {
                base = "D";
            } else {
                base = "d";
            }
        } else {
            if (match(indel_bases,/[ACGT]/)) {
                base = "I";
            } else {
                base = "i";
            }            
        }

        read = pre base post;
    }      
# print all fields
    print read;
}'
