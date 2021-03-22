#!/bin/sh

### cleans samtools mpileup output
# first any end-of-read markers like ^I or ^[ will be removed
# cleanpileup takes extra care not to change any of the quality fields 
# .. which might accidentally contain such traces
# next, the A+12T for inserts and the T-12A for deletions will be converted to I/i and D/d

mawk '
NR == 1 {
    samples = (NF - 3) / 3;
    printf("Chr\tStart\tRef");
    for (s=0;s++<samples;) {
        printf("\tDepth%s\tRead%s",s,s);
    }
    printf("\n");
}
{   
    # print generic fields
    printf($1);
    for (c=1;c++<3;) {
        printf("\t%s",$c);
    }

    # cycle through the reads
    for (i = 0; i++ < samples;) {
        col = (i*3) + 2;
        # remove position traces from all read fields
        gsub(/\^[^\t]|\$/,"",$col);
        # 
        while (match($col,/[+-][0-9]+/)) {
            pre = substr($col,1,RSTART-2);
            indel = substr($col,RSTART,1); # + or -
            base = substr($col,RSTART-1,1); # base before the deletion
            l = substr($col,RSTART+1,RLENGTH-1);
            indel_bases = substr($col,RSTART+RLENGTH,l);
            post = substr($col,RSTART+RLENGTH+l);
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
            $col = pre base post;
        }
        printf("\t%s\t%s",$(col-1),$col);
    }
    printf("\n");
}   '