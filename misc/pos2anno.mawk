#!/bin/sh
# modify start and end positions including the left1-step with indels compatible with annovar

mawk ' NR == 1{ print "Chr\tStart\tEnd\tRef\tAlt\tMut_ID\tMut_ID_old\tcount\ttype\MutStatus\n" } NR > 1 { start= 2;𝑒𝑛𝑑= 2; RL=length( 3);𝐴𝐿=𝑙𝑒𝑛𝑔𝑡ℎ( 4);

# variant is an insertion
if (AL > 1) {
    if (RL == 1) {
        start=start-1;
        end=start-1;
        $3="-";
            $4=substr($4,2)
        }
    } else if (RL > 1) {
        if (AL == 1) {
            start=start+1;
            end=start + RL -2;
            $3=substr($3,2);
            $4="-";
        }
    }
    # print all fields
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,start,end,$3,$4,$5,$6,$7,$8);
}'