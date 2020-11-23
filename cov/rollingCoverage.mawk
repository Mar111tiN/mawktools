#!/bin/sh

# ################ GENOMIC ROLLING COVERAGE ##############
# takes the output from chromCoverage and outputs the rolling average coverage per 100 (default) bases
# takes frame with rows Chr Pos Coverage

mawk '

BEGIN {
  # print Header
    printf("%s\t%s\t%s\n","Chr","Pos","Coverage");
    width=int("'${1:-100}'" / 2) * 2; # parameter for the rolling scope (default = 100)
    half=width / 2;
    first_line=1;
}

$0 ~ /Chr\t/ {
    next;
}

# if chrom changes 
$1 != currentChrom {
    currentChrom = $1;
    # READ LINE
    lastPos=$2;
    lastCov=$3;
    # posCenter is the middle between sumL and sumR
    posCenter=int(lastPos / half) * half;
    next;
}
{
    # READ LINE
    chrom=$1;
    pos=$2;
    cov=$3;
    # first check, whether a posCenter border has been crossed
    # sumL is left of posCenter [posCenter-half..posCenter)
    # sumR is right of posCenter including posCenter [posCenter..posCenter + half)
    # print("posCenter", posCenter);
    if (pos >= posCenter + half) { # jumped over a border 
        # fill up sumB
        ##########
        sumR += lastCov * (half - lastPos + posCenter);
        covMean = (sumL + sumR) / width;
        printf("%s\t%s\t%s\n",chrom,posCenter, covMean);

        if (pos >= posCenter + width) { # stepped over two borders
        # spill the right side coverage
        sumL = lastCov * half;
        covMean = (sumR + sumL) / width;
        printf("%s\t%s\t%s\n",chrom,posCenter + half,covMean);
        posCenter=int(pos / half) * half;
    
        } else { # stepped over one border
        # step up the posCenter and fill sumR
        posCenter = posCenter + half;
        sumL = sumR;
        }
        sumR = lastCov * (pos - posCenter + 1);
    } else { # did not jump a border
        sumR += lastCov * (pos - lastPos);
    }
    lastCov = cov;
    lastPos = pos;
    # printf("Pos: %s\tCov: %s\tsumL: %s\tsumR: %s\n",pos, cov, sumL, sumR);
 }
 END {
   print(posCenter,sumL, sumR)
   printf("%s\t%s\t%s\n",chrom,posCenter,(sumL + sumR) / width);
   printf("%s\t%s\t%s\n",chrom,posCenter + half, sumR / width);
 }'