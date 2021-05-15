#!/bin/sh

# takes genmap txt file as input and returns file with genomic positions


####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # set precision
        -p|--precision)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            precision=$2
            shift 2
        else
            echo "<convertGenmap> Error: Precision argument is missing\n[-p|--precision (default=8]" >&2
            exit 1
        fi
        ;;
        -u|--unique-only)
        uniqueOnly=1
        shift
        ;;
        -*|--*=) # unsupported flags
        echo "<chromSplit> Error: Unsupported flag $1" >&2
        exit 1
        ;;

        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done

# I have no positional args
# # set positional arguments in their proper place
eval set -- "$PARAMS"

mawk '
BEGIN {
    uniqueOnly="'${uniqueOnly-0}'";
    prc='${precision-8}'
    # set the record separator to space
    RS=" ";
    printf("Chr\tPos\tmap\n");
}
NR==1 {
    # find the chrom marker
    sub(">", "", $0);
    # split of the first base position (\n is not the RS anymore)
    split($0, S, "\n");
    chrom=S[1];
    pos=1;
    map=S[2];
    printf("%s\t%s\t%s\n", chrom, pos, map);
    lastMap=map;
    next;
}
/>chr.+$/ { # at chrom change
    # split of the first base position (\n is not the RS anymore)
    # output is "M\n>chrX\nM
    split($0, S, "\n");
    # print the last map of current chrom
    printf("%s\t%s\t%s\n", chrom, ++pos, S[1])
    chrom = S[2];
    # find the chrom marker
    sub(">", "", chrom);
    # split of the first base position (\n is not the RS anymore)
    pos=1;
    map=S[3];
    printf("%s\t%s\t%s\n", chrom, pos, map);
    lastMap=map;
    next;
}
{
    pos++;
    map = (prc > 7) ? $0 : (int((($0 + 0.5) * 10^prc)) / 10^prc) -0.5;
    if ((map==lastMap) && (uniqueOnly==1)) next;
    printf("%s\t%s\t%s\n", chrom, pos, map);
    lastMap=map;
}
'
