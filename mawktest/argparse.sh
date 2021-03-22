#!/bin/bash

PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -x|--use-exon-pos)
        useExonicCoords=1
        shift
        ;;
        # snp_file output
        -c|--chrom)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            filterChrom=$2
            shift 2
        else
            echo "Error: chromosome argument is missing\n[-c|--chrom (default=all)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "Error: Unsupported flag $1" >&2
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

# echo $1;
# echo $2;

useExonicCoords=${useExonicCoords-0};
filterChrom=${filterChrom-""};
bedFile=$1;

echo "bedFile:" $bedFile;
echo "useExonicCoords:" $useExonicCoords;
echo "filterChrom:" $filterChrom;