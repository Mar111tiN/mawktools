#!/bin/bash

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # coverage Window
        -w|--window-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            WINDOW=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide size for coverage window [-w|--coverage-window-size (default=100)]" >&2
            exit 1
        fi
        ;;
        # step Size
        -s|--step-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            stepSize=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide size for step size [-s|--step-size (default=10)]" >&2
            exit 1
        fi
        ;;
        # line length of input genome file
        -l|--genome-line-length)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            LINE=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide minimum VAF for heteroSNP output [-v|--min-vaf (default=0)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<pile2CNV> Error: Unsupported flag $1" >&2
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

WINDOW=${WINDOW-100};
stepSize=${stepSize-10};
LINE=${LINE-50};

echo "stepSize:" $stepSize;
echo "WINDOW:" $WINDOW;
echo "LINE:" $LINE;