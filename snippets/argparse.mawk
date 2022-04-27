#!/bin/sh

# Snippet for loading arguments via command line into mawk 
# with short and long options 
# and arguments
# arguments are loaded in shell wrapper and transfered to mawk 

# USAGE: 
# toolName [options] arg1 arg2
# [     -f | --flag                 <Flag=False>                toggle function in tool           ]
# [     -o | --option1              <INT=0>                      set parameter in tool            ]
# [     -O | --option2              <Int=100>                   set parameter in tool             ]
# [     -p | --path                 <path to file>              describe usage in tool            ]


####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # flag
        -f|--flag1)
        flag=1
        shift
        ;;
        # option1
        -o|--option1)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            option1=$2
            shift 2
        else
            echo "<toolName> Error: param for option1 is missing\n[-o|--option1]" >&2
            exit 1
        fi
        ;;
        # option2
        -O|--option2)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            option2=$2
            shift 2
        else
            echo "<toolName> Error: param for option2 is missing\n[-O|--option2]" >&2
            exit 1
        fi
        ;;
        # file path
        -p|--path)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            path=$2
            shift 2
        else
            echo "<toolName> Error: Path to file is missing\n[-p|--path]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<toolName> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done
# # set positional arguments in their proper place
eval set -- "$PARAMS"

mawk '
#############################################
############# BEGIN #########################
BEGIN {  ### GET/WRITE HEADER
    ##### params and args (with defaults)
    # params
    path="'${path-"test.txt"}'";

    option1='${option1-0}';
    option2='${option2-100}';
    flag='${flag-0}';
    # args
    arg1="'${1-"arg1"}'";
    arg2="'${2-"arg2"}'";
    printf("<toolName> Reading file %s [option1=%s|option2=%s|flag=%s|arg1=%s|arg2=%s\n", path, option1, option2, flag, arg1, arg2) >> "/dev/stderr";
}'