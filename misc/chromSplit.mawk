#!/bin/sh

# file with genomic pos in Pos (!MUST HAVE HEADER!!) | filterBed <bedfile> [-x ..useExonicCoords] [-c chrom ..which chrom to use; default all]
# filters any position-based file to positions included in a bed file
# takes bedfile and as parameter
# works only on stdin in a pipe


####### ARGPARSE ##################
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
        -f|--split-folder)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            folder=$2
            shift 2
        else
            echo "<chromSplit> Error: split-folder argument is missing\n[-f|--split-chrom (default=./split]" >&2
            exit 1
        fi
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
    file="'$1'";
    splitFolder="'${folder-"./split"}'";
    # add the slash if missing
    gsub(/([^/])$/, "&/", splitFolder);
    printf("<chromSplit> Splitting file %s by chrom into %s\n", file, splitFolder) > "/dev/stderr";
    if (file ~ /.gz$/) {
        readcmd = "cat " file " 2>/dev/null | gunzip ";
        n = split(file, S, ".");
        # get the proper extension to be added to file
        if (n > 2) {
            ext = "." S[n-1] "." S[n];
        } else {
            ext = ".gz"
        }

    } else {
        readcmd = "cat " file " 2>/dev/null ";
        n = split(file, S, ".");
        ext = "." S[n];
    }
    # get the basename of the file
    n = split(file, S, "/");
    # strip the path
    base = S[n];
    # strip the extension
    gsub(ext, "", base);
    if ((readcmd | getline) == 0) {
        printf("<chromSplit> file %s not found!\n", file) > "/dev/stderr";
            exit;
    } else {                   
        printf("<chromSplit> Reading file %s.\n", file) > "/dev/stderr";
    }
    if ($0 ~ "Chr\t") {
        printf("<chromSplit> Header detected\n") > "/dev/stderr";
        hasHeader=1;
        header=$0;
    }

    ##### Reading #################
    while ((readcmd | getline) > 0) {
        chrom = $1;
        if (chrom !~ "chr") {
            chrom = "chr" chrom;
        };
        if (currentChrom != chrom) {
            currentChrom = chrom;
            splitFile = splitFolder base "." chrom ext;
            print(header)>splitFile;
            printf("<chromSplit> Writing %s to %s\n", currentChrom, splitFile) > "/dev/stderr";
        }
        print($0)>>splitFile;
    }
}'
