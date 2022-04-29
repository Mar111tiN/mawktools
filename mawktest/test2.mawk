#!/bin/sh

# >>>extractPBAdapters<<<
# v0.9
# takes a fastq(gz) file and extracts the oligo sequences given in file



# USAGE: 
# gunzip < fastq.gz | extractPBadapters
# [     -m | --allowMismatch       <Flag=False>                 search for oligos with hamming dist 1   ]
# [     -o | --oligoFile           <path to file>               oligo-file containing name and seq      ]
# [     -h | --Nhits               <Int=3>                      stop searching after N hits in seq      ]

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # mismatch
        -m|--allowMismatch)
        allowMismatch=1
        shift
        ;;
        # oligoFile
        -o|--oligoFile)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            oligoFile=$
            shift 2
        else
            echo "<extractPBAdapters> Error: Path to oligoFile is missing\n[-o|--oligoFile]" >&2
            exit 1
        fi
        ;;
        # number of hits
        -h|--Nhits)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            Nhits=$2
            shift 2
        else
            echo "<extractPBAdapters> Error: Value for param Nhits is missing\n[-h|--Nhits]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<extractPBAdapters> Error: Unsupported flag $1" >&2
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
############## functions ######################
function revcomp(seq) {
    revseq="";
    len=split(seq,SEQ,"");
    for (l=0;l++<len;) {
        revseq=revseq REV[SEQ[len-l+1]];
    }
    return revseq
}

############# BEGIN #########################
BEGIN {  ### GET/WRITE HEADER
    ##### params and args (with defaults)
    # params
    oligoFile="'${oligoFile-""}'";
    if (oligoFile == "") {
        printf("<extractPBAdapters> You have to provide an oligo file!\n") >> "/dev/stderr";
        exit;
    }

    # set read command and check for file existence
    readCMD = "cat " oligoFile " 2>/dev/null "
    # first line is header and is not needed 
    if ((readCMD | getline) == 0) { 
        printf("<extractPBAdapters> File %s not found!\n", file) > "/dev/stderr";
        exit;
    }
    Nhits='${Nhits-3}';
    allowMismatch='${allowMismatch-1}';
    printf("<extractPBAdapters> Reading oligo file %s [Nhits=%s|allowMismatch=%s]\n", oligoFile, Nhits, allowMismatch) >> "/dev/stderr";
    
    
    ### reading the oligo file into array along with reverse complements
    # set reverser as global
    REV["A"]="T";
    REV["T"]="A";
    REV["G"]="C";
    REV["C"]="G";

    while ((readCMD | getline) > 0) {
        OLIGOS[$1] = $2;
        REVOLIGOS[$1] = revcomp($2)
        # get the tags for the oligos for the replacement
        # get length of oligo name and seq
        l = length($2);
        n = length($1);
        # forward tags
        tag = "[" $1 "]";
        fwdtag=tag;
        # first write to DELTAG (TAG in case of deletion), then add the last <
        for (i=0;i++<l-n-3;) {
            fwdtag = fwdtag ">";
        }
        DELTAGS[$1] = fwdtag;
        TAGS[$1] = fwdtag ">";

        # reverse tags
        revtag=tag;
        for (i=0;i++<l-n-3;) {
            revtag = "<" revtag;
        }
        DELREVTAGS[$1] = revtag;
        REVTAGS[$1] = "<" revtag;

        ### get the hamming1 oligos (MMOLIS)for each oligo and revoligo
        # DEL version is representing the deletion
        if (allowMismatch == 1) {
            for (i=0;i++<l;) {
                MMOLIGOS[$1 "-" i] = substr($2,1,i-1) "." substr($2,i+1)
                DELOLIGOS[$1 "-" i] = substr($2,1,i-1) "" substr($2,i+1)
                MMREVOLIGOS[$1 "-" i] = substr(REVOLIGOS[$1],1,i-1) "." substr(REVOLIGOS[$1],i+1)
                DELREVOLIGOS[$1 "-" i] = substr(REVOLIGOS[$1],1,i-1) "" substr(REVOLIGOS[$1],i+1) 
            }
        }
    }

    # for (o in OLIGOS) {
    #     print(o, "\n", OLIGOS[o], REVOLIGOS[o], "\n", MMOLIGOS[o "-" 3], MMREVOLIGOS[o "-" 3], "\n", TAGS[o], REVTAGS[o], "\n", DELTAGS[o], DELREVTAGS[o]) ;
    # }
    RS="\n@m64316_";
    FS="\n";
    OFS="\t";
}

# remove 
NR==1{
    $1=substr($1,9);
}
function search_replace(pattern,tag) {
    while (match(seq, pattern)) { 
        ls=substr(seq,1,RSTART-1);
        rs=substr(seq,RSTART+RLENGTH);
        seq=ls tag rs;
        ls=substr(qual,1,RSTART-1);
        rs=substr(qual,RSTART+RLENGTH);
        qual=ls tag rs;
        hits++;
        #### 
        # print("H",pattern,hits);
        ####
    }
}

{   
    seq=$2;
    qual=$4;
    # cycle through the oligos and replace oligos
    for (o in OLIGOS) {
        # init the hits
        hits = 0;

        # first the actual oligos
        search_replace(OLIGOS[o],TAGS[o]);
        search_replace(REVOLIGOS[o],REVTAGS[o]);

        # here Nhits could also be oligo-specific
        if (hits >= Nhits) continue;

        # get the oligo length for all
        l=length(OLIGOS[o]);
        if (allowMismatch == 1) {
            for (i=0;i++<l;) {
                # print(DELOLIGOS[o "-" i]);
                if (hits < Nhits) {
                    #####
                    # print(i, MMOLIGOS[o "-" i],MMREVOLIGOS[o "-" i],DELOLIGOS[o "-" i],DELREVOLIGOS[o "-" i]);
                    #####
                    search_replace(MMOLIGOS[o "-" i],TAGS[o]);
                    search_replace(MMREVOLIGOS[o "-" i],REVTAGS[o]);
                    search_replace(DELOLIGOS[o "-" i],DELTAGS[o]);
                    search_replace(DELREVOLIGOS[o "-" i],DELREVTAGS[o]);                    
                };
            }
        }
    }

    ### retrieve the short string
    info=seq;
    gsub(">+",">",info);
    gsub("<+","<", info);
    gsub("[ACTG][ACTG]+","N",info);
    print($1,info,seq,qual);
    
}
'