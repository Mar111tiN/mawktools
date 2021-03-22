#!/bin/sh

# INPUT
# takes clean pileup from samtools mpileup with ref (...,,,AT.Ii..)
# and base qualities removed
# can contain ExonPos (auto-detected)
# 1 <= samples < N
# 	Chr	Start	Ref	Cov1	Read1	Cov2	Read2	...	CovN	ReadN	[ ExonPos ]	

# USAGE: 
# samtools mpileup -f $HG38 -q 20 -Q 25 -l $BED -r chr? | cleanpileup | filterBed $BED 1 chr? |
# pile2CNV snp_file [minCoverage=0] []

# OUTPUT
# Chr Start [ExonPos] Ref Depth1 Read1 Depth2	Read2...
# last sample should be the reference sample
# ARGUMENTS:
# 	ARG1: output file for redirected SNP-file
#   ARG2: minCoverage for output
#   ARG3: coverage Window
#   ARG4: step size for new Window

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    case "$1" in
        -x|--keep-exon-pos)
        keepExon=1
        shift
        ;;
        # snp_file output
        -o|--output-snp-file)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            snpFile=$2
            shift 2
        else
            echo "<pile2CNV> Error: output file for heteroSNP is missing\n[-o|--output-snp-file]" >&2
            exit 1
        fi
        ;;
        # minimum coverage
        -c|--min-coverage)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minCov=$2
            shift 2
        else
            echo "<pile2CNV> Error: Argument for $1 is missing" >&2
            exit 1
        fi
        ;;
        # coverage Window
        -w|--coverage-window-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            covWindow=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide size for coverage window [-w|--coverage-window-size (default=100)]" >&2
            exit 1
        fi
        ;;
        # step Size
        -s|--step-size)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            covWindow=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide size for step size [-s|--step-size (default=10)]" >&2
            exit 1
        fi
        ;;
        # minVAF
        -v|--min-vaf|min-VAF)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minVAF=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide minimum VAF for heteroSNP output [-v|--min-vaf (default=0)]" >&2
            exit 1
        fi
        ;;
        # minDepth
        -d|--min-depth)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            minDepth=$2
            shift 2
        else
            echo "<pile2CNV> Error: Provide minimum VAF for heteroSNP output [-d|--min-depth (default=10)]" >&2
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

snpFile=${snpFile-"test.snp"};
covWindow=${covWindow-100};
minCov=${minCov-0};
stepSize=${stepSize-10};
minVAF=${minVAF-0};
minDepth=${minDepth-10};

### DEBUG ######
# echo "snpFile:" $snpFile;
# echo "minCoverage: " $minCov;
# echo "windowSize: " $covWindow;
# echo "stepSize:" $stepSize;
# echo "minVAF:" $minVAF;
# echo "minDepth:" $minDepth;
### DEBUG ######

mawk '
#### BEGIN
NR == 1 {
    ##### ARGS
    snpFile="'$snpFile'";
    if (snpFile !~ ".snp") snpFile = snpFile ".snp"
    # setup SNP file and write its Header
    printf("Writing snpData to %s\n", snpFile) >> "/dev/stderr";
    covWindow='$covWindow';
    stepSize='$stepSize';
    printf("Rolling coverage: WindowSize=%s | stepSize=%s | minCoverage=%s\n", covWindow, stepSize, '$minCov') >> "/dev/stderr";

    #### INIT Length L of DATA ARRAY 
    L=int(covWindow / stepSize);

    ###### HEADER
    baseHeader = "Chr\tStart"
    # detect XPos
    if ($NF == "ExonPos") hasXPos = 1;
    if (hasXPos) { 
        # print("ExonPos detected") > "/dev/stderr";
        baseHeader = baseHeader "\tExonPos";
    }
    printf(baseHeader);
    printf(baseHeader) > snpFile;

    # detect samples
    samples = (NF-3-hasXPos)/2;
    print(samples, "samples detected") > "/dev/stderr"

    # assign the reference column (last read column)
    refCol = NF - hasXPos;

    for (s=0; s++ < samples;) {
        printf("\tCov%s", s);
        printf("\tVAF%s", s) >> snpFile;
        # COL stores the col number for each sample
        COL[s] = 2 + (2 * s);


        # !!!!!!!!! INIT NOT necessary !!!
        # init the ARRAYS
        COVSUM[s] = 0;
        for (bin=-1; ++bin<L;) {
            COVBIN[bin "-" s] = 0;
        }

    }
    printf("\n");
    printf("\n") >> snpFile;
    next;
}

######## VAF detection ###############
$refCol ~ /[.,]*[AaCcTtGgDdI]+[.,]*/ {
    print($0) >> snpFile;
    a=1;
}


##### COV detection first row
NR == 2 {
    ## get the exonShift 
    if (hasXPos) {
        exonShift=$2-$NF;
        # print("ExonShift", exonShift)
    }
    # determine the bin for bin-wise COVSUM accumulation
    # bins in [0..L-1]
    # set the lastStep and store the coverage
    lastStep=int(($2-1) / stepSize);
    bin=lastStep%L;
    
    ####### DEBUG ########
    # print($2,$4, lastStep, bin)
    ####### DEBUG ########

    for (s in COL) {
        COVBIN[bin "-" s] = $COL[s];
    }
    next;
}

function output(step, doReset) {

    ### helper to output and update the coverage
    ### get position of previous windows center
    centerPos =  (step - (L / 2) + 1) * stepSize;
    # output centerPos
    printf("%s\t%s",$1,centerPos);
    # print exonCoords
    if (exonShift) printf("\t%s", centerPos-exonShift);

    thisBin=step%L;
    nextBin=(step+1)%L;

    for (s=0;s++<samples;) {
        # output the COVSUM array
        printf("\t%s", COVSUM[s] / covWindow);

        ### DEBUG #####
        # if (s==1) {
        #     print("reset=", doReset * 1);
        #     printf("\nS=%s: -bin%s=%s + bin%s=%s --> S=%s\n", COVSUM[s], nextBin, COVBIN[nextBin "-" s], thisBin, COVBIN[thisBin "-" s],COVSUM[s] - COVBIN[nextBin "-" s] + COVBIN[thisBin "-" s]);
        # }
        ### DEBUG #####

        # update Tip and Tail of COVSUM array
        COVSUM[s] = COVSUM[s] - COVBIN[nextBin "-" s] + COVBIN[thisBin "-" s]
        # reset bin with new coverage / zero
        COVBIN[nextBin "-" s] = (doReset) ? 0 : $COL[s];

        ########### DEBUG #############
        # if (s==1) {
        #     for (c=-1; ++c< L;) {
        #         print(c, COVBIN[c "-" s]);
        #     }
        # }
        ########### DEBUG #############
    }
    printf("\n");
}

# ######## COV detection ###############
{
    ### create an ARRAY for covWindow / stepSize bins
    pos=$2;

    currentStep=int((pos -1) / stepSize);
    stepDist=currentStep - lastStep;

    bin=currentStep%L;

    ####### DEBUG ########
    # print($2,$4, currentStep, bin)
    ####### DEBUG ########

    ## SIMPLE CASE: same bin
    if (stepDist == 0) {
        # add coverage to respective bin
        for (s in COL) {
            COVBIN[bin "-" s] += $COL[s];
        }
        lastStep=currentStep;
        next;
    }

    ## SIMPLE CASE: move into next bin
    if (stepDist == 1) {

        output(lastStep);

        # update step
        lastStep=currentStep;
        next;
    }

    ## COMPLEX CASE: jump over bins
    if (stepDist > 1) {
        # print("JUMP");
        #### loop from lastStep to currentStep
        # lastStep --> lastStep + 2
        for (step=lastStep; step++<currentStep;) {
            # output
            output(step-1, 1);
            if (step-lastStep>3) break       
        }
        output(currentStep - 1)
        # update
        lastStep=currentStep;
    }
}' | mawk 'BEGIN{minCov='$minCov'} NR==1 || ($NF>minCov) || ($(NF-1)>minCov)'
