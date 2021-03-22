#!/bin/sh

# USAGE:
# cat chrom.fa | base2GC [primarywindow=10bp] []-->

# ARGS

# OUTPUT:
# HEADER    Start   GC/AT
#           51      0.3


WINDOW=${1-100}; # Window Size for rolling GC

stepSize=${2-10};  # length of primary split window
LINE=${3-50};  ######## line length in genome file
THREADS=8;
################################
# SPLIT LINES
################################
mawk '
NR == 1 {


    stepSize='$stepSize';
    # make the split pattern from the desired length
    dot=".";
    for (i=0;i++<stepSize;) {
        pattern = pattern dot;
    }
    next;
}
/[AGTCagtc]/{
    gsub(pattern,"&\n", $0);
    gsub(/[AaTt]/,"",$0);
    gsub(/[GgCc]/,"C",$0);
    printf($0);
    next;
}
{
    print("-N-");
}' | mawk '

###############################
COMPUTE PRIMARY WINDOW
###############################

NR ==1 {
    # get the arguments from bash
    stepSize='$stepSize';
    LINE='$LINE';
    cols="Start,GC/AT";
    split(cols,COLS,",");
    for (i = 0; i++ < length(COLS)-1;) {
        printf("%s\t", COLS[i]);
    }
    printf("%s\n",COLS[length(COLS)]);
    basePos=1;
}
/-N-/ {
    basePos+=LINE;
    next;
}
{
    GC = gsub("C", "", $0);
    N = gsub("N", "", $0);
    ratio = (GC + (N/2)) / stepSize;
    printf("%s\t%s\n",basePos, ratio);
    basePos+=stepSize;
}' | mawk '

################################
# COMPUTE ROLLING WINDOW
################################

NR==1 {
    # get the arguments from bash args
    WIN='$WINDOW';  # window size
    stepSize='$stepSize'; # base interval of GC input

    L=WIN/stepSize; # length of HOLD ARRAY
    # the shift for getting the center of the window corresp to coords
    shift=(WIN / 2) - stepSize + 1;
    # printf("WIN:%s,stepSize:%s,L:%s,shift:%s\n",WIN,stepSize,L,shift)
    ########## SHOW ARRAY
    # for (i=0; i++<L;) {
    #     print(i,HOLD[i]);
    # }
    # INIT THE QC-SUM
    # is updated constantly for removing TIP and adding TAIL of HOLD
    # start with a half full neutral hold sum
    # in case the first row already has bases
    # INIT THE HOLD pointer as full
    pointer=10;
    lastPos=-10000;

    ## print the header
    print;
    next;
}


function update(minus, plus) {



    # remove TAIL
    SUM -= minus;
    # add new TIP
    HOLD[pointer] = plus;
    SUM += plus;
    printf("%s\t%.2f\n",lastPos-shift, SUM/L);

    ####=== DEBUG ===####
        # for (i=0;i++<L;){
        #     print(i, HOLD[i]);
        # }
    ####^^^ DEBUG ^^^####

    return SUM
}

NR==2 {

    # action for first row
    # find the pointer from the base position
    pointer=(($1-1) % WIN) / stepSize + 1;
    # print("NewPointerPos",newPointer);
    
    # fill up the HOLD 
    # re-init a neutral HOLD and its SUM
    SUM=L*0.5;
    for (i=0; i++ < L;){
        p=(pointer + i-1)%L + 1
        HOLD[p]=0.5;
    }
    # update lastPos
    lastPos=$1;   
    # add the new position
    update(0.5, $2)
    next;
}

# if pos difference is greater than stepSize we have a jump --> 
# find the pointer position based on lastPos
# fill up the array and flush

$1 - lastPos > stepSize {
    # cycle through pointer until original pointer position

    # keep the last pointer position
    lastPointer = pointer;
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;
    ###
    # print("JUMP", lastPos, "-->", $1);
    ###
    # increase current coords until new position is reached OR 
    # one round through pointers is finished

    while (pointer != lastPointer) {
        # increase coord with every step
        lastPos += stepSize;
        update(HOLD[pointer], 0.5)
        # new position is reached
        if (lastPos >= $1) {
            # update lastPos and exit 
            lastPos=$1;
            next;
        }
        # else, up the pointer (pointer > L --> L=1)
        pointer = (pointer)%L + 1;
    }

    ######### GAP TOO LARGE
    # find the pointer from the base position
    pointer=(($1-1) % WIN) / stepSize + 1;
    # print("NewPointerPos",newPointer);
    
    # fill up the HOLD 
    # re-init a neutral HOLD and its SUM
    SUM=L*0.5;
    for (i=0; i++ < L;){
        p=(pointer + i-1)%L + 1
        HOLD[p]=0.5;
    }
    # update lastPos
    lastPos=$1;   
    # add the new position
    update(0.5, $2)
    next;
}

# just the next in the interval
{   
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;
    lastPos=$1;

    update(HOLD[pointer], $2)

    # move pointer up by one and flush if at StepBorder
    # update lastPos
    lastPos=$1;
}

END {
    ## final flush
    # keep the last pointer position
    lastPointer = pointer;
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;
        # one round through pointers is finished

    while (pointer != lastPointer) {
        # increase lastPos with every step
        lastPos += stepSize;
        update(HOLD[pointer], 0.5)
        # else, up the pointer (pointer > L --> L=1)
        pointer = (pointer)%L + 1;
    }
}' | gzip