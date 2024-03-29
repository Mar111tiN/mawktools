#!/bin/sh

# ASSUMES empty rows to be NNN



WINDOW=${1-100}; # Window Size for rolling GC

# ASSUMES 10-base intervals as input
# can be set as second arg
LEN=${2-10};

mawk '
NR==1 {
    # get the arguments from bash args
    WIN='$WINDOW';  # window size
    LEN='$LEN'; # base interval of GC input

    L=WIN/LEN; # length of HOLD ARRAY

    # the shift for getting the center of the window corresp to coords
    shift=(WIN / 2) - LEN + 1;
    printf("WIN:%s,LEN:%s,L:%s,shift:%s\n",WIN,LEN,L,shift)
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
    coord=-10000;

    ## print the header
    print;
    next;
}

function update( HOLD, SUM, coords, pointer, minus, plus, L,  i) {

    # remove TAIL
    SUM -= minus;
    # add new TIP
    HOLD[pointer] = plus;
    SUM += plus;

    if (pointer%(L/2)==0) {
    ####=== DEBUG ===####
        # for (i=0;i++<L;){
        #     print(i, HOLD[i]);
        # }
    ####^^^ DEBUG ^^^####

    printf("%s\t%.2f\n",coords, SUM/L);
    }
    return SUM
}

# if coord difference is greater than 10
# we have a jump --> 
# find the pointer position based on coords
# fill up the array and flush
$1 - coord > LEN {
    # cycle through pointer until original pointer position

    # keep the last pointer position
    lastPointer = pointer;
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;
    ###
    # print("JUMP", coord, "-->", $1);
    ###
    # increase current coords until new position is reached OR 
    # one round through pointers is finished
    # break while if coords < 0 (start)
    if (coord > 0) {
        while (pointer != lastPointer) {
            # increase coord with every step
            coord += LEN;
            SUM = update(HOLD, SUM, coord-shift, pointer, HOLD[pointer], 0.5, L)
            # new position is reached
            if (coord >= $1) {
                # update coords and exit
                coord=$1;
                next;
            }
            # else, up the pointer (pointer > L --> L=1)
            pointer = (pointer)%L + 1;
        }
    }

    ######### GAP TOO LARGE
    # find the pointer from the base position
    pointer=(($1-1) % WIN) / LEN + 1;
    # print("NewPointerPos",newPointer);
    
    # fill up the HOLD 
    # re-init a neutral HOLD and its SUM
    SUM=L*0.5;
    for (i=0; i++ < L;){
        p=(pointer + i-1)%L + 1
        HOLD[p]=0.5;
    }
    # update coord
    coord=$1;   
    # add the new position
    SUM = update(HOLD, SUM, coord-shift, pointer, 0.5, $2, L)

    ####
    # print($1, pointer);
    ####

 
    # print($0, pointer);
    next;
}

# just the next in the interval
{   
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;

    ####
    # print($1, pointer);
    ####

    SUM = update(HOLD, SUM, $1- shift, pointer, HOLD[pointer], $2, L)

    # move pointer up by one and flush if at 50
    # update coord
    coord=$1;
}
END {
    ## final flush
    # keep the last pointer position
    lastPointer = pointer;
    # up the pointer (pointer > L --> L=1)
    pointer = (pointer)%L + 1;
        # one round through pointers is finished
    # break while if coords < 0 (start)
    while (pointer != lastPointer) {
        # increase coord with every step
        coord += LEN;
        SUM = update(HOLD, SUM, coord-shift, pointer, HOLD[pointer], 0.5, L)
        # else, up the pointer (pointer > L --> L=1)
        pointer = (pointer)%L + 1;
    }
}
'