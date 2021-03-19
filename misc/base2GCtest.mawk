#!/bin/sh

# ASSUMES 10-base intervals as input
# ASSUMES empty rows to be >=5N

mawk '
BEGIN {
    WIN=100;
    L=WIN/10;
    shift=L-9;
    ### INIT THE HOLD ARRAY
    for (i=0; i++<L;){
        HOLD[i]=0.5;
    }

    ########## SHOW ARRAY
    # for (i=0; i++<L;) {
    #     print(i,HOLD[i]);
    # }
    # INIT THE QC-SUM
    # is updated constantly for removing TIP and adding TAIL of HOLD
    SUM=L*0.5;
    # INIT THE HOLD pointer
    pos=1;
    coord=0;
}
# special case for first encountered row
NR==1{

}
# if coord difference is greater than 10
# we have a jump --> 
# find the pointer position based on coords
# fill up the array and flush
$1 - coord > 10 {
    coord=$1;
    print("JUMP");
    print($0);
    next;
}
{   
    # move pointer up by one and flush if at 50
    coord=$1;
    print($0);
}
'