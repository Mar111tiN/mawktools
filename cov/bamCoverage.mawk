#!/fast/users/szyskam_c/work/miniconda/bin/mawk -f

BEGIN {
    # get the minCoverage
    minCov = (ARGV[1] == "") ? 0 : ARGV[1];
    # ARGV has to be reset to "" or ARGV[1] will be read as file
    ARGV[1] = "";
    cigPat = "^[0-9]+[NMDIS]";  
    # init pointer and coverage
    cov = 0;
    }
    
function flush( COVscan, cov, minCov, start,  i, sorted, len) {
    ########## OUTPUT PREVIOUS ##############
    # get intermediate array COValt for summing up and down 
    # and delete consumed up and down
    # if position is not passed, it will be 0 and a total flush will be performed
    for (p in COVscan) {
        # for some reason 99911600 is not smaller than 100000000
        # bypassing that problem
        if (p / 10 < start / 10 || start == 0) {
            COVsort[++i]=p;
            #########
            # printf("indexing: i=%s COVsort[%s]=%s p=%s\n", i, i, COVsort[i], p) >> "/dev/stderr";
            #########
        }
    }
    ######## SORT COValt #####
    len=length(COVsort);

    # keep sorting until sorted
    while (sorted != 1 && len>1) {
        # presume it is sorted
        sorted = 1;
        for (i=1;i++<len;) {
            #########
            # printf("Sorting: i=%s %s i=%s %s\n", i-1, COVsort[i-1], i, COVsort[i]) >> "/dev/stderr";
            #########

            # if first is smaller then last, swap
            if (COVsort[i-1] > COVsort[i]) {
                tmp = COVsort[i];
                COVsort[i] = COVsort[i-1];
                COVsort[i-1] = tmp;
                # set unsorted
                sorted = 0;
            }
        }
    }
    ########## PRINT OUT DATA #############
    for (i=0;i++<len;) {
        p = COVsort[i];
        # add the sorted coverage change to the running cov
        cov += COVscan[p];

        ########
        # print("COValt", p, COValt[p], "cov", cov)
        ########

        delete COVscan[p];
        if (cov >= minCov) {  # this makes it sparse output
            printf("%s\t%s\t%s\n", currentChrom, p,cov);     
        } 
    }
    delete COVsort;
    return cov; 
}
 
NR == 1 {
    currentChrom = $3;
    # print Header
    printf("%s\t%s\t%s\n","Chr","Pos","Coverage");
}

####### DATA FLUSH @ CHROM JUMP ########
$3 != currentChrom {
    print("chromFLUSH");
    flush(COVscan,cov, minCov);
    cov = 0;
    currentChrom = $3;
}

# look at all hits
$6 !~ /\*/ {
    # get the data
    chrom = $3;
    # I want chr everywhere
    if (chrom !~"chr") {
        chrom = "chr" chrom;
    }

    start = $4; # start of current view window
    cigar = $6;

    ###########
    # print("READ", start, cigar);
    ###########
    cov = flush(COVscan, cov, minCov, start);
    ########## CIGAR SCAN ###############
    readPos = start # readPos is current read position during cigar parse
    # get the start end end coords for the coverage from the cigar
    while (match(cigar, cigPat)) {
        # length of cigar block
        l = substr(cigar,RSTART,RLENGTH-1)
        # type of cigar block
        t = substr(cigar,RSTART+RLENGTH-1,1)
        # if match block, push coords into coverage array
        if (t == "M") {
        # with every read, total coverage increases by one between start and end
        COVscan[readPos] ++;
        readPos += l; # readPos is moved up by the length of the Mcigar

        #####
        # print(start, " +1 - ", readPos, " -1");
        #####

        COVscan[readPos] --; # coverage goes down at new start pos
        # jump down over intron gaps
        } else if (t ~ /[ND]/) {
        readPos += l; # N and D result in no coverage
        }
        # reduce cigar string
        cigar = substr(cigar, RSTART+RLENGTH);
    }
  # get the current window's end position
}

END { # spill out all the remaining data arrays
    # print("END FLUSH") >> "/dev/stderr";
    flush(COVscan, cov, minCov);
}