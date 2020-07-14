#!/bin/sh

mawk '
BEGIN {
    ################ VCF2CSV TRANSLATOR ###################

    ########## SPECS ###############
    # here go the standard 8 tab-separated fields as described in the VCF-specs
    vcfspex="Chr,Pos,ID,Ref,Alt,Qual,Filter,Info";
    # get the SPEX array to get the data from
    # SPEX[1] = "Chr"
    # SPEX[2] = "Pos"
    # --
    # SPEX[8] = "Info"

    spexCount = split(vcfspex,SPEX,",");
    # here come the actual output fields. omitting some specs cols
    # can be flexibly changed
    spexOut="Chr,Pos,ID,Ref,Alt";
    # split the outputSpex into the FIELDS array
    spexOutCount = split(spexOut,FIELDS,",");
    # FIELDS[1] = "Chr"
    # ..
    # FIELDS[5] = "Alt"

    ######### TAGS #####################
    # get extra tags from command argument $1 (with default for dbSNP)
    # provide in format Tag:Header
    # Header is used in HeaderColumns
    tags="'${1:-FREQ:AlleleFreq,VC:VariantClass}'";
    tagsCount = split(tags, TAGS, ",");

    fieldCount = spexOutCount + tagsCount;

    for (i=0; i++ < tagsCount;) {
        # split TagName and Header translation
        split(TAGS[i], TAG, ":");
        FIELDS[spexoutCount + i] = TAG[1];
        TAGHEADER[i] = TAG[2];
    }
    for (i = 0; i++ < fieldCount;) {
        print("test", TAGHEADER[i], FIELDS[i]);
    }

    # reverse the array for easy check
    for (i = 0; i++ < fieldCount;) {
      COLNUM[FIELDS[i]] = i;
    }
    # COLNUM["Chr"] = 1   ....

    ######## HEADER OUTPUT #############
    for (col = 0; col++ < spexOutCount;) {
      printf("%s\t",FIELDS[col]);
    }
    for (col = 0; col++ < tagsCount - 1;) {
        printf("%s\t",TAGHEADER[col]);
    }
    printf("%s\n",TAGHEADER[tagsCount]);
}



##### PROCESS LINES #################
/^[^#]/ {
    lines++;
    start=$2;
    R=$4
    A=$5
    RL=length($4);
    AL=length($5);

    # convert VCF tags to array Data
    Info=$8;
    In = split(Info,IF,";");
    for (i = 0; ++i <= In;) {
        split(IF[i],s,"=");
        # convert the SS field value 1,2,3 to Germline, Somatic or LOH
        if (s[1] == "SS") {
            if (s[2] == 1) {
                Data[s[1]] = "Germline";
            } else if (s[2] == 2) {
                Data[s[1]] = "Somatic";
            } else {
                Data[s[1]] = "LOH";
            }
        } else {
            Data[s[1]] = s[2];
        }
        
    }
    # IF is Info tag array
    # normal values
    Format=$9;
    n = split(Format,FFields,":");
    N=$10;   
    split(N,NFields,":");
    # tumor values
    T=$11;
    split(T,TFields,":");
    # get the data values into an array
    # N-DP4 value
    # T-DP4 value etc.
    for (i = 0; ++i <= n;) {
        if (FFields[i] == "DP4") {
            split(TFields[i],TDP,",");
            split(NFields[i],NDP,",");
            for (j=0; j++ < 5;) {
                Data["T-" DPF[j]] = TDP[j];
                Data["N-" DPF[j]] = NDP[j];
            }
        } else {
            Data["T-" FFields[i]]=TFields[i];
            Data["N-" FFields[i]]=NFields[i]; 
        }

    }
    # ######### DEBUG THE DATA FIELDS ########
    # for (i in Data) {
    #     printf("%s:%s;",i,Data[i])
    # }
    # print "\n";
    # ########################################

    ######### CHANGE INDEL COORDS ############
    # variant is an insertion/complex substitution 
    if (AL > 1) {
        if (RL == 1) { # an insertion
            end=start;
            # if first base is the same 
            # remove the redundant base in Ref and Alt
            if (R == substr(A,1)) {             
                R="-";
                A=substr(A,2);
            }
        } else { # is a complex substitution AC -> TG  
            end=start+RL-1;
        }
    # variant is a deletion
    } else if (RL > 1) {
        # is a deletion
        # make the pos of deleted bases to start and end
        start=start+1;
        end=start + RL -2;
        # remove the redundant base in Ref and Alt
        R=substr(R,2);
        A="-";
    } else { # it is a simple SNP
        end = start;
    }

    ######## OUTPUT #################
    printf("%s\t%s\t%s\t%s\t%s",$1,start,end,R,A);
    if (ALL == 1) {
        for (i in Data) {
            printf("\t%s",Data[i]);
        }
    } else {
        for (i = 4; i++ < 15;) {
        printf("\t%s",Data[VCF[HEADER[i]]]);
        }
    }
    printf("\n");

}'