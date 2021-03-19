#!/bin/sh

# only use the header
mawk '
NR == 1 { in header
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
    spexOut="Chr,Pos,Ref,Alt,Filter";
    # split the outputSpex into the FIELDS array
    spexOutCount = split(spexOut,FIELDS,",");
    # FIELDS[1] = "Chr"
    # ..
    # FIELDS[5] = "Alt"

    ######### TAGS #####################
    # get extra tags from command argument $1 (with default for dbSNP)
    # provide in format Tag:Header,Tag[:Header],..)
    # optional Header is used in HeaderColumns
    common_all="G5A:MinorAlleleAllPop>5,G5:MinorAllele+1";
    dbSNP="FREQ:AlleleFreq,VC:VariantClass";
    varscan="SS:somaticStatus,GPV:variantP,SPV:somaticP";

    tags=varscan;
    
    # provide the tags and translations for the format cols
    format="RD:R1,AD:R2,DP4:readCounts";

    tagsCount = split(tags, TAGS, ",");

    fieldCount = spexOutCount + tagsCount;

    for (i=0; i++ < tagsCount;) {
        # split TagName and Header translation
        split(TAGS[i], TAG, ":");
        FIELDS[spexOutCount + i] = TAG[1];
        # if no header translation is given, use the field name itself
        if (TAG[2]) {
            TAGHEADER[i] = TAG[2];
        } else {
            TAGHEADER[i] = TAG[1];
        }
        
    }

    # invert the array for easy check
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
{
    # get normal fields with COLNUM check
    # SPEX stores the original location acc. to specs
    # 
    for (col = 0; col++ < spexCount;) {
      if (SPEX[col] in COLNUM) {
        Data[SPEX[col]]=$col;
      }
    }
    # get TAGS 
    for (col = spexOutCount; col++ < fieldCount;) {
      # create the pattern from the DataField and the value format field  (like FREQ=0.2342)
      tag=FIELDS[col];
      tagLen=length(tag)+1;
      pattern=tag "=[^\\t;$]+";
      if (match($0, pattern)) {
        Data[tag]=substr($0,RSTART+tagLen,RLENGTH-tagLen);
      } else {
        Data[tag]=".";
      }
    }

    ######## OUTPUT #################
    for (col = 0; col++ < fieldCount - 1;) {
      if (FIELDS[col] in COLNUM) {
        printf("%s\t",Data[FIELDS[col]]);
      }
    }
    printf("%s\n",Data[FIELDS[fieldCount]]);

}'