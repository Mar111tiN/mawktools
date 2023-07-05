#!/bin/sh
####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # standard field output
        -s|--specs)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            spexOut=$2
            shift 2
        else
            echo "<vcf2csv> Error: standard specs output argument is missing\n[-s|--specs (default=Chr,Pos,Ref,Alt,Filter)]" >&2
            exit 1
        fi
        ;;
        # info tags
        -i|--info)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            itags=$2
            shift 2
        else
            echo "<vcf2csv> Error: info tag argument is missing\n[-i|--info (default="")]" >&2
            exit 1
        fi
        ;;
        # format tags
        -f|--format)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            ftags=$2
            shift 2
        else
            echo "<vcf2csv> Error: format tag argument is missing\n[-f|--format (default=DP:Depth,AF:VAF)]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<bed2index> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done

# only use the header
mawk '
BEGIN {
  ################ VCF2CSV TRANSLATOR ###################
  ########## SPECS ###############
  # here go the standard 8 tab-separated fields as described in the VCF-specs
  vcfspex="Chr,Pos,ID,Ref,Alt,Qual,Filter,Info";
  split(vcfspex,SPEX,",");
  # get the col number for the Spex
  for (col in SPEX) {
    SPEXCOL[SPEX[col]] = col;
    ###
    # SPEXCOL["Chr"] == 1
    ###
  }

  ### get input from args
  # here come the actual output fields. omitting some specs cols
  # can be flexibly changed
  spexOut="'${spexOut-"Chr,Pos,Ref,Alt,Filter"}'";
  # split the outputSpex into the FIELDS array
  spexCount = split(spexOut,FIELDS,",");

  ######### TAGS #####################
  # get extra tags from command argument $1 (with default) for format tags
  # and $2 for info tags

  # FTAGS ############
  # Tags to extract data from FORMAT column
  # provide in format Tag:Header,Tag[:Header],..)
  # optional Header is used in HeaderColumns
  # different defaults are: 
  common_all="G5A:MinorAlleleAllPop>5,G5:MinorAllele+1";
  dbSNP="FREQ:AlleleFreq,VC:VariantClass";
  varscan="SS:somaticStatus,GPV:variantP,SPV:somaticP,DP:Depth,RD:R1,AD:R2,DP4:readCounts";

  # from args
  ftags="'${ftags-"CDP:Depth,AF:VAF"}'";
  # split the ftags into arrays
  ftagsCount = split(ftags, FTAGS, ",");

  # split the format tags (ftags) into the FTAG array and translations into FTAGHEADER
  for (i=0; i++<ftagsCount;) {
    # split TagName and Header translation
    split(FTAGS[i], FTAG, ":");
    FTAGS[i]=FTAG[1];
    if (FTAG[2]) {
        FTAGNAME[FTAG[1]] = FTAG[2];
    } else {
        FTAGNAME[FTAG[1]] = FTAG[1];
    }

    ###
    # FTAGS[1] == "Freq"
    # FTAGHEADER["Freq"] == "AlleleFreq"]
    ###
  }

  # ITAGS #############
  # tags to extract from INFO column
  itags="'${itags:-""}'";
  # split the itags into arrays
  itagsCount = split(itags, ITAGS, ",");


  # split the info tags (itags) into the ITAG array and translations into ITAGHEADER
  for (i=0; i++<ftagsCount;) {
    # split TagName and Header translation
    split(ITAGS[i], ITAG, ":");
    ITAGS[i]=ITAG[1];
    if (ITAG[2]) {
        ITAGNAME[ITAG[1]] = ITAG[2];
    } else {
        ITAGNAME[ITAG[1]] = ITAG[1];
    }

    ###
    # ITAGS[1] == "Freq"
    # ITAGHEADER["Freq"] == "AlleleFreq"]
    ###
  }

}

readData { # only becomes active after the header scan
  # SPEX
  # first line extra
  printf("%s", $SPEXCOL[FIELDS[1]]);
  for (i=1; i++<spexCount;) {
    printf("\t%s",$SPEXCOL[FIELDS[i]]);
  }

  ######## INFO FIELDS ##########
  for (i=spexCount; i++< spexCount+IL;) {
    tag=FIELDS[i];
    tagLen=length(tag)+1;
    pattern="[\t;]" tag "(=[^\\t;$]+)?[\t;]";
    if (match($8, pattern)) {
      value=substr($8,RSTART+tagLen+1,RLENGTH-tagLen-2);
      if (value == "") {
        value="+";
      }
      printf("\t%s", value);
    } else {
      printf("\t-");
    }
  }
  ####### FORMAT FIELDS #########
  if (FL) { # FORMAT IS USED
    if (firstLine) { # find the order of the format tags
      fCount=split($9,FORMAT,":");
      for (i=0;i++<fCount;) {
        FCOL[FORMAT[i]]=i;
      }
      firstLine=0;
    }
    for (s=0; s++<sampleCount;) {
      for (i=spexCount+IL; i++<spexCount+IL+FL;) {
        split($(9+s),FDATA,":");
        printf("\t%s",FDATA[FCOL[FIELDS[i]]]);
      }
    }
  }
  printf("\n");
  next;
}


/^##/ { HEADER RUN
  # get INFO on FIELDS and alocate to info and format arrays

  if ((length(ftags) == 0) && (length(itags) == 0)) next; # skip if there are no tags to use
  for (i=0; i++<ftagsCount;) {
    tag=FTAGS[i];

    ftagLen=length(ftag);  ### ????
    # make regex pattern
    pattern="##FORMAT=<ID=" tag ",";
    if (match($0, pattern)) {
      FORMAT[++FL]=tag;
    } 
  }

  for (i=0; i++<itagsCount;) {
    tag=ITAGS[i];
    itagLen=length(tag);  ### ????
    # make regex pattern
    pattern="##INFO=<ID=" tag ",";
    if (match($0, pattern)) {
      INFO[++IL]=tag;
    } 
  }
  next;
}
/^#/ {
  # get the number and names of the samples defined in format
  for (col=9; col++<NF;) {
    sampleCount++;
    SAMPLES[++s]=$col;
  }

  # FL is switch for processing of FORMAT data in the read part
  # allocate the found tags to the FIELD array in order INFO --> FORMAT
  
  # INFO TAGS
  for (t=0;t++<IL;) {
    FIELDS[spexCount+t] = INFO[t]
  }
    
  for (f=0;f++<FL;) {
    for (s=0;s++<sampleCount;){
      fieldPos=spexCount+IL+f+(s-1)*FL;
      FIELDS[fieldPos]=FORMAT[f];
    }
  }
  ############ PRINT HEADER ##########################
  printf("%s",FIELDS[1]);
  for (i=1; i++<spexCount;) {
    printf("\t%s",FIELDS[i]);
  }

  if (IL) { # only if Info tags are requested
    for (i=spexCount; i++< spexCount+IL;) {
      printf("\t%s",ITAGNAME[FIELDS[i]]);
    }
  }

  if (FL) { # only if Format tags are requested
    for (s=0; s++<sampleCount;) {
      for (i=spexCount+IL; i++<spexCount+IL+FL;) {
        line=SAMPLES[s] "-" FTAGNAME[FIELDS[i]];
        printf("\t%s",line);
      }
    }
  }
  printf("\n");
  readData=1;
  firstLine=1;
}'