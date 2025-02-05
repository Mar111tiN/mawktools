#!/bin/bash

#v2.1

# tool to transpose vcf files into tab-separated tsv files
# dynamically takes arguments to shape output parameters
# output from FORMAT fields is dynamically split per sample (sample name is taken from header)
# output from INFO field is auto-detected as Flag and converted to +/- output
# v2.1>> output from FUNC field (FUNC=[{<'FIELD':'VALUE'}]) is extracted with -F/--func option
# v2.2>> add -S flag for simple sample names (1_<FIELD>, .... 2_<FIELD>)

# USAGE: 
# cat file.vcf | vcf2csv 
#    -s Chr,Pos,Ref,Alt 
#    -f DP:Depth,AD:Alleles
#    -i SB:Strand,MMQ,MBQ
#    -F gene:Gene,coding:Nchange,protein:AAchange
#    -S
# [ -s | --specs]     <FIELD,[FIELD,...]>                             these fields of the official vcf specs are output unchanged. Use snakecase colnames (Chr,Pos,ID,Ref,Alt,Qual,Filter)                                                                          ]
# [ -i | --info       <FIELD[:FieldName],[FIELD[:FieldName],...]>     these tags of the info field are extracted and output in column FIELD (opt. FieldName)                                                ]
# [ -F | --FUNC       <FIELD[:FieldName],[FIELD[:FieldName],...]>     these tags of the FUNC field are extracted from the sample-value field and output in column FIELD (opt. FieldName)                    ]
# [ -f | --format     <FIELD[:FieldName],[FIELD[:FieldName],...]>     these tags of the format field are extracted from the sample-value field and output in column Sample_FIELD (opt. Sample_FieldName)    ]
# [ -S | --sss | simple_sample_names]     FLAG                        if set, the sample names in the data header are not used for naming the format columns (for >1 samples, <FIELD>_N is used)            ]

####### ARGPARSE ##################
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        ### FLAG for simple sample names
        -S|--simple_sample_names|--sss)
        simpleNames=1;
        shift
        ;;
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
            echo "<vcf2csv> Error: info tag argument is missing\n[-i|--info (default=\x27\x27)]" >&2
            exit 1
        fi
        ;;
        # format tags
        -f|--format)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            ftags=$2
            shift 2
        else
          shift 2;
            echo "<vcf2csv> Error: format tag argument is missing\n[-f|--format (default=\x27\x27)]" >&2
            exit 1
        fi
        ;;
        # FUNC tags  # for FUNC=[{'FIELD':'VALUE'}] annotation
        -F|--FUNC)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            FUNCtags=$2
            shift 2
        else
          shift 2;
            echo "<vcf2csv> Error: FUNC tag argument is missing\n[-F|--FUNC (default=\x27\x27)]" >&2
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
    # SPEXCOL["Pos"] == 2
    # ...
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
  # different presets are: 
  common_all="G5A:MinorAlleleAllPop>5,G5:MinorAllele+1";
  dbSNP="FREQ:AlleleFreq,VC:VariantClass";
  varscan="SS:somaticStatus,GPV:variantP,SPV:somaticP,DP:Depth,RD:R1,AD:R2,DP4:readCounts";
  Mutect2="DP:Depth,AD:Alleles,F1R2:Fwd,F2R1:Rwd,SB:Strand";
  Georg="AF:VAF,AO:TDepth,DP:Depth,GQ:GenoQual,SAF:TR1,SAR:TR2,SRF:NR1,SRR:NR2";
  other="CDP:Depth,AF:VAF"

  # from args
  ftags="'${ftags-""}'";

  # split the ftags into arrays
  FCount = split(ftags, FTAGS, ",");
  # split the format tags (ftags) into the FTAG array and translations into FTAGNAME
  for (i=0; i++<FCount;) {
    # split TagName and Header translation
    split(FTAGS[i], FTAG, ":");
    FTAGS[i]=FTAG[1];
    if (FTAG[2]) {
        FTAGNAME[FTAG[1]] = FTAG[2];
    } else {
        FTAGNAME[FTAG[1]] = FTAG[1];
    }
    
    ###
    # FTAGS[1] == "DP"   FTAGNAME["DP"] == "Depth"
    # FTAGS[2] == "AD"   FTAGNAME["AD"] == "Alleles"
    ###
  }

  # ITAGS #############
  # tags to extract from INFO column
  # presets:
  # Mutect2: "PON:PoN,SB:Strand,MMQ,MBQ"

  itags="'${itags:-""}'";
  # split the itags into arrays
  ICount = split(itags, ITAGS, ",");
  # split the info tags (itags) into the ITAG array and translations into ITAGNAME
  for (i=0; i++<ICount;) {
    # split TagName and Header translation
    split(ITAGS[i], ITAG, ":");
    ITAGS[i]=ITAG[1];
    if (ITAG[2]) {
        ITAGNAME[ITAG[1]] = ITAG[2];
    } else {
        ITAGNAME[ITAG[1]] = ITAG[1];
    }

    ###
    # ITAGS[1] == "PON"   ITAGHEADER["PON"] == "PoN"
    # ITAGS[2] == "SB"   ITAGHEADER["SB"] == "Strand"
    ###
  }

  # FUNCTAGS #############
  # tags to extract from FUNC column
  # presets:
  # Georg: "gene:Gene,coding:Nchange,protein:AAchange,transcript:TX"

  FUNCtags="'${FUNCtags:-""}'";
  # split the FUNCtags into arrays and set FUNCount variable
  FUNCount = split(FUNCtags, FUNCTAGS, ",");
  # split the FUNC tags (FUNCtags) into the FUNCTAG array and translations into FUNCTAGNAME
  for (i=0; i++<FUNCount;) {
    # split TagName and Header translation
    split(FUNCTAGS[i], FUNCTAG, ":");
    FUNCTAGS[i]=FUNCTAG[1];
    if (FUNCTAG[2]) {
        FUNCTAGNAME[FUNCTAG[1]] = FUNCTAG[2];
    } else {
        FUNCTAGNAME[FUNCTAG[1]] = FUNCTAG[1];
    }

    ###
    # FUNCTAGS[1] == "gene"   FUNCTAGNAME["gene"] == "Gene"
    # FUNCTAGS[2] == "coding"   FUNCTAGNAME["coding"] == "Nchange"
    ###
  }
}

readData { # only becomes active after the header scan
  # SPEX
  # first line extra
  printf("%s", $SPEXCOL[FIELDS[1]]);
  fieldPos=1;
  for (i=1; i++<spexCount;) {
    printf("\t%s",$SPEXCOL[FIELDS[++fieldPos]]);
  }
  ####### Extracting FUNC ###############

  tag="FUNC";
  tagLen=length(tag)+1;
  pattern="[\t;]" tag "=\[\{[^\\t;$\\]}]+\}\]";

  if (match($8, pattern)) {
    # storing FUNC key values in func variable and stripping quotes
    # + 3 comes from =[{  , -5 comes from removing the 3 again and another }]
    func=substr($8,RSTART+tagLen+3,RLENGTH-tagLen-5);
    gsub("\x27", "", func);
    # removing FUNC from INFO variable (necessary?)
    $8=substr($8,1,RSTART-1) substr($8,RSTART+RLENGTH);
  }


  ######## INFO FIELDS ##########
  for (i=0; i++< IL;) {
    tag=FIELDS[++fieldPos];
    pattern="(;|^)" tag "(=[^\\t;]+)?(;|$)";
    if (match($8, pattern)) {
      startSkip=length(tag)+1;
      # skip the ; if at start
      if (substr($8,RSTART,1) == ";") startSkip++;
      endSkip=startSkip;
      # skip the ; at the end
      if (substr($8,RSTART+RLENGTH-1,1) == ";") endSkip++;
      value=substr($8,RSTART+startSkip,RLENGTH-endSkip);
      if (value == "") {
        value="+";
      }
      printf("\t%s", value);
    } else {
      printf("\t-");
    }
  }

  ######## FUNC FIELDS ##########
  for (i=0; i++< FUNCount;) {
    tag=FIELDS[++fieldPos];

    pattern=tag ":[^:]+[,$]";
    # print(pattern);
    # print($8);
    if (match(func, pattern)) {
      startSkip=length(tag)+1;
      value=substr(func,RSTART+startSkip,RLENGTH-startSkip-1);
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
      split($(9+s),FDATA,":");
      for (i=0; i++<FL;) {
        printf("\t%s",FDATA[FCOL[FIELDS[++fieldPos]]]);
      }
    }
  }

  printf("\n");
  next;
}

/^##/ { HEADER RUN

  ### get tumor and normal sample into variables
    # these are stored in array eg: TN[234]="N_" TN[234]="T_"
    if ($0 ~ /^##normal_sample=/) {
      split($0, a, "=")
      N = a[2]
      # 
      TN[N]="N_"
  } else if ($0 ~ /^##tumor_sample=/) {
      split($0, a, "=")
      T = a[2]
      TN[T]="T_"
  }

  ### get INFO on FIELDS and alocate to info and format arrays
  if (!(FCount || ICount)) next; # skip if there are no tags to use

  # check if FORMAT TAGS are defined and if so, add them to the FORMAT array
  for (i=0; i++<FCount;) {
    tag=FTAGS[i];
    # make regex pattern
    pattern="##FORMAT=<ID=" tag ",";
    if (match($0, pattern)) {
      FORMAT[++FL]=tag;
    } 
  }
  # check if INFO TAGS are defined and if so, add them to the INFO array
  for (i=0; i++<ICount;) {
    tag=ITAGS[i];
    # make regex pattern
    pattern="##INFO=<ID=" tag ",";
    if (match($0, pattern)) {
      INFO[++IL]=tag;
    } 
  }
  next;
}
/^#/ { ###### reached the DATA column header of the VCF file
  # get the number and names of the samples defined in format

  # get Sample Flag from args
  simpleNames='${simpleNames-0}'
  for (col=9; col++<NF;) {  # count columns after FORMAT field
    sampleCount++;
    if (simpleNames) {
      SAMPLES[++s]=s "_";
    } else {
      SAMPLES[++s]=$col "_";  # get the sample names
    }
  }
  # remove name for single sample in simpleSample mode
  if (simpleNames && (sampleCount ==1)) SAMPLES[1]="";

  if (simpleNames && (sampleCount ==2) && N != "") {
    # exchange the tumor/normal numerics with T and N
    SAMPLES[1]=TN[$10];
    SAMPLES[2]=TN[$11];
    }

  # FL is switch for processing of FORMAT data in the read part
  # allocate the found tags to the FIELD array in order INFO --> FUNC --> FORMAT
  
  # INFO TAGS
  # set the fieldPos variable
  fieldPos = spexCount
  for (t=0;t++<IL;) {
    FIELDS[++fieldPos] = INFO[t]
  }
  
  for (f=0;f++<FUNCount;) {
    FIELDS[++fieldPos] = FUNCTAGS[f];
  }

  # FORMAT TAGS
  for (s=0;s++<sampleCount;){
    for (f=0;f++<FL;) {
      FIELDS[++fieldPos]=FORMAT[f];
    }
  }
  ### DEBUG
  # for (f=0;f++<fieldPos;) print(FIELDS[f])
  ####

  ############ PRINT HEADER ##########################
  printf("%s",FIELDS[1]);
  printPos=1;
  for (i=1; i++<spexCount;) {
    printf("\t%s",FIELDS[++printPos]);
  }

  for (i=0; i++<IL;) {
    printf("\t%s",ITAGNAME[FIELDS[++printPos]]);
  }

  for (i=0; i++<FUNCount;) {
    printf("\t%s",FUNCTAGNAME[FIELDS[++printPos]]);
  }

  for (s=0; s++<sampleCount;) {
    for (i=0; i++<FL;) {
      line=SAMPLES[s] FTAGNAME[FIELDS[++printPos]]
      printf("\t%s",line);
    }
  }

  printf("\n");
  readData=1;
  firstLine=1;
}'
