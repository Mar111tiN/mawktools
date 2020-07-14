#!/bin/sh

# prints out the bam files with the base at a certain position
# used for MutDetect in 10x data
# outputs single reads that cover that mutation along with base and base quality
# USAGE
# $ samtools view bam_file [region] | bam2csv | bamPileup region

mawk '
NR == 1 {
  param="'$1'"; # takes chr:start-end as parameter for the mut position
  paramCount = split(param,params,":");
  chrom=params[1];
  pos=params[2];
  split(pos, startend, "-");
  mutPos=startend[2];
  cigPat = "^[0-9]+[NMDIS]";

  # get the field number to check as a criteria in the main data rows
  # some rows do not have CB value and thus NF < standard NF
  fields=NF;

  # get the coords for the fields
  # for flexible data structures
  n = split("SEQ,QUAL,CIGAR,START",FIELDS,",");
  for (col in FIELDS) {
    for (i = 0; i++ < NF;) {
      if ($i == FIELDS[col]) {
        COORD[FIELDS[col]]=i;
      }
    }
  } 
  printf "Piling up bam at %s:%s\n", chrom,mutPos > "/dev/stderr";
  # print Header
  printf("%s\t%s\t%s\n",$0, "QueryBase", "QueryQual");
}
# check for chromosome and existence of CB (with NF17)
($3 == chrom) && (NF == fields) {
  # those field numbers should be ok hard-coded because they are the standard bam fields
    seq=$COORD["SEQ"];
    qual=$COORD["QUAL"];
    cigar=$COORD["CIGAR"];
    # genStart is the genomic start position of the block
    genStart=$COORD["START"];

    # blockStart is the absolute start position in the sequence string
    blockStart = 1;
    # extract cigar info
    while (match(cigar, cigPat)) {
      # length of cigar block
      l = substr(cigar,RSTART,RLENGTH-1)
      # type of cigar block
      t = substr(cigar,RSTART+RLENGTH-1,1)

      # if match block, check if coords are in 
      if (t == "M") {
        # stepped over mutPos
        if (genStart > mutPos) {
          cigar = "";
          # -1 :: no coverage at that position
          hit = "-1";
          q = "-1";
        } else {
          # contained in block
          if (genStart + l >= mutPos) {
            seqPos = mutPos - genStart + blockStart;
            cigar = "";
            hit = substr(seq,seqPos,1);
            q = substr(qual,seqPos,1);
            # ####### DEBUG ###############
            # before = substr(seq,blockStart,seqPos-blockStart);

            # after = substr(seq,seqPos+1,l - seqPos + blockStart);
            # printf("(%s<<%s>>%s)",before, hit, after);
            #############################

          } else {
            # block is upstream of mutation
          # ####### DEBUG ###############
          # string=substr(seq,blockStart,l);
          # printf("(%s)", string)
          # #############################
            genStart += l;
            blockStart += l;   
          }
        }
      # jump down over intron gaps
      } else if (t ~ /[ND]/) {
        genStart += l;
      }  else if (t ~ /[SI]/) {
        blockStart += l;
      }
      
      # reduce cigar string
      cigar = substr(cigar, RSTART+RLENGTH);
    }
    if (hit != "-1") {
      printf("%s\t%s\t%s\n",$0,hit,q)
    }
 }'