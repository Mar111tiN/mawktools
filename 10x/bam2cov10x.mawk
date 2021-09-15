#!/bin/sh

# extracts the CB UB reads covering a given transcript via ID as $1
# outputs CB UB TransID Start End and strand of transcript
 
# USAGE:
# $ samtools view bam_file region | 10xbam2cov TransID + 
# stdout can be piped into 10xcov2pileup

mawk '
BEGI {
    ################ BAM2CSV TRANSLATOR ###################

    # get strand orientation of transcript for adjusting Cigar location

    transOrient = "'$2'";
    # define the bam Tags for the output:
    cols = "CB,UB"

    # split the cols string into the DATAFIELDS array
    fieldCount = split(cols,FIELDS,",");
    # FIELDS[1] = "CB"
    # FIELDS[2] = "UB
    # FIELDS[3] = "TX"
    # FIELDS[4] = "AN"
    # set the pattern for the coverage string in the TX or AN tag
    ENSTpattern="'$1'" ",[^;]+[;\\t]";
    ######## HEADER OUTPUT #############
    for (col = 0; col++ < fieldCount;) {
      # now I can check for cols
      printf("%s\t",FIELDS[col]);
    }
    printf("TransID\tStart\tEnd\tStrand\n");
}

##### PROCESS LINES #################
 $0 ~ "CB:Z:" {  # process only lines containing the TransID
    # get TAGS 
    for (col=0; col++ < fieldCount;) {
      # create the pattern from the DataField and the value format field  (like CB:Z:ATTGCCTG)
      pattern=FIELDS[col] ":[ZiB]:[^\\t]+\\t?[$]?";
      if (match($0, pattern)) {
        Data[col]=substr($0,RSTART+5,RLENGTH-6);
      } else {
        Data[col]=".";
      }
    }

    ######## OUTPUT #################
    for (col = 0; col++ < fieldCount;) {
      printf("%s\t",Data[col]);
    }
    if (match($0, ENSTpattern)) {
      coverage=substr($0,RSTART,RLENGTH-1);
      split(coverage,COV,",");
      enst=COV[1];
      strand=substr(COV[2],1,1);
      pos=substr(COV[2],2);

      # get the Mcigar from the cigar
      match(COV[3], "[0-9]+M");
      len=substr(COV[3],RSTART,RLENGTH-1);

      # !!!!!! HAVE TO CHECK THE POSITIONS IN transcript with transorient "+"
      # works!

      if (strand == "+") {
        start=pos;
        end=pos+len;
      } else {
        start=(pos-len > 1) ? pos-len : 1;
        end=pos;
      }
    }
    printf("%s\t%s\t%s\t%s\n",enst, start, end, strand);
}' | sort -k4n

