#!/bin/sh

# 4 CORES

# quick tool to separate mutation database lists for multiple alleles and DB entries
# this here is hard-coded for the output of DBSNP database files
# these files have to be preprocessed in the following steps:

# 1. recompress into bgzip file for processing with bcftools
# gunzip < dbSNPxxx.vcf.gz | ## unzip
# sed -e 's/NC_00000?([1-9XYM][0-9T])\.[1-9]+/chr\1/'   ## convert NC_00001.11 to chr1
#     -e 's/chr23/chrX/' -e 's/chr24/chrY/' |   ## convert chr23 to chrX etc
# bgzip > dbSNPxxx.vcf.gz ## rezip the file with bgzip

# 2. left-normalize the mutations (probably only works for single ALTs)
# [ tabix -p vcf dbSNPxxx.vcf.gz  ## index the file ] only neccessary for stepwise processing
# bcftools norm -f $HG38 < dbSNPxxx.vcf.gz | gzip > dbSNPxxx.ln.vcf.gz

# 3. convert the left-normalized vcf via vcf2csv into a csv file
# gunzip < dbSNPxxx.ln.vcf.gz | vcf2csv | expandDBSNP | pigz -p -4 > dbSNPxxx.ln.sep.vcf.gz


mawk '
######## STEP1 ########
# filter reads that have a freq
$6 != "." ' | mawk '

######## STEP2 ########
# split by DB
$6 ~ /\|/ { #FREQ field has several DBs
  dbCount = split($6,FREQ,"|");
  for (db=0;db++<dbCount;) {
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,FREQ[db],$7);
  }
  next;
}
{ # print all other lines
  print $0;
}' | mawk '

######### STEP3 #########
# separate the Freq Field
# $6 FREQ: "DB:RefFreq,AltFreq1,AltFreq2,..
# --> DB = DB1; RefFreq = RefFreq; AltFreq = AltFreq1,AltFreq2,..
NR == 1 { # HEADER ######
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,"RefFreq","AltFreq","DB",$7);
  next;
}
{
  split($6,AF1,":");
  db = AF1[1];
  freqs = AF1[2];
  # split the Alts into ref and alt
  split(freqs,AF2,",");
  ref = AF2[1];
  # remove the reffreq from the freqs
  refcomma = ref ",";
  sub(refcomma,"",freqs);
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,ref,freqs,db,$7);
}' | mawk '

# STEP4 ########## separate the alternative alleles
# 
$4 ~ /,/ { # Alt has separate alleles
  altCount = split($4,ALT,",");
  split($7,AF,",");
  for (alt=0;alt++<altCount;) {
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,ALT[alt],$5,$6,AF[alt],$8,$9);
  }
  next;
}
{ # print all other lines (including HEADER)
  print $0;
}' | mawk ' $7 != "."'
# # STEP5 ################### only print alleles with freq