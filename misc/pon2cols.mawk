#!/bin/sh
# creates input for a cut command depending on the length of a sample file input as stream
# only the Chr Start and respective reads of the samples are returned by the cut step
# used in a samtools mpileup workflow as follows:
# samtools mpileup bam ... | cut -f $(pon2cols < sample_list) | cleanpileupReadsOnly

mawk '
END {
    printf("%s,%s",1,2);
    for (i=0; i++<NR;) {
        printf(",%s",2+3*i)
    }
}'