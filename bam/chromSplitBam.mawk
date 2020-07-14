#!/bin/sh

# takes one bam file as input and splits the bam file into chromosome-split bams into provided folder
# DEPENDENCIES:
# -mawk
# -samtools

input=$1;
folder=$2;
mkdir -p $folder;
for i in {1..22}
    do samtools view -H $input > $folder/chr${i}.sam;
done;

for i in {X,Y}
    do samtools view -H $input > $folder/chr${i}.sam;
done;

samtools view $input | mawk '
BEGIN {
    current=1;
    for (x = 0; x++ < 24;) {
        if (x < 23) {
            Chrom[x] = x;
    } else {
        if (x == 23) {
            Chrom[x] = "X";
        } else {
            Chrom[x] = "Y";
        }
    }
}
}
{   
    if ($3 ~ Chrom[current]) {
        {print >> "chr"Chrom[current]".sam"};
    } else {
        if ($3 ~ Chrom[current+1]) {
            {print >> "chr"Chrom[current+1]".sam"};
            current+=1;
        }
        else {
            {print >> "chrOther.sam"};
        }
    }
}
END {
    for (x = 0; x++ < 24;) {
        printf("%s--", Chrom[x]);
        }
    printf("\n");
}' 

for i in {1..22}
    do samtools view -Sb chr${i}.sam > ${folder}/chr${i}.bam;
    samtools index ${folder}/chr${i}.bam;
    rm chr${i}.sam;
done;

for i in {X,Y}
    do samtools view -Sb chr${i}.sam > ${folder}/chr${i}.bam;
    samtools index ${folder}/chr${i}.bam;
    rm chr${i}.sam;
done;