#!/bin/bash

inp=$1
samout="${inp}_temp_samsorted.bam"
samtools sort -@ 4 $inp -o $samout
ARGout="${inp}__temp_sort_addgroup.bam"
gatk AddOrReplaceReadGroups --INPUT $samout  --OUTPUT $ARGout --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM 20 --SORT_ORDER coordinate
MDout="${inp}_temp_sort_addgroup_rawdedupped.bam"
MDmetrix="${inp}_temp_sort_addgroup_rawdedupped_metrics.bam"
gatk MarkDuplicates --INPUT $ARGout --OUTPUT $MDout --METRICS_FILE $MDmetrix
rm $samout
rm $ARGout
SNCout="${inp}_temp_sort_addgroup_dedupped_splitN.bam"
gatk SplitNCigarReads -R ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I $MDout -O $SNCout
rm $MDmetrix
rm $MDout
varfile="vars/${inp}_gatk_variants.vcf"
gatk HaplotypeCaller -R ~/gary/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa -I $SNCout --dont-use-soft-clipped-bases  --standard-min-confidence-threshold-for-calling 20.0 -O $varfile
rm $SNCout
rm *_splitN.*
