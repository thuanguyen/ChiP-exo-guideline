#!/bin/sh
WHERE_BOWTIE='bowtie'
WHERE_MAKEGFF='makegff.py'
echo '#bowtie in' $WHERE_BOWTIE
echo '#makegff.py in' $WHERE_MAKEGFF

FILE_PE_LEFT='../fastq/yiep1.fastq'
FILE_NAME='yiep1_align'

FILE_OUTPUT='../alignment/'$FILE_NAME
FILE_GFF='../gff/'$FILE_NAME'.gff'

BOWTIE_OPTION='-p 8 -S'
BOWTIE_REFERENCE='../reference/ecoli_mg1655'
BOWTIE_OUTPUT=$FILE_OUTPUT'.sam'
BOWTIE_UNALIGNED=$FILE_OUTPUT'_unaligned.fastq'

echo '#command line: ' $WHERE_BOWTIE $BOWTIE_OPTION $BOWTIE_REFERENCE $FILE_PE_LEFT '--un' $BOWTIE_UNALIGNED '>' $BOWTIE_OUTPUT
$WHERE_BOWTIE $BOWTIE_OPTION $BOWTIE_REFERENCE $FILE_PE_LEFT --un $BOWTIE_UNALIGNED > $BOWTIE_OUTPUT

echo '#bowtie run is done.'

SAM_BAM_UNSORTED=$FILE_OUTPUT'.unsorted.bam'
SAM_BAM_SORTED=$FILE_OUTPUT

echo '#command line: ' samtools view '-bS' $BOWTIE_OUTPUT '-o' $SAM_BAM_UNSORTED
samtools view -bS $BOWTIE_OUTPUT -o $SAM_BAM_UNSORTED

echo '#command line: ' samtools sort $SAM_BAM_UNSORTED $SAM_BAM_SORTED
samtools sort $SAM_BAM_UNSORTED $SAM_BAM_SORTED
rm $SAM_BAM_UNSORTED

echo '#command line: ' samtools index $SAM_BAM_SORTED'.bam'
samtools index $SAM_BAM_SORTED'.bam'

echo '#command line: ' 'samtools index' $SAM_BAM_SORTED'.bam'
samtools index $SAM_BAM_SORTED'.bam'

echo '#command line: ' python makegff.py --separate_strand $SAM_BAM_SORTED'.bam' $FILE_GFF
python makegff.py --separate_strand $SAM_BAM_SORTED'.bam' $FILE_GFF
