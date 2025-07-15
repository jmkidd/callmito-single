# callmito-single

Pipeline for calling canine mitochondrial variation from single-end Illumina data.

This is a modification of the [callmito](https://github.com/jmkidd/callmito) pipeline that
works with single-end Illumina data, as in ancient DNA applications.  It uses bwa aln/samse
to process the reads.

# Summary

This pipeline calls mitochondrial variation from Illumina data. 

The process starts with aligning single-end reads to a genome file that contains
the standard and rotated mitochondrial reference genomes.  This accounts for reads
that align across the end of the circular mitochondria. The resulting bam file is then 
given as input to callmito-single/process-sample.py which employs the same basic steps
as [callmito](https://github.com/jmkidd/callmito) but for single end reads. To accomodate
lower coverage data, regions with a depth less than 3 are masked to N.

# Example workflow

## Step 1: Align reads to mitochondrial genomes

```
bwa aln -l 1024 -n 0.01 -o 2 callmito-single/refs/mito-both.fa READ.fastq.gz > READ.sai 
bwa samse callmito-single/refs/mito-both.fa READ.sai  READ.fastq.gz | samtools view -F 4 -h -b - | samtools sort - > SAMPLE.bam 
```

## Step 2: Run callmito-single pipeline

```

python callmito-single/process-sample.py \
--ref callmito-single/refs/mito-both.fa  \
--finaldir OUTPUT-DIR/  \
--name SAMPLE \
--cram SAMPLE.bam \
--coords callmito-single/to-extract.bed \
--mitoFa callmito-single/refs/NC_002008.4.fa  \
--mitoFaRotated callmito-single/refs/NC_002008.4.rotate8k.fa \
--chainfile callmito-single/refs/rotatedToOriginal.liftOver \
--diagnosticTable callmito-single/fregel-haplogroups.txt 

```

# Software required

The following software and versions are used
```
bwa aln/samse (we used version 0.7.17)
gatk version 4.2.5.0
samtools version >= 1.9
liftOver
bgzip
tabix
bcftools version >= 1.9.
```