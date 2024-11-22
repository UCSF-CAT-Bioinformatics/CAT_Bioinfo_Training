#!/bin/bash

## assumes star version 2.7.11b
## assumes STAR is available on the Path

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
FASTA="../GRCh38.primary_assembly.genome.fa.gz"


wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz
gunzip gencode.v47.basic.annotation.gtf.gz
GTF="../gencode.v47.basic.annotation.gtf"

mkdir star.overlap100.gencode.v47
cd star.overlap100.gencode.v47

call="STAR
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 100"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
