#!/bin/bash

## runs HTStream across a number all samples in a project
## assumes htstream is available on the Path
start=`date +%s`
echo $HOSTNAME

RRNA="References/human_rrna.fasta"
inpath="00-RawData"
outpath="01-HTS_Preproc"
[[ -d ${outpath} ]] || mkdir ${outpath}

for sample in `cat samples.txt`
do
  echo "SAMPLE: ${sample}"

  call="hts_Stats -L ${outpath}/${sample}.json -N 'initial stats' \
            -1 ${inpath}/${sample}_L008_R1_001.fastq.gz \
            -2 ${inpath}/${sample}_L008_R1_001.fastq.gz | \
        hts_SeqScreener -A ${outpath}/${sample}.json -N 'screen phix' | \
        hts_SeqScreener -A ${outpath}/${sample}.json -N 'count the number of rRNA reads'\
            -r -s $RRNA | \
        hts_SuperDeduper -A ${outpath}/${sample}.json -N 'remove PCR duplicates' | \
        hts_Overlapper -A ${outpath}/${sample}.json -N 'trim adapters' | \
        hts_PolyATTrim --no-left --skip_polyT  -A ${outpath}/${sample}.json -N 'remove polyAT tails' | \
        hts_NTrimmer -A ${outpath}/${sample}.json -N 'remove any remaining N characters' | \
        hts_QWindowTrim -A ${outpath}/${sample}.json -N 'quality trim the ends of reads' | \
        hts_LengthFilter -A ${outpath}/${sample}.json -N 'remove reads < 50bp' \
            -n -m 50 | \
        hts_Stats -A ${outpath}/${sample}.json -N 'final stats' \
            -f ${outpath}/${sample}"

  echo $call
  eval $call
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
