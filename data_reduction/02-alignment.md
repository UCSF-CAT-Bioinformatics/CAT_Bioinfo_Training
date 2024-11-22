# Alignment to Read Counts & Visualization in IGV

## **Learning Objectives**
1. **Distinguish Key Concepts in Mapping and Assembly:**
   - Explain the differences between mapping (alignment to a reference) and assembly (constructing sequences without a reference).
   - Identify the challenges and solutions associated with mapping RNA-Seq data, such as repetitive regions, splicing, and genetic variation.

2. **Comprehend Essential Alignment Metrics:**
   - Define and interpret multimappers, duplicates, soft-clipping, and splicing.
   - Evaluate mapping quality and understand its impact on downstream analyses.

3. **Utilize RNA-Seq Alignment Tools:**
   - Run RNA-Seq alignment using STAR, including setting appropriate parameters for counting reads at the gene level.
   - Generate and interpret mapping statistics using custom scripts (e.g., `star_stats.sh`).

4. **Understand Reference Genome Requirements:**
   - Select appropriate genome FASTA and GTF files for RNA-Seq mapping.
   - Identify the importance of matching chromosome names and including organelles or additional sequences if needed.

5. **Apply Tools for Data Quality Assessment:**
   - Use mapping statistics to ensure alignment quality across samples.
   - Troubleshoot common errors in alignment outputs.

6. **Analyze and Visualize Alignment Outputs:**
   - Load and explore BAM files in IGV to inspect gene-level alignments.
   - Identify patterns in aligned reads and relate them to experimental outcomes (e.g., high expression regions).

7. **Prepare Data for Differential Expression Analysis:**
   - Generate gene count matrices from RNA-Seq alignments.
   - Compare different counting strategies (e.g., union, intersection-strict) and their implications for downstream analyses.

8. **Collaborate in Problem-Solving Exercises:**
   - Interpret STAR alignment results and mapping summaries collaboratively.
   - Investigate specific genomic regions using IGV to link experimental observations to biological insights.


---
## Mapping vs Assembly

**Assembly** seeks to put together the puzzle without knowing what the picture is.

- The focus is on the pieces, how they fit together.

**Mapping** (or alignment to a reference) tries to put the puzzle pieces directly onto an image of the picture._
- The focus is on the puzzle, regions of the puzzle that contain certain characteristics (ex. what background) that will help you place the piece onto the puzzle.  
- In mapping the question is more, given a small chunk of sequence, where in the genome did this sequence most likely come from.
- The goal then is to find the match(es) with either the “best” edit distance (smallest difference), or all matches with edit distance less than max edit distance. Main issues are:
* Large search space
* Regions of similarity (aka repeats)
* Gaps (INDELS)
* Complexity (RNA, splicing, transcripts)

### Alignment concepts

* Multimappers:
  * Reads that align equally well to more than one reference location.
  * Generally, multimappers are discounted in variant detection, and are often discounted in counting applications (like RNA-Seq ... would “cancel” out anyway).
  * Note: multimapper “rescue” in some algorithms (RSEM, eXpress).
* Duplicates:
  * Reads or read pairs arising from the same original library fragment, either during library preparation (PCR duplicates) or during sequencing (optical duplicates).
  * Generally, duplicates can only be detected reliably with paired-end sequencing. If PE, they’re discounted in variant detection, and discounted in counting applications (like RNA-Seq).
* Clipping vs Splicing  
  * soft-clipped: bases in 5' and 3' of the read are NOT part of the alignment.
  * hard-clipped: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been
  removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
  * [From biostars](https://www.biostars.org/p/119537/)

<img style="padding-left:100px" src="alignment_figures/alignment_figure3.png" alt="alignment_figure3" width="60%"/>  

* Inner length, insert size, fragment length  

<img src="alignment_figures/alignment_figure4.jpg" alt="alignment_figure4" width="60%"/>  
*From [This Biostars answer](https://www.biostars.org/p/106291/)*


#### Considerations when mapping
* Placing reads in regions that do not exist in the reference genome (reads extend off the end of linearized fragments) [ mitochondrial, plasmids, structural variants, etc.].
* Sequencing errors and genetic variation: alignment between read and true source in genome may have more differences than alignment with some other copy of repeat.
* What if the closest fully sequenced genome is too divergent?
* Placing reads in repetitive regions: Some algorithms only return 1 mapping; If multiple: map quality = 0
* Algorithms that use paired-end information => might prefer correct distance over correct alignment.

In RNAseq data, you must also consider effect of splice junctions, reads may span an intron.

<img src="alignment_figures/alignment_figure1.png" alt="alignment_figure1" width="80%"/>

<iframe width="80%" src="https://www.youtube.com/embed/_asGjfCTLNE" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

---
## Aligners/Mappers

Many [alignment algorithms](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software
) to choose from. Examples include:
* Spliced Aligners
  * STAR
  * HiSAT2 (formerly Tophat [Bowtie2])
  * GMAP/GSNAP
  * SOAPsplice
  * MapSplice
* Aligners that can ’clip’
  * bwa-mem
  * Bowtie2 in local mode

#### Pseudo-aligners (salmon and kalisto)

* Quasi-mapping
* Probabilistic
* Map to transcripts, not genome
* Does transcript quantifications (or gene)
* Blazing FAST and can run on most laptops
* Experience suggests differences between “traditional” mappers are in the low abundance genes.

---

## Mapping against the genome vs transcriptome

May seem intuitive to map RNAseq data to transcriptome, but it is not that simple.
  * Transcriptomes are rarely complete.
  * Which transcript of a gene should you map to? canonical transcript (which is that)?
  * Shouldn’t map to all splice variants as these would show up as ‘multi-mappers’.

More so, a aligner will try to map every read, somewhere, provided the alignment meets its minimum requirements.
  * Need to provide a mapper with all possible places the read could have arisen from, which is best represented by the genome. Otherwise, you get mis-mapping because its close enough.

### Genome and Genome Annotation

Genome sequence fasta files and annotation (gff, gtf) files go together! These should be identified at the beginning of analysis.

* Genome fasta files should include all primary chromosomes, unplaced sequences and un-localized sequences, as well as any organelles. Should not contain any contigs that represent patches, or alternative haplotypes.
* If you expect contamination, or the presence of additional sequence/genome, add the sequence(s) to the genome fasta file.
* Annotation file should be GTF (preferred), and should be the most comprehensive you can find.
  * Chromosome names in the GTF must match those in the fasta file (they don’t always do).

---
## Counting reads as a proxy for gene expression

The more you can count (and HTS sequencing systems can count a lot) the better the measure of copy number for even rare transcripts in a population.
* Most RNA-seq techniques deal with count data. Reads are mapped to a reference genome, transcripts are detected, and the number of reads that map to a transcript (or gene) are counted.
* Read counts for a transcript are roughly proportional to the gene’s length and transcript abundance (in whole transcript methods).

Technical artifacts should be considered during counting
* Mapping quality
* Map-ability (uniqueness), the read is not ambiguous

Options are (STAR, HTSEQ, featureCounts)

### The HTSEQ way

* Given a sam/bam file with aligned sequence reads and a list of genomic feature (genes locations), we wish to count the number of reads (fragments) than overlap each feature.
  * Features are defined by intervals, they have a start and stop position on a chromosome.
  * For this workshop and analysis, features are genes which are the union of all its exons. You could consider each exon as a feature, for alternative splicing.
* Htseq-count has three overlapping modes
  * union
  * intersection-strict
  * intersection-nonempty

<img src="alignment_figures/alignment_figure2.png" alt="alignment_figure2" width="60%"/>
*from the [HTSeq Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4287950/)*

#### Star Implementation
Counts coincide with Htseq-counts under default parameters (union and tests all orientations). Need to specify GTF file at genome generation step or during mapping.
* Output, 4 columns
  * GeneID
  * Counts for unstranded
  * Counts for first read strand
  * Counts for second read strand

<div class="output">

N_unmapped	213761	213761	213761
N_multimapping	132491	132491	132491
N_noFeature	309643	2619257	322976
N_ambiguous	116750	892	47031
ENSMUSG00000102693	0	0	0
ENSMUSG00000064842	0	0	0
ENSMUSG00000051951	0	0	0
ENSMUSG00000102851	0	0	0
ENSMUSG00000103377	0	0	0
ENSMUSG00000104017	0	0	0
</div>

Choose the appropriate column given the library preparation characteristics and generate a matrix expression table, columns are samples, rows are genes.

What does stranded and unstranded mean? Which is better and why? [Stranded vs Unstranded](https://www.ecseq.com/support/ngs/how-do-strand-specific-sequencing-protocols-work)

---
## Initial Setup

*This document assumes [preproc htstream](./02-alignment-indexref.md) has been completed.*
To catch up to where we are:

```
cd /mnt/analysis/cat_users/$USER/rnaseq_example

refcheck=$(egrep "DONE: Genome generation" ./References/star.overlap100.gencode.v47/Log.out)
if [[ ! -z $refcheck ]]
then
  ln -s /mnt/analysis/workshop/CAT_Training/References/star.overlap100.gencode.v47 /mnt/analysis/cat_users/$USER/rnaseq_example/References/.
fi
```

---
## Alignments

1. Then run the star commands

```bash
mkdir tmp
STAR \
--runThreadN 4 \
    --genomeDir /mnt/analysis/cat_users/$USER/rnaseq_example/References/star.overlap100.gencode.v47 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix ./tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_ \
    --readFilesCommand zcat \
    --readFilesIn  ./01-HTS_Preproc/NEB_Mixed-10-ng-1_1M_S383_R1.fastq.gz  ./01-HTS_Preproc/NEB_Mixed-10-ng-1_1M_S383_R1.fastq.gz
```

In the command, we are telling star to count reads on a gene level ('--quantMode GeneCounts'), the prefix for all the output files will be ./tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_, the command to unzip the files (zcat), and finally, the input file pair.


###  Now let's take a look at an alignment in IGV.

1.  We first need to index the bam file, will use 'samtools' for this step, which is a program to manipulate SAM/BAM files. Take a look at the options for samtools and 'samtools index'.

```bash
samtools
samtools index
```

We need to index the BAM file:

```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example/tmp
samtools index NEB_Mixed-10-ng-1_1M_S383.htstream_Aligned.sortedByCoord.out.bam
```

**IF for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy.**

```bash
cp /mnt/analysis/workshop/CAT_Training/tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_Aligned.sortedByCoord.out.bam* /mnt/analysis/cat_users/$USER/rnaseq_example/tmp
```

2. Transfer ./tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_Aligned.sortedByCoord.out.bam and ./tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_Aligned.sortedByCoord.out.bam.bai (the index file) to your computer using scp,FileZilla or winSCP.

Mac/Linux users can use scp. In a new shell session on my laptop. **NOT logged into tadpole. Replace [your_username] with your username.**

**ON YOUR PERSONAL COMPUTER**

```bash
mkdir ~/rnaseq_workshop
cd ~/rnaseq_workshop
scp [your_username]@cat-bergamo1.ucsf.edu:/mnt/analysis/cat_users/[your_username]/rnaseq_example/tmp/NEB_Mixed-10-ng-1_1M_S383.htstream_Aligned.sortedByCoord.out.bam* .
```

Its ok if the mkdir command fails ("File exists") because we aleady created the directory earlier.


1. Now we are ready to use IGV.

Go to the [IGV page at the Broad Institute](http://software.broadinstitute.org/software/igv/).

<img src="alignment_figures/index_igv1.png" alt="index_igv1" width="80%" style="border:5px solid #ADD8E6;"/>

And then navigate to the download page, [IGV download](http://software.broadinstitute.org/software/igv/download)

<img src="alignment_figures/index_igv2.png" alt="index_igv2" width="80%" style="border:5px solid #ADD8E6;"/>

Here you can download IGV for your respective platform (Window, Mac OSX, Linux), but we are going to use the web application they supply, [IGV web app](https://igv.org/app).

<img src="alignment_figures/index_igv3.png" alt="index_igv3" width="80%" style="border:5px solid #ADD8E6;"/>

1. The first thing we want to do is load the Mouse genome. Click on "Genomes" in the menu and choose "Mouse (GRCm39/mm39)".

<img src="alignment_figures/index_igv4.png" alt="index_igv4" width="80%" style="border:5px solid #ADD8E6;"/>

1. Now let's load the alignment bam and index files. Click on "Tracks" and choose "Local File ...".

<img src="alignment_figures/index_igv5.png" alt="index_igv5" width="80%" style="border:5px solid #ADD8E6;"/>

Navigate to where you transferred the bam and index file and select them **both**.

<img src="alignment_figures/index_igv6.png" alt="index_igv6" width="80%" style="border:5px solid #ADD8E6;"/>

Now your alignment is loaded. Any loaded bam file aligned to a genome is called a "track".

<img src="alignment_figures/index_igv7.png" alt="index_igv7" width="80%" style="border:5px solid #ADD8E6;"/>

1. Lets take a look at the alignment associated with the gene __Fn1__, and if for some reason it doesn't find HBB (web IGV can be fickle) go to position __chr1:71,610,633-71,628,073__. If you don't see any reads, this likely means your in the wrong genome, double check that it says **mm39** in the top left.

<img src="alignment_figures/index_igv8.png" alt="index_igv8" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv9.png" alt="index_igv9" width="80%" style="border:5px solid #ADD8E6;"/>

You can zoom in by clicking on the plus sign (top right) or zoom out by clicking the negative sign (click it twice). You also may have to move around by clicking and dragging in the BAM track window.

You can also zoom in by clicking and dragging across the number line at the top. That section will highlight, and when you release the button, it will zoom into that section.

<img src="alignment_figures/index_igv10.png" alt="index_igv10" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv11.png" alt="index_igv11" width="80%" style="border:5px solid #ADD8E6;"/>

Reset the window by searching for Fn1 again, **and** then zoom in 2 steps.

<img src="alignment_figures/index_igv12.png" alt="index_igv12" width="80%" style="border:5px solid #ADD8E6;"/>

1. See that the reads should be aligning within the exons in the gene. This makes sense, since RNA-Seq reads are from exons. Play with the settings on the right hand side a bit and selecting reads.

<img src="alignment_figures/index_igv13.png" alt="index_igv13" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv14.png" alt="index_igv14" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv15.png" alt="index_igv15" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv16.png" alt="index_igv16" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="alignment_figures/index_igv17.png" alt="index_igv17" width="80%" style="border:5px solid #ADD8E6;"/>


---
## Running STAR on the experiment

- work through Running STAR on all samples
- work through QA/QC of the experiment
- complete the questions at the end

1. We can now run STAR across all samples on the real data using a script, [star.sh](../scripts/star.sh), that we should take a look at now.

```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example  # We'll run this from the main directory
wget https://raw.githubusercontent.com/ucsf-cat-bioinformatics/CAT_BIOINFO_TRAINING/master/software_scripts/scripts/star.sh
less star.sh
```

<pre class="prettyprint"><code class="language-py" style="background-color:333333">#!/bin/bash

## assumes star version 2.7.11b
## assumes STAR is available on the Path

start=`date +%s`
echo $HOSTNAME

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}

REF="./References/star.overlap100.gencode.v47"

for sample in `cat samples.txt`
do
  [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
  echo "SAMPLE: ${sample}"

  call="STAR
       --runThreadN 8 \
       --genomeDir $REF \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesCommand zcat \
       --readFilesIn 01-HTS_Preproc/${sample}_R1.fastq.gz 01-HTS_Preproc/${sample}_R2.fastq.gz \
       --quantMode GeneCounts \
       --outFileNamePrefix ${outpath}/${sample}/${sample}_ \
       > ${outpath}/${sample}/${sample}-STAR.stdout 2> ${outpath}/${sample}/${sample}-STAR.stderr"

  echo $call
  eval $call
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
</code></pre>


When you are done, type "q" to exit.

{:start="2"}
1. After looking at the script, lets run it.

```bash
bash star.sh  # moment of truth!
```

---
## Quality Assurance - Mapping statistics as QA/QC.

1. Once your jobs have finished successfully.

Use a script of ours, [star_stats.sh](../software_scripts/scripts/star_stats.sh) to collect the alignment stats. Don't worry about the script's contents at the moment; you'll use very similar commands to create a counts table in the next section. For now:

```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example # We'll run this from the main directory
wget https://raw.githubusercontent.com/ucsf-cat-bioinformatics/CAT_BIOINFO_TRAINING/master/software_scripts/scripts/star_stats.sh
bash star_stats.sh
```

<pre class="prettyprint"><code class="language-py" style="background-color:333333">#!/bin/bash

echo -en "sample_names" > names.txt
echo -en "total_in_feature\t" > totals.txt
cat 02-STAR_alignment/*/*ReadsPerGene.out.tab | head -4 | cut -f1 > stats.txt
cat samples.txt | while read sample; do
    echo ${sample}
    echo -en "\t${sample}" >> names.txt
    head -4 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | cut -f4 > temp1
    paste stats.txt temp1 > temp2
    mv temp2 stats.txt
    tail -n +5 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | cut -f4 | \
        perl -ne '$tot+=$_ }{ print "$tot\t"' >> totals.txt
done
echo -en "\n" >> names.txt
cat names.txt stats.txt totals.txt > temp1
mv temp1 summary_star_alignments.txt
rm stats.txt
rm names.txt
rm totals.txt
</code></pre>


Transfer `summary_star_alignments.txt` to your computer.


Or in case of emergency, download this copy: [summary_star_alignments.txt](../datasets/summary_star_alignments.txt)



**Questions:**
1. Look through the files in an output directory and check out what is present and discuss what each of them mean. (for example: `cd /mnt/analysis/cat_users/$USER/rnaseq_example/02-STAR_alignment/mouse_110_WT_C` )
2. Come up with a brief command you might use to check that all of the sample alignments using STAR have a reasonable output and/or did not produce any errors.
3. Open `summary_star_alignments.txt` in excel (or excel like application), and review. The table that this script creates ("summary_star_alignments.txt"). Are all samples behaving similarly? Discuss ...
4. If time, find some other regions/genes with high expression using IGV with your group. (Looking at genes the paper references is a great place to start)
