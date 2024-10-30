# DASHit

## References 

- Gu, W. et al. [Depletion of Abundant Sequences by Hybridization (DASH): using Cas9 to remove unwanted high-abundance species in sequencing libraries and molecular counting applications.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0904-5) Genome Biology 17, 41 (2016).
- Dynerman, D and Lyden, A. et al. DASHit. Preprint soon! ?????
- DASHit [Github page(https://github.com/czbiohub-sf/dashit/tree/master)].
- Lyden, A. et al. [DASH Protocol V.4](https://www.protocols.io/view/dash-protocol-yxmvm7z99v3p/v4).

## Learning Objectives

*Run the DashIt pipeline to produce a set of 20-mer Guide RNAs*

### The Basic Proceedure

*Prepare the Reads*
1. From a datasets fastq file(s) (or bam file of mapped reads).
2. Extract a set of unwanted reads, so over-abundant reads/regions.
3. Clean the reads, remove things like primers, adapters, etc.
4. Convert the fastq file to a fasta file.
5. Reduce the data?
*Produce candidate guides*
6. Run DASHit crispr_sites to produce an initial candidate set of 20-mer guides
*Filter candidates on/off target* (Optional)
7. Create a ontarget and offtarget set of 'guides'
8. Use dashit_filter to remove all guides RNAs present in offtarget, and those not present in ontarget_sites.txt
*Optimize Guides*
9. Use optimize_guides to return a predefined number of reads meeting a certain requirement.
10. Score final guides against the input reads.
    
    

