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
    
### Starting up DASHit

SO code seems to be pretty old and I can't get it to install on the server, probably could get it to work on my laptop, but thats not the point.

Next the code does have a docker container that works, 

* Docker is a platform for creating, deploying, and managing lightweight, isolated environments called containers, which bundle applications with their dependencies to ensure consistent performance across different systems. It allows developers to build once and run anywhere, improving portability and scalability.

However, docker requires you to run it under admin (sudo)

* Singularity, is a container platform designed for high-performance computing (HPC) and scientific applications, enabling users to run applications in isolated environments without needing root privileges. It’s optimized for compatibility with shared computing resources, making it ideal for securely running complex applications and reproducible research across different systems.

* The main difference between Docker and Singularity is their focus on use cases and security: Docker is designed primarily for software development, microservices, and scalable web applications, often requiring root privileges to manage containers. Singularity, on the other hand, is tailored for high-performance computing (HPC) environments, allowing users to run containers securely without root access, making it ideal for shared computing resources in academic and research settings.

I've created a singularity docker sandbox that we can run on ther server.

Lets create our environment
```bash
cd /mnt/analysis/cat_users/$USER
mkdir -p dashit
cd dashit
```

lets copy over some example data

```bash
cp -r /mnt/analysis/software/dashit/misc/tests/data .
```

NOW lets launch DASHit

```bash
singularity exec --bind /mnt/analysis/software/dashit/misc/tests/data:/data /mnt/analysis/software/dashit/dashit_sandbox bash
```

This should produce 
```
Singularity>
```

Now cd into data
```bash
cd /data
```

Run crispr_sites -r to find candidate Cas9-gRNAs in your input.

cat test_reads.fasta | crispr_sites -r > input_sites_to_reads.txt

cat test_expected_reads_dashed.fasta | crispr_sites > ontarget_sites.txt

cat test_expected_reads_undashed.fasta | crispr_sites > offtarget_sites.txt

#doesn't work, offtarget app fails
#ontarget and offtarget filtering require port 8080 to be available on your computer.
dashit_filter --gc_freq_min 5 --gc_freq_max 15 --ontarget ontarget.txt --offtarget offtarget.txt input_sites_to_reads.txt > input_sites_to_reads_filtered.txt

dashit_filter --gc_freq_min 5 --gc_freq_max 15 input_sites_to_reads.txt > input_sites_to_reads_filtered.txt

#Find 25 guides that hit the largest number of reads
optimize_guides input_sites_to_reads_filtered.txt 25 1 > guides.csv

The number of times to cover each read option, here set to 1, is how many guides need to hit a read in input.fasta before that read is considered covered. In principle you could use this for additional redundancy, but in practice we’ve never needed anything other than 1

score_guides guides.csv test_reads.fasta
