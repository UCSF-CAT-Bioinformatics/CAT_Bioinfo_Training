
<script>
function buildQuiz(myq, qc){
  // variable to store the HTML output
  const output = [];

  // for each question...
  myq.forEach(
    (currentQuestion, questionNumber) => {

      // variable to store the list of possible answers
      const answers = [];

      // and for each available answer...
      for(letter in currentQuestion.answers){

        // ...add an HTML radio button
        answers.push(
          `<label>
            <input type="radio" name="question${questionNumber}" value="${letter}">
            ${letter} :
            ${currentQuestion.answers[letter]}
          </label><br/>`
        );
      }

      // add this question and its answers to the output
      output.push(
        `<div class="question"> ${currentQuestion.question} </div>
        <div class="answers"> ${answers.join('')} </div><br/>`
      );
    }
  );

  // finally combine our output list into one string of HTML and put it on the page
  qc.innerHTML = output.join('');
}

function showResults(myq, qc, rc){

  // gather answer containers from our quiz
  const answerContainers = qc.querySelectorAll('.answers');

  // keep track of user's answers
  let numCorrect = 0;

  // for each question...
  myq.forEach( (currentQuestion, questionNumber) => {

    // find selected answer
    const answerContainer = answerContainers[questionNumber];
    const selector = `input[name=question${questionNumber}]:checked`;
    const userAnswer = (answerContainer.querySelector(selector) || {}).value;

    // if answer is correct
    if(userAnswer === currentQuestion.correctAnswer){
      // add to the number of correct answers
      numCorrect++;

      // color the answers green
      answerContainers[questionNumber].style.color = 'lightgreen';
    }
    // if answer is wrong or blank
    else{
      // color the answers red
      answerContainers[questionNumber].style.color = 'red';
    }
  });

  // show number of correct answers out of total
  rc.innerHTML = `${numCorrect} out of ${myq.length}`;
}
</script>

# Lesson Project Setup

# Project Setup

##  Creating a Project Directory

First, create a directory for you and the example project in the Lesson share directory:

```bash
cd
mkdir -p /mnt/analysis/cat_users/$USER/rnaseq_example
```

## Link raw fastq files

Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the sample fastq files that contains the raw data.


```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example
mkdir 00-RawData
cd 00-RawData/
ln -s /mnt/analysis/workshop/CAT_Training/htstream/smdata/* .
```


This directory now contains all of the sample fastq files for each sample.

## Create a Sample Sheet

Let's create a sample sheet for the project and store sample names in a file called samples.txt

```bash
ls *R1_001.fastq.gz | sed 's/\.R1\.fastq\.gz$//' > ../samples.txt
cat ../samples.txt
```

## Getting to know your data

Now, take a look at the raw data directory.

```bash
ls /mnt/analysis/cat_users/$USER/rnaseq_example/00-RawData
```


Lets get a better look at all the files.


```bash
ls -lah *
```

Pick a directory and go into it. View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):


```bash
zless mouse_110_WT_C.R1.fastq.gz
```


Make sure you can identify which lines correspond to a read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen.


Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:


```bash
zcat mouse_110_WT_C.R1.fastq.gz | wc -l
```


Divide this number by **4** and you have the number of reads in this file.


One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

```bash
zcat mouse_110_WT_C.R1.fastq.gz  | head -2 | tail -1
```


Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block.


Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):


```bash
echo -n [sequence] | wc -c
```


This will give you the length of the read. Also can do the bash one liner:


```bash
echo -n $(zcat mouse_110_WT_C.R1.fastq.gz  | head -2 | tail -1) | wc -c
```


See if you can figure out how this command works.

This will give you the read count without doing any division. See if you can figure out how this command works:


```bash
zcat mouse_110_WT_C.R1.fastq.gz | grep -c "^@A00461:28"
```

## Prepare our experiment folder for analysis

Now go back to your 'rnaseq_example' directory and create two directories called 'slurmout' and '01-HTS_Preproc':


```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example
mkdir References
mkdir 01-HTS_Preproc
```


### Learning Objectives

1. **Understand the Lesson Sample Dataset**
   - Familiarize yourself with the key findings and objectives of the study by Dorothée Selimoglu-Buet et al. on the miR-150/TET3 pathway and its role in monocyte differentiation.

2. **Perform Data Organization and Management**
   - Set up a well-structured project directory for RNA-seq data analysis, including the creation of directories for raw data, reference files, and preprocessing outputs.

   - Organize your experiment folder with appropriate subdirectories to streamline the analysis process, ensuring all steps are well-documented and reproducible.

   - Learn best practices for maintaining a clear and organized directory structure, enabling easy replication of the analysis process by others.

3. **Link and Organize Raw Data Files**
   - Create symbolic links to raw FASTQ files and generate a sample sheet to keep track of all sample names.

4. **Review and Preview Raw Data**
   - Explore the contents of FASTQ files, including identifying read structures, counting the number of reads, and determining read lengths.

5. **Evaluate and Troubleshoot Sequencing Data**
   - Assess the integrity and quality of sequencing data, and understand the sequencing run details, including the sequencer used, run number, and lane information.


## The Dataset

Dorothée Selimoglu-Buet, et al. ["A miR-150/TET3 pathway regulates the generation of mouse and human non-classical monocyte subset."](https://www.nature.com/articles/s41467-018-07801-x) Nature Communications volume 9, Article number: 5455 (2018)

The project on the European Nucleotide Archive (ENA), [PRJEB29201](https://www.ebi.ac.uk/ena/browser/view/PRJEB29201).

Relevant RNA sections of the paper
### Abstract
Non-classical monocyte subsets may derive from classical monocyte differentiation and the proportion of each subset is tightly controlled. Deregulation of this repartition is observed in diverse human diseases, including chronic myelomonocytic leukemia (CMML) in which non-classical monocyte numbers are significantly decreased relative to healthy controls. Here, we identify a down-regulation of hsa-miR-150 through methylation of a lineage-specific promoter in CMML monocytes. Mir150 knock-out mice demonstrate a cell-autonomous defect in non-classical monocytes. Our pulldown experiments point to Ten-Eleven-Translocation-3 (TET3) mRNA as a hsa-miR-150 target in classical human monocytes. We show that Tet3 knockout mice generate an increased number of non-classical monocytes. Our results identify the miR-150/TET3 axis as being involved in the generation of non-classical monocytes.

### Mouse models
Wild-type CD45.1 and CD45.2, C57BL/6 animals expressing the GFP under human ubiquitin C promoter38 and miR-150−/− mice were purchased from Jackson Laboratories (Charles River France, L’Arbresle, France). Mice harboring Tet3 allele with the coding sequences of exon 11 flanked by two loxP, a strategy similar that described with the Tet2 allele69, were generated by the Plateforme Recombinaison homolog (Institut Cochin, Paris, France) and were intercrossed with mice expressing tamoxifen-inducible Cre (Cre-ERT) transgene under control of the Scl/Tal1 promoter/enhancer. To delete Tet3-floxed alleles, tamoxifen was solubilized at 20 mg/ml in sunflower oil (Sigma-Aldrich, Saint-Quentin Fallavier, France) and 8 mg tamoxifen were administered to mice once per day for 2 days via oral gavage. For competitive and rescue experiments, recipient mice were housed in a barrier facility under pathogen-free conditions after transplantation. Cell transfer experiments were performed in 8- to 12-week-old female mice.

### RNA-sequencing
RNA integrity (RNA integrity score ≥7.0) was checked on the Agilent 2100 Bioanalyzer (Agilent) and quantity was determined using Qubit (Invitrogen). SureSelect Automated Strand Specific RNA Library Preparation Kit was used according to the manufacturer’s instructions with the Bravo Platform. Briefly, 100 ng of total RNA sample was used for poly-A mRNA selection using oligo(dT) beads and subjected to thermal mRNA fragmentation. The fragmented mRNA samples were subjected to cDNA synthesis and were further converted into double-stranded DNA using the reagents supplied in the kit, and the resulting double-stranded DNA was used for library preparation. The final libraries were sequenced on a NovaSeq 6000 for mice samples (Illumina) in paired-end 100 bp mode in order to reach at least 30 millions reads per sample at Gustave Roussy.

For mouse sample analysis, Fastq files quality have been analyzed with FastQC (v0.11.7) and aggregated with MultiQC (v1.5). The quantification was performed on Gencode mouse M18 (GRCm38p6) transcriptome and comprehensive gene annotation, with Salmon (v0.10.2). The index was build with the default k-mer length of 31, with the genecode flag on, the perfect hash option, and all 1569 sequence duplicates within the genome were kept. The quantification was done with the default expected maximization algorithm, verified through 100 boostrap rounds, sequence-specific bias correction, fragment GC-bias correction, and automatic library detection parameter. Clustered heatmaps were performs from normalized counts with pheatmap, an R package, using Pearson’s coefficient as distance metric for rows and column and the Ward.D2 method for the clustering. Volcano plot were built using ggplot package. The differential analysis was performed with Sleuth (v0.29.0), on data converted by wasabi (v0.2).

### Statistical analysis
Student's t test were performed using the Prism software.

### Results
RNA-sequencing of classical and nonclassical monocyte subsets of four healthy donors identified 2176 differentially expressed genes with an adjusted P value <0.01 (Supplementary Table 3 and Supplementary Figure 10C), including TET3 whose expression was down-regulated in nonclassical monocytes (−1.6-fold, P = 0.001) (Fig. 8d). Finally, an abnormal repartition of monocyte subsets, with an increase in Ly6Clow monocytes at the expanse of Ly6Chigh cells, was detected in the blood of mice carrying inactivated Tet3 alleles compared to wild-type littermates (Fig. 8e, f and Supplementary Figure 10D–F), without any change in the mean fluorescence intensity of CD115 and CX3CR1 at the surface of Tet3−/− mouse monocytes (Supplementary Figure 10G). Hence, these experiments argued for TET3 as a target of miR-150 whose down-regulation is required for the differentiation of classical into nonclassical monocytes, in mice and in humans.

# Project Setup

Let's set up a project directory for the week, and talk a bit about project philosophy..

##  Creating a Project Directory

First, create a directory for you and the example project in the Lesson share directory:

```bash
cd
mkdir -p /mnt/analysis/cat_users/$USER/rnaseq_example
```

## Link raw fastq files

Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the sample fastq files that contains the raw data.


```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example
mkdir 00-RawData
cd 00-RawData/
ln -s ln -s /mnt/analysis/cat_users/original_dataset/00-RawData/* .
```


This directory now contains all of the sample fastq files for each sample.

## Create a Sample Sheet

Let's create a sample sheet for the project and store sample names in a file called samples.txt

```bash
ls *.R1.fastq.gz | sed 's/\.R1\.fastq\.gz$//' > ../samples.txt
cat ../samples.txt
```

## Getting to know your data

Now, take a look at the raw data directory.

```bash
ls /mnt/analysis/cat_users/$USER/rnaseq_example/00-RawData
```


Lets get a better look at all the files.


```bash
ls -lah *
```

Pick a directory and go into it. View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):


```bash
zless mouse_110_WT_C.R1.fastq.gz
```


Make sure you can identify which lines correspond to a read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen.


Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:


```bash
zcat mouse_110_WT_C.R1.fastq.gz | wc -l
```


Divide this number by **4** and you have the number of reads in this file.


One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

```bash
zcat mouse_110_WT_C.R1.fastq.gz  | head -2 | tail -1
```


Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block.


Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):


```bash
echo -n [sequence] | wc -c
```


This will give you the length of the read. Also can do the bash one liner:


```bash
echo -n $(zcat mouse_110_WT_C.R1.fastq.gz  | head -2 | tail -1) | wc -c
```


See if you can figure out how this command works.

This will give you the read count without doing any division. See if you can figure out how this command works:


```bash
zcat mouse_110_WT_C.R1.fastq.gz | grep -c "^@A00461:28"
```

## Prepare our experiment folder for analysis

Now go back to your 'rnaseq_example' directory and create two directories called 'slurmout' and '01-HTS_Preproc':


```bash
cd /mnt/analysis/cat_users/$USER/rnaseq_example
mkdir References
mkdir 01-HTS_Preproc
```


We'll put reference sequence, genome, etc. in the References directory. The results of all our slurm script will output .out and .err files into the slurmout folder. The results of our preprocessing steps will be put into the 01-HTS_Preproc directory. The next step after that will go into a "02-..." directory, etc. You can collect scripts that perform each step, and notes and metadata relevant for each step, in the directory for that step. This way anyone looking to replicate your analysis has limited places to search for the commands you used. In addition, you may want to change the permissions on your original 00-RawData directory to "read only", so that you can never accidentally corrupt (or delete) your raw data. We won't worry about this here, because we've linked in sample folders.

Your directory should then look like the below:


<div class="output">
$ ls
00-RawData  01-HTS_Preproc  References  samples.txt
</div>


### Questions you should now be able to answer.

<div id="quiz_setup1" class="quiz"></div>
<button id="submit_setup1">Submit Quiz</button>
<div id="results_setup1" class="output"></div>
<script>
quizContainer1 = document.getElementById('quiz_setup1');
resultsContainer1 = document.getElementById('results_setup1');
submitButton1 = document.getElementById('submit_setup1');

myQuestions1 = [
  {
    question: "How many samples are in this experiment?",
    answers: {
      a: "1",
      b: "6",
      c: "22",
      d: "44"
    },
    correctAnswer: "c"
  },
  {
    question: "How many reads are in the sample you checked?",
    answers: {
      a: "25000000",
      b: "4000000",
      c: "1256832",
      d: "1"
    },
    correctAnswer: "b"
  },
  {
    question: "How many basepairs is R1, how many is R2?",
    answers: {
      a: "100 and 100",
      b: "100 and 50",
      c: "101 and 101",
      d: "50"
    },
    correctAnswer: "c"
  },
  {
    question: "What is the name of the sequencer this dataset was run on?",
    answers: {
      a: "A00461",
      b: "UCSFCAT",
      c: "@A00461:28",
      d: "HF2WYDMXX"
    },
    correctAnswer: "a"
  },
  {
    question: "Which run number is this for that sequencer?",
    answers: {
      a: "HF2WYDMXX",
      b: "1",
      c: "1101",
      d: "28"
    },
    correctAnswer: "d"
  },
  {
    question: "What lane was this ran on?",
    answers: {
      a: "HF2WYDMXX",
      b: "1",
      c: "1101",
      d: "28"
    },
    correctAnswer: "b"
  },
  {
    question: "Randomly check a few samples, were they all run on the same sequencer, run, and lane?, Can you check them all quckly",
    answers: {
      a: "TRUE",
      b: "FALSE",
    },
    correctAnswer: "a"
  }
];

buildQuiz(myQuestions1, quizContainer1);
submitButton1.addEventListener('click', function() {showResults(myQuestions1, quizContainer1, resultsContainer1);});
</script>
