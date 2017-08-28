# TECtool (Terminal exon characterization Tool)


## Features

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

## Input and output files

Input files
* A file containing all chromosomes in fasta format
* A file with the corresponding annotation in GTF format
* A file with genome coordinates of 3’ end processing sites in BED format
* A file containing spliced alignments of mRNA-seq reads to the corresponding genome (in BAM format, sorted by coordinates and indexed)

Output files
* An augmented annotation file (in GTF format)
* A file containing the novel terminal exons

## Installation instructions


### Step 1: Download Miniconda3


On Linux:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

On MacOS X:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Step 2: Create a new environment with all required packages/software

```bash
conda create -n TECtool python=2.7 -c bioconda --file requirements.txt
```


### Step 3: Activate the virtual environment

```bash
source activate TECtool
```

To exit the environment (after finishing the usage of the pipeline), just execute
```bash
source deactivate
```

### Step 4: Clone the git repository

```bash
git clone https://git@git.scicore.unibas.ch:2222/zavolan_public/TECtool.git
```

or

```bash
git clone ssh://git@git.scicore.unibas.ch:2222/zavolan_public/TECtool.git
```

## Run TECtool

In order see all the option that are provided to TECtool then just go to the 
TECtool directory

```bash
cd TECtool/tectool
```

and run (after you activate the virtual environment that you previously generated)

```bash
python tectool.py --help
```

### Options

The following options are available and should be set by the user:
* *--annotation* FILE: The annotation file GTF format (as provided in EMSEMBL)
* *--polyasites* FILE: Bed (bed6) file that contains the polya sites
* *--bam* FILE: The BAM file that should be analysed. Note that the BAM file should be shorted by coordinates. An index file should be also present in the same directory.
* *--sequencing_direction* SEQUENCING_DIRECTION: Are the reads annotated with 'First' in the BAM file mapped forward (on the right strand) or reverse (on the opposite strand)? Please choose 'forward' or 'reverse', respectively. If the data are unstranded please select 'unstranded'. Currently, only 'forward' or 'unstranded' options are recommended.
* *--genome* FILE: Genome in fasta format
* *--minimum_spliced_reads_for_cryptic_exon_start_site* MINIMUM_SPLICED_READS_FOR_CRYPTIC_EXON_START_SITE: Minimum number of spliced reads required to characterize the start site of a cryptic exon. (default=5)
* *--min_region_overlap* MIN_REGION_OVERLAP: min_region_overlap (default=1) (It will be suppressed in future versions)
* *--max_splice_fuzziness* MAX_SPLICE_FUZZINESS: Maximum splice fuzziness (default=3) (It will be suppressed in future versions)
* *--feature_length_threshold* FEATURE_LENGTH_THRESHOLD: Some poly(A) sites are located a few bases away from annotated exons. It is difficult for TECtool to predict if the reads correspond to splicing of an annotated exon or if there is actual terminal exons. For this reason TECtool reports such cases as computationally predicted. The flag feature_length_threshold specifies the distance of the poly(A) site to the annotated exons to allow such predictions. (default=10) (It will be suppressed in future versions)
* *--drop_manually_selected_features*: Flag to not use the manually selected features (['ReadsOUTvsIN_all', 'entropy_efficiency']). When this flag is used all the features will be selected by greedy.
* *--output_dir* OUTPUT_DIR: The path to the output directory.

### Output files

The output of TECtool:
* An augmented annotation file in gtf format named enriched_annotation.gtf. The gtf file contains genes, transcripts, exons, CDS, START and STOP sites.
* A file containing the novel terminal exons named classified_as_terminal_with_probabilities.tsv: The table contains the terminal exon region, the gene id, the features that were used, the probability that this region is terminal (terminal_probability), the probability that this region is intermediate (intermediate_probability), the probability that the region is background (background_probability), the type that was selected (terminal/intermediate/background) and the genomic coordinates of the region (chromosome, start, end, strand).


## Licence

* Free software: MIT license