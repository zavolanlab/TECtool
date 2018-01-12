# TECtool (Terminal exon characterization Tool)

## Features

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

## Input and output files

Input files

* A file containing all chromosomes in fasta format. **Important note:** The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 annotated as "1", then the fasta should have a header called ">1". No white spaces or trailing text should be included.
* A file with the corresponding annotation in GTF format. **Note:** Currently only gtf files in ENSEMBL (tested with ENSEMBL v87), or GENCODE are supported. 
* A file with genome coordinates of 3’ end processing sites in BED format.
* A file containing spliced alignments of mRNA-seq reads to the corresponding genome (in BAM format, sorted by coordinates and indexed) (tested with STAR aligner).

Output files
* An augmented annotation file (in GTF format)
* A file containing the novel terminal exons

## INSTALLATION

TECtool as of version 0.2 is written in Python 3. You can use conda package manager to install it or you can install it via pip.

### Conda installation

#### Step 1: Download Miniconda3 (if not already installed)

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
### Step 2

Create a new conda environment

```bash
conda create --name TECtool --channel bioconda python=3.6.2
```

Activate the virtual environment

```bash
source activate TECtool
```

Install dependencies

```bash
conda install --channel bioconda htseq==0.9.1
conda install --channel bioconda pybedtools==0.7.10
conda config --add channels conda-forge
conda install --channel conda-forge bzip2
conda install --channel bioconda pyfasta==0.5.2
conda install --channel coda-forge scikit-learn==0.19.0
```

### Pip intallation

#### Create a virtual env

```bash
python3 -m virtualenv venvTECtool
```

#### Install requirements

```bash
pip install -r requirements.txt
```

## Run TECtool

Clone the repository

```bash
git clone https://git.scicore.unibas.ch/zavolan_public/TECtool.git
```
In order see all the option that are provided to TECtool then just go to the TECtool directory

```bash
cd scripts
```

and run (after you activate the virtual environment that you previously generated)

```bash
tectool --help
```

### Options

The following options are available and should be set by the user:

* *--annotation* FILE: Annotation file gtf format (tested with ENSEMBL gtf v87) [REQUIRED].

* *--polyasites* FILE: Bed file (bed6) that contains polya sites [REQUIRED].

* *--bam* FILE: The BAM file that should be analysed. [REQUIRED] Note that the BAM file should be shorted by coordinates. An index file should be also present in the same directory [REQUIRED].

* *--sequencing_direction* SEQUENCING_DIRECTION: Are the reads annotated with 'First' in the BAM file mapped forward (on the right strand) Please choose 'forward'. If the data are unstranded please select 'unstranded' [REQUIRED].

* *--genome* FILE: Genome in fasta format. The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 annotated as "1", then the fasta should have a header called ">1". No white spaces or trailing text should be included. [REQUIRED].

* *--minimum_spliced_reads_for_cryptic_exon_start_site* Minimum number of spliced reads required to characterize the start site of a cryptic exon. [default=5]

* *--min_region_overlap* MIN_REGION_OVERLAP: min_region_overlap (default=10) (It will be suppressed in future versions)

* *--max_splice_fuzziness* MAX_SPLICE_FUZZINESS: Maximum splice fuzziness (default=0) (It will be suppressed in future versions)

* *--drop_manually_selected_features*: Flag to not use the manually selected features (['ReadsOUTvsIN_all', 'entropy_efficiency']). When this flag is used all the features will be selected by greedy (default=False).

* *--drop_intronic_polya_sites_of_overlapping_genes*: Flag to ignore intronic polya sites shared among overlapping genes (default=False).

* *--output_dir* OUTPUT_DIR: The path to the output directory.

### Output files

The output of TECtool:
* An augmented annotation file in gtf format named enriched_annotation.gtf. The gtf file contains genes, transcripts, exons, CDS, START and STOP lines.
* A file containing the novel terminal exons named classified_as_terminal_with_probabilities.tsv: The table contains the terminal exon region, the gene id, the features that were used, the probability that this region is terminal (terminal_probability), the probability that this region is intermediate (intermediate_probability), the probability that the region is background (background_probability), the type that was selected (terminal/intermediate/background) and the genomic coordinates of the region (chromosome, start, end, strand).

## Recommended files for testing

Coming soon

## Licence and documentation

* Free software: MIT license
* Official website: http://tectool.unibas.ch/




