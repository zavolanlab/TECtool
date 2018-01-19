# TECtool (Terminal Exon Characterization tool)

## Features

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

## INSTALLATION

** TECtool uses genome sequence, annotation and RNA-seq data. Therefore, ~10 GB of disk space are needed for installation and testing. **

TECtool as of version 0.2 is written in Python 3. The recommended way to install TECtool is via the conda package manager, because it can install non Python dependencies (for example bedtools). 

If you do not want to use conda to install TECtool, other options are described below. 

### Conda installation

#### Step 1: Download miniconda 3 installation file (if not already installed)

You can do this by:
    1. filling in the URL for the appropriate file in a browser window and saving the file
    for Linux:
    ```
    https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    ```   
    for Mac OSX:
    ```
    https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    ```
    
    2. using wget:
    for Linux:
    ```
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    ```
     for Mac OSX:
    ```
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    ```
   
    3. using curl:
    for Linux:
    ```
    curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3-latest-Linux-x86_64.sh
    ```
     for Mac OSX:
    ```
    curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    ```

#### Step 2: Install miniconda 3

Make sure that you run the 'bash' shell and execute:
```
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Step 2: Create a new conda environment with conda create

Create a new conda environment
```bash
conda create --name TECtool --channel bioconda --channel conda-forge --channel fgypas tectool
```

Activate the virtual environment
```bash
source activate TECtool
```

Check if tectool returns the necessary argparse options by typing
```bash
tectool --help
```

### Step 2 (alternative): Create a new conda environment and install the dependencies manually

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
conda install --channel bioconda bedtools==2.26.0
conda install --channel bioconda pybedtools==0.7.10
conda install --channel conda-forge bzip2
conda install --channel bioconda pyfasta==0.5.2
conda install --channel coda-forge scikit-learn==0.19.0
conda install --channel fgypas tectool
```

Check if tectool returns the necessary argparse options by typing
```bash
tectool --help
```

### Step 2 (alternative): Install tectool in an existing environment or globally

Install tectool in an existing conda environment or globally
```bash
conda install --channel conda-forge --channel bioconda --channel fgypas tectool
```

Check if tectool returns the necessary argparse options by typing
```bash
tectool --help
```

### Non conda installation

For users that do not want to use conda

Create a virtual environment with virtualenv
```bash
virtualenv TECtool
```

Clone the TECtool repository
```bash
git clone https://git.scicore.unibas.ch/zavolan_public/TECtool.git
```

Enter the directory
```bash
cd TECtool
```

Install dependencies with
```bash
pip install .
```

or

```bash
pip install -r requirements.txt
```

**Important Note:** The requirements that will be installed include only Python modules. Users should additionally install bedtools=2.26 in their system. TECtool is not checking if the correct version of bedtools is installed and this might lead to run errors.

Install TECtool
```bash
python setup.py install
```

Check if tectool returns the necessary argparse options by typing
```bash
tectool --help
```

## Available options

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

## Input files

Input files

* A file containing all chromosomes in fasta format. **Important note:** The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 annotated as "1", then the fasta should have a header called ">1". No white spaces or trailing text should be included.
* A file with the corresponding annotation in GTF format. **Important note:** Currently only gtf files in ENSEMBL (tested with ENSEMBL v87).
* A file with genome coordinates of 3’ end processing sites in BED format.
* A file containing spliced alignments of mRNA-seq reads to the corresponding genome (in BAM format, sorted by coordinates and indexed) (tested with STAR aligner).

## Output files

The output of TECtool:
* An augmented annotation file in gtf format named enriched_annotation.gtf. The gtf file contains genes, transcripts, exons, CDS, START and STOP lines.
* A file containing the novel terminal exons named classified_as_terminal_with_probabilities.tsv: The table contains the terminal exon region, the gene id, the features that were used, the probability that this region is terminal (terminal_probability), the probability that this region is intermediate (intermediate_probability), the probability that the region is background (background_probability), the type that was selected (terminal/intermediate/background) and the genomic coordinates of the region (chromosome, start, end, strand).

## Recommended files for testing

Download the test data
```bash
wget http://tectool.unibas.ch/data/test_data.tar.gz
```

Uncompress the files
```bash
tar xzvf test_data.tar.gz
```

Enter the directory
```bash
cd test_data
```

Run TECtool with the following options
```
tectool \
--annotation Homo_sapiens.GRCh38.87.chr.support_level_5.correct_gene_coordinates.chr1.14.22.X.16.gtf \
--polyasites polya_sites.merged.anno.hg38.ENSEMBL.chr1.14.22.X.16.bed \
--bam GSM1502499_RNA_seq_control_rep2.chr22.bam \
--genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.fixed.fa \
--sequencing_direction forward \
--minimum_spliced_reads_for_cryptic_exon_start_site 5 \
--min_region_overlap 10 \
--max_splice_fuzziness 0 \
--output_dir results \
--verbose
```

## Licence and documentation

* Free software: MIT license
* Official website: http://tectool.unibas.ch/
