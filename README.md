# TECtool (Terminal Exon Characterization tool)

## Features

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

## INSTALLATION

**TECtool uses genome sequence, annotation and RNA-seq data. Therefore, ~10 GB of disk space are needed for installation and testing.**

TECtool as of version 0.2 is written in Python 3. Current instructions are written for Python 3 (>=3.4) ONLY. Installation with Python 2 will not work. The recommended way to install TECtool is via the conda package manager, because it can install non Python dependencies (for example bedtools).

If you do not want to use conda to install TECtool, other options are described below. 

### Installation of TECtool using conda

#### Step 1: Download miniconda 3 installation file (if not already installed)

You can do this with one of the following options:

* filling in the **URL** for the appropriate file in a browser window and saving the file

for Linux:
```
https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
for Mac OSX:
```
https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
    
* using **wget**:

for Linux:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
for Mac OSX:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
   
* using **curl**:

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

for Linux:
```
bash Miniconda3-latest-Linux-x86_64.sh
```
for Mac OSX:
```
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Step 3: Create a new conda environment

Create a new conda environment
```
conda create --name TECtool --channel bioconda --channel conda-forge --channel fgypas tectool
```

Activate the virtual environment
```
source activate TECtool
```

Check the installation and options for the tool
```
tectool --help
```

### Step 3 (alternative 1): Create a new conda environment and install the dependencies manually

Create a new conda environment that only contains python 3
```
conda create --name TECtool --channel bioconda python=3.6.2
```

Activate the virtual environment
```
source activate TECtool
```

Install dependencies one by one
```
conda install --channel bioconda htseq==0.9.1
conda install --channel bioconda bedtools==2.26.0
conda install --channel bioconda pybedtools==0.7.10
conda install --channel conda-forge bzip2
conda install --channel bioconda pyfasta==0.5.2
conda install --channel coda-forge scikit-learn==0.19.0
conda install --channel fgypas tectool
```

Check the installation and options for the tool
```
tectool --help
```

### Step 3 (alternative 2): Install tectool in an existing environment or globally

Install tectool in an existing conda environment or globally
```
conda install --channel bioconda --channel conda-forge --channel fgypas tectool
```

Check the installation and options for the tool
```
tectool --help
```

### Installation of TECtool without conda

For users that do not want to use conda, but have a python installation (3.4 and above).
**Note:** Normally venv (https://docs.python.org/3/library/venv.html) should be pre-installed in Python 3.
Some distributions separate it to a different package, such as python3-venv on Ubuntu/Debian, so users need to additionally install it before running the following command.


Create a virtual environment with venv named envTECtool.
```
python -m venv envTECtool
```

Activate the virtual environment
```
source envTECtool/bin/activate
```

In the activated virtual environment you might need to upgrade pip
```
pip install --upgrade pip
```

Clone the TECtool repository
```
git clone https://git.scicore.unibas.ch/zavolan_public/TECtool.git
```

Now you should see a direcroty for the virtual environment (envTECtool) and one directory for the TECtool package (TECtool)

Enter the cloned directory
```
cd TECtool
```

**Important Note:** The requirements that will be installed include only Python modules.
Users should additionally install **bedtools version 2.26** in their system.
For installation instruction please see here: http://bedtools.readthedocs.io/en/latest/content/quick-start.html#install-bedtools.
TECtool is not checking if the correct version of bedtools is installed and this might lead to run-time errors.

Install dependencies with
```
pip install -r requirements.txt
```

You can check the version of bedtools by typing
```
bedtools --version
```

Install TECtool
```
python setup.py install
```

Check the installation and options for the tool
```
tectool --help
```

**Note:** TECtool will be available only when the virtual environment is active.

## TECtool options

The following options are available and should be set by the user:

* *--annotation* FILE: Annotation file gtf format (tested with ENSEMBL gtf v87) [REQUIRED].

* *--polyasites* FILE: Bed file (bed6) that contains polya sites [REQUIRED].

* *--bam* FILE: The BAM file that should be analysed. [REQUIRED] Note that the BAM file should be sorted by coordinates. An index file should be also present in the same directory [REQUIRED].

* *--sequencing_direction* SEQUENCING_DIRECTION: If the data is generated by a stranded protocol, reads should be provided in the forward orientation and set the option to 'forward'. If the reads were generated by an unstranded protocol select 'unstranded' [REQUIRED].

* *--genome* FILE: Genome sequence in fasta format. The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 denoted by "1", then the sequence in the fasta file should have a description line of the form ">1". No white spaces or trailing text should be included. [REQUIRED].

* *--minimum_spliced_reads_for_cryptic_exon_start_site* Minimum number of spliced reads required to consider the start site of a novel exon. [default=5]

* *--min_region_overlap* MIN_REGION_OVERLAP: min_region_overlap (default=10) Minimum number of bases that a read should share with an exon to be considered towards feature calculations for the exon (It will be suppressed in future versions)

* *--max_splice_fuzziness* MAX_SPLICE_FUZZINESS: Maximum splice fuzziness (default=0) Maximum number of bases by which the 5' boundary of a candidate novel exon is allowed to vary, considering all reads that splice into that candidate exon (It will be suppressed in future versions)

* *--drop_manually_selected_features*: Select features for the exon classification from scratch, using a greedy approach (default=False). Default is to include a minimum set of features (['ReadsOUTvsIN_all', 'entropy_efficiency']), which we found most informative. 

* *--drop_intronic_polya_sites_of_overlapping_genes*: Ignore intronic polya sites that occur in regions of overlap between genes (default=False).

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
```
wget http://tectool.unibas.ch/data/test_data.tar.gz
```

Uncompress the files
```
tar xzvf test_data.tar.gz
```

Enter the directory
```
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
