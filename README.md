# TECtool (Terminal Exon Characterization tool)

## Features

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

## INSTALLATION

**TECtool uses genome sequence, annotation and RNA-seq data. Therefore, ~10 GB of disk space are needed for installation and testing. The installation should take about 5 to 10 minutes.**

TECtool as of version 0.2 is written in Python 3. Current instructions are written for Python 3 ONLY (>=3.4 should work, but extensively tested with Python 3.6). Installation with Python 2 will not work. The recommended way to install TECtool is via the conda package manager, because it can install non Python dependencies (for example bedtools).

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

If installation was sucessfull, you can start testing (see section [Recommended files for testing](README.md#recommended-files-for-testing))

### Step 3 (alternative 1): Install tectool in an existing environment or globally

Install tectool in an existing conda environment or globally
```
conda install --channel bioconda --channel conda-forge --channel fgypas tectool
```

Check the installation and options for the tool
```
tectool --help
```

If installation was sucessfull, you can start testing (see section [Recommended files for testing](README.md#recommended-files-for-testing))

### Step 3 (alternative 2): Create a new conda environment and install the dependencies manually

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
conda install --channel conda-forge progress==1.3
conda install --channel fgypas tectool
```

Check the installation and options for the tool
```
tectool --help
```

If installation was sucessfull, you can start testing (see section [Recommended files for testing](README.md#recommended-files-for-testing))

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
git clone https://github.com/zavolanlab/TECtool.git
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

If installation was sucessfull, you can start testing (see section [Recommended files for testing](README.md#recommended-files-for-testing))

### Docker

A docker image is also available with all dependencies installed. Users should have docker installed in their system with root access.

Pull the docker image
```
sudo docker pull fgypas/tectool:0.4
```

Start and enter the container
```
sudo docker run -it fgypas/tectool:0.4 bash
```

If installation was sucessfull, you can start testing (see section [Recommended files for testing](README.md#recommended-files-for-testing))

## TECtool options

The following options are available and should be set by the user:

* *--annotation* FILE: Annotation file in ENSEMBL GTF format  (tested with ENSEMBL v87) [REQUIRED]. Other GTF formats are currently not supported.

* *--polyasites* FILE: Bed file that contains polya sites [REQUIRED].

* *--genome* FILE: Genome in fasta format [REQUIRED]. Note for versions <0.3: The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 denoted by "1", then the sequence in the fasta file should have a description line of the form ">1". No white spaces or trailing text should be included.

* *--bam* FILE: The BAM file that should be analysed. [REQUIRED] Note that the BAM file should be sorted by coordinates. An index file should be also present in the same directory.

* *--sequencing_direction* SEQUENCING_DIRECTION: If the data is generated by a stranded protocol, reads should be provided in the forward orientation and set the option to 'forward'. If the reads were generated by an unstranded protocol select 'unstranded'. [default=unstranded]

* *--minimum_spliced_reads_for_cryptic_exon_start_site* MINIMUM_SPLICED_READS_FOR_CRYPTIC_EXON_START_SITE: Minimum number of spliced reads required to consider the start site of a novel exon. [default=5]

* *--min_region_overlap* MIN_REGION_OVERLAP: Minimum number of bases that a read should share with an exon to be considered towards feature calculations for the exon (It will be suppressed in future versions). [default=10]

* *--max_splice_fuzziness* MAX_SPLICE_FUZZINESS: Maximum number of bases by which the 5' boundary of a candidate novel exon is allowed to vary, considering all reads that splice into that candidate exon (It will be suppressed in future versions). [default=0]

* *--drop_manually_selected_features*: Select features for the exon classification from scratch, using a greedy approach. Default is to include a minimum set of features (['ReadsOUTvsIN_all', 'entropy_efficiency']), which we found most informative. [default=False]

* *--drop_intronic_polya_sites_of_overlapping_genes*: Ignore intronic polya sites that occur in regions of overlap between genes. [default=False]

* *--use_precalculated_training_set*: Use precalculated training set (skips the generation of training and validation sets). This option should be provided in combination with the option: --training_set_directory. [default=False]

* *--training_set_directory* DIRECTORY: Training set directory created by another TECtool run using the same annotation files. This option should be provided in combination with the option: --use_precalculated_training_set, otherwise it is ignored. [No default option provided]

* *--output_dir* DIRECTORY: The path to the output directory. [default="." (current working directory)]

* *--verbose* Be Verbose

* *--version* show program's version number and exit

* *--help* show help message and exit

## Input files

Input files

* A file containing all chromosomes in fasta format. **Important note for versions <0.3::** The file should have the same chromosome names (header lines) as the ones specified in the gtf file. For example if the gtf file has chromosome 1 annotated as "1", then the fasta should have a header called ">1". No white spaces or trailing text should be included.
* A file with the corresponding annotation in GTF format. **Important note:** Currently only gtf files in ENSEMBL (tested with ENSEMBL v87).
* A file with genome coordinates of 3’ end processing sites in BED format.
* A file containing spliced alignments of mRNA-seq reads to the corresponding genome (in BAM format, sorted by coordinates and indexed) (tested with STAR aligner).

## Output files

The output of TECtool:
* An augmented annotation file in gtf format named enriched_annotation.gtf. The gtf file contains genes, transcripts, exons, CDS, START and STOP lines.
* A file containing the novel terminal exons named classified_as_terminal_with_probabilities.tsv: The table contains the terminal exon region, the gene id, the features that were used, the probability that this region is terminal (terminal_probability), the probability that this region is intermediate (intermediate_probability), the probability that the region is background (background_probability), the type that was selected (terminal/intermediate/background) and the genomic coordinates of the region (chromosome, start, end, strand).
* A file containing the root and the novel transcript ids
* When TECtool is run without the options --use_precalculated_training_set --training_set_directory a directory called training_data is generated. This can be used as input when the options the options --use_precalculated_training_set --training_set_directory are provided.

## Plot novel exons

A supplementary script (written in R) is also provided that uses one of the outputs of TECtool and visualizes the novel terminal exons. The script is called plot_novel_exons.R and is available in the scripts directory of TECtool.
In order to run it users should have R installed (>=3.4) (tested with R 3.4.1 on CentOS 7.3) with the following packages: optparse, rtracklayer, Gviz, biomaRt and GenomicFeatures.

**Note for users that installed tectool via conda**: The default environment for running tectool does not contain any R installation. In order to run the plotting script please create a new conda environment that contains both TECtool and the R dependencies. You can do this as following:

```
conda create --name TECtool_plot_novel_exons --channel bioconda --channel conda-forge --channel r --channel fgypas r-base=3.4.1 bioconductor-gviz r-optparse openblas tectool
```

**Note:** that this is a time consuming step and many packages are installed.

Acivate the virtual environment
```
source activate TECtool_plot_novel_exons
```

The following options are available in the script:
* *--gtf*: GTF file with annotated transcripts
* *--polyasites*: BED file with polya sites
* *--bam*: Alignment file in BAM format
* *--tectool_exons*: TECtool exons file in tsv format. Output of TECtool (classified_as_terminal_with_probabilities.tsv).
* *--output_dir*: Output directory
* *--help*: Show help message
* *--verbose*: Be Verbose

## Recommended files for testing

*Using the recommended files testing TECtool takes around 10 minutes run time on a desktop computer (using one core).*

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
--output_dir results
```

**Note:** Some warnings will appear because the dataset that we use is small.

In order to test the vizualization script please first check the section [Plot novel exons](README.md#plot-novel-exons) (for proper installation of R dependencies) and then run the following example:
```
plot_novel_exons.R \
--gtf Homo_sapiens.GRCh38.87.chr.support_level_5.correct_gene_coordinates.chr1.14.22.X.16.gtf \
--polyasites polya_sites.merged.anno.hg38.ENSEMBL.chr1.14.22.X.16.bed \
--bam GSM1502499_RNA_seq_control_rep2.chr22.bam \
--tectool_exons results/classified_as_terminal_with_probabilities.tsv \
--output_dir plots
```

In the output directory a pdf file is generated with the identified novel terminal exons. An example can be found here with two novel terminal exons: http://tectool.unibas.ch/data/example_plots.pdf . The novel exons are marked with red boxes. The closest upstream and downstream exons are shown in the plot.

**Note for docker**: You can copy files from a docker container to your host OS using the command.

```
docker cp <container id>:/path/to/file/in/container /path/to/file/in/host/os
```

## Annotation files

In the following links you can find genome files for GRCh38 and GRCm38, their corresponding poly(A) sites and ENSEMBL annotation files wiht support level 1 and support level 5.
* Humnan, GRCh38, Ensembl version 87, support level 1: http://tectool.unibas.ch/data/GRCh38.87_support_level_1.tar.gz
* Humnan, GRCh38, Ensembl version 87, support level 5: http://tectool.unibas.ch/data/GRCh38.87_support_level_5.tar.gz
* Mouse, GRCm38, Ensembl version 87, support level 1: http://tectool.unibas.ch/data/GRCm38.87_support_level_1.tar.gz
* Mouse, GRCm38, Ensembl version 87, support level 5: http://tectool.unibas.ch/data/GRCm38.87_support_level_5.tar.gz

**Note**: The script tectool_filter_gtf_by_transcript_support_level in the scritps directory of TECtool takes as input a gtf file and keeps only the transcripts with the user specified support level. You can run it as following:
```
tectool_filter_gtf_by_transcript_support_level \
--gtf <ORIGINAL GTF FILE> \
--out <FILTERED GTF FILE> \
--support_level <Transcript support level to choose [1,2,3,4,5] \
--fix_gene_coordinates
```


## Licence and documentation

* Free software: MIT license
* Official website: http://tectool.unibas.ch/
