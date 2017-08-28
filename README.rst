===============================
TECtool
===============================

.. image:: https://img.shields.io/pypi/v/tectool.svg
        :target: https://pypi.python.org/pypi/tectool

.. image:: https://img.shields.io/travis/fgypas/tectool.svg
        :target: https://travis-ci.org/fgypas/tectool

.. image:: https://readthedocs.org/projects/tectool/badge/?version=latest
        :target: https://tectool.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/fgypas/cookiecutter-django/shield.svg
     :target: https://pyup.io/repos/github/fgypas/tectool/
     :alt: Updates

Terminal exon characterization Tool
-----------------------------------

* Free software: MIT license
* Documentation: https://tectool.readthedocs.io (coming soon...)

Features
--------

TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.
* Identify novel terminal exons
* Infer novel transcripts
* Annotate CDS for novel transcripts

Input and output files
--------

Input files
* A file containing all chromosomes in fasta format
* A file with the corresponding annotation in GTF format
* A file with genome coordinates of 3’ end processing sites in BED format
* A file containing spliced alignments of mRNA-seq reads to the corresponding genome (in BAM format, sorted by coordinates and indexed)

Output files
* An augmented annotation file (in GTF format)
* A file containing the novel terminal exons with the corresponding probability and the features that were calculated

Instructions
--------

#### Step 1: Download Miniconda3


On Linux:

.. code:: bash

  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh

On MacOS X:


.. code:: bash

  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
