import sys
from setuptools import setup

if sys.version_info < (3, 4):
    sys.exit('Sorry, TECtool requires Python >= 3.4 (Tested with Python 3.6.2)')

requirements = [
    "numpy>=0.13",
    "scipy==0.19",
    "htseq==0.9.1",
    "pybedtools==0.7.10",
    "pyfasta==0.5.2",
    "scikit-learn==0.19.0",
    "pandas==0.20.3",
    "matplotlib==2.0.2",
    "progress==1.3"
]

setup(
    name='tectool',
    version='0.4',
    description="TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.",
    author="Foivos Gypas",
    author_email='foivos.gypas@unibas.ch',
    url='https://github.com/zavolanlab/TECtool.git',
    packages=['tectool'],
    package_dir={'tectool': 'tectool'},
    include_package_data=True,
    scripts=['scripts/tectool',
             'scripts/tectool_add_novel_transcripts_to_gtf_file',
             'scripts/tectool_create_mapping_lists',
             'scripts/tectool_extract_terminal_exons',
             'scripts/tectool_filter_gtf_by_transcript_support_level',
             'scripts/tectool_quantify_split_reads',
             'scripts/plot_novel_exons.R'],
    install_requires=requirements,
    keywords='tectool',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
