import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 4):
    sys.exit('Sorry, TECtool requires Python >= 3.4 (Tested with Python 3.6.2)')

requirements = [
    "scipy==0.19",
    "htseq==0.9.1",
    "pybedtools==0.7.10",
    "pyfasta==0.5.2",
    "scikit-learn==0.19.0",
    "pandas==0.20.3",
    "matplotlib==2.0.2",
]

setup(
    name='tectool',
    version='0.2',
    description="TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.",
    author="Foivos Gypas",
    author_email='foivos.gypas@unibas.ch',
    url='https://git.scicore.unibas.ch/zavolan_public/TECtool.git',
    packages=['tectool'],
    package_dir={'tectool':'tectool'},
    include_package_data=True,
    scripts=['scripts/tectool'],
    install_requires=requirements,
    keywords='tectool',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
