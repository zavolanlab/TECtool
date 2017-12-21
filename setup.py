import sys
from setuptools import setup

if sys.version_info < (3, 6):
    sys.exit('Sorry, TECtool requirs Python >= 3.6 (Tested with Python 3.6.2)')

requirements = [
    "htseq>=0.9.1",
    "pybedtools>=0.7.10",
    "bzip2",
    "pyfasta>=0.5.2",
    "scikit-learn>=0.19.0",
]

setup(
    name='TECtool',
    version='0.2',
    description="TECtool is a method that uses mRNA and 3’ end sequencing data to identify novel terminal exons.",
    author="Foivos Gypas",
    author_email='foivos.gypas@unibas.ch',
    url='https://git.scicore.unibas.ch/zavolan_public/TECtool.git',
    packages=['tectool'],
    entry_points={
        'console_scripts': [
            'tectool=tectool.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='TECtool',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
