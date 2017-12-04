#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import os

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='tectool',
    version='0.1.0',
    description="Terminal exon characterization",
    long_description=readme + os.linesep + os.linesep + history,
    author="Foivos Gypas",
    author_email='foivos.gypas@unibas.ch',
      url='https://git.scicore.unibas.ch/zavolan_public/TECtool.git',
    packages=[
        'tectool',
    ],
    package_dir={'tectool':
                 'tectool'},
    entry_points={
        'console_scripts': [
            'tectool=tectool.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='tectool',
    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
