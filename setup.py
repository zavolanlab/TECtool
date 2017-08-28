#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

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
    long_description=readme + '\n\n' + history,
    author="Andreas Gruber",
    author_email='andreas.j.gruber@gmail.com',
      url='https://git@git.scicore.unibas.ch:2222/zavolan_public/TECtool',
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
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
