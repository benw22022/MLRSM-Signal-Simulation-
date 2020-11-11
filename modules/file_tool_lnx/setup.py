# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:09:13 2020

@author: benwi
"""

# -*- coding: utf-8 -*-

import os
from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='file_tools_lnx',
    version='1.0.0',
    author='B.Wilson',
    author_email='benjamin.james.wilson@cern.ch',
    license='',
    packages=['d:\\PhD\Work\\Majorana\\signal_simulation\\modules\\file_tools_lnx'],
    description='Some useful tools for dealing with files',
    long_description='Useful tools for creating root files from csv',
    keywords='physics HEP, computing particle',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Physicists',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: ',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.5.5',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=requirements,
    zip_safe=False
)