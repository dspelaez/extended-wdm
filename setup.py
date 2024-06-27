#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

"""
@author: Daniel Peláez-Zapata
@github: http://github.com/dspelaez
@created: 2024-03-09
"""

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="ewdm",
    version="0.1",
    author="Daniel Peláez-Zapata",
    author_email="daniel.pelaez-zapata@ucdconnect.ie",
    description="EWDM: A package for a wavelet-based directional wave spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dspelaez/ewdm",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "xarray",
        "matplotlib",
        "tqdm"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)

