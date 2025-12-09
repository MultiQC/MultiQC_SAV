#!/usr/bin/env python
"""
MultiQC_SAV is a plugin for MultiQC that adds InterOp-based visualizations
to the core SAV module, providing detailed sequencing metrics similar to
Illumina's Sequencing Analysis Viewer application.
"""

from setuptools import setup, find_packages

version = "0.1.0"

setup(
    name="multiqc_sav",
    version=version,
    author="Matthias De Smet",
    author_email="11850640+matthdsm@users.noreply.github.com",
    description="MultiQC plugin for Illumina SAV InterOp visualizations",
    long_description=__doc__,
    keywords="bioinformatics illumina sequencing multiqc",
    url="https://github.com/MultiQC/MultiQC_SAV",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=[
        "interop>=1.1.23",
        "multiqc>=1.25",
        "pandas",
        "numpy",
    ],
    entry_points={
        "multiqc.hooks.v1": [
            "sav_extra = multiqc_sav.multiqc_sav:sav_extra_hook",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
