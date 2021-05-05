#!/usr/bin/env python
"""
MultiQC_SAV is a plugin for MultiQC, leveraging the [InterOp python API](https://github.com/Illumina/interop)
to generate the most used tables and graphs from the Illumina SAV
"""

from setuptools import setup, find_packages

version = "0.1"

setup(
    name="multiqc_sav",
    version=version,
    author="Matthias De Smet",
    author_email="11850640+matthdsm@users.noreply.github.com",
    description="MultiQC plugin for Illumina SAV Graphs and tables",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/CenterForMedicalGeneticsGhent/MultiQC_SAV",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["interop>=1.1.23", "multiqc>=1.10"],
    entry_points={
        "multiqc.hooks.v1": ["after_modules = multiqc_sav.multiqc_sav:MultiQC_SAV",],
        "multiqc.cli_options.v1": ["illumina_dir = multiqc_sav.cli:illumina_dir",],
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
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)

