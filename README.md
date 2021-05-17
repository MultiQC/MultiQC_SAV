# [<img src="docs/images/MultiQC_logo.png" width="250" title="MultiQC">](https://github.com/ewels/MultiQC)

**MultiQC_SAV is a plugin for MultiQC, leveraging the [InterOp python API](https://github.com/Illumina/interop) to generate the most used tables and graphs from the Illumina SAV**

For more information about MultiQC, see [http://multiqc.info](http://multiqc.info)

For more information about Illumina SAV, see [the Illumina support page](https://support.illumina.com/sequencing/sequencing_software/sequencing_analysis_viewer_sav/downloads.html)

## Description

The MultiQC_SAV plugin parses the InterOp files in an Illumina Sequencer output directory and generates tables and graphs of the most important metrics

## Installation

This plugin can be installed using the following methods

- using `pip`:

```bash
pip install --upgrade --force-reinstall git+https://github.com/CenterForMedicalGeneticsGhent/MultiQC_SAV.git
```

- using `conda`:

```bash
conda install -c bioconda multiqc_sav
```

- using `setuptools`:

```bash
git clone https://github.com/CenterForMedicalGeneticsGhent/MultiQC_SAV
cd MultiQC_SAV
python setup.py install
```

## Usage

This plugin adds a QC module and searches for the `RunInfo.xml` and `RunParameters.xml` files. It tries to infer if all required files are present. No special params are needed.

### Required files

The illumina directory should contain the following files the directory structure as dictated by the sequencer

```bash
illumina_dir
├── InterOp
│   ├── *.bin
├── RunInfo.xml
├── RunParameters.xml

```
