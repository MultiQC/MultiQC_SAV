# CLAUDE.md

This file provides guidance for Claude Code (claude.ai/code) when working with the MultiQC_SAV plugin.

## Project Overview

MultiQC_SAV is a plugin for [MultiQC](https://multiqc.info/) that parses InterOp files from Illumina sequencers and generates tables and graphs for quality control metrics. It leverages the [Illumina InterOp Python API](https://github.com/Illumina/interop) to read binary metric files.

## Project Structure

```
MultiQC_SAV/
├── multiqc_sav/                    # Main plugin package
│   ├── __init__.py                 # Package metadata
│   ├── multiqc_sav.py              # Plugin hook (before_config)
│   └── modules/
│       ├── __init__.py
│       └── sav/
│           ├── __init__.py
│           └── sav.py              # Main SAV module implementation
├── test_data/                      # Test datasets for various sequencers
│   ├── HiSeq3000/
│   ├── MiSeq/
│   ├── MiSeqI100/
│   ├── NextSeq500/
│   ├── NextSeq2000/
│   ├── NovaSeq6000/
│   └── NovaSeqX/
├── .devcontainer/                  # VS Code devcontainer configuration
├── .github/workflows/              # CI/CD workflows
│   ├── linux.yaml                  # Build and test workflow
│   ├── lint.yaml                   # Linting workflow (pre-commit)
│   └── publish.yaml                # PyPI publish workflow
├── pyproject.toml                  # Project configuration and dependencies
├── .pre-commit-config.yaml         # Pre-commit hooks configuration
├── .prettierrc.js                  # Prettier formatter config
└── .prettierignore                 # Prettier exclusions
```

## Development Setup

```bash
# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install

# Or use pixi (if available)
pixi install
```

## Common Commands

```bash
# Run linting
ruff check .
ruff format --check .

# Auto-fix linting issues
ruff check --fix .
ruff format .

# Run prettier on markdown/yaml
prettier --check "**/*.{md,yaml,yml,json}"

# Run tests with test data
multiqc --strict -v --no-version-check -m SAV test_data/MiSeq
multiqc --strict -v --no-version-check -m SAV test_data/HiSeq3000
multiqc --strict -v --no-version-check -m SAV test_data/NextSeq500
multiqc --strict -v --no-version-check -m SAV test_data/NextSeq2000
multiqc --strict -v --no-version-check -m SAV test_data/NovaSeq6000
multiqc --strict -v --no-version-check -m SAV test_data/NovaSeqX

# Run all pre-commit hooks
pre-commit run --all-files
```

## Code Style Guidelines

- **Line length**: 120 characters
- **Python version**: 3.9+ (tested on 3.11, 3.12, 3.13)
- **Formatting**: Ruff (format + lint)
- **Lint rules**: E, F, W, I (pycodestyle, pyflakes, isort)
- **Imports**: Sorted by isort (via ruff)

## MultiQC Module Architecture

### Plugin Registration

The plugin uses entry points in `pyproject.toml`:

```toml
[project.entry-points."multiqc.hooks.v1"]
before_config = "multiqc_sav.multiqc_sav:sav_execution_start"

[project.entry-points."multiqc.modules.v1"]
SAV = "multiqc_sav.modules.sav.sav:SAVModule"
```

### Hook System

`multiqc_sav.py` contains the `sav_execution_start()` hook which runs before config loading:

- Logs the plugin version
- Registers the SAV module in the module order with tags (DNA, RNA, BCL, Demultiplex)
- Disables the built-in InterOp module to avoid duplicate data
- Configures search patterns:
  - `SAV/RunInfo`: `RunInfo.xml`
  - `SAV/RunParameters`: `RunParameters.xml`
  - `SAV/InterOp`: `InterOp/*.bin`
  - `bclconvert/runinfo`: shared `RunInfo.xml`

### Main Module

`modules/sav/sav.py` contains the `SAVModule` class which extends `BaseMultiqcModule`:

1. **Initialization**: Minimal class that sets up module name and anchor
2. **InterOp Processing**: The `add_interop_sections()` function processes InterOp data
3. **Metrics Loading**: Uses InterOp API to read binary metric files
4. **Data Processing**: Parses metrics into pandas DataFrames
5. **Visualization**: Generates MultiQC plots (tables, bargraphs, heatmaps, linegraphs, scatter plots)

### Key Functions

- `add_interop_sections(module)`: Main entry point for adding InterOp visualizations
- `_add_summary_sections()`: Read and lane summary tables
- `_add_qscore_sections()`: Q-score heatmaps and histograms
- `_add_imaging_sections()`: Imaging metrics visualizations

## Supported Sequencers

- MiSeq
- MiSeq (Illumina Connected)
- HiSeq 3000/4000
- NextSeq 500/550
- NextSeq 1000/2000
- NovaSeq 6000
- NovaSeq X/X Plus

## Key Dependencies

- `interop>=1.7.0,<2` - Illumina InterOp Python API for reading binary metrics
- `multiqc>=1.25` - MultiQC framework
- `pandas` - Data manipulation
- `numpy` - Numerical operations

## Testing

Tests are run via GitHub Actions on Python 3.11, 3.12, and 3.13. Each test:
1. Runs MultiQC with the SAV module on test data
2. Verifies that `multiqc_report.html` is generated
3. Checks that the SAV module appears in the log

## CI/CD

- **Build workflow** (`linux.yaml`): Tests installation and module execution with report validation
- **Lint workflow** (`lint.yaml`): Runs pre-commit hooks (ruff, prettier)
- **Publish workflow** (`publish.yaml`): Publishes to PyPI on release using trusted publishing
