# CLAUDE.md

This file provides guidance for Claude Code (claude.ai/code) when working with the MultiQC_SAV plugin.

## Project Overview

MultiQC_SAV is a plugin for [MultiQC](https://multiqc.info/) that parses InterOp files from Illumina sequencers and generates tables and graphs for quality control metrics. It leverages the [Illumina InterOp Python API](https://github.com/Illumina/interop) to read binary metric files.

## Project Structure

```
MultiQC_SAV/
├── multiqc_sav/                    # Main plugin package
│   ├── __init__.py                 # Package identifier
│   ├── multiqc_sav.py              # Plugin hook/config registration
│   └── modules/
│       └── sav.py                  # Main SAV module implementation
├── test_data/                      # Test datasets for various sequencers
│   ├── HiSeq/
│   ├── MiSeq/
│   ├── NextSeq500/
│   ├── NextSeq2000/
│   └── NovaSeq/
├── docs/                           # Documentation and example reports
├── .github/workflows/              # CI/CD workflows
│   ├── linux.yaml                  # Build and test workflow
│   └── lint.yaml                   # Linting workflow
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
multiqc -m SAV test_data/MiSeq
multiqc -m SAV test_data/HiSeq
multiqc -m SAV test_data/NextSeq500
multiqc -m SAV test_data/NextSeq2000
multiqc -m SAV test_data/NovaSeq

# Run all pre-commit hooks
pre-commit run --all-files
```

## Code Style Guidelines

- **Line length**: 120 characters
- **Python version**: 3.9+
- **Formatting**: Ruff (format + lint)
- **Type hints**: Required for all function signatures
- **Docstrings**: Required for all public functions
- **Imports**: Sorted by isort (via ruff)

## MultiQC Module Architecture

### Plugin Registration

The plugin uses entry points in `pyproject.toml`:

```toml
[project.entry-points."multiqc.hooks.v1"]
config_loaded = "multiqc_sav.multiqc_sav:update_config"

[project.entry-points."multiqc.modules.v1"]
SAV = "multiqc_sav.modules.sav:SAV"
```

### Hook System

`multiqc_sav.py` contains the `update_config()` hook which:

- Registers the SAV module in the module order
- Sets module tags (DNA, RNA, BCL, Demultiplex)
- Disables the built-in InterOp module to avoid duplicate data
- Configures search patterns for RunInfo.xml and RunParameters.xml

### Main Module

`modules/sav.py` contains the `SAV` class which extends `BaseMultiqcModule`:

1. **File Discovery**: Uses `find_log_files("SAV/xml")` to locate XML files
2. **Run Info Parsing**: Extracts metadata from RunInfo.xml
3. **Metrics Loading**: Uses InterOp API to read binary metric files
4. **Data Processing**: Parses metrics into pandas DataFrames
5. **Visualization**: Generates MultiQC plots (tables, bargraphs, heatmaps, linegraphs, scatter plots)

## Supported Sequencers

- MiSeq
- HiSeq 3000/4000
- NextSeq 500/2000
- NovaSeq 6000

## Key Dependencies

- `interop>=1.1.23` - Illumina InterOp Python API for reading binary metrics
- `multiqc>=1.10` - MultiQC framework
- `pandas` - Data manipulation
- `numpy` - Numerical operations

## Testing

Tests are run via GitHub Actions on Python 3.9-3.12. Each test verifies that the module can process test data without errors by running MultiQC with the SAV module flag.

## CI/CD

- **Build workflow** (`linux.yaml`): Tests installation and module execution
- **Lint workflow** (`lint.yaml`): Runs ruff format check, ruff lint check, and prettier
