# MultiQC_SAV Changelog

## [0.2.0](https://github.com/MultiQC/MultiQC_SAV/releases/tag/v0.2.0) - 2024-12-18

This release includes a major refactor to follow MultiQC best practices and modernize the codebase.

### Major Changes

- Converted to hook-based plugin architecture ([#7](https://github.com/MultiQC/MultiQC_SAV/pull/7))
- Refactored SAV module to follow MultiQC best practices ([#11](https://github.com/MultiQC/MultiQC_SAV/pull/11))

### Added

- Comprehensive type hints throughout the codebase
- VS Code devcontainer configuration for easier development
- CI/CD tooling with GitHub Actions workflows
- Developer documentation (CLAUDE.md)
- mypy configuration for static type checking
- New test data for additional sequencer platforms (NextSeq 2000, NovaSeq X)

### Changed

- Migrated from `setup.py` to modern `pyproject.toml` packaging
- Applied ruff UP and SIM rules for cleaner code
- Removed Python 3.13 upper bound constraint
- Updated MultiQC logo

### Fixed

- Search patterns configuration
- Sample name ignore list now correctly prevents report generation when all samples are filtered

## [0.0.3](https://github.com/MultiQC/MultiQC_SAV/releases/tag/v0.0.3) - 2022-02-22

Fix parsing of date notation in NSQ2K `RunInfo.xml`

## [0.0.2](https://github.com/MultiQC/MultiQC_SAV/releases/tag/v0.0.2) - 2021-06-09

This release improves the filename matching and updates the search patterns for the SAV module.

## [0.0.1](https://github.com/MultiQC/MultiQC_SAV/releases/tag/v0.0.1) - 2021-05-19

Initial release of the MultiQC SAV plugin.

### Added

- InterOp-based visualizations for Illumina sequencing metrics
- Read and lane summary tables
- Q-score heatmap and histogram
- Clusters per lane plot
- Occupancy vs PF (pass filter) plot
- By-cycle metrics visualization
- Support for MiSeq, HiSeq, NextSeq, and NovaSeq platforms
- Test data for multiple sequencer types
- GitHub Actions CI workflow
