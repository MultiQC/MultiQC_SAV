"""
MultiQC SAV Plugin - Adds InterOp-based visualizations to the core SAV module.

This plugin extends the core SAV module with advanced visualizations by parsing
Illumina InterOp binary files.
"""

import logging

from multiqc_sav.sav_interop import add_interop_sections

log = logging.getLogger(__name__)


def sav(module):
    """
    Plugin hook called by the core SAV module.

    This function adds InterOp-based visualizations to the SAV module.
    It receives the module instance and can:
    - Add new sections with plots
    - Augment data_by_sample with additional metrics
    - Add general stats columns

    Args:
        module: The SAV MultiqcModule instance
    """

    add_interop_sections(module)
