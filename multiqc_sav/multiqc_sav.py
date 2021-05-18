#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_sav_version = get_distribution("multiqc_sav").version
log.info("Running MultiQC SAV Plugin v{}".format(config.multiqc_sav_version))


def update_config() -> None:
    """
    Update MultiQC config object
    * Update module order
    * Disable unnecessary modules to avoid duplicate data
    * Update search patterns
    """

    log.debug("Updating config")
    # Add module to module order
    config.module_order.append({"SAV": {"module_tag": ["DNA", "RNA", "BCL", "Demultiplex"]}})

    # Move module to the top
    config.top_modules.append("SAV")

    # Disable InterOp module to avoid duplicate data
    disabled_modules = ["interop"]
    for module in disabled_modules:
        del config.avail_modules[module]

    # Update search patterns
    if "SAV/xml" not in config.sp:
        config.update_dict(config.sp, {"SAV/xml": {"fn_re": ".*(RunInfo|RunParameters)\.xml"}})
