#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_sav_version = get_distribution("multiqc_sav").version
log.info("Running MultiQC SAV Plugin v{}".format(config.multiqc_sav_version))


def update_defaults():
    log.info("SAV: Updating search patterns")
    # Update search patterns
    if "sav/runinfo" not in config.sp:
        config.update_dict(config.sp, {"SAV/runinfo": {"fn": "RunInfo.xml"}})
    if "sav/runparameters" not in config.sp:
        config.update_dict(
            config.sp, {"SAV/runparameters": {"fn": "RunParameters.xml"}}
        )
