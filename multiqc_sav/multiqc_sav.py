#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_sav_version = get_distribution("multiqc_sav").version
log.info("Running MultiQC SAV Plugin v{}".format(config.multiqc_sav_version))


def update_config():
    log.info("SAV: Updating search patterns")
    # Update search patterns
    if "SAV/xml" not in config.sp:
        config.update_dict(config.sp, {"SAV/xml": {"fn_re": ".*(RunInfo|RunParameters)\.xml"}})
