#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_sav_version = get_distribution("multiqc_sav").version
log.info("Running MultiQC SAV Plugin v{}".format(config.multiqc_sav_version))
