import logging
import re
from collections import OrderedDict

import interop
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils import config, report


class MultiQC_SAV(BaseMultiqcModule):
    def __init__(self):

        self.log = logging.getLogger("multiqc")

        self.log.debug("Running MultiQC SAV plugin")

        # Check if the plugin has an input directory
        if not config.kwargs.get("illumina_dir", None):
            self.log.info("Skipping MultiQC_SAV")
            return None

        super(MultiQC_SAV, self).__init__(name="Illumina SAV", anchor="sav")

        self.intro = """<p>The <a href="https://github.com/CenterForMedicalGeneticsGhent/MultiQC_SAV" target="_blank">MultiQC SAV</a> 
                        plugin replicates the most used tables and graphs from SAV</p>"""

