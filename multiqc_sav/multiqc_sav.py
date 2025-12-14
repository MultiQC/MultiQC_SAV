"""MultiQC SAV plugin configuration and hooks."""

import logging
from importlib.metadata import version

from multiqc.utils import config

log = logging.getLogger("multiqc")

# Save this plugin's version number to the MultiQC config
config.multiqc_sav_version = version("multiqc_sav")
log.info(f"Running MultiQC SAV Plugin v{config.multiqc_sav_version}")


def update_config() -> None:
    """
    Update MultiQC config object.

    - Update module order
    - Disable unnecessary modules to avoid duplicate data
    - Update search patterns
    """
    log.debug("SAV - Updating config")

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
        config.update_dict(
            config.sp, {"SAV/xml": {"fn_re": ".*([Rr]un[Ii]nfo|[Rr]un[Pp]arameters)\\.xml", "shared": True}}
        )
    config.update_dict(config.sp, {"bclconvert/runinfo": {"fn": "RunInfo.xml", "shared": True}})
