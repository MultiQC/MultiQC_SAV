import logging

import importlib_metadata
from multiqc import config
from multiqc.utils.util_functions import update_dict

log = logging.getLogger("multiqc")


def sav_execution_start():
    # Plugin's version number defined in pyproject.toml:
    version = importlib_metadata.version("multiqc_sav")
    log.debug(f"Running MultiQC SAV Plugin v{version}")

    log.debug("SAV - Updating config")
    # Add module to module order
    config.module_order.append({"SAV": {"module_tag": ["DNA", "RNA", "BCL", "Demultiplex"]}})

    # Move module to the top
    config.top_modules.append("SAV")

    # Disable InterOp module to avoid duplicate data
    disabled_modules = ["interop"]
    for module in disabled_modules:
        del config.avail_modules[module]
    log.debug("SAV - Disabled modules: {}".format(", ".join(disabled_modules)))

    # Update search patterns
    # Set RunInfo to shared for the bclconvert module
    update_dict(config.sp, {"bclconvert/runinfo": {"fn": "RunInfo.xml", "shared": True}})
    # Set SAV file search patterns
    update_dict(
        config.sp,
        {
            "SAV/RunInfo": {"fn": "RunInfo.xml", "shared": True},
            "SAV/RunParameters": {"fn": "RunParameters.xml", "shared": True},
            "SAV/InterOp": {"fn_re": "InterOp/.*\\.bin"},
        },
    )
