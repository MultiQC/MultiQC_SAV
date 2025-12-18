"""
InterOp-based visualizations for the SAV module.

This module provides functions to parse Illumina InterOp binary files
and generate advanced visualizations.
"""

import contextlib
import glob
import logging
import os
import re
from typing import Any, Optional

import interop
import numpy as np
import pandas as pd
from interop import py_interop_plot
from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap, linegraph, scatter, table
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

# Table headers for summary metrics
HEADERS: dict[str, dict] = {
    "Error Rate": {
        "title": "Error Rate (%)",
        "description": "The calculated error rate, as determined by a PhiX spike-in",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "OrRd",
    },
    "Error Rate 35": {
        "title": "Error Rate 35 (%)",
        "description": "The calculated error rate for cycles 1-35.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "OrRd",
        "hidden": True,
    },
    "Error Rate 50": {
        "title": "Error Rate 50 (%)",
        "description": "The calculated error rate for cycles 1-50.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "OrRd",
        "hidden": True,
    },
    "Error Rate 75": {
        "title": "Error Rate 75 (%)",
        "description": "The calculated error rate for cycles 1-75.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "OrRd",
        "hidden": True,
    },
    "Error Rate 100": {
        "title": "Error Rate 100 (%)",
        "description": "The calculated error rate for cycles 1-100.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "OrRd",
        "hidden": True,
    },
    "First Cycle Intensity": {
        "title": "Intensity Cycle 1",
        "description": "The average of the A channel intensity measured at the first cycle",
        "scale": "Blues",
    },
    "% Aligned": {
        "title": "Aligned (%)",
        "description": "Percentage of reads that aligned to the PhiX genome",
        "suffix": "%",
        "min": 0,
        "max": 100,
        "format": "{:,.2f}",
        "scale": "RdYlGn",
    },
    "% >= Q30": {
        "title": "% >= Q30",
        "description": "Percentage of reads with quality phred score of 30 or above",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "RdYlGn",
    },
    "% Occupancy Proxy": {
        "title": "Occupancy Proxy (%)",
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "Blues",
    },
    "% Occupied": {
        "title": "Occupied (%)",
        "description": "The percentage of nanowells occupied by clusters, +/- 1 standard deviation.",
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "Blues",
    },
    "Projected Yield G": {
        "title": f"Projected Yield ({config.base_count_prefix})",
        "description": f"The expected number of bases sequenced ({config.base_count_desc} base pairs over all 'usable cycles'",
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0) * config.base_count_multiplier,
        "scale": "Greens",
        "hidden": True,
    },
    "Yield G": {
        "title": f"Yield ({config.base_count_prefix})",
        "description": f"The number of bases sequenced ({config.base_count_desc} base pairs over all 'usable cycles'",
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0) * config.base_count_multiplier,
        "scale": "Greens",
    },
    "Cluster Count": {
        "title": f"Clusters ({config.read_count_prefix})",
        "description": f"Number of clusters for each tile ({config.read_count_desc})",
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
        "scale": "Blues",
    },
    "Cluster Count Pf": {
        "title": f"Clusters PF ({config.read_count_prefix})",
        "description": f"Number of clusters PF for each tile ({config.read_count_desc})",
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
        "scale": "Blues",
    },
    "% Pf": {
        "title": "Reads PF (%)",
        "description": "Percentage of clusters Passing Filter",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "scale": "RdYlGn",
    },
    "Density": {
        "title": "Density",
        "description": "The density of clusters (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.",
        "scale": "Purples",
        "hidden": True,
    },
    "Density Pf": {
        "title": "Density PF",
        "description": "The density of clusters PF (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.",
        "scale": "Purples",
        "hidden": True,
    },
    "Phasing": {
        "title": "Phasing",
        "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) the current cycle within a read.",
        "scale": "OrRd",
    },
    "Phasing Offset": {
        "title": "Phasing Offset",
        "description": "The best-fit offset of the phasing corrections, calculated from the entire read.",
        "scale": "Greys",
        "hidden": True,
    },
    "Phasing Slope": {
        "title": "Phasing Slope",
        "description": "The best-fit slope of the phasing corrections, calculated from the entire read.",
        "scale": "Greys",
        "hidden": True,
    },
    "Prephasing": {
        "title": "Prephasing",
        "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing jumps ahead (prephasing) the current cycle within a read.",
        "scale": "OrRd",
    },
    "Prephasing Offset": {
        "title": "Prephasing Offset",
        "description": "The best-fit offset of the prephasing corrections, calculated from the entire read.",
        "scale": "Greys",
        "hidden": True,
    },
    "Prephasing Slope": {
        "title": "Prephasing Slope",
        "description": "The best-fit slope of the prephasing corrections, calculated from the entire read.",
        "scale": "Greys",
        "hidden": True,
    },
    "Reads": {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"The number of reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
        "scale": "Blues",
    },
    "Reads Pf": {
        "title": f"{config.read_count_prefix} PF Reads",
        "description": f"The number of passing filter reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
        "scale": "Blues",
    },
    "Tile Count": {
        "title": "Tiles",
        "description": "The number of tiles per lane.",
        "scale": "Greys",
        "hidden": True,
    },
    "Total Pf Reads": {
        "title": f"{config.read_count_prefix} PF Reads",
        "description": f"The total number of passing filter reads for this lane ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
        "scale": "Blues",
    },
    "Total Reads": {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"The total number of reads for this lane ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
        "scale": "Blues",
    },
    "Mapped Reads Cv": {
        "title": "CV",
        "description": "The coefficient of variation for the number of counts across all indexes.",
        "format": "{:,.2f}",
        "scale": "Oranges",
    },
    "Max Mapped Reads": {
        "title": f"{config.read_count_prefix} Max Mapped Reads",
        "description": f"The highest representation for any index ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
        "scale": "Blues",
    },
    "Min Mapped Reads": {
        "title": f"{config.read_count_prefix} Min Mapped Reads",
        "description": f"The lowest representation for any index ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
        "scale": "Blues",
    },
    "Total Fraction Mapped Reads": {"hidden": True},
    "Fraction Mapped": {"hidden": True},
    "Index1": {
        "title": "Index 1 (I7)",
        "description": "The sequence for the first Index Read.",
    },
    "Index2": {
        "title": "Index 2 (I5)",
        "description": "The sequence for the second Index Read",
    },
    "Project Name": {
        "title": "Project Name",
        "description": "Sample Project Name",
    },
    "Sample Id": {
        "title": "Sample ID",
        "description": "The Sample ID given in the SampleSheet",
    },
}


class SAVModule(BaseMultiqcModule):
    """
    SAV (Sequencing Analysis Viewer) module for Illumina InterOp metrics.

    This module parses InterOp binary files from Illumina sequencing runs
    and generates quality metrics similar to Illumina's SAV application.

    Supported sequencers:

    - MiSeq / MiSeq (Illumina Connected)
    - HiSeq 3000/4000
    - NextSeq 500/550 / NextSeq 1000/2000
    - NovaSeq 6000 / NovaSeq X/X Plus

    Required files:

    - `RunInfo.xml`
    - `InterOp/*.bin` files
    """

    def __init__(self) -> None:
        super().__init__(
            name="SAV",
            anchor="sav",
            href="https://github.com/Illumina/interop",
            info="Parses Illumina InterOp binary files for sequencing run metrics.",
            # No DOI for the InterOp library
        )

        # Find run directories with RunInfo.xml
        self.sav_data: dict = {}
        for f in self.find_log_files("SAV/RunInfo"):
            run_dir = os.path.dirname(f["root"])
            run_name = os.path.basename(run_dir) if run_dir else f["s_name"]

            # Check for InterOp directory
            interop_dir = os.path.join(f["root"], "InterOp")
            if not os.path.isdir(interop_dir):
                log.debug(f"No InterOp directory found for {run_name}")
                continue

            interop_files = glob.glob(os.path.join(interop_dir, "*.bin"))
            if not interop_files:
                log.debug(f"No InterOp .bin files found for {run_name}")
                continue

            if run_name in self.sav_data:
                log.debug(f"Duplicate run name found! Overwriting: {run_name}")

            self.add_data_source(f, s_name=run_name)
            self.add_software_version(None, sample=run_name)
            self.sav_data[run_name] = {"_run_dir": f["root"]}

        # Filter ignored samples
        self.sav_data = self.ignore_samples(self.sav_data)

        if len(self.sav_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.sav_data)} sequencing run(s)")

        # Add InterOp sections
        add_interop_sections(self)

        # Write data file at the END
        self.write_data_file(self.sav_data, "multiqc_sav")


def add_interop_sections(module: SAVModule) -> None:
    """
    Add InterOp-based sections to the SAV module.

    Args:
        module: The SAV MultiqcModule instance
    """
    # Process each run found
    for run_name, data in module.sav_data.items():
        run_dir = data.get("_run_dir")
        if not run_dir:
            continue

        log.info(f"Loading InterOp metrics from {run_dir}")

        try:
            run_metrics = interop.read(
                run=run_dir,
                valid_to_load=interop.load_imaging_metrics(),
                finalize=True,
            )
        except OSError as e:
            log.warning(f"Failed to load InterOp metrics from {run_dir}: {e}")
            continue

        # Add summary sections
        _add_summary_sections(module, run_metrics, run_name)

        # Add Q-score sections
        _add_qscore_sections(module, run_metrics, run_name)

        # Add imaging sections
        _add_imaging_sections(module, run_metrics, run_name)

        # Store summary data for general stats
        _store_summary_data(module, run_metrics, run_name)

    # Add general stats after processing all runs
    _add_general_stats(module)


def _add_summary_sections(module: SAVModule, run_metrics: Any, run_name: str) -> None:
    """Add read and lane summary table sections."""
    log.info(f"Gathering summary metrics for {run_name}")

    try:
        summary_read = pd.DataFrame(interop.summary(run_metrics, level="Read"))
        summary_nonindex = pd.DataFrame(interop.summary(run_metrics, level="NonIndex"))
        summary_total = pd.DataFrame(interop.summary(run_metrics, level="Total"))

        read_data = _parse_read_summary(summary_read, summary_nonindex, summary_total)

        if read_data:
            headers = {h: HEADERS[h] for h in interop.summary_columns(level="Lane") if h in HEADERS}
            module.add_section(
                name="Summary Read Metrics",
                anchor="sav-read-summary",
                description="Summary metrics per Read",
                plot=table.plot(
                    read_data,
                    headers if isinstance(headers, dict) else None,  # type: ignore[arg-type]
                    pconfig={
                        "title": "SAV: Read Metrics Summary",
                        "namespace": "SAV",
                        "id": "sav-read-metrics-summary-table",
                        "col1_header": "Read",
                    },
                ),
            )
    except (ValueError, TypeError) as e:
        log.debug("Could not generate read summary: %s", e)

    try:
        summary_lane = pd.DataFrame(interop.summary(run_metrics, level="Lane"))
        lane_data = _parse_lane_summary(summary_lane)

        if lane_data:
            headers = {h: HEADERS[h] for h in interop.summary_columns(level="Lane") if h in HEADERS}
            module.add_section(
                name="Summary Lane Metrics",
                anchor="sav-lane-summary",
                description="Summary metrics per Lane per Read",
                plot=table.plot(
                    lane_data,
                    headers if isinstance(headers, dict) else None,  # type: ignore[arg-type]
                    pconfig={
                        "title": "SAV: Lane Metrics Summary",
                        "namespace": "SAV",
                        "id": "sav-lane-metrics-summary-table",
                        "col1_header": "Lane - Read",
                    },
                ),
            )

            # Add clusters/reads per lane plot
            module.add_section(
                name="Clusters/Reads per Lane",
                anchor="sav-clusters-lane",
                description="Total Cluster/Read count per Lane",
                plot=_clusters_lane_plot(lane_data),
            )
    except (ValueError, TypeError) as e:
        log.debug("Could not generate lane summary: %s", e)


def _add_qscore_sections(module: SAVModule, run_metrics: Any, run_name: str) -> None:
    """Add Q-score heatmap and histogram sections."""
    log.info(f"Generating Q-score plots for {run_name}")

    # Q-score heatmap
    try:
        options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
        rows = py_interop_plot.count_rows_for_heatmap(run_metrics)
        cols = py_interop_plot.count_columns_for_heatmap(run_metrics)

        if rows > 0 and cols > 0:
            data_buffer = np.zeros((rows, cols), dtype=np.float32)
            data = py_interop_plot.heatmap_data()

            with contextlib.suppress(py_interop_plot.invalid_filter_option):
                py_interop_plot.plot_qscore_heatmap(run_metrics, options, data, data_buffer.ravel())

            # data_buffer shape is (rows=qscores, cols=cycles)
            # heatmap expects: rows match ycats, cols match xcats
            plot_data = data_buffer.tolist()
            x_cats = list(range(cols))  # cycles (x-axis)
            y_cats = list(range(rows))  # qscores (y-axis)

            module.add_section(
                name="Qscore Heatmap",
                anchor="sav-qscore-heatmap",
                description="The Qscore Heat Map provides an overview of quality scores across cycles.",
                plot=heatmap.plot(
                    plot_data,
                    xcats=x_cats,
                    ycats=y_cats,
                    pconfig={
                        "id": "sav-qscore-heatmap-plot",
                        "title": "SAV: Qscore Heatmap",
                        "xlab": "Cycle",
                        "ylab": "Qscore",
                        "square": False,
                        "colstops": [
                            [0, "#FFFFFF"],
                            [0.1, "#1a9850"],
                            [0.2, "#66bd63"],
                            [0.3, "#a6d96a"],
                            [0.4, "#d9ef8b"],
                            [0.5, "#ffffbf"],
                            [0.6, "#fee08b"],
                            [0.7, "#fdae61"],
                            [0.8, "#f46d43"],
                            [0.9, "#d73027"],
                            [1, "#a50026"],
                        ],
                    },
                ),
            )
    except (ValueError, TypeError) as e:
        log.debug("Could not generate Q-score heatmap: %s", e)

    # Q-score histogram
    try:
        bar_data = py_interop_plot.bar_plot_data()
        options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
        py_interop_plot.plot_qscore_histogram(run_metrics, options, bar_data)

        hist: dict = {}
        qscore: list = []
        reads: list = []
        binsize: list = []

        for i in range(bar_data.size()):
            qscore = [bar_data.at(i).at(j).x() for j in range(bar_data.at(i).size())]
            reads = [bar_data.at(i).at(j).y() for j in range(bar_data.at(i).size())]
            binsize = [bar_data.at(i).at(j).width() for j in range(bar_data.at(i).size())]

        i = 0
        while i < len(qscore):
            j = 0
            while j < binsize[i]:
                hist[qscore[i] + j] = reads[i]
                j += 1
            i += 1

        if hist:
            plot_data = {bar_data.title(): hist}
            module.add_section(
                name="Qscore Histogram",
                anchor="sav-qscore-histogram",
                description="Qscore Histogram graphs the number of bases by quality score. The quality score is cumulative for the current cycle. Only bases from reads that pass the quality filter are included.",
                plot=linegraph.plot(
                    plot_data,
                    pconfig={
                        "id": "sav-qscore-histogram-plot",
                        "title": "SAV: Qscore Histogram",
                        "xlab": "Qscore",
                        "ylab": "Reads (Billion)",
                    },
                ),
            )
    except (ValueError, TypeError) as e:
        log.debug("Could not generate Q-score histogram: %s", e)


def _add_imaging_sections(module: SAVModule, run_metrics: Any, run_name: str) -> None:
    """Add imaging-related sections (intensity per cycle, % PF vs % occupied)."""
    log.info(f"Gathering imaging metrics for {run_name}")

    try:
        imaging = pd.DataFrame(interop.imaging(run_metrics))
        plot_data = _parse_imaging_table(imaging)

        # Intensity per cycle plot
        if plot_data.get("intensity_cycle"):
            module.add_section(
                name="Intensity per Cycle",
                anchor="sav-intensity-cycle",
                description="Intensity by color and cycle of the 90% percentile of the data for each tile",
                plot=_intensity_cycle_plot(plot_data["intensity_cycle"]),
            )

        # % PF vs % Occupied plot
        if plot_data.get("occ_vs_pf"):
            module.add_section(
                name="% PF vs % Occupied",
                anchor="sav-imaging-pf-vs-occ",
                description="% Clusters passing filter vs % Wells Occupied",
                plot=_occ_vs_pf_plot(plot_data["occ_vs_pf"]),
            )
    except (ValueError, TypeError) as e:
        log.debug("Could not generate imaging plots: %s", e)  # noqa: PD011


def _parse_read_summary(
    read_metrics: pd.DataFrame,
    non_index_metrics: pd.DataFrame,
    total_metrics: pd.DataFrame,
) -> dict:
    """Parse read summary DataFrames into dict format."""
    table_data: dict = _parse_reads(read_metrics)

    for _, data in non_index_metrics.iterrows():
        table_data["Non-Indexed"] = data.to_dict()

    for _, data in total_metrics.iterrows():
        table_data["Total"] = data.to_dict()

    return table_data


def _parse_lane_summary(data: pd.DataFrame) -> dict:
    """Parse lane summary DataFrame into dict format."""
    lanes = data.groupby("Lane")
    table_data: dict = {}

    for lane, reads in lanes:
        reads_dict = _parse_reads(reads, key_prefix=f"Lane {lane}")
        table_data.update(reads_dict)

    return table_data


def _parse_reads(reads_df: pd.DataFrame, key_prefix: Optional[str] = None) -> dict:
    """Utility function to parse a reads DataFrame to dict."""
    reads_dict: dict = {}
    reads_df = reads_df.set_index("ReadNumber")

    for read, data in reads_df.iterrows():
        key = f"Read {read}" + " (I)" if data["IsIndex"] == 89 else f"Read {read}"
        if key_prefix:
            key = f"{key_prefix} - {key}"
        reads_dict[key] = data.drop("IsIndex").to_dict()

    return reads_dict


def _parse_imaging_table(data: pd.DataFrame) -> dict:
    """Parse imaging table DataFrame for intensity and occupancy plots."""
    cscale = mqc_colour.mqc_colour_scale()
    colors = cscale.get_colours("Dark2")

    per_lane = data.groupby("Lane")
    occ_pf: dict = {}
    intensity_cycle: dict = {}

    for lane, lane_data in per_lane:
        lane_int = None
        try:
            lane_int = int(float(lane))  # type: ignore[arg-type]
        except (ValueError, TypeError):
            lane_int = None

        # Determine channel set
        CHANNEL_SETS = [
            {"P90/RED", "P90/GREEN"},
            {"P90/Red", "P90/Green"},
            {"P90/G", "P90/A", "P90/T", "P90/C"},
        ]
        channels: set = set()
        for channel_set in CHANNEL_SETS:
            if channel_set.issubset(lane_data.columns):
                channels = channel_set

        # Initialize occupancy data for this lane
        if lane_int is not None and f"Lane {lane_int}" not in occ_pf:
            occ_pf[f"Lane {lane_int}"] = []
        prev_occ = 0
        prev_pf = 0

        for _, row in lane_data.iterrows():
            # Intensity per cycle
            cycle = int(row["Cycle"])
            for channel in channels:
                intensity = float(row[channel])
                if channel not in intensity_cycle:
                    intensity_cycle[channel] = {}
                if cycle not in intensity_cycle[channel]:
                    intensity_cycle[channel][cycle] = 0
                intensity_cycle[channel][cycle] += intensity

            # Occupancy vs PF
            if {"% Occupied", "% Pass Filter"}.issubset(lane_data.columns):
                occ = float(row["% Occupied"])
                pf = float(row["% Pass Filter"])

                if occ != prev_occ or pf != prev_pf:
                    prev_occ = occ
                    prev_pf = pf
                    if lane_int is not None and isinstance(lane_int, int):
                        # Use modulo to cycle through colors for any number of lanes
                        color_idx = (lane_int - 1) % len(colors)
                        occ_pf[f"Lane {lane_int}"].append({"x": occ, "y": pf, "color": colors[color_idx]})
            else:
                occ_pf = {}

    return {"intensity_cycle": intensity_cycle, "occ_vs_pf": occ_pf}


def _clusters_lane_plot(data: dict) -> Any:
    """Generate clusters/reads per lane bar plot."""
    cluster_data: dict = {}
    read_data: dict = {}

    for value in data.values():
        lane = value.get("Lane")
        if lane is None:
            continue
        lane = int(lane)
        lane_key = f"Lane {lane}"

        cluster_count = value.get("Cluster Count", 0)
        cluster_count_pf = value.get("Cluster Count Pf", 0)
        reads = value.get("Reads", 0)
        reads_pf = value.get("Reads Pf", 0)

        if lane_key not in cluster_data:
            cluster_data[lane_key] = {
                "clusters_pf": cluster_count_pf,
                "clusters_diff": cluster_count - cluster_count_pf,
            }
            read_data[lane_key] = {
                "reads_pf": reads_pf,
                "reads_diff": reads - reads_pf,
            }
        else:
            cluster_data[lane_key]["clusters_pf"] += cluster_count_pf
            cluster_data[lane_key]["clusters_diff"] += cluster_count - cluster_count_pf
            read_data[lane_key]["reads_pf"] += reads_pf
            read_data[lane_key]["reads_diff"] += reads - reads_pf

    cats = [
        {"clusters_pf": {"name": "Clusters PF"}, "clusters_diff": {"name": "Clusters not PF"}},
        {"reads_pf": {"name": "Reads PF"}, "reads_diff": {"name": "Reads not PF"}},
    ]

    return bargraph.plot(
        [cluster_data, read_data],
        cats,
        pconfig={
            "id": "sav-summary-clusters-reads-lane-plot",
            "title": "SAV: Cluster/Reads per Lane",
            "data_labels": ["Clusters", "Reads"],
            "ylab": "Count",
        },
    )


def _intensity_cycle_plot(data: dict) -> Any:
    """Generate intensity per cycle line plot."""
    key_color_dict: dict = {}

    for key in data:
        if re.match(r"\w+/red", key, re.IGNORECASE):
            key_color_dict[key] = "red"
        elif re.match(r"\w+/green", key, re.IGNORECASE):
            key_color_dict[key] = "green"
        elif re.match(r"\w+/G", key):
            key_color_dict[key] = "blue"
        elif re.match(r"\w+/A", key):
            key_color_dict[key] = "black"
        elif re.match(r"\w+/T", key):
            key_color_dict[key] = "green"
        elif re.match(r"\w+/C", key):
            key_color_dict[key] = "red"

    return linegraph.plot(
        data,
        pconfig={
            "id": "sav-intensity-vs-cycle-plot",
            "title": "SAV: Intensity per cycle",
            "xlab": "Cycle",
            "ylab": "Intensity",
            "colors": key_color_dict,
        },
    )


def _occ_vs_pf_plot(data: dict) -> Any:
    """Generate % PF vs % Occupied scatter plot."""
    return scatter.plot(
        data,
        pconfig={
            "id": "sav-pf-vs-occ-plot",
            "title": "SAV: % Passing Filter vs % Occupied",
            "xlab": "% Occupied",
            "ylab": "% Passing Filter",
            "xmin": 0,
            "xmax": 100,
            "ymin": 0,
            "ymax": 100,
        },
    )


def _store_summary_data(module: SAVModule, run_metrics: Any, run_name: str) -> None:
    """Store summary data in module.sav_data for general stats and data export."""
    try:
        summary_total = pd.DataFrame(interop.summary(run_metrics, level="Total"))
        if not summary_total.empty:
            # Get total row data
            total_data = summary_total.iloc[0].to_dict()
            # Store relevant metrics
            module.sav_data[run_name].update(
                {
                    "% >= Q30": total_data.get("% >= Q30"),
                    "Yield G": total_data.get("Yield G"),
                    "% Pf": total_data.get("% Pf"),
                    "Error Rate": total_data.get("Error Rate"),
                    "Cluster Count": total_data.get("Cluster Count"),
                    "Cluster Count Pf": total_data.get("Cluster Count Pf"),
                }
            )
    except (ValueError, TypeError) as e:
        log.debug(f"Could not store summary data for {run_name}: {e}")


def _add_general_stats(module: SAVModule) -> None:
    """Add key metrics to the general statistics table."""
    # Prepare data for general stats (exclude internal keys starting with _)
    gs_data = {}
    for run_name, data in module.sav_data.items():
        gs_data[run_name] = {k: v for k, v in data.items() if not k.startswith("_") and v is not None}

    if not gs_data:
        return

    headers = {
        "% >= Q30": {
            "title": "% >= Q30",
            "description": "Percentage of reads with quality score >= 30",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        },
        "Yield G": {
            "title": f"Yield ({config.base_count_prefix})",
            "description": f"Total yield ({config.base_count_desc})",
            "shared_key": "base_count",
            "scale": "Blues",
            "modify": lambda x: (x * 1e9) * config.base_count_multiplier,
        },
        "% Pf": {
            "title": "% PF",
            "description": "Percentage of clusters passing filter",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
            "hidden": True,
        },
        "Error Rate": {
            "title": "Error %",
            "description": "Error rate as determined by PhiX spike-in",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "OrRd",
            "format": "{:,.2f}",
            "hidden": True,
        },
    }

    headers = module.get_general_stats_headers(all_headers=headers)
    if headers:
        module.general_stats_addcols(gs_data, headers, namespace="sav")
