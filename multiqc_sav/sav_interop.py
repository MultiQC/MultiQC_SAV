"""
InterOp-based visualizations for the SAV module.

This module provides functions to parse Illumina InterOp binary files
and generate advanced visualizations.
"""

import glob
import logging
import os
import re
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from multiqc import config
from multiqc.plots import bargraph, heatmap, linegraph, scatter, table
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

# Table headers for summary metrics
HEADERS: Dict[str, Dict] = {
    "Error Rate": {
        "title": "Error Rate (%)",
        "description": "The calculated error rate, as determined by a PhiX spike-in",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
    },
    "Error Rate 35": {
        "title": "Error Rate 35 (%)",
        "description": "The calculated error rate for cycles 1-35.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "hidden": True,
    },
    "Error Rate 50": {
        "title": "Error Rate 50 (%)",
        "description": "The calculated error rate for cycles 1-50.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "hidden": True,
    },
    "Error Rate 75": {
        "title": "Error Rate 75 (%)",
        "description": "The calculated error rate for cycles 1-75.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "hidden": True,
    },
    "Error Rate 100": {
        "title": "Error Rate 100 (%)",
        "description": "The calculated error rate for cycles 1-100.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
        "hidden": True,
    },
    "First Cycle Intensity": {
        "title": "Intensity Cycle 1",
        "description": "The average of the A channel intensity measured at the first cycle",
    },
    "% Aligned": {
        "title": "Aligned (%)",
        "description": "Percentage of reads that aligned to the PhiX genome",
        "suffix": "%",
        "min": 0,
        "max": 100,
        "format": "{:,.2f}",
    },
    "% >= Q30": {
        "title": "% >= Q30",
        "description": "Percentage of reads with quality phred score of 30 or above",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
    },
    "% Occupancy Proxy": {
        "title": "Occupancy Proxy (%)",
        "suffix": "%",
        "format": "{:,.2f}",
    },
    "% Occupied": {
        "title": "Occupied (%)",
        "description": "The percentage of nanowells occupied by clusters, +/- 1 standard deviation.",
        "suffix": "%",
        "format": "{:,.2f}",
    },
    "Projected Yield G": {
        "title": f"Projected Yield ({config.base_count_prefix})",
        "description": f"The expected number of bases sequenced ({config.base_count_desc} base pairs over all 'usable cycles'",
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0) * config.base_count_multiplier,
        "hidden": True,
    },
    "Yield G": {
        "title": f"Yield ({config.base_count_prefix})",
        "description": f"The number of bases sequenced ({config.base_count_desc} base pairs over all 'usable cycles'",
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0) * config.base_count_multiplier,
    },
    "Cluster Count": {
        "title": f"Clusters ({config.read_count_prefix})",
        "description": f"Number of clusters for each tile ({config.read_count_desc})",
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Cluster Count Pf": {
        "title": f"Clusters PF ({config.read_count_prefix})",
        "description": f"Number of clusters PF for each tile ({config.read_count_desc})",
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "% Pf": {
        "title": "Reads PF (%)",
        "description": "Percentage of clusters Passing Filter",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.2f}",
    },
    "Density": {
        "title": "Density",
        "description": "The density of clusters (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.",
        "hidden": True,
    },
    "Density Pf": {
        "title": "Density PF",
        "description": "The density of clusters PF (in thousands per mm2) detected by image analysis, +/- 1 standard deviation.",
        "hidden": True,
    },
    "Phasing": {
        "title": "Phasing",
        "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing falls behind (phasing) the current cycle within a read.",
    },
    "Phasing Offset": {
        "title": "Phasing Offset",
        "description": "The best-fit offset of the phasing corrections, calculated from the entire read.",
        "hidden": True,
    },
    "Phasing Slope": {
        "title": "Phasing Slope",
        "description": "The best-fit slope of the phasing corrections, calculated from the entire read.",
        "hidden": True,
    },
    "Prephasing": {
        "title": "Prephasing",
        "description": "The value used by RTA for the percentage of molecules in a cluster for which sequencing jumps ahead (prephasing) the current cycle within a read.",
    },
    "Prephasing Offset": {
        "title": "Prephasing Offset",
        "description": "The best-fit offset of the prephasing corrections, calculated from the entire read.",
        "hidden": True,
    },
    "Prephasing Slope": {
        "title": "Prephasing Slope",
        "description": "The best-fit slope of the prephasing corrections, calculated from the entire read.",
        "hidden": True,
    },
    "Reads": {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"The number of reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Reads Pf": {
        "title": f"{config.read_count_prefix} PF Reads",
        "description": f"The number of passing filter reads ({config.read_count_desc})",
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Tile Count": {
        "title": "Tiles",
        "description": "The number of tiles per lane.",
        "hidden": True,
    },
    "Total Pf Reads": {
        "title": f"{config.read_count_prefix} PF Reads",
        "description": f"The total number of passing filter reads for this lane ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Total Reads": {
        "title": f"{config.read_count_prefix} Reads",
        "description": f"The total number of reads for this lane ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Mapped Reads Cv": {
        "title": "CV",
        "description": "The coefficient of variation for the number of counts across all indexes.",
        "format": "{:,.2f}",
    },
    "Max Mapped Reads": {
        "title": f"{config.read_count_prefix} Max Mapped Reads",
        "description": f"The highest representation for any index ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Min Mapped Reads": {
        "title": f"{config.read_count_prefix} Min Mapped Reads",
        "description": f"The lowest representation for any index ({config.read_count_desc})",
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
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


def add_interop_sections(module) -> None:
    """
    Add InterOp-based sections to the SAV module.

    Args:
        module: The SAV MultiqcModule instance
    """
    try:
        import interop
        from interop import py_interop_plot
    except ImportError:
        log.warning("InterOp library not installed. Run: pip install interop")
        return

    # Find the run directory from the module's data
    for s_name, data in module.data_by_sample.items():
        run_dir = data.get("_run_dir")
        if not run_dir:
            continue

        # Check for InterOp directory
        interop_files = glob.glob(os.path.join(run_dir, "InterOp", "*.bin"))
        if not interop_files:
            log.debug(f"No InterOp files found in {run_dir}")
            continue

        log.info(f"Loading InterOp metrics from {run_dir}")

        try:
            run_metrics = interop.read(
                run=run_dir,
                valid_to_load=interop.load_imaging_metrics(),
                finalize=True,
            )
        except Exception as e:
            log.warning(f"Failed to load InterOp metrics from {run_dir}: {e}")
            continue

        # Add summary sections
        _add_summary_sections(module, run_metrics, interop)

        # Add Q-score sections
        _add_qscore_sections(module, run_metrics, py_interop_plot)

        # Add imaging sections
        _add_imaging_sections(module, run_metrics, interop)

        # Only process first run for now
        break


def _add_summary_sections(module, run_metrics, interop) -> None:
    """Add read and lane summary table sections."""
    log.info("Gathering summary metrics")

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
                    headers,
                    pconfig={
                        "title": "SAV: Read Metrics Summary",
                        "namespace": "SAV",
                        "id": "sav-read-metrics-summary-table",
                        "col1_header": "Read",
                    },
                ),
            )
    except Exception as e:
        log.debug(f"Could not generate read summary: {e}")

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
                    headers,
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
    except Exception as e:
        log.debug(f"Could not generate lane summary: {e}")


def _add_qscore_sections(module, run_metrics, py_interop_plot) -> None:
    """Add Q-score heatmap and histogram sections."""
    log.info("Generating Q-score plots")

    # Q-score heatmap
    try:
        options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
        rows = py_interop_plot.count_rows_for_heatmap(run_metrics)
        cols = py_interop_plot.count_columns_for_heatmap(run_metrics)

        if rows > 0 and cols > 0:
            data_buffer = np.zeros((rows, cols), dtype=np.float32)
            data = py_interop_plot.heatmap_data()

            try:
                py_interop_plot.plot_qscore_heatmap(run_metrics, options, data, data_buffer.ravel())
            except py_interop_plot.invalid_filter_option:
                pass

            # data_buffer shape is (rows=qscores, cols=cycles)
            # heatmap expects: rows match ycats, cols match xcats
            plot_data = data_buffer.tolist()
            x_cats = list(range(cols))   # cycles (x-axis)
            y_cats = list(range(rows))   # qscores (y-axis)

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
    except Exception as e:
        log.debug(f"Could not generate Q-score heatmap: {e}")

    # Q-score histogram
    try:
        bar_data = py_interop_plot.bar_plot_data()
        options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
        py_interop_plot.plot_qscore_histogram(run_metrics, options, bar_data)

        hist: Dict = {}
        qscore: List = []
        reads: List = []
        binsize: List = []

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
    except Exception as e:
        log.debug(f"Could not generate Q-score histogram: {e}")


def _add_imaging_sections(module, run_metrics, interop) -> None:
    """Add imaging-related sections (intensity per cycle, % PF vs % occupied)."""
    log.info("Gathering imaging metrics")

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
    except Exception as e:
        log.debug(f"Could not generate imaging plots: {e}")


def _parse_read_summary(
    read_metrics: pd.DataFrame,
    non_index_metrics: pd.DataFrame,
    total_metrics: pd.DataFrame,
) -> Dict:
    """Parse read summary DataFrames into dict format."""
    table_data: Dict = _parse_reads(read_metrics)

    for _, data in non_index_metrics.iterrows():
        table_data["Non-Indexed"] = data.to_dict()

    for _, data in total_metrics.iterrows():
        table_data["Total"] = data.to_dict()

    return table_data


def _parse_lane_summary(data: pd.DataFrame) -> Dict:
    """Parse lane summary DataFrame into dict format."""
    lanes = data.groupby("Lane")
    table_data: Dict = {}

    for lane, reads in lanes:
        reads_dict = _parse_reads(reads, key_prefix=f"Lane {lane}")
        table_data.update(reads_dict)

    return table_data


def _parse_reads(reads_df: pd.DataFrame, key_prefix: Optional[str] = None) -> Dict:
    """Utility function to parse a reads DataFrame to dict."""
    reads_dict: Dict = {}
    reads_df = reads_df.set_index("ReadNumber")

    for read, data in reads_df.iterrows():
        key = f"Read {read}" + " (I)" if data["IsIndex"] == 89 else f"Read {read}"
        if key_prefix:
            key = f"{key_prefix} - {key}"
        reads_dict[key] = data.drop("IsIndex").to_dict()

    return reads_dict


def _parse_imaging_table(data: pd.DataFrame) -> Dict:
    """Parse imaging table DataFrame for intensity and occupancy plots."""
    cscale = mqc_colour.mqc_colour_scale()
    colors = cscale.get_colours("Dark2")

    per_lane = data.groupby("Lane")
    occ_pf: Dict = {}
    intensity_cycle: Dict = {}

    for lane, lane_data in per_lane:
        lane = int(lane)

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
        if f"Lane {lane}" not in occ_pf:
            occ_pf[f"Lane {lane}"] = []
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
                    occ_pf[f"Lane {lane}"].append({"x": occ, "y": pf, "color": colors[lane]})
            else:
                occ_pf = {}

    return {"intensity_cycle": intensity_cycle, "occ_vs_pf": occ_pf}


def _clusters_lane_plot(data: Dict):
    """Generate clusters/reads per lane bar plot."""
    cluster_data: Dict = {}
    read_data: Dict = {}

    for key, value in data.items():
        lane = int(value["Lane"])
        lane_key = f"Lane {lane}"

        if lane_key not in cluster_data:
            cluster_data[lane_key] = {
                "clusters_pf": value["Cluster Count Pf"],
                "clusters_diff": value["Cluster Count"] - value["Cluster Count Pf"],
            }
            read_data[lane_key] = {
                "reads_pf": value["Reads Pf"],
                "reads_diff": value["Reads"] - value["Reads Pf"],
            }
        else:
            cluster_data[lane_key]["clusters_pf"] += value["Cluster Count Pf"]
            cluster_data[lane_key]["clusters_diff"] += value["Cluster Count"] - value["Cluster Count Pf"]
            read_data[lane_key]["reads_pf"] += value["Reads Pf"]
            read_data[lane_key]["reads_diff"] += value["Reads"] - value["Reads Pf"]

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


def _intensity_cycle_plot(data: Dict):
    """Generate intensity per cycle line plot."""
    key_color_dict: Dict = {}

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


def _occ_vs_pf_plot(data: Dict):
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
