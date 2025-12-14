"""MultiQC SAV module for parsing Illumina InterOp data."""

from __future__ import annotations

import contextlib
import glob
import logging
import os
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Any

import interop
import numpy as np
import pandas as pd
from interop import py_interop_plot
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, heatmap, linegraph, scatter, table
from multiqc.utils import mqc_colour

log = logging.getLogger("multiqc")

HEADERS: dict[str, dict[str, Any]] = {
    "Error Rate": {
        "title": "Error Rate (%)",
        "description": "The calculated error rate, as determined by a PhiX spike-in",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
    },
    "Error Rate 35": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-35.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
        "hidden": True,
    },
    "Error Rate 50": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-50.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
        "hidden": True,
    },
    "Error Rate 75": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-75.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
        "hidden": True,
    },
    "Error Rate 100": {
        "title": "Error Rate 100 Cycles (%)",
        "description": "The calculated error rate for cycles 1-100.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
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
        "format": "{:,.0f}",
    },
    "% >= Q30": {
        "title": "% >= Q30",
        "description": "Percentage of reads with quality phred score of 30 or above",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",
    },
    "% Occupancy Proxy": {
        "title": "Occupancy Proxy (%)",
        "suffix": "%",
        "format": "{:,.0f}",
    },
    "% Occupied": {
        "title": "Occupied (%)",
        "description": "The percentage of nanowells occupied by clusters, +/- 1 standard deviation.",
        "suffix": "%",
        "format": "{:,.0f}",
    },
    "Projected Yield G": {
        "title": f"Projected Yield ({config.base_count_prefix})",
        "description": f"The expected number of bases sequenced ({config.base_count_desc} base pairs over all 'usable cycles'",
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0) * config.base_count_multiplier,
        "hidden": True,
    },
    "Yield G": {
        "title": f"Yield ({config.read_count_prefix})",
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
        "format": "{:,.0f}",
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


class SAV(BaseMultiqcModule):
    """
    Generate SAV tables and graphs including:
    - GRAPH: Intensity/Cycle/Channel
    - GRAPH: Clusters/Lane
    - GRAPH: Qscore Heatmap
    - GRAPH: Qscore Histogram
    - GRAPH: %Occ/%PF
    - TABLE: Run Summary
    """

    def __init__(self) -> None:
        super().__init__(
            name="Illumina SAV",
            anchor="sav",
            info=" - Sequencing Metrics from Illumina sequencers",
        )

        # Set variables
        run_info_xml = ""
        run_parameters_xml = ""
        illumina_dir = ""

        # Check if required files are found
        for f in self.find_log_files("SAV/xml"):
            if re.match(r".*[Rr]un[Ii]nfo\.xml", f["fn"]):
                run_info_xml = os.path.join(f["root"], f["fn"])
            if re.match(r".*[Rr]un[Pp]arameters\.xml", f["fn"]):
                run_parameters_xml = os.path.join(f["root"], f["fn"])

        # Assume single run for now
        if (os.path.dirname(run_info_xml) == os.path.dirname(run_parameters_xml)) and len(
            glob.glob(os.path.join(os.path.dirname(run_info_xml), "InterOp/*.bin"))
        ) > 0:
            illumina_dir = os.path.dirname(run_info_xml)
        else:
            log.debug("Skipping MultiQC_SAV, required files were not found or not in the right structure.")
            return

        self.run_metrics: Any = None
        self.set_run_info(run_info_xml)
        self.load_metrics(illumina_dir)
        self.summary_qc()
        self.q_summary()
        self.imaging_qc()

    def load_metrics(self, illumina_dir: str) -> None:
        """Load run metrics from InterOp directory."""
        log.info("Loading Run Metrics")
        self.run_metrics = interop.read(
            run=illumina_dir,
            valid_to_load=interop.load_imaging_metrics(),
            finalize=True,
        )

    #############
    # RUN INFO
    #############

    def set_run_info(self, run_info_xml: str) -> None:
        """Parse and display run information from RunInfo.xml."""
        log.info("Loading Run Info")
        tree = ET.parse(run_info_xml)
        root = tree.getroot()

        for run in root:
            run_number = run.attrib["Number"]
            flowcell = next(fc.text for fc in run.iter("Flowcell"))
            instrument_id = next(fc.text for fc in run.iter("Instrument"))
            run_date = next(fc.text for fc in run.iter("Date"))

            parsed_run_date = self._parse_run_date(run_date)

            read_info = ""
            for read in run.iter("Read"):
                key = (
                    f"Read {read.attrib['Number']} (I)"
                    if read.attrib["IsIndexedRead"] == "Y"
                    else f"Read {read.attrib['Number']}"
                )
                read_info += f"<li><b>{key}</b>: {read.attrib['NumCycles']} Cycles</li>"

        self.add_section(
            name="Run Info",
            anchor="sav-run-info",
            content=f"""
                <div class="container-fluid">
                  <div class="row">
                    <div class="col-sm-4">
                      <h4>Instrument</h4>
                      <ul>
                        <li><b>Instrument ID:</b> {instrument_id}</li>
                        <li><b>Flowcell:</b> {flowcell}</li>
                        <li><b>Run Number:</b> {run_number}</li>
                        <li><b>Run Date:</b> {parsed_run_date}</li>
                      </ul>
                    </div>
                    <div class="col-sm-4">
                      <h4>Settings</h4>
                      <ul>
                        {read_info}
                      </ul>
                    </div>
                  </div>
                </div>
            """,
        )

    def _parse_run_date(self, run_date: str) -> str:
        """
        Parse run date from various Illumina sequencer formats.

        Supported formats:
        - MiSeq/NextSeq500/HiSeq3000: YYMMDD
        - NovaSeq6000: M/D/YYYY H:MM:SS AM/PM
        - NextSeq2000: YYYY-MM-DDTHH:MM:SSZ
        """
        date_formats = [
            ("%y%m%d", "%d-%m-%Y"),  # MiSeq/NextSeq500/HiSeq3000
            ("%m/%d/%Y %I:%M:%S %p", "%d-%m-%Y"),  # NovaSeq6000
            ("%Y-%m-%dT%H:%M:%SZ", "%d-%m-%Y"),  # NextSeq2000
        ]

        for input_fmt, output_fmt in date_formats:
            try:
                return datetime.strftime(datetime.strptime(run_date, input_fmt), output_fmt)
            except ValueError:
                continue

        # Return raw date if no format matched
        return run_date

    #############
    # SUMMARY QC
    #############

    def summary_qc(self) -> None:
        """Generate MultiQC sections related to Summary tables."""
        log.info("Gathering Read summary metrics")
        summary_read = pd.DataFrame(interop.summary(self.run_metrics, level="Read"))
        summary_nonindex = pd.DataFrame(interop.summary(self.run_metrics, level="NonIndex"))
        summary_total = pd.DataFrame(interop.summary(self.run_metrics, level="Total"))

        self.add_section(
            name="Summary Read Metrics",
            anchor="sav-read-summary",
            description="Summary metrics per Read",
            plot=self.read_summary_table(self.parse_read_summary(summary_read, summary_nonindex, summary_total)),
        )

        log.info("Gathering Lane summary metrics")
        summary_lane = pd.DataFrame(interop.summary(self.run_metrics, level="Lane"))

        self.add_section(
            name="Summary Lane Metrics",
            anchor="sav-lane-summary",
            description="Summary metrics per Lane per Read",
            plot=self.lane_summary_table(self.parse_lane_summary(summary_lane)),
        )

        # - GRAPH: Clusters/Lane
        log.info("Generating 'Clusters/Lane' plot")
        self.add_section(
            name="Clusters/Reads per Lane",
            anchor="sav-clusters-lane",
            description="Total Cluster/Read count per Lane",
            plot=self.clusters_lane_plot(self.parse_lane_summary(summary_lane)),
        )

    def parse_read_summary(
        self,
        read_metrics: pd.DataFrame,
        non_index_metrics: pd.DataFrame,
        total_metrics: pd.DataFrame,
    ) -> dict[str, dict[str, Any]]:
        """Parse "Read Summary" table DataFrame."""
        table_data: dict[str, dict[str, Any]] = self._parse_reads(read_metrics)

        for _read, data in non_index_metrics.iterrows():
            table_data["Non-Indexed"] = data.to_dict()

        for _read, data in total_metrics.iterrows():
            table_data["Total"] = data.to_dict()

        return table_data

    def read_summary_table(self, data: dict[str, dict[str, Any]]) -> table.plot:
        """Format "Read Summary" data dict and add plot config."""
        headers = {header: HEADERS[header] for header in interop.summary_columns(level="Lane")}

        table_config = {
            "namespace": "SAV",
            "id": "sav-read-metrics-summary-table",
            "col1_header": "Read",
        }

        return table.plot(data, headers, table_config)

    def parse_lane_summary(self, data: pd.DataFrame) -> dict[str, dict[str, Any]]:
        """Parse "Lane Summary" table DataFrame."""
        lanes = data.groupby("Lane")
        table_data: dict[str, dict[str, Any]] = {}
        for lane, reads in lanes:
            reads_dict = self._parse_reads(reads, key_prefix=f"Lane {lane}")
            table_data.update(reads_dict)

        return table_data

    def lane_summary_table(self, data: dict[str, dict[str, Any]]) -> table.plot:
        """Format "Lane Summary" data dict and add plot config."""
        headers = {header: HEADERS[header] for header in interop.summary_columns(level="Lane")}
        table_config = {
            "namespace": "SAV",
            "id": "sav-lane-metrics-summary-table",
            "col1_header": "Lane - Read",
        }

        return table.plot(data, headers, table_config)

    def clusters_lane_plot(self, data: dict[str, dict[str, Any]]) -> bargraph.plot:
        """Format "Clusters per Lane" data dict and add plot config."""
        cluster_data: dict[str, dict[str, float]] = {}
        read_data: dict[str, dict[str, float]] = {}

        for value in data.values():
            lane = int(value["Lane"])
            lane_key = f"Lane {lane}"

            if lane_key not in cluster_data:
                cluster_data[lane_key] = {
                    "clusters": value["Cluster Count"],
                    "clusters_pf": value["Cluster Count Pf"],
                    "clusters_diff": value["Cluster Count"] - value["Cluster Count Pf"],
                }
                read_data[lane_key] = {
                    "reads": value["Reads"],
                    "reads_pf": value["Reads Pf"],
                    "reads_diff": value["Reads"] - value["Reads Pf"],
                }
            else:
                cluster_data[lane_key]["clusters"] += value["Cluster Count"]
                cluster_data[lane_key]["clusters_pf"] += value["Cluster Count Pf"]
                cluster_data[lane_key]["clusters_diff"] += value["Cluster Count"] - value["Cluster Count Pf"]
                read_data[lane_key]["reads"] += value["Reads"]
                read_data[lane_key]["reads_pf"] += value["Reads Pf"]
                read_data[lane_key]["reads_diff"] += value["Reads"] - value["Reads Pf"]

        cats = [
            {"clusters_pf": {"name": "Clusters PF"}, "clusters_diff": {"name": "Clusters not PF"}},
            {"reads_pf": {"name": "Reads PF"}, "reads_diff": {"name": "Reads not PF"}},
        ]

        plot_config = {
            "id": "sav-summary-clusters-reads-lane-plot",
            "title": "SAV: Cluster/Reads per Lane",
            "data_labels": ["Clusters", "Reads"],
            "ylab": "Lane",
        }

        return bargraph.plot([cluster_data, read_data], cats, plot_config)

    def _parse_reads(
        self,
        reads_df: pd.DataFrame,
        key_prefix: str | None = None,
    ) -> dict[str, dict[str, Any]]:
        """Parse a "Reads" dataframe to dict."""
        reads_dict: dict[str, dict[str, Any]] = {}
        reads_df = reads_df.set_index("ReadNumber")

        for read, data in reads_df.iterrows():
            key = f"Read {read}" + " (I)" if data["IsIndex"] == 89 else f"Read {read}"
            if key_prefix:
                key = f"{key_prefix} - {key}"
            reads_dict[key] = data.drop("IsIndex").to_dict()

        return reads_dict

    #############
    # Q SUMMARY
    #############

    def q_summary(self) -> None:
        """Generate MultiQC sections related to Qscore."""
        # - GRAPH: Qscore Heatmap
        log.info("Generating 'Qscore Heatmap' plot")
        self.add_section(
            name="Qscore Heatmap",
            anchor="sav-qscore-heatmap",
            description="The Qscore Heat Map provides an overview of quality scores across cycles.",
            plot=self.qscore_heatmap_plot(),
        )

        # - GRAPH: Qscore Histogram
        log.info("Generating 'Qscore Histogram' plot")
        self.add_section(
            name="Qscore Histogram",
            anchor="sav-qscore-histogram",
            description=(
                "Qscore Histogram graphs the number of bases by quality score. "
                "The quality score is cumulative for the current cycle. "
                "Only bases from reads that pass the quality filter are included."
            ),
            plot=self.qscore_histogram_plot(),
        )

    def qscore_heatmap_plot(self) -> heatmap.plot:
        """
        Get heatmap data from run_metrics object.

        Note: this function has room for improvement, but we need to wait
        for further developments in the InterOp library.
        """
        options = py_interop_plot.filter_options(self.run_metrics.run_info().flowcell().naming_method())
        rows = py_interop_plot.count_rows_for_heatmap(self.run_metrics)
        cols = py_interop_plot.count_columns_for_heatmap(self.run_metrics)
        data_buffer = np.zeros((rows, cols), dtype=np.float32)
        data = py_interop_plot.heatmap_data()

        with contextlib.suppress(py_interop_plot.invalid_filter_option):
            py_interop_plot.plot_qscore_heatmap(self.run_metrics, options, data, data_buffer.ravel())

        plot_data = data_buffer.transpose().tolist()
        # cycles
        x_cats = list(range(0, cols))
        # qscore
        y_cats = list(range(0, rows))

        plot_config = {
            "id": "sav-qscore-heatmap-plot",
            "title": "SAV: Qscore Heatmap",
            "xTitle": "Cycle",
            "yTitle": "Qscore",
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
        }
        return heatmap.plot(plot_data, x_cats, y_cats, plot_config)

    def qscore_histogram_plot(self) -> linegraph.plot:
        """
        Get histogram data from run_metrics object.

        Note: this function has room for improvement, but we need to wait
        for further developments in the InterOp library.
        """
        bar_data = py_interop_plot.bar_plot_data()
        options = py_interop_plot.filter_options(self.run_metrics.run_info().flowcell().naming_method())
        py_interop_plot.plot_qscore_histogram(self.run_metrics, options, bar_data)

        hist: dict[float, float] = {}
        qscore: list[float] = []
        reads: list[float] = []
        binsize: list[float] = []

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

        plot_data = {bar_data.title(): hist}

        plot_config = {
            "id": "sav-qscore-histogram-plot",
            "title": "SAV: Qscore Histogram",
            "xlab": "Qscore",
            "ylab": "Reads (Billion)",
        }
        return linegraph.plot(plot_data, plot_config)

    #############
    # IMAGING QC
    #############

    def imaging_qc(self) -> None:
        """
        Generate MultiQC sections related to Imaging.

        This includes:
            - Plot: Intensity/Cycle/Channel
            - Plot: %Occ/%PF
        """
        log.info("Gathering Imaging metrics")
        imaging = pd.DataFrame(interop.imaging(self.run_metrics))

        plot_data = self.parse_imaging_table(imaging)

        # - GRAPH: Intensity/Cycle/Channel
        if len(plot_data.get("intensity_cycle", [])) > 0:
            log.info("Generating 'Intensity per Cycle' plot")
            self.add_section(
                name="Intensity per Cycle",
                anchor="sav-intensity-cycle",
                description="Intensity by color and cycle of the 90% percentile of the data for each tile",
                plot=self.intensity_cycle_plot(plot_data.get("intensity_cycle", {})),
            )

        # - GRAPH: %Occ/%PF
        log.info("Generating '% PF vs % Occupied' plot")
        if len(plot_data.get("occ_vs_pf", [])) > 0:
            self.add_section(
                name="% PF vs % Occupied",
                anchor="sav-imaging-pf-vs-occ",
                description="% Clusters passing filter vs % Wells Occupied",
                plot=self.occ_vs_pf_plot(plot_data.get("occ_vs_pf", {})),
            )

    def parse_imaging_table(self, data: pd.DataFrame) -> dict[str, Any]:
        """
        Parse full imaging table DataFrame.

        Returns dict containing data for intensity per cycle plot (key: "intensity_cycle")
        and %occ vs %pf plot (key: "occ_vs_pf").
        """
        # set color scale for occ_pf
        cscale = mqc_colour.mqc_colour_scale()
        colors = cscale.get_colours("Dark2")

        per_lane = data.groupby("Lane")
        occ_pf: dict[str, list[dict[str, Any]]] = {}
        intensity_cycle: dict[str, dict[int, float]] = {}

        for lane_num, lane_data in per_lane:
            lane = int(lane_num)

            # prep intensity_cycle
            channel_sets = [
                {"P90/RED", "P90/GREEN"},
                {"P90/Red", "P90/Green"},
                {"P90/G", "P90/A", "P90/T", "P90/C"},
            ]
            channels: set[str] = set()
            for channel_set in channel_sets:
                if channel_set.issubset(lane_data.columns):
                    channels = channel_set

            # prep occ_pf
            lane_key = f"Lane {lane}"
            if lane_key not in occ_pf:
                occ_pf[lane_key] = []
                prev_occ = 0.0
                prev_pf = 0.0

            # parse imaging table lane
            for _, row in lane_data.iterrows():
                # intensity_cycle
                cycle = int(row["Cycle"])
                for channel in channels:
                    intensity = float(row[channel])
                    if channel not in intensity_cycle:
                        intensity_cycle[channel] = {}
                    if cycle not in intensity_cycle[channel]:
                        intensity_cycle[channel][cycle] = 0
                    intensity_cycle[channel][cycle] += intensity

                # occ_pf
                if {"% Occupied", "% Pass Filter"}.issubset(lane_data.columns):
                    occ = float(row["% Occupied"])
                    pf = float(row["% Pass Filter"])

                    if occ != prev_occ or pf != prev_pf:
                        prev_occ = occ
                        prev_pf = pf
                        occ_pf[lane_key].append({"x": occ, "y": pf, "color": colors[lane]})
                else:
                    occ_pf = {}

        return {"intensity_cycle": intensity_cycle, "occ_vs_pf": occ_pf}

    def intensity_cycle_plot(self, data: dict[str, dict[int, float]]) -> linegraph.plot:
        """Format Intensity per Cycle data dict and add plot config."""
        # get keys from data
        key_color_dict: dict[str, str] = {}
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

        plot_config = {
            "id": "sav-intensity-vs-cycle-plot",
            "title": "SAV: Intensity per cycle",
            "xlab": "Cycle",
            "ylab": "Intensity",
            "colors": key_color_dict,
        }

        return linegraph.plot(data, plot_config)

    def occ_vs_pf_plot(self, data: dict[str, list[dict[str, Any]]]) -> scatter.plot:
        """Format %Occ vs %PF data dict and add plot config."""
        plot_config = {
            "id": "sav-pf-vs-occ-plot",
            "title": "SAV: % Passing Filter vs % Occupied",
            "xlab": "% Occupied",
            "ylab": "% Passing Filter",
            "xmin": 0,
            "xmax": 100,
            "ymin": 0,
            "ymax": 100,
        }
        return scatter.plot(data, plot_config)
