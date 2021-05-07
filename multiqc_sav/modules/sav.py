#!/usr/bin/env python

import logging
from collections import OrderedDict

import interop
import pandas as pd
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

HEADERS = {
    "Error Rate": {
        "title": "Error Rate (%)",
        "description": "The calculated error rate, as determined by a PhiX spike-in",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
    },
    "Error Rate 35": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-35.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
        "hidden": True,
    },
    "Error Rate 50": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-50.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
        "hidden": True,
    },
    "Error Rate 75": {
        "title": "Error Rate 35 Cycles (%)",
        "description": "The calculated error rate for cycles 1-75.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
        "hidden": True,
    },
    "Error Rate 100": {
        "title": "Error Rate 100 Cycles (%)",
        "description": "The calculated error rate for cycles 1-100.",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
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
        "format": "{:,.0f}",  # No decimal places please
    },
    "% >= Q30": {
        "title": "% >= Q30",
        "description": "Percentage of reads with quality phred score of 30 or above",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
    },
    "% Occupancy Proxy": {
        "title": "Occupancy Proxy (%)",
        # "description": "",
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
    },
    "% Occupied": {
        "title": "Occupied (%)",
        "description": "The percentage of nanowells occupied by clusters, +/- 1 standard deviation.",
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
    },
    "Projected Yield G": {
        "title": "Projected Yield ({})".format(config.base_count_prefix),
        "description": "The expected number of bases sequenced ({} base pairs over all 'usable cycles'".format(
            config.base_count_desc
        ),
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0)
        * config.base_count_multiplier,  # number is already in gigabases
    },
    "Yield G": {
        "title": "Yield ({})".format(config.base_count_prefix),
        "description": "The number of bases sequenced ({} base pairs over all 'usable cycles'".format(
            config.base_count_desc
        ),
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0)
        * config.base_count_multiplier,  # number is already in gigabases
    },
    "Cluster Count": {
        "title": "Cluster Count",
        "description": "Number of clusters for each tile",
        "format": "{:,.0f}",  # No decimal places please
    },
    "Cluster Count Pf": {
        "title": "Cluster Count PF",
        "description": "Number of clusters PF for each tile",
        "format": "{:,.0f}",  # No decimal places please
    },
    "% Pf": {
        "title": "Reads PF (%)",
        "description": "Percentage of clusters Passing Filter",
        "min": 0,
        "max": 100,
        "suffix": "%",
        "format": "{:,.0f}",  # No decimal places please
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
        "title": "{} Reads".format(config.read_count_prefix),
        "description": "The number of reads ({})".format(config.read_count_desc),
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Reads Pf": {
        "title": "{} PF Reads".format(config.read_count_prefix),
        "description": "The number of passing filter reads ({})".format(
            config.read_count_desc
        ),
        "shared_key": "read_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Tile Count": {
        "title": "Tiles",
        "description": "The number of tiles per lane.",
        "hidden": True,
    },
    "Mapped Reads Cv": {},
    "Max Mapped Reads": {},
    "Min Mapped Reads": {},
    "Total Fraction Mapped Reads": {},
    "Total Pf Reads": {},
    "Total Reads": {},
    "Cluster Count": {},
    "Fraction Mapped": {},
    "Id": {},
    "Index1": {},
    "Index2": {},
    "Project Name": {},
    "Sample Id": {},
}


class SAV(BaseMultiqcModule):
    """
    Generate SAV tables an Graphs including:
    - GRAPH: Intensity/Cycle/Channel
    - GRAPH: Clusters/Lane
    - GRAPH: Qscore Heatmap
    - GRAPH: Qscore Histogram
    - GRAPH: %Occ/%PF
    - TABLE: Run Summary
    - TABLE: Indexing Summary
    """

    def __init__(self):

        # Check if the plugin has an input directory
        if not config.kwargs.get("illumina_dir", None):
            log.info("Skipping MultiQC_SAV")
            return None
        else:
            log.debug("Running MultiQC SAV plugin")
            self.illumina_dir = config.kwargs.get("illumina_dir")

        super(SAV, self).__init__(
            name="Illumina SAV",
            anchor="sav",
            info=" - Sequencing Metrics from Illumina sequencers",
        )

        self.load_metrics()
        self.summary_qc()
        self.indexing_qc()
        # self.imaging_qc()

    def load_metrics(self):
        log.info("Loading run metrics from {}".format(self.illumina_dir))
        self.run_metrics = interop.read(run=self.illumina_dir, finalize=True)

    def summary_qc(self):
        log.info("Gathering Read summary metrics")
        summary_read = pd.DataFrame(interop.summary(self.run_metrics, level="Read"))
        summary_nonindex = pd.DataFrame(
            interop.summary(self.run_metrics, level="NonIndex")
        )
        summary_total = pd.DataFrame(interop.summary(self.run_metrics, level="Total"))

        self.add_section(
            name="Summary Read Metrics",
            anchor="sav-read-summary",
            description="Summary metrics per Read",
            plot=self.read_summary_table(
                self.parse_read_summary(summary_read, summary_nonindex, summary_total)
            ),
        )

        log.info("Gathering Lane summary metrics")
        summary_lane = pd.DataFrame(interop.summary(self.run_metrics, level="Lane"))

        self.add_section(
            name="Summary Lane Metrics",
            anchor="sav-lane-summary",
            description="Summary metrics per Lane per Read",
            plot=self.lane_summary_table(self.parse_lane_summary(summary_lane)),
        )

    def parse_read_summary(self, read_metrics, non_index_metrics, total_metrics):
        table_data: dict = self._parse_reads(read_metrics)

        for read, data in non_index_metrics.iterrows():
            table_data["Non-Indexed"] = data.to_dict()

        for read, data in total_metrics.iterrows():
            table_data["Total"] = data.to_dict()

        return table_data

    def read_summary_table(self, data):
        headers = {
            header: HEADERS[header] for header in interop.summary_columns(level="Lane")
        }

        table_config = {
            "namespace": "SAV",
            "id": "sav-read-metrics-summary-table",
            "col1_header": "Read",
        }

        return table.plot(data, headers, table_config)

    def parse_lane_summary(self, data):
        lanes = data.groupby("Lane")
        table_data: dict = {}
        for lane, reads in lanes:
            lane_data = {}
            reads_dict = self._parse_reads(reads, key_prefix=f"Lane {lane}")
            table_data.update(reads_dict)

        return table_data

    def lane_summary_table(self, data):
        headers = {
            header: HEADERS[header] for header in interop.summary_columns(level="Lane")
        }
        table_config = {
            "namespace": "SAV",
            "id": "sav-lane-metrics-summary-table",
            "col1_header": "Lane - Read",
        }

        return table.plot(data, headers, table_config,)

    def indexing_qc(self):
        log.info("Gathering Lane Indexing metrics")
        index_summary_lane = pd.DataFrame(
            interop.index_summary(self.run_metrics, level="Lane")
        )
        self.add_section(
            name="Indexing Lane Metrics",
            anchor="sav-lane-index",
            description="Indexing metrics per Lane",
            plot=self.lane_index_summary_table(
                self.parse_lane_index_summary(index_summary_lane)
            ),
        )

        log.info("Gathering Barcode Indexing metrics")
        indes_summar_barcode = pd.DataFrame(
            interop.indexing(self.run_metrics, level="Barcode")
        )

        # self.add_section(
        #     name="Indexing Barcode Metrics",
        #     anchor="sav-barcode-index",
        #     description="Indexing metrics per Barcode",
        #     plot=lane_index_summary_table(
        #         parse_lane_index_summary(indexing_lane_summary)
        #     ),
        # )

    #def parse_lane_index_summary(self, data):
        data = data.set_index("Lane")
        table_data = {}
        for lane, lane_data in data.iterrows():
            table_data[f"Lane {lane}"] = lane_data.to_dict()
        return table_data


    #def lane_index_summary_table(self, data):
        headers = {
            header: HEADERS[header]
            for header in interop.index_summary_columns(level="Lane")
        }
        table_config = {
            "namespace": "SAV",
            "id": "sav-lane-index-metrics-summary-table",
            "col1_header": "Lane",
        }

        return table.plot(data, headers, table_config,)


    def parse_barcode_index_summary(self, data):

    def barcode_index_summary_table(self, data):

    def imaging_qc(self):
        log.debug("SAV: Gathering Imaging metrics")
        imaging = pd.DataFrame(interop.imaging(self.run_metrics))

        #####
        # GRAPH: %Occ/%PF
        #####
        self.add_section(
            name="Imaging Metrics", anchor="sav-imaging", description="",
        )

    def _parse_reads(self, reads_df, key_prefix: str = None):
        reads_dict = {}
        reads_df = reads_df.set_index("ReadNumber")
        for read, data in reads_df.iterrows():
            key = f"Read {read}" + " (I)" if data["IsIndex"] == 89 else f"Read {read}"
            if key_prefix:
                key = f"{key_prefix} - {key}"
            reads_dict[key] = data.drop("IsIndex").to_dict()
        return reads_dict
