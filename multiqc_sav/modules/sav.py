#!/usr/bin/env python

import glob
import logging
import os
from collections import OrderedDict

import interop
import pandas as pd
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, beeswarm, heatmap, linegraph, scatter, table
from multiqc.utils import mqc_colour

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
        "hidden": True,
    },
    "Yield G": {
        "title": "Yield ({})".format(config.read_count_prefix),
        "description": "The number of bases sequenced ({} base pairs over all 'usable cycles'".format(
            config.base_count_desc
        ),
        "shared_key": "base_count",
        "modify": lambda x: (x * 1000000000.0)
        * config.base_count_multiplier,  # number is already in gigabases
    },
    "Cluster Count": {
        "title": "Clusters ({})".format(config.read_count_prefix),
        "description": "Number of clusters for each tile ({})".format(
            config.read_count_desc
        ),
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
    },
    "Cluster Count Pf": {
        "title": "Clusters PF ({})".format(config.read_count_prefix),
        "description": "Number of clusters PF for each tile ({})".format(
            config.read_count_desc
        ),
        "shared_key": "cluster_count",
        "modify": lambda x: x * config.read_count_multiplier,
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
    "Total Pf Reads": {
        "title": "{} PF Reads".format(config.read_count_prefix),
        "description": "The total number of passing filter reads for this lane ({})".format(
            config.read_count_desc
        ),
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Total Reads": {
        "title": "{} Reads".format(config.read_count_prefix),
        "description": "The total number of reads for this lane ({})".format(
            config.read_count_desc
        ),
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Mapped Reads Cv": {
        "title": "CV",
        "description": "The coefficient of variation for the number of counts across all indexes.",
        "format": "{:,.2f}",  # 2 decimal places please
    },
    "Max Mapped Reads": {
        "title": "{} Max Mapped Reads".format(config.read_count_prefix),
        "description": "The highest representation for any index ({})".format(
            config.read_count_desc
        ),
        "modify": lambda x: float(x) * config.read_count_multiplier,
        "format": "{:,.2f}",
        "shared_key": "read_count",
    },
    "Min Mapped Reads": {
        "title": "{} Min Mapped Reads".format(config.read_count_prefix),
        "description": "The lowest representation for any index ({})".format(
            config.read_count_desc
        ),
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
    "Project Name": {"title": "Project Name", "description": "Sample Project Name",},
    "Sample Id": {
        "title": "Sample ID",
        "description": "The Sample ID given in the SampleSheet",
    },
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

        super(SAV, self).__init__(
            name="Illumina SAV",
            anchor="sav",
            info=" - Sequencing Metrics from Illumina sequencers",
        )

        # Check if required files are found
        runinfo_files = []
        runparam_files = []
        for f in self.find_log_files("SAV/runinfo"):
            runinfo_files.append(os.path.join(f["root"], f["fn"]))
        for f in self.find_log_files("SAV/runparameters"):
            runparam_files.append(os.path.join(f["root"], f["fn"]))

        # Assume single run for now
        if (
            (os.path.dirname(runinfo_files[0]) == os.path.dirname(runparam_files[0]))
            and len(
                glob.glob(
                    os.path.join(os.path.dirname(runinfo_files[0]), "InterOp/*.bin")
                )
            )
            > 0
        ):
            illumina_dir = os.path.dirname(runinfo_files[0])
        else:
            log.warning(
                "Skipping MultiQC_SAV, required files were not found or not in the right structure."
            )
            return None

        self.load_metrics(illumina_dir)
        self.summary_qc()
        # self.indexing_qc()
        self.imaging_qc()

    def load_metrics(self, illumina_dir):
        log.info("Loading run metrics from {}".format(illumina_dir))
        self.run_metrics = interop.read(run=illumina_dir, finalize=True)

    #############
    # SUMMARY QC
    #############

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

        # - GRAPH: Clusters/Lane
        log.info("Generating 'Clusters/Lane' plot")
        self.add_section(
            name="Summary Metrics - Clusters/Lane",
            anchor="sav-imaging-clusters-lane",
            description="",
            plot=self.clusters_lane_plot(self.parse_lane_summary(summary_lane)),
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

    def clusters_lane_plot(self, data):
        cluster_data = {}
        read_data = {}
        for key, value in data.items():
            lane = int(value["Lane"])
            if f"Lane {lane}" not in cluster_data:
                cluster_data[f"Lane {lane}"] = {
                    "clusters": value["Cluster Count"],
                    "clusters_pf": value["Cluster Count Pf"],
                    "clusters_diff": value["Cluster Count"] - value["Cluster Count Pf"],
                }
                read_data[f"Lane {lane}"] = {
                    "reads": value["Reads"],
                    "reads_pf": value["Reads Pf"],
                    "reads_diff": value["Reads"] - value["Reads Pf"],
                }
            else:
                cluster_data[f"Lane {lane}"]["clusters"] += value["Cluster Count"]
                cluster_data[f"Lane {lane}"]["clusters_pf"] += value["Cluster Count Pf"]
                cluster_data[f"Lane {lane}"]["clusters_diff"] += (
                    value["Cluster Count"] - value["Cluster Count Pf"]
                )
                read_data[f"Lane {lane}"]["reads"] += value["Reads"]
                read_data[f"Lane {lane}"]["reads_pf"] += value["Reads Pf"]
                read_data[f"Lane {lane}"]["reads_diff"] += (
                    value["Reads"] - value["Reads Pf"]
                )

        cats = [OrderedDict(), OrderedDict()]
        cats[0]["clusters_pf"] = {"name": "Clusters PF"}
        cats[0]["clusters_diff"] = {"name": "Clusters not PF"}
        cats[1]["reads_pf"] = {"name": "Reads PF"}
        cats[1]["reads_diff"] = {"name": "Reads not PF"}

        plot_config = {
            "id": "sav-summary-clusters-reads-lane-plot",
            "title": "SAV - Cluster/Reads per Lane",
            "data_labels": ["Clusters", "Reads"],
        }

        return bargraph.plot([cluster_data, read_data], cats, plot_config)

    def _parse_reads(self, reads_df, key_prefix: str = None):
        reads_dict = {}
        reads_df = reads_df.set_index("ReadNumber")
        for read, data in reads_df.iterrows():
            key = f"Read {read}" + " (I)" if data["IsIndex"] == 89 else f"Read {read}"
            if key_prefix:
                key = f"{key_prefix} - {key}"
            reads_dict[key] = data.drop("IsIndex").to_dict()
        return reads_dict

    #############
    # WIP: INDEXING QC
    #############
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
        try:
            index_summary_barcode = pd.DataFrame(
                interop.indexing(self.run_metrics, level="Barcode")
            )

            self.add_section(
                name="Indexing Barcode Metrics",
                anchor="sav-barcode-index",
                description="Indexing metrics per Barcode",
                plot=self.lane_index_summary_table(
                    self.parse_lane_index_summary(index_summary_barcode)
                ),
            )
        except IndexError as e:
            log.warning(f"Unable to parse Barcode Indexing Metrics: {e}")

    def parse_lane_index_summary(self, data):
        data = data.set_index("Lane")
        table_data = {}
        for lane, lane_data in data.iterrows():
            table_data[f"Lane {lane}"] = lane_data.to_dict()
        return table_data

    def lane_index_summary_table(self, data):
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
        data.set_index("Id")
        table_data = {}
        lanes = index_summary_barcode.groupby("Lane")
        for lane, lane_data in lanes:
            lane_data = lane_data.set_index("Sample Id")
            for sample, sample_data in lane_data.iterrows():
                table_data[f"{sample} - Lane {lane}"] = sample_data.drop(
                    "Lane"
                ).to_dict()
        return table_data

    def barcode_index_summary_table(self, data):
        headers = {
            header: HEADERS[header]
            for header in interop.index_summary_columns(level="Lane")
        }
        table_config = {
            "namespace": "SAV",
            "id": "sav-barcode-index-metrics-summary-table",
            "col1_header": "Sample - Lane",
        }

    #############
    # IMAGING QC
    #############
    def imaging_qc(self):
        log.info("Gathering Imaging metrics")
        imaging = pd.DataFrame(interop.imaging(self.run_metrics))

        # # - GRAPH: Intensity/Cycle/Channel
        # log.info("Generating 'Intensity/Cycle' plot")
        # self.add_section(
        #     name="Imaging Metrics - Intensity/Cycle",
        #     anchor="sav-imaging-intensity-cycle",
        #     description="",
        #     plot=self.intensity_cycle_plot(imaging),
        # )

        # # - GRAPH: Qscore Heatmap
        # log.info("Generating 'Qscore Heatmap' plot")
        # self.add_section(
        #     name="Imaging Metrics - Qscore Heatmap",
        #     anchor="sav-imaging-qscore-heatmap",
        #     description="",
        #     plot=self.qscore_heatmap_plot(imaging),
        # )

        # # - GRAPH: Qscore Histogram
        # log.info("Generating 'Qscore Histogram' plot")
        # self.add_section(
        #     name="Imaging Metrics - Qscore Histogram",
        #     anchor="sav-imaging-qscore-histogram",
        #     description="",
        #     plot=self.qscore_histogram_plot(imaging),
        # )

        # - GRAPH: %Occ/%PF
        log.info("Generating '% PF vs % Occupied' plot")
        self.add_section(
            name="Imaging Metrics - % PF vs % Occupied",
            anchor="sav-imaging-pf-vs-occ",
            description="",
            plot=self.occ_vs_pf_plot(imaging),
        )

    def intensity_cycle_plot(self, data):
        plot_data = {}
        return linegraph.plot(plot_data)

    def qscore_heatmap_plot(self, data):
        plot_data = {}
        return heatmap.plot(plot_data)

    def qscore_histogram_plot(self, data):
        plot_data = {}
        return bargraph.plot(plot_data)

    def occ_vs_pf_plot(self, data):
        # set color scale
        cscale = mqc_colour.mqc_colour_scale()
        colors = cscale.get_colours("Dark2")

        # filter relevant colums
        data = data[["Lane", "% Pass Filter", "% Occupied"]]
        # split by lane
        per_lane = data.groupby("Lane")
        plot_data = {}
        for lane, lane_data in per_lane:
            lane = int(lane)
            if not f"Lane {lane}" in plot_data:
                plot_data[f"Lane {lane}"] = []
                prev_x = 0
                prev_y = 0
                for idx, tile in lane_data.iterrows():
                    x = float(tile["% Occupied"])
                    y = float(tile["% Pass Filter"])
                    if x != prev_x or y != prev_y:
                        prev_x = x
                        prev_y = y
                        plot_data[f"Lane {lane}"].append(
                            {"x": x, "y": y, "color": colors[lane]}
                        )

        plot_config = {
            "id": "sav-imaging-pf-vs-occ-plot",
            "title": "SAV - %Occupied/%Passing Filter",
            "xlab": "% Occupied",
            "ylab": "% Passing Filter",
            "xmin": 0,
            "xmax": 100,
            "ymin": 0,
            "ymax": 100,
        }

        return scatter.plot(plot_data, plot_config)

