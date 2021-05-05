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
        self.imaging_qc()

    def load_metrics(self):
        log.debug("SAV: Loading run metrics from {}".format(self.illumina_dir))
        self.run_metrics = interop.read(run=self.illumina_dir)

    def summary_qc(self):
        log.info("SAV: Gathering Read summary metrics")
        summary_read = pd.DataFrame(interop.summary(self.run_metrics, level="Read"))
        self.parse_read_summary(summary_read)

        self.add_section(
            name="Summary Read Metrics",
            anchor="sav-read-summary",
            description="Summary metrics per Read",
        )

        log.info("SAV: Gathering Lane summary metrics")
        summary_lane = (
            pd.DataFrame(interop.summary(self.run_metrics, level="Lane"))
            .transpose()
            .to_dict()
        )
        headers = {key: {} for key in interop.summary_columns(level="Lane")}

        self.add_section(
            name="Summary Lane Metrics",
            anchor="sav-lane-summary",
            description="Summary metrics per Lane per Read",
            plot=table.plot(summary_lane, headers),
        )

    def indexing_qc(self):
        try:
            log.info("SAV: Gathering Lane Indexing metrics")
            indexing_lane_summary = (
                pd.DataFrame(interop.indexing(self.run_metrics, level="Lane"))
                .transpose()
                .to_dict()
            )

            self.add_section(
                name="Indexing Lane Metrics",
                anchor="sav-lane-index",
                description="Indexing metrics per Lane",
                plot=table.plot(indexing_lane_summary),
            )

            log.info("SAV: Gathering Barcode Indexing metrics")
            indexing_barcode_summary = (
                pd.DataFrame(interop.indexing(self.run_metrics, level="Barcode"))
                .transpose()
                .to_dict()
            )

            self.add_section(
                name="Indexing Barcode Metrics",
                anchor="sav-barcode-index",
                description="Indexing metrics per Barcode",
                plot=table.plot(indexing_barcode_summary),
            )

        except ValueError:
            log.debug("SAV: No indexing metrics available")

    def imaging_qc(self):
        log.debug("SAV: Gathering Imaging metrics")
        imaging = pd.DataFrame(interop.imaging(self.run_metrics))

        self.add_section(
            name="Imaging Metrics", anchor="sav-imaging", description="",
        )

    def parse_read_summary(self, metrics):
        # format pandas dataframe
        metrics.set_index("ReadNumber").transpose()
        parsed_metrics = {}
        for read_number, read_data in metrics.iterrows():
            metrics.at[read_number, "IsIndex"] = (
                False if read_data["IsIndex"] == 89 else True
            )
        print(metrics)

    # def parse_lane_summary(self, metrics):
    #     print("test")

    # def parse_lane_indexing(self, metrics):
    #     print("test")

    # def parse_barcode_indexing(self, metrics):
    #     print("test")

    def summary_read_table(self, data):
        headers = interop.summary_columns()
        headers = OrderedDict()
        headers["Yield"] = {
            "rid": "summary_Yield",
            "title": "{}p Yield".format(config.base_count_prefix),
            "description": 'The number of bases sequenced ({} base pairs over all "usable cycles"'.format(
                config.base_count_desc
            ),
            "scale": "PuOr",
            "shared_key": "base_count",
            "modify": lambda x: (x * 1000000000.0)
            * config.base_count_multiplier,  # number is already in gigabases
        }
        headers["Aligned"] = {
            "rid": "summary_Aligned",
            "title": "Aligned (%)",
            "description": "The percentage of the sample that aligned to the PhiX genome",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "PiYG",
        }
        headers["Error Rate"] = {
            "title": "Error Rate (%)",
            "description": "",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "OrRd",
        }
        headers["Intensity C1"] = {
            "rid": "summary_Intensity_C1",
            "title": "Intensity Cycle 1",
            "description": "The intensity statistic at cycle 1.",
        }
        headers["%>=Q30"] = {
            "rid": "summary_Q30",
            "title": "% >= Q30",
            "description": "Percentage of reads with quality phred score of 30 or above",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        }
        table_config = {
            "namespace": "interop",
            "id": "interop-runmetrics-summary-table",
            "table_title": "Read metrics summary",
            "col1_header": "Run - Read",
        }

        tdata = {}
        for s_name in data:
            for key in data[s_name]["summary"]:
                tdata["{} - {}".format(s_name, key)] = data[s_name]["summary"][key]

        return table.plot(tdata, headers, table_config)
