#!/usr/bin/env python
""" MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters. """

import click

illumina_dir = click.option(
    "--illumina_dir",
    type=click.Path(exists=True, file_okay=False, readable=True),
    help="MultiQC_SAV plugin: Illumina sequencer output directory",
)

