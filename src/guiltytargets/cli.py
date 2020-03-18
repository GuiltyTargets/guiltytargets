# -*- coding: utf-8 -*-

"""Command line interface for GuiltyTargets."""

import logging
import os
import warnings

import click
from easy_config.contrib.click import args_from_config
from sklearn.exceptions import UndefinedMetricWarning

from .constants import EMOJI, GuiltyTargetsConfig
from .pipeline import run

__all__ = [
    'main',
]

logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore', category=UndefinedMetricWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)


@click.command()
@args_from_config(GuiltyTargetsConfig)
def main(
    input_directory,
    output_directory,
    targets_file_name,
    ppi_graph_file_name,
    dge_file_name,
    max_adj_p,
    max_log2_fold_change,
    min_log2_fold_change,
    entrez_id_header,
    log2_fold_change_header,
    adj_p_header,
    base_mean_header,
    entrez_delimiter,
    ppi_edge_min_confidence,
    auc_output_file_name,
    ranked_targets_output_file_name,
) -> None:
    """Run the GuiltyTargets pipeline."""
    if not os.path.exists(input_directory):
        raise FileNotFoundError(input_directory)

    targets_path = os.path.join(input_directory, targets_file_name)
    ppi_graph_path = os.path.join(input_directory, ppi_graph_file_name)
    dge_path = os.path.join(input_directory, dge_file_name)

    os.makedirs(output_directory, exist_ok=True)
    auc_output_path = os.path.join(output_directory, auc_output_file_name)
    probs_output_path = os.path.join(output_directory, ranked_targets_output_file_name)

    click.echo(f'{EMOJI} starting GuiltyTargets')
    run(
        input_directory,
        targets_path,
        ppi_graph_path,
        dge_path,
        auc_output_path,
        probs_output_path,
        max_adj_p,
        max_log2_fold_change,
        min_log2_fold_change,
        entrez_id_header,
        log2_fold_change_header,
        adj_p_header,
        base_mean_header,
        entrez_delimiter,
        ppi_edge_min_confidence,
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
