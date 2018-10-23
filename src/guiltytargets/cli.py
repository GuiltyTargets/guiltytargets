import click
import logging
import os

from guiltytargets.pipeline import rank_targets
from guiltytargets.config import get_config
from ppi_network_annotation.pipeline import generate_ppi_network
from ppi_network_annotation.parsers import parse_gene_list

logger = logging.getLogger(__name__)

@click.command()
@click.option('--config-path',
              prompt='Please enter the path to the config file',
              help='Path to config file.')

def main(config_path: str):
    """Run GuiltyTargets.

    :param str config_path: The path to the configuration file.
    """
    logging.basicConfig(level=logging.INFO)

    click.secho('getting config', color='cyan')
    cfp = get_config(config_path)

    ppi_path = cfp['paths']['ppi_path']
    assert os.path.exists(ppi_path)

    data_path = cfp['paths']['data_path']
    assert os.path.exists(data_path)

    targets_path = cfp['paths']['targets_path']
    assert os.path.exists(targets_path)

    hippie_min_edge_weight = cfp.getfloat('default', 'interaction_confidence_cutoff')
    current_disease_ids_path = cfp.get('paths', 'disease_ids')
    disease_associations_path = cfp.get('paths', 'disease_associations')

    maximum_adjusted_p_value = cfp.getfloat('default', 'maximum_adjusted_p_value')
    maximum_log2_fold_change = cfp.getfloat('default', 'maximum_log2_fold_change')
    minimum_log2_fold_change = cfp.getfloat('default', 'minimum_log2_fold_change')

    entrez_id_header = cfp['dge']['entrez_id']
    log_fold_change_header = cfp['dge']['log2_fold_change']
    adjusted_p_value_header = cfp['dge']['adjusted_p_value']
    split_char = cfp['dge']['split_character']
    base_mean_header = cfp.get('dge', 'base_mean')

    gat2vec_input_directory = cfp['gat2vec']['home']
    os.makedirs(gat2vec_input_directory, exist_ok=True)
    gat2vec_result_output = cfp['gat2vec']['auc_g2v']
    gat2vec_predictions_path = cfp['gat2vec']['predictions_path']

    adjacency_list_path = cfp['gat2vec']['adjacency_list']
    attribute_adjacency_list_path = cfp['gat2vec']['attribute_adjacency_list']

    mapped_labels_path = cfp['gat2vec']['mapped_labels']

    del cfp

    click.secho('generating PPI network', color='cyan')
    network = generate_ppi_network(
        ppi_graph_path=ppi_path,
        gene_expression_file_path=data_path,
        maximum_adjusted_p_value=maximum_adjusted_p_value,
        maximum_log2_fold_change=maximum_log2_fold_change,
        minimum_log2_fold_change=minimum_log2_fold_change,
        entrez_id_header=entrez_id_header,
        log_fold_change_header=log_fold_change_header,
        adjusted_p_value_header=adjusted_p_value_header,
        base_mean_header=base_mean_header,
        split_char=split_char,
        hippie_min_edge_weight=hippie_min_edge_weight,
        current_disease_ids_path=current_disease_ids_path,
        disease_associations_path=disease_associations_path,
    )

    click.secho('ranking targets', color='cyan')
    targets = parse_gene_list(targets_path, network.graph)

    rank_targets(
        home_dir=gat2vec_input_directory,
        targets=targets,
        ranked_targets_path=gat2vec_predictions_path,
        network=network,
        auc_path=gat2vec_result_output
    )
