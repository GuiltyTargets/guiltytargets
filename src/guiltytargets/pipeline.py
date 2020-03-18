# -*- coding: utf-8 -*-

"""Pipeline for GuiltyTargets."""

from typing import List, Tuple

import pandas as pd

from .constants import gat2vec_config
from .gat2vec import Classification, Gat2Vec, gat2vec_paths
from .ppi_network_annotation import AttributeNetwork, LabeledNetwork, Network, generate_ppi_network, parse_dge
from .ppi_network_annotation.parsers import parse_gene_list

__all__ = [
    'run',
    'rank_targets',
]


def run(
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
) -> None:
    """Run the GuiltyTargets pipeline."""
    gene_list = parse_dge(
        dge_path=dge_path,
        entrez_id_header=entrez_id_header,
        log2_fold_change_header=log2_fold_change_header,
        adj_p_header=adj_p_header,
        entrez_delimiter=entrez_delimiter,
        base_mean_header=base_mean_header,
    )
    network = generate_ppi_network(
        ppi_graph_path=ppi_graph_path,
        dge_list=gene_list,
        max_adj_p=max_adj_p,
        max_log2_fold_change=max_log2_fold_change,
        min_log2_fold_change=min_log2_fold_change,
        ppi_edge_min_confidence=ppi_edge_min_confidence,
    )

    targets = parse_gene_list(targets_path, network.graph)

    auc_df, probs_df = rank_targets(
        directory=input_directory,
        targets=targets,
        network=network,
    )

    probs_df.to_csv(
        probs_output_path,
        sep="\t",
    )

    auc_df.to_csv(
        auc_output_path,
        encoding="utf-8",
        sep="\t",
        index=False,
    )


def write_gat2vec_input_files(network: Network, targets: List[str], home_dir: str) -> None:
    """Write the input files for gat2vec tool.

    :param network: Network object with attributes overlayed on it.
    :param targets:
    :param home_dir:
    """
    network.write_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "graph"))

    attribute_network = AttributeNetwork(network)
    attribute_network.write_attribute_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "na"))

    labeled_network = LabeledNetwork(network)
    labeled_network.write_index_labels(targets, gat2vec_paths.get_labels_path(home_dir))


def rank_targets(
    network: Network,
    targets: List[str],
    directory: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Rank proteins based on their likelihood of being targets.

    :param network: The PPI network annotated with differential gene expression data.
    :param targets: A list of targets.
    :param directory: Home directory for Gat2Vec.
    :return: A 2-tuple of the auc dataframe and the probabilities dataframe?
    """
    write_gat2vec_input_files(network=network, targets=targets, home_dir=directory)

    g2v = Gat2Vec(directory, directory, label=False, tr=gat2vec_config.training_ratio)
    model = g2v.train_gat2vec(
        gat2vec_config.num_walks,
        gat2vec_config.walk_length,
        gat2vec_config.dimension,
        gat2vec_config.window_size,
        output=True,
    )
    classifier = Classification(directory, directory, tr=gat2vec_config.training_ratio)

    auc_df = classifier.evaluate(model, label=False, evaluation_scheme="cv")
    probs_df = get_rankings(classifier, model, network)

    return auc_df, probs_df


def get_rankings(
    classifier: Classification,
    embedding: pd.DataFrame,
    network: Network,
) -> pd.DataFrame:
    """Save the predicted rankings to a file.

    :param classifier: Classification model.
    :param embedding: Embedding model
    :param network: PPI network with annotations
    """
    probs_df = pd.DataFrame(classifier.get_prediction_probs_for_entire_set(embedding))
    probs_df['Entrez'] = network.get_attribute_from_indices(
        probs_df.index.values,
        attribute_name='name',
    )
    return probs_df
