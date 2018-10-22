#!/usr/bin/env python3

"""Pipeline for emig-reimplementation.

This can be run as a module with ``python -m gene_prioritization.cli``
or through the installed command ``gene_prioritization.
"""

import logging
import pandas as pd
import os

from GAT2VEC.evaluation.classification import Classification
from GAT2VEC.gat2vec import Gat2Vec
from ppi_network_annotation.model.network import Network
from ppi_network_annotation.model.attribute_network import AttributeNetwork
from ppi_network_annotation.model.labeled_network import LabeledNetwork
from guiltytargets.constants import *

__all__ = [
    'rank_targets'
]

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))


def write_gat2vec_input_files(
        network: Network,
        targets: list,
        adjacency_list_path: str,
        attribute_adjacency_list_path: str,
        mapped_labels_path: str):
    """Write the input files for gat2vec tool.

    :param Network network: Network object with attributes overlayed on it.
    """
    network.write_adj_list(adjacency_list_path)

    attribute_network = AttributeNetwork(network)
    attribute_network.write_attribute_adj_list(attribute_adjacency_list_path)

    labeled_network = LabeledNetwork(network)
    labeled_network.write_index_labels(targets, mapped_labels_path)


def rank_targets(
        network: Network,
        targets: list,
        input_directory: str,
        adjacency_list_path: str,
        attribute_adjacency_list_path: str,
        mapped_labels_path: str,
        ranked_targets_path: str,
        auc_path: str
) -> pd.DataFrame:
    """
    Rank proteins based on their likelihood of being targets

    :param network: The PPI network annotated with differential gene expression data.
    :param targets: A list of targets.
    :param input_directory: Input directory for Gat2Vec.
    :param adjacency_list_path: The path for writing the structural network as an adjacency list.
    :param attribute_adjacency_list_path: The path for writing the attribute network as an adjacency list.
    :param mapped_labels_path: The path for writing the labels of node indices.
    :param ranked_targets_path: The path for writing the ranked proteins
    :param auc_path: The path for writing the cross validation results
    :return: The classification model
    """
    write_gat2vec_input_files(
        network=network,
        adjacency_list_path=adjacency_list_path,
        attribute_adjacency_list_path=attribute_adjacency_list_path,
        targets=targets,
        mapped_labels_path=mapped_labels_path
    )

    g2v = Gat2Vec(input_directory, input_directory, label=False, tr=TR)
    model = g2v.train_gat2vec(
        NUM_WALKS,
        WALK_LENGTH,
        DIMENSION,
        WINDOW_SIZE,
        output=True,
    )
    clf_model = Classification(input_directory, input_directory, tr=TR)

    if ranked_targets_path and network:
        probs_df = pd.DataFrame(clf_model.get_prediction_probs_for_entire_set(model))
        entrez_ids = network.get_attribute_from_indices(probs_df.index.values,
                                                        attribute_name="name")
        probs_df["Entrez"] = entrez_ids
        probs_df.to_csv(ranked_targets_path, sep="\t")

    results_model = clf_model.evaluate(model, label=False, evaluation_scheme="cv")

    results_model.to_csv(
        auc_path,
        encoding="utf-8",
        sep="\t",
        index=False,
    )

    return results_model
