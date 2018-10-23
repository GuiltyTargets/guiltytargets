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
from GAT2VEC import paths as gat2vec_paths
from ppi_network_annotation.model.network import Network
from ppi_network_annotation.model.attribute_network import AttributeNetwork
from ppi_network_annotation.model.labeled_network import LabeledNetwork
from guiltytargets.constants import *

__all__ = [
    'rank_targets'
]

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))


def write_gat2vec_input_files(network: Network, targets: list, home_dir: str):
    """Write the input files for gat2vec tool.

    :param Network network: Network object with attributes overlayed on it.
    """
    network.write_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "graph"))

    attribute_network = AttributeNetwork(network)
    attribute_network.write_attribute_adj_list(gat2vec_paths.get_adjlist_path(home_dir, "na"))

    labeled_network = LabeledNetwork(network)
    labeled_network.write_index_labels(targets, gat2vec_paths.get_labels_path(home_dir))


def rank_targets(network: Network, targets: list, home_dir: str) -> pd.DataFrame:
    """
    Rank proteins based on their likelihood of being targets

    :param network: The PPI network annotated with differential gene expression data.
    :param targets: A list of targets.
    :param home_dir: Home directory for Gat2Vec.
    :param ranked_targets_path: The path for writing the ranked proteins
    :param auc_path: The path for writing the cross validation results
    :return: The classification model
    """
    write_gat2vec_input_files(network=network, targets=targets, home_dir=home_dir)

    g2v = Gat2Vec(home_dir, home_dir, label=False, tr=TR)
    model = g2v.train_gat2vec(
        NUM_WALKS,
        WALK_LENGTH,
        DIMENSION,
        WINDOW_SIZE,
        output=True,
    )
    clf_model = Classification(home_dir, home_dir, tr=TR)

    results_model = clf_model.evaluate(model, label=False, evaluation_scheme="cv")

    save_rankings(clf_model, home_dir, model, network)

    results_model.to_csv(
        os.path.join(home_dir, AUC_FILE_NAME),
        encoding="utf-8",
        sep="\t",
        index=False,
    )

    return results_model


def save_rankings(clf_model, home_dir, emb, network):
    """Save the predicted rankings to a file.

    :param clf_model: Classification model.
    :param home_dir: Home directory
    :param emb: Embedding model
    :param network: PPI network with annotations
    """
    probs_df = pd.DataFrame(clf_model.get_prediction_probs_for_entire_set(emb))
    entrez_ids = network.get_attribute_from_indices(
        probs_df.index.values,
        attribute_name="name"
    )
    probs_df["Entrez"] = entrez_ids
    probs_df.to_csv(os.path.join(home_dir, RANKED_TARGETS_FILE), sep="\t")
