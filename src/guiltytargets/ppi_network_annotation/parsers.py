# -*- coding: utf-8 -*-

"""Parser methods to read the input files."""

import logging
import os
from collections import defaultdict
from typing import List, Set

import igraph
import pandas as pd

from .model.gene import Gene

__all__ = [
    'parse_csv',
    'parse_disease_associations',
    'parse_disease_ids',
    'parse_excel',
    'parse_ppi_graph',
    'parse_gene_list',
]

logger = logging.getLogger(__name__)


def parse_ppi_graph(path: str, min_edge_weight: float = 0.0) -> igraph.Graph:
    """Build an undirected graph of gene interactions from edgelist file.

    :param path: The path to the edgelist file. This file has two columns, and is tab separated.
    :param min_edge_weight: Cutoff to keep/remove the edges, default is 0, but could also be 0.63.
    :return: Protein-protein interaction graph
    """
    logger.info("In parse_ppi_graph()")
    path = os.path.expanduser(path)
    try:
        graph = igraph.read(path, format="ncol", directed=False, names=True)
    except Exception:
        logger.warning(f'Could not read {path}')
        raise
    graph.delete_edges(graph.es.select(weight_lt=min_edge_weight))
    graph.delete_vertices(graph.vs.select(_degree=0))
    logger.info(f"Loaded PPI network.\n"
                f"Number of proteins: {len(graph.vs)}\n"
                f"Number of interactions: {len(graph.es)}\n")
    return graph


def parse_excel(
    file_path: str,
    entrez_id_header,
    log_fold_change_header,
    adjusted_p_value_header,
    entrez_delimiter,
    base_mean_header=None,
) -> List[Gene]:
    """Read an excel file on differential expression values as Gene objects.

    :param file_path: The path to the differential expression file to be parsed.
    :param config.Params params: An object that includes paths, cutoffs and other information.
    :return list: A list of Gene objects.
    """
    logger.info("In parse_excel()")

    df = pd.read_excel(file_path)

    return _handle_dataframe(
        df,
        entrez_id_name=entrez_id_header,
        log2_fold_change_name=log_fold_change_header,
        adjusted_p_value_name=adjusted_p_value_header,
        entrez_delimiter=entrez_delimiter,
        base_mean=base_mean_header,
    )


def parse_csv(
    file_path: str,
    entrez_id_header,
    log_fold_change_header,
    adjusted_p_value_header,
    entrez_delimiter,
    base_mean_header=None,
    sep=",",
) -> List[Gene]:
    """Read a csv file on differential expression values as Gene objects.

    :param str file_path: The path to the differential expression file to be parsed.
    :param config.Params params: An object that includes paths, cutoffs and other information.
    :return list: A list of Gene objects.
    """
    logger.info("In parse_csv()")

    df = pd.read_csv(file_path, sep=sep)

    return _handle_dataframe(
        df,
        entrez_id_name=entrez_id_header,
        log2_fold_change_name=log_fold_change_header,
        adjusted_p_value_name=adjusted_p_value_header,
        entrez_delimiter=entrez_delimiter,
        base_mean=base_mean_header,
    )


def _handle_dataframe(
    df: pd.DataFrame,
    entrez_id_name,
    log2_fold_change_name,
    adjusted_p_value_name,
    entrez_delimiter,
    base_mean=None,
) -> List[Gene]:
    """Convert data frame on differential expression values as Gene objects.

    :param df: Data frame with columns showing values on differential
    expression.
    :param cfp: An object that includes paths, cutoffs and other information.
    :return list: A list of Gene objects.
    """
    logger.info("In _handle_df()")

    if base_mean is not None and base_mean in df.columns:
        df = df[pd.notnull(df[base_mean])]

    df = df[pd.notnull(df[entrez_id_name])]
    df = df[pd.notnull(df[log2_fold_change_name])]
    df = df[pd.notnull(df[adjusted_p_value_name])]

    # try:
    #     import bio2bel_hgnc
    # except ImportError:
    #     logger.debug('skipping mapping')
    # else:
    #     manager = bio2bel_hgnc.Manager()
    #     # TODO @cthoyt

    return [
        Gene(
            entrez_id=entrez_id,
            log2_fold_change=data[log2_fold_change_name],
            padj=data[adjusted_p_value_name],
        )
        for _, data in df.iterrows()
        for entrez_id in str(data[entrez_id_name]).split(entrez_delimiter)
    ]


def parse_gene_list(path: str, graph: igraph.Graph, anno_type: str = "name") -> List:
    """Parse a list of genes and return them if they are in the network.

    :param path: The path of input file.
    :param graph: The graph with genes as nodes.
    :param anno_type: The type of annotation with two options:name-Entrez ID, symbol-HGNC symbol.
    :return: A list of genes, all of which are in the network.
    """
    # read the file
    genes = pd.read_csv(path, header=None)[0].tolist()

    # get those genes which are in the network
    if anno_type == "name":
        genes = [str(int(gene)) for gene in genes]
        ind = graph.vs.select(name_in=genes).indices
    elif anno_type == "symbol":
        ind = graph.vs.select(symbol_in=genes).indices
    else:
        raise Exception(f"The type can either be name or symbol, {anno_type} is not supported")

    genes = graph.vs[ind][anno_type]
    return genes


def parse_disease_ids(path: str) -> Set[str]:
    """Parse the disease identifier file.

    :param path: Path to the disease identifier file.
    :return: List of disease identifiers.
    """
    if os.path.isdir(path) or not os.path.exists(path):
        logger.info("Couldn't find the disease identifiers file. Returning empty list.")
        return set()

    df = pd.read_csv(path, names=["ID"])
    return set(df["ID"].tolist())


def parse_disease_associations(path: str, excluded_disease_ids: set):
    """Parse the disease-drug target associations file.

    :param str path: Path to the disease-drug target associations file.
    :param list excluded_disease_ids: Identifiers of the disease for which drug targets are being predicted.
    :return: Dictionary of drug target-disease mappings.
    """
    if os.path.isdir(path) or not os.path.exists(path):
        logger.info("Couldn't find the disease associations file. Returning empty list.")
        return {}

    disease_associations = defaultdict(list)

    with open(path) as input_file:
        for line in input_file:
            target_id, disease_id = line.strip().split(" ")
            if disease_id not in excluded_disease_ids:
                disease_associations[target_id].append(disease_id)

    return dict(disease_associations)
