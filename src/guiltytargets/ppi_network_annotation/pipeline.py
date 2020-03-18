# -*- coding: utf-8 -*-

"""Functions to easily set up the network."""

import logging
from typing import List, Optional

from .model.gene import Gene
from .model.network import Network
from .parsers import parse_csv, parse_disease_associations, parse_disease_ids, parse_excel, parse_ppi_graph

__all__ = [
    'generate_ppi_network',
    'parse_dge',
]

logger = logging.getLogger(__name__)


def generate_ppi_network(
    ppi_graph_path: str,
    dge_list: List[Gene],
    max_adj_p: float,
    max_log2_fold_change: float,
    min_log2_fold_change: float,
    ppi_edge_min_confidence: Optional[float] = None,
    current_disease_ids_path: Optional[str] = None,
    disease_associations_path: Optional[str] = None,
) -> Network:
    """Generate the protein-protein interaction network.

    :return Network: Protein-protein interaction network with information on differential expression.
    """
    # Compilation of a protein-protein interaction (PPI) graph (HIPPIE)
    protein_interactions = parse_ppi_graph(ppi_graph_path, ppi_edge_min_confidence)
    protein_interactions = protein_interactions.simplify()

    if disease_associations_path is not None and current_disease_ids_path is not None:
        current_disease_ids = parse_disease_ids(current_disease_ids_path)
        disease_associations = parse_disease_associations(disease_associations_path,
                                                          current_disease_ids)
    else:
        disease_associations = None

    # Build an undirected weighted graph with the remaining interactions based on Entrez gene IDs
    network = Network(
        protein_interactions,
        max_adj_p=max_adj_p,
        max_l2fc=max_log2_fold_change,
        min_l2fc=min_log2_fold_change,
    )
    network.set_up_network(dge_list, disease_associations=disease_associations)

    return network


def parse_dge(
    dge_path: str,
    entrez_id_header: str,
    log2_fold_change_header: str,
    adj_p_header: str,
    entrez_delimiter: str,
    base_mean_header: Optional[str] = None,
) -> List[Gene]:
    """Parse a differential expression file.

    :param dge_path: Path to the file.
    :param entrez_id_header: Header for the Entrez identifier column
    :param log2_fold_change_header: Header for the log2 fold change column
    :param adj_p_header: Header for the adjusted p-value column
    :param entrez_delimiter: Delimiter between Entrez ids.
    :param base_mean_header: Header for the base mean column.
    :return: A list of genes.
    """
    if dge_path.endswith('.xlsx'):
        return parse_excel(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            entrez_delimiter=entrez_delimiter,
            base_mean_header=base_mean_header,
        )

    if dge_path.endswith('.csv'):
        return parse_csv(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            entrez_delimiter=entrez_delimiter,
            base_mean_header=base_mean_header,
        )

    if dge_path.endswith('.tsv'):
        return parse_csv(
            dge_path,
            entrez_id_header=entrez_id_header,
            log_fold_change_header=log2_fold_change_header,
            adjusted_p_value_header=adj_p_header,
            entrez_delimiter=entrez_delimiter,
            base_mean_header=base_mean_header,
            sep="\t",
        )

    raise ValueError(f'Unsupported extension: {dge_path}')
