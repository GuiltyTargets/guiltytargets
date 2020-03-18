# -*- coding: utf-8 -*-

"""This module contains the class FilteredNetwork."""

import logging

import numpy as np
from igraph import Graph

from .network import Network

__all__ = [
    'FilteredNetwork',
]

logger = logging.getLogger(__name__)


class FilteredNetwork:
    """Mimic encapsulation of filtered PPI networks."""

    def __init__(self, network: Network):
        """Initialize the network object.

        :param network: A PPI network annotated with differential gene expression
        """
        self.graph = network.graph

    def get_upregulated_genes_network(self) -> Graph:
        """Get the graph of up-regulated genes.

        :return Graph: Graph of up-regulated genes.
        """
        logger.info("In get_upregulated_genes_network()")

        deg_graph = self.graph.copy()  # deep copy graph
        not_diff_expr = self.graph.vs(up_regulated_eq=False)

        # delete genes which are not differentially expressed or have no connections to others
        deg_graph.delete_vertices(not_diff_expr.indices)
        deg_graph.delete_vertices(deg_graph.vs.select(_degree_eq=0))

        return deg_graph

    def get_downregulated_genes_network(self) -> Graph:
        """Get the graph of down-regulated genes.

        :return Graph: Graph of down-regulated genes.
        """
        logger.info("In get_downregulated_genes_network()")

        deg_graph = self.graph.copy()  # deep copy graph
        not_diff_expr = self.graph.vs(down_regulated_eq=False)

        # delete genes which are not differentially expressed or have no connections to others
        deg_graph.delete_vertices(not_diff_expr.indices)
        deg_graph.delete_vertices(deg_graph.vs.select(_degree_eq=0))

        return deg_graph

    def get_shortest_paths_graph(
        self,
        genes_to_keep: list = None,
        keep_isolated_nodes: bool = False,
    ):
        """Get the shortest paths graph between differentially expressed + special genes.

        :genes_to_keep list: A list of special genes.
        :keep_isolated_nodes bool: Removes the vertices with no neighbors when False.
        :return Graph: The shortest paths graph between special genes.
        """
        logger.info("In get_shortest_paths_graph()")
        sp_graph = self.graph.copy()
        weights = list(1 - np.array(sp_graph.es['weight']))
        sp_graph.es['weight'] = weights

        # Get the indices of all genes which are not to be deleted
        # (genes_to_keep + diff.expr.)
        relevant_gene_ind = self.graph.vs.select(diff_expressed=True).indices
        if genes_to_keep is not None:
            genes_to_keep_ind = self.graph.vs.select(name_in=genes_to_keep).indices
            relevant_gene_ind = set(relevant_gene_ind).union(set(genes_to_keep_ind))

        # Calculate the shortest paths between relevant genes and save the edges
        # that reside in shortest paths to a set
        shortest_path_edges = set()
        for ind in relevant_gene_ind:
            shortest_paths = sp_graph.get_shortest_paths(
                ind,
                to=relevant_gene_ind,
                weights='weight',
            )
            for path in shortest_paths:
                for i in range(len(path) - 1):
                    eid = sp_graph.get_eid(path[i], path[i + 1])
                    shortest_path_edges.add(eid)

        # get and remove irrelevant edges
        irrelevant_edges = set(sp_graph.es.indices) - shortest_path_edges
        sp_graph.delete_edges(irrelevant_edges)
        if not keep_isolated_nodes:
            sp_graph.delete_vertices(sp_graph.vs.select(_degree_eq=0))

        return sp_graph
