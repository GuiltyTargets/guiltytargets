# -*- coding: utf-8 -*-

"""Module to test network module under model package."""

import unittest

import numpy as np
from igraph import Graph

from guiltytargets.ppi_network_annotation.model.filtered_network import FilteredNetwork
from guiltytargets.ppi_network_annotation.model.gene import Gene
from guiltytargets.ppi_network_annotation.model.network import Network


class NetworkTest(unittest.TestCase):
    """Class to test network module."""

    def setUp(self):
        """Initialize two dummy graphs and a list of genes for the test."""
        self.interact_network = Graph()
        self.interact_network.add_vertices(11)
        self.interact_network.vs["name"] = [
            str(0), str(1), str(2), str(3),
            str(4), str(5), str(6), str(7),
            str(8), str(9), str(10),
        ]
        self.interact_network.add_edges([
            (0, 0),
            (0, 1),
            (0, 3),
            (1, 2),
            (1, 3),
            (2, 6),
            (3, 4),
            (4, 4),
            (4, 5),
            (9, 10),
        ])
        self.interact_network.es["weight"] = [0.9, 0.9, 0.7, 0.7, 0.9, 0.7, 0.9,
                                              0.9, 0.7, 0.9]

        self.symbols = [str(0), str(1), str(2), str(3), str(4), str(5), str(6),
                        str(7), str(8)]
        self.l2fcs = [2, 1, -1, 2, 1, -2, 1, -1, 2]
        self.padjs = [0.01, 0.01, 0.1, 0.01, 0.01, 0.049, 0.01, 0.01, 0.05]
        self.diff_expressed = [True, False, False, True, False, True, False,
                               False, False]
        self.up_regulated = [True, False, False, True, False, False, False,
                             False, False]
        self.down_regulated = [False, False, False, False, False, True, False,
                               False, False]

        self.mapped_network = Graph()
        self.mapped_network.add_vertices(9)
        self.mapped_network.add_edges(
            [(0, 0), (0, 1), (0, 3), (1, 2), (1, 3), (2, 6), (3, 4), (4, 4),
             (4, 5)])
        self.mapped_network.vs["name"] = self.symbols
        self.mapped_network.vs["l2fc"] = self.l2fcs
        self.mapped_network.vs["padj"] = self.padjs
        self.mapped_network.vs["diff_expressed"] = self.diff_expressed
        self.mapped_network.vs["up_regulated"] = self.up_regulated
        self.mapped_network.vs["down_regulated"] = self.down_regulated
        self.mapped_network.es["weight"] = [0.9, 0.9, 0.7, 0.7, 0.9, 0.7, 0.9,
                                            0.9, 0.7]

        self.protein_list = [
            Gene(
                entrez_id=str(i),
                log2_fold_change=self.l2fcs[i],
                symbol=self.symbols[i],
                padj=self.padjs[i],
            )
            for i in range(9)
        ]

    def test_init(self):
        """Test the constructor and set_up methods."""
        n = Network(
            self.interact_network,
            max_adj_p=0.05,
            max_l2fc=-1,
            min_l2fc=1,
        )
        n.set_up_network(self.protein_list, gene_filter=True)
        self.__check_for_graph_eq(n.graph, self.mapped_network)

    def test_get_upregulated_genes_network(self):
        """Test the method to get upregulated genes network."""
        de_up = self.mapped_network.copy()
        de_up.delete_vertices([1, 2, 4, 5, 6, 7, 8])
        de_up.delete_vertices(de_up.vs.select(_degree_eq=0))

        n = Network(
            self.interact_network,
            max_adj_p=0.05,
            max_l2fc=-1,
            min_l2fc=1,
        )
        n.set_up_network(self.protein_list)

        fn = FilteredNetwork(n)
        de_up_to_test = fn.get_upregulated_genes_network()
        self.__check_for_graph_eq(de_up, de_up_to_test)

    def test_get_downregulated_genes_network(self):
        """Test the method to get downregulated genes network."""
        de_down = self.mapped_network.copy()
        de_down.delete_vertices([0, 1, 2, 3, 4, 6, 7, 8])
        de_down.delete_vertices(de_down.vs.select(_degree_eq=0))

        n = Network(
            self.interact_network,
            max_adj_p=0.05,
            max_l2fc=-1,
            min_l2fc=1,
        )
        n.set_up_network(self.protein_list)

        fn = FilteredNetwork(n)
        de_down_to_test = fn.get_downregulated_genes_network()
        self.__check_for_graph_eq(de_down, de_down_to_test)

    def test_get_shortest_paths_graph(self):
        """Test method to get shortest paths graph."""
        shortest_path_graph = self.mapped_network.copy()
        shortest_path_graph.delete_vertices([2, 6, 7, 8])
        shortest_path_graph.simplify(combine_edges=max)
        shortest_path_graph.delete_vertices(
            shortest_path_graph.vs.select(_degree_eq=0))
        eid = shortest_path_graph.get_eid("0", "3")
        shortest_path_graph.delete_edges(eid)
        weights = list(1 - np.array(shortest_path_graph.es['weight']))
        shortest_path_graph.es['weight'] = weights

        n = Network(
            self.interact_network,
            max_adj_p=0.05,
            max_l2fc=-1,
            min_l2fc=1,
        )
        n.set_up_network(self.protein_list)

        fn = FilteredNetwork(n)
        shortest_path_graph_to_test = fn.get_shortest_paths_graph()

        self.__check_for_graph_eq(shortest_path_graph,
                                  shortest_path_graph_to_test)

    def __check_for_graph_eq(self, g1, g2):
        """Check if two graphs are the same."""
        self.assertEqual(g1.vs["name"], g2.vs["name"])
        self.assertEqual(g1.vs["l2fc"], g2.vs["l2fc"])
        self.assertEqual(g1.vs["padj"], g2.vs["padj"])
        self.assertEqual(g1.vs["diff_expressed"], g2.vs["diff_expressed"])
        self.assertEqual(g1.vs["up_regulated"], g2.vs["up_regulated"])
        self.assertEqual(g1.vs["down_regulated"], g2.vs["down_regulated"])

        self.assertEqual(g1.es["weight"], g2.es["weight"])
        for e1 in g1.es:
            has_match = any(
                e1.source == e2.source and e1.target == e2.target and e1["weight"] == e2["weight"]
                for e2 in g2.es
            )
            self.assertTrue(has_match)

        for e1 in g2.es:
            has_match = any(
                e1.source == e2.source and e1.target == e2.target and e1["weight"] == e2["weight"]
                for e2 in g1.es
            )
            self.assertTrue(has_match)
