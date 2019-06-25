# -*- coding: utf-8 -*-

"""This module contains the class LabeledNetwork."""

import logging

from .network import Network

__all__ = [
    'LabeledNetwork',
]

logger = logging.getLogger(__name__)


class LabeledNetwork:
    """Mimic encapsulation of a labeled and annotated PPI network for Gat2Vec."""

    def __init__(self, network: Network):
        """Initialize the network object.

        :param network: A PPI network annotated with differential gene expression and disease association.
        """
        self.graph = network.graph

    def write_index_labels(self, targets, output_path):
        """Write the mappings between vertex indices and labels(target vs. not) to a file.

        :param list targets: List of known targets.
        :param str output_path: Path to the output file.
        """
        label_mappings = self.get_index_labels(targets)

        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                print(k, v, sep='\t', file=file)

    def get_index_labels(self, targets):
        """Get the labels(known target/not) mapped to indices.

        :param targets: List of known targets
        :return: Dictionary of index-label mappings
        """
        target_ind = self.graph.vs.select(name_in=targets).indices
        rest_ind = self.graph.vs.select(name_notin=targets).indices
        label_mappings = {i: 1 for i in target_ind}
        label_mappings.update({i: 0 for i in rest_ind})
        return label_mappings
