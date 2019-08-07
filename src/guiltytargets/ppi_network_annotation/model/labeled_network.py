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

    def write_index_labels(self, targets, output_path, sample_scores: dict = None):
        """Write the mappings between vertex indices and labels(target vs. not) to a file.

        :param list targets: List of known targets.
        :param str output_path: Path to the output file.
        :param str sample_scores: Sample scores from OpenTarget.
        """
        label_mappings = self.get_index_labels(targets)
        print('labeled_network.write_index_labels')
        print('labeled_network._convert_score_to_weightz')
        print('known targets have weight fixed to 1.')

        with open(output_path, "w") as file:
            for k, v in label_mappings.items():
                if sample_scores:
                    if self.graph.vs[k]["name"] in sample_scores:
                        score = self._convert_score_to_weight(v, sample_scores[self.graph.vs[k]["name"]])
                        print(k, v, score, sep='\t', file=file)
                    else:
                        print(k, v, '1.', sep='\t', file=file)
                else:
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

    @staticmethod
    def _convert_score_to_weight(label: int, score: float) -> float:
        """Convert the association score into a weight for the weighted classification. If the
        label is positive the weight is the score. If negative, the weight is 1 - score. This means,
        a high score for a negative label will imply some uncertainty about it being a target.

        :param label: 1 for positive, 0 for negative.
        :param score: The association score.
        :return: The weight.
        """
        if label:
            # return score
            return 1.  # Fix positive labels to weight 100% always.
        else:
            return 1 - score
