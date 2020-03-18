# -*- coding: utf-8 -*-

"""Constants for gene-prioritization."""

import os
from typing import Tuple

from easy_config import EasyConfig

EMOJI = '🦑'

_CONFIG_DIRECTORY = os.path.join(os.path.expanduser('~'), '.config')
_GUILTY_TARGETS_CONFIG_DIRECTORY = os.path.join(_CONFIG_DIRECTORY, 'guiltytargets')

CONFIG_FILE_PATHS = [
    os.path.join(_CONFIG_DIRECTORY, 'guiltytargets.cfg'),
    os.path.join(_CONFIG_DIRECTORY, 'guiltytargets.ini'),
    os.path.join(_GUILTY_TARGETS_CONFIG_DIRECTORY, 'cfg.ini'),
    os.path.join(_GUILTY_TARGETS_CONFIG_DIRECTORY, 'config.ini'),
]


class GuiltyTargetsConfig(EasyConfig):
    """GuiltyTargets configuration."""

    NAME = 'guiltytargets'
    FILES = CONFIG_FILE_PATHS

    """Required"""

    input_directory: str
    output_directory: str

    # These defaults match the GEO2R files
    entrez_id_header: str = 'Gene.ID'
    log2_fold_change_header: str = 'logFC'
    adj_p_header: str = 'adj.P.Val'
    base_mean_header: str = None

    #: Delimiter between Entrez Gene identifiers
    entrez_delimiter: str = '///'

    """Optional"""

    ppi_graph_file_name: str = 'string.edgelist'
    dge_file_name: str = 'DifferentialExpression.tsv'
    targets_file_name: str = 'targets.txt'

    ppi_edge_min_confidence: float = 0.0

    max_adj_p: float = 0.05

    max_log2_fold_change: float = -1.0
    min_log2_fold_change: float = +1.0

    """Output configuration"""

    #:
    auc_output_file_name: str = 'auc_g2v.tsv'

    #:
    ranked_targets_output_file_name: str = 'rankings.tsv'

    """Derived configuration properties"""

    @property
    def auc_output_path(self) -> str:  # noqa: D102
        return os.path.join(self.output_directory, self.auc_output_file_name)

    @property
    def ranked_targets_output_path(self) -> str:  # noqa: D102
        return os.path.join(self.output_directory, self.ranked_targets_output_file_name)

    @property
    def ppi_graph_path(self) -> str:  # noqa: D102
        return os.path.join(self.input_directory, self.ppi_graph_file_name)

    @property
    def dge_path(self) -> str:  # noqa: D102
        return os.path.join(self.input_directory, self.dge_file_name)

    @property
    def targets_path(self) -> str:  # noqa: D102
        return os.path.join(self.input_directory, self.targets_file_name)


class Gat2VecConfig(EasyConfig):
    """Gat2Vec configuration."""

    NAME = 'gat2vec'
    FILES = CONFIG_FILE_PATHS

    #:
    num_walks: int = 10

    #:
    walk_length: int = 80

    #:
    save_output: bool = True

    #:
    dimension: int = 128

    #:
    window_size: int = 5

    #:
    multilabel: bool = False

    #:
    training_ratio: Tuple[float] = (0.1, 0.3, 0.5)


gat2vec_config = Gat2VecConfig.load()
