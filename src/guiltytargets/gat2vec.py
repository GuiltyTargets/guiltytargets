# -*- coding: utf-8 -*-

"""Wrappers around GAT2VEC functions."""

from GAT2VEC import paths as gat2vec_paths
from GAT2VEC.evaluation.classification import Classification
from GAT2VEC.gat2vec import Gat2Vec

__all__ = [
    'gat2vec_paths',
    'Classification',
    'Gat2Vec',
]
