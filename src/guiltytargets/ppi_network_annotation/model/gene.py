# -*- coding: utf-8 -*-

"""This module contains the class Gene."""

from dataclasses import dataclass

from dataclasses_json import dataclass_json

__all__ = [
    'Gene',
]


@dataclass_json
@dataclass
class Gene:
    """Encapsulate a gene and its attributes."""

    #: Entrez Gene identifier
    entrez_id: str

    #: log2 fold change of the gene
    log2_fold_change: float = 0.0

    #: Adjusted p-value
    padj: float = 1.0

    #: Gene symbol (from its species-nomenclature consortium)
    symbol: str = ''
