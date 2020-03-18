# -*- coding: utf-8 -*-

"""For annotating a protein protein interaction network with differential gene expression."""

from .model import AttributeNetwork, FilteredNetwork, Gene, LabeledNetwork, Network  # noqa: F401
from .pipeline import generate_ppi_network, parse_dge  # noqa: F401
