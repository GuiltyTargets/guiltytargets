# -*- coding: utf-8 -*-

"""Utilities for GuiltyTargets-Results."""

import os
from typing import Optional, TextIO

import pandas as pd
from mygene import MyGeneInfo
from opentargets import OpenTargetsClient

__all__ = [
    'download_hippie',
    'download_targets_for_disease',
]


def download_hippie(*, url: str, path: str) -> str:
    """Download HIPPIE data and return the path where it is."""
    if os.path.exists(path):
        return path
    cols = ['symbol1', 'entrez1', 'symbol2', 'entrez2', 'confidence', 'description']
    df = pd.read_csv(url, sep='\t', header=None, names=cols)
    df[['entrez1', 'entrez2', 'confidence']].to_csv(path, sep='\t', header=False, index=False)
    return path


def download_targets_for_disease(
    disease_efo_id: str,
    open_targets_client: Optional[OpenTargetsClient] = None,
    my_gene_info: Optional[MyGeneInfo] = None,
    file: Optional[TextIO] = None,
) -> None:
    """Download targets for a given disease.

    :param disease_efo_id: A disease's EFO identifier
    :param open_targets_client: An OpenTargetsClient
    :param my_gene_info: A MyGeneInfo client
    :param file: Place to output targets for disease
    """
    if open_targets_client is None:
        open_targets_client = OpenTargetsClient()
    associations = open_targets_client.get_associations_for_disease(
        disease_efo_id,
        fields=[
            'associationscore.datatypes',
            'target.id',
        ],
    ).filter(
        datatype='known_drug',
    )
    ensembl_list = [
        association['target']['id']
        for association in associations
    ]

    if my_gene_info is None:
        my_gene_info = MyGeneInfo()

    id_mappings = my_gene_info.getgenes(ensembl_list, fields="entrezgene")

    print('efo', 'ncbigene', file=file, sep='\t')
    for mapping in id_mappings:
        entrez_gene_id = mapping.get('entrezgene')
        if entrez_gene_id is not None:
            print(disease_efo_id, entrez_gene_id, file=file, sep='\t')
