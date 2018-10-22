# -*- coding: utf-8 -*-

"""Configuration for gene-prioritization."""

import os
from configparser import ConfigParser

__all__ = [
    'get_config',
]


def get_config(path) -> ConfigParser:
    assert os.path.exists(path), f'path does not exist: {path}'

    cfp = ConfigParser()
    cfp.read(path)

    # Read input file names and join the name with the input directory path
    # TODO: add disease associations
    input_directory = cfp['paths']['input_directory']
    ppi_path = os.path.join(input_directory, cfp.get('paths', "protein_protein_interaction_graph"))
    data_path = os.path.join(input_directory, cfp.get('paths', 'differential_gene_expression'))
    targets_path = os.path.join(input_directory, cfp.get('paths', "drug_targets"))
    cfp.set('paths', 'ppi_path', ppi_path)
    cfp.set('paths', 'data_path', data_path)
    cfp.set('paths', 'targets_path', targets_path)

    # Set up the intermediary and output file names
    output_directory = cfp.get('paths', 'output_directory')
    dataset = os.path.split(output_directory)[1]
    #Gat2Vec
    cfp.add_section('gat2vec')
    cfp.set('gat2vec', 'home', output_directory)
    cfp.set('gat2vec', 'adjacency_list', os.path.join(output_directory, f'{dataset}_graph.adjlist'))
    cfp.set('gat2vec', 'attribute_adjacency_list', os.path.join(output_directory, f'{dataset}_na.adjlist'))
    cfp.set('gat2vec', 'mapped_labels', os.path.join(output_directory, "labels_maped.txt"))
    cfp.set('gat2vec', 'auc_g2v', os.path.join(output_directory, "auc_g2v.tsv"))
    cfp.set('gat2vec', 'predictions_path', os.path.join(output_directory, "probs.tsv"))

    return cfp
