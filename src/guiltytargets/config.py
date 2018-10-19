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

    input_directory = cfp['paths']['input_directory']
    ppi_path = os.path.join(
        input_directory,
        cfp.get('paths', "protein_protein_interaction_graph")
    )
    cfp.set('paths', 'ppi_path', ppi_path)
    # set values for differentiating differentially expressed genes

    data_path = os.path.join(
        input_directory,
        cfp.get('paths', 'differential_gene_expression')
    )
    cfp.set('paths', 'data_path', data_path)

    targets_path = os.path.join(input_directory, cfp.get('paths', "drug_targets"))
    cfp.set('paths', 'targets_path', targets_path)

    output_directory = cfp.get('paths', 'output_directory')

    cfp.add_section('gat2vec')

    # GAT2VEC
    dataset_name = cfp.get('options', 'dataset')
    gat2vec_home = os.path.join(output_directory, dataset_name)
    cfp.set('gat2vec', 'home', gat2vec_home)
    cfp.set('gat2vec', 'adjacency_list',
            os.path.join(gat2vec_home, f'{dataset_name}_graph.adjlist'))
    cfp.set('gat2vec', 'attribute_adjacency_list',
            os.path.join(gat2vec_home, f'{dataset_name}_na.adjlist'))
    cfp.set('gat2vec', 'mapped_labels', os.path.join(gat2vec_home, "labels_maped.txt"))

    cfp.set('gat2vec', 'auc_g2v', os.path.join(gat2vec_home, "auc_g2v.tsv"))
    cfp.set('gat2vec', 'predictions_path', os.path.join(gat2vec_home, "probs.tsv"))

    # TODO: add disease associations
    # self.DISEASE_ASSOCIATIONS_PATH = os.path.join(self.INPUT_DIR, paths["disease_associations"])
    # self.CURRENT_DISEASE_IDS_PATH = os.path.join(self.INPUT_DIR, paths["disease_ids"])

    return cfp
