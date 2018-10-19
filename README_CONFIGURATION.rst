DEFAULT
-------
**minimum_log2_fold_change**: Minimum value of log2 fold change for determining up-regulated genes

**maximum_log2_fold_change**: Maximum value of log2 fold change for determining down-regulated genes

**maximum_adjusted_p_value**: Maximum value of adjusted p-value, used as additional criteria for determining differentially expressed genes


PATHS
-----
**input_directory**: Home directory for input files. If you prefer not to have such a directory, leave empty and then give absolute paths for input files

**output_directory**: Home directory for output files

**protein_protein_interaction_graph**: Path to protein-protein interaction graph w.r.t. input_directory. Details can be found in README.md

**experiment_file**: Path to experiment file w.r.t. input_directory. Details can be found in README.md

**drug_targets_file**: Path to drug targets file w.r.t. input_directory. Details can be found in README.md

**transcription_factors_file**: Path to transcription factors file w.r.t. input_directory. Details can be found in README.md. This file is not mandatory, leave blank after = sign, if you prefer not to give this information.


EXP_FILE
-------
**base_mean_name**: Column name for base mean

**log2_fold_change_name**: Column name for log2 fold change

**adjusted_p_value_name**: Column name for adjusted p value

**entrez_id_name**: Column name for Entrez id

**split_character**: If there are multiple Entrez ids and symbols in a cell, the character that separates them

**sheet_name**: The name of the sheet. Can be left empty if the experiment_file is a csv file
