GuiltyTargets
=============
OVERVIEW
--------
This is a tool for therapeutic target prioritization using network representation learning. 

REQUIREMENTS
------------
Python 3.6 or above

INSTALLATION
------------
Download this repository, go to the directory it resides and run:

pip3 install -e .


After that, you can use it as a library in Python

.. code-block:: python
  from guiltytargets.pipeline import run
  run(input_directory,
      targets_path,
    ppi_graph_path,
    dge_path,
    auc_output_path,
    probs_output_path,
    max_adj_p=max_padj,
    max_log2_fold_change=lfc_cutoff * -1,
    min_log2_fold_change=lfc_cutoff,
    entrez_id_header=entrez_id_name,
    log2_fold_change_header=log_fold_change_name,
    adj_p_header=adjusted_p_value_name,
    base_mean_header=base_mean_name,
    entrez_delimiter=split_char,
    ppi_edge_min_confidence=confidence_cutoff)

This will create files in paths auc_output_path and probs_output_path, where the former shows the AUC values of cross validation and the latter shows the predicted targets.

The parameters are explained below.

INPUT FILES
-----------
There are 3 files which are necessary to run this program. All input files should be found under input_directory 

1. ppi_graph_path: A protein-protein interaction network in the format of:

    **EntrezID** **EntrezID** **CONFIDENCE**
    
    
    Such as:
    
    216 216 0.76
    
    3679 1134 0.73
    
    55607 71 0.65
    
    5552 960 0.63
    
    2886 2064 0.9
    
    5058 2064 0.73
    
    1742 2064 0.87
    
    An example of such a network can be found [here](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)


2. dge_path: An experiment file, in tsv format. Rows show individual entries, columns are the values of the following properties:
  - **Base mean**
  - **Log fold change**
  - **Adjusted p value**
  - **Entrez id**

  The file may contain other columns too, but the indices and names of the above columns must be entered to the configuration file.

3. targets_path: A list of Entrez ids of known targets, in the format of

    EntrezID1
    
    EntrezID2
    
    ...
    
    
    Such as:
    
    1742
    
    3996
    
    150
    
    152
    
    151





OPTIONS
-------
The options that should be set are:
max_adj_p: Maximum value for adjusted p-value for a gene to be considered differentially expressed.
max_log2_fold_change: Maximum value for log2 fold change for a gene to be considered differentially expressed
min_log2_fold_change: Minimum value for log2 fold change for a gene to be considered differentially expressed
ppi_edge_min_confidence: Minimum confidence score for the edges in PPI network.
entrez_id_header: The column name for the Entrez id in the differential expression file.
log2_fold_change_header: The column name for the log2 fold change in the differential expression file.
adj_p_header: The column name for the adjusted p-value in the differential expression file.
base_mean_header: The column name for the base mean in the differential expression file.
entrez_delimiter: If there is more than one Entrez id per row in the diff. expr. file, the separator betweem them.
