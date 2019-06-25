Protein-Protein Interaction (PPI) Network Annotation |build|
============================================================
A library for annotating a protein-protein interaction network using differentially expressed gene
information

.. |build| image:: https://travis-ci.com/GuiltyTargets/ppi-network-annotation.svg?branch=master
    :target: https://travis-ci.com/GuiltyTargets/ppi-network-annotation

INPUT FILES
-----------
There are 4 files which are necessary to run this program:

1. A protein-protein interaction network in the format of:

    **Entrez ID** **Entrez ID** **Confidence score**

    Such as:

    216 216 0.76

    3679 1134 0.73

    55607 71 0.65

    5552 960 0.63

    2886 2064 0.9

    5058 2064 0.73

    1742 2064 0.87


    An example of such a network can be found [here](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)


2. An experiment file, in Excel format. Rows show individual entries, columns are the values of the following properties:

    - **Log2 fold change**
    - **Adjusted p value**
    - **Entrez id**

    The file may contain other columns too, but the indices and names of the above columns must
    be entered as parameters
