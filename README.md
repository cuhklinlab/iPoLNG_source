This is the repo to store the processed data used in the iPoLNG manuscript. We also include a code example ``preprocess_example_code`` to preprocess the original data.

Each data folder consists of the following files:

``<dataname>_barcodes.csv``: the cell barcodes in the data.

``<dataname>_DNA20k.mtx``: the feature by cell matrix for ATAC / histone modification data.

``<dataname>_DNA20kbins.csv``: the selected features for ATAC / histone modification data.

``<dataname>_RNA5k.mtx``: the feature by cell matrix for RNA data.

``<dataname>_RNA5kgenes.csv``: the selected features for RNA data.

For dataset with ground truth labels (Paired-Tag and SHARE-seq), an additional file ``<dataname>_celltype.csv`` is attached.

For 10xPBMC10k data, the ATAC data is divided into two parts for the data storage issue. ``10xPBMC10k_DNA20k_part1.mtx`` is the feature by cell matrix with the first 10,000 features, while ``10xPBMC10k_DNA20k_part2.mtx`` is the feature by cell matrix with the last 10,000 features. Interested users should combine the features together to get the full matrix.
