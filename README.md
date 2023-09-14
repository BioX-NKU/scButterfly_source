# The source code for scButterfly

The source code for the reproducibility of "A versatile single-cell cross-modality translation method via dual-aligned variational autoencoders".

We also provide [a software of scButterfly](https://github.com/BioX-NKU/scButterfly) and [a detailed tutorial and document](http://scbutterfly.readthedocs.io/). 

## Usage instructions

```
|   experiment                            # experiemnt for scButterfly
|
|---RNA_ATAC                              # translation between RNA and ATAC profiles
|   |
|   |---cross_validation_by_cell          # five-fold cross-validation by cell
|   |---defense_noise                     # anti-noise experiment on MCC dataset
|   |---cross_validation_by_celltype      # three-fold cross-validation by cell type
|   |---cross_validation_by_batch         # four-fold cross-validation by batch
|   |---GO_analysis                       # gene ontology analysis on MDS dataset
|   |---integrative_analysis              # integrative analysis in BMMC dataset
|   |---data_enhancement                  # data enhancement on MK and MB datasets
|   |---unpaired_data                     # five-fold cross-validation with unpaired data
|
|---RNA_ADT                               # translation between RNA and ADT profiles
|   |
|   |---cross_validation_by_cell          # five-fold cross-validation by cell
|   |---cross_validation_by_batch         # cross-validation by batch
|   |---consecutive_translation           # consecutive translation from ATAC to RNA to ADT
|
|---Perturbational_responses              # translation between single-cell perturbation responses

```
We provide reproducible methods through `.sh` and `.ipynb` files in each folder.