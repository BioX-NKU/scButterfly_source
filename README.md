# The source code for scButterfly

The source code for the reproducibility of "A versatile single-cell cross-modality translation method via dual-aligned variational autoencoders".

We also provide [a software of scButterfly](https://github.com/BioX-NKU/scButterfly) and [a detailed tutorial and document](http://scbutterfly.readthedocs.io/). 

## Usage instructions

```
experiment                                # experiemnt for scButterfly
|
|---RNA_ATAC                              # translation between RNA and ATAC profiles
|   |
|   |---batch_correction_step             # add batch correction step before training model
|   |   |
|   |   |---paired_data                   # add batch correction step for paired data
|   |   |---unpaired_data                 # add batch correction step for unpaired data
|   |   |---simulated_data                # simulating data with higher batch effects based on the CL dataset
|   |
|   |---CCA_anchors                       # data augmentation incorporating with CCA anchors
|   |---cross_validation_by_cell          # five-fold cross-validation by cell
|   |---cross_validation_by_celltype      # three-fold cross-validation by cell type
|   |---cross_validation_by_batch         # four-fold cross-validation by batch
|   |---data_enhancement                  # data enhancement on the MK and MB datasets
|   |---defense_noise                     # anti-noise experiment on the MCC dataset
|   |
|   |---down_sampling                     # down sampling for dataset
|   |   |
|   |   |---cells                         # down sampling for cells on the MCC dataset
|   |   |---features                      # down sampling for features on the CL dataset
|   |
|   |---GO_analysis                       # gene ontology analysis on the MDS dataset
|   |---integrative_analysis              # integrative analysis in the BMMC dataset
|   |---SNPsea_analysis                   # SNPsea analysis for translated profiles on the PBMC dataset
|   |
|   |---sensitive_analyses                # sensitive analyses for multiple hyperparameters
|   |   |
|   |   |---augmented_count               # cell count of data augmentation
|   |   |---epochs_patience               # pretraining epochs, integrative training epochs, and patience for early-stop 
|   |   |---loss_weights                  # four weights for loss calculation
|   |
|   |---TF_regulatory_network_inference   # TF regulatory network inference with DeepTFni on the PBMC dataset
|   |---unpaired_data                     # five-fold cross-validation with unpaired data
|
|---RNA_ADT                               # translation between RNA and ADT profiles
|   |
|   |---consecutive_translation           # consecutive translation from ATAC to RNA to ADT
|   |---cross_validation_by_cell          # five-fold cross-validation by cell
|   |---cross_validation_by_batch         # cross-validation by batch
|
|---Perturbational_responses              # translation between single-cell perturbation responses

```
We provide reproducible methods through `.sh` and `.ipynb` files in each folder.
