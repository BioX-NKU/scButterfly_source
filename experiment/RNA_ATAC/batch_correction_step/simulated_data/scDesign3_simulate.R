library(SingleCellExperiment)
library(anndata)
library(scDesign3)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(umap)

rna_sce = readRDS('CL_rna.rds')
BATCH_data <- construct_data(
  sce = rna_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("batch"),
  corr_by = "1",
  ncell=3000,
)

BATCH_marginal <- fit_marginal(
  data = BATCH_data,
  predictor = "gene",
  mu_formula = "cell_type + batch",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 8,
  usebam = FALSE
)

BATCH_copula <- fit_copula(
  sce = rna_sce,
  assay_use = "counts",
  marginal_list = BATCH_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 8,
  input_data = BATCH_data$dat
)

BATCH_marginal_alter <- lapply(BATCH_marginal, function(x) {
  x$fit$coefficients[length(x$fit$coefficients)] <- rnorm(1, mean = 6, sd = 1)
  x$fit$coefficients[length(x$fit$coefficients)-1] <- rnorm(1, mean = 4, sd = 1)
  x
})

BATCH_marginal_alter <- lapply(BATCH_marginal, function(x) {
  indices <- grepl("^batch", names(x$fit$coefficients))  
  x$fit$coefficients[indices] <- rnorm(sum(indices), mean = 2, sd = 4)
  x
})


BATCH_para_alter <- extract_para(
  sce = rna_sce,
  marginal_list = BATCH_marginal_alter,
  n_cores = 8,
  family_use = "nb",
  new_covariate = BATCH_data$newCovariate,
  data = BATCH_data$dat
)

BATCH_newcount_alter <- simu_new(
  sce = rna_sce,
  mean_mat = BATCH_para_alter$mean_mat,
  sigma_mat = BATCH_para_alter$sigma_mat,
  zero_mat = BATCH_para_alter$zero_mat,
  quantile_mat = NULL,
  copula_list = BATCH_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = BATCH_data$dat,
  new_covariate = BATCH_data$newCovariate,
  important_feature = BATCH_copula$important_feature,
  filtered_gene = BATCH_data$filtered_gene
)

rna_new <- AnnData(X=t(BATCH_newcount_alter))
rna_new$obs$cell_type = BATCH_data$newCovariate$cell_type
rna_new$obs$batch = BATCH_data$newCovariate$batch
write_h5ad(rna_new, "CL_simulated_rna_3000.h5ad")


