#### BOX 1 ####
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj,
  receiver = "CD8 T",
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
  condition_colname = "aggregate",
  condition_oi = "LCMV", condition_reference = "SS",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

#### BOX 2 ####
# Define cross-validation folds and number of iterations
k <- 3; n <- 10

# Build random forest model and obtain prediction values
predictions_list <- lapply(1:n, assess_rf_class_probabilities,
                           folds = k, geneset = geneset_oi,
                           background_expressed_genes = background_expressed_genes,
                           ligands_oi = best_upstream_ligands,
                           ligand_target_matrix = ligand_target_matrix)

# Get classification metrics of the models, then calculate mean across all rounds
performances_cv <- bind_rows(lapply(predictions_list,
                                    classification_evaluation_continuous_pred_wrapper))
colMeans(performances_cv)

# Calculate fraction of target genes and non-target genes that are among the top 5% predicted targets
fraction_cv <- bind_rows(lapply(predictions_list,
                                calculate_fraction_top_predicted, quantile_cutoff = 0.95), .id = "round")

# In this case, ~30% of target genes are in the top targets compared to ~1% of the non-target genes
mean(filter(fraction_cv, true_target)$fraction_positive_predicted)
mean(filter(fraction_cv, !true_target)$fraction_positive_predicted)

# Perform Fischer's exact test
lapply(predictions_list, calculate_fraction_top_predicted_fisher,
       quantile_cutoff = 0.95)

# Get which genes had the highest prediction values
top_predicted_genes <- lapply(1:n, get_top_predicted_genes, predictions_list)
top_predicted_genes <- reduce(top_predicted_genes, full_join,
                              by = c("gene","true_target"))

#### BOX 3 ####
# Download each part of the network
zenodo_path <- "https://zenodo.org/record/7074291/files/"
lr_network <- readRDS(url(paste0(zenodo_path, "lr_network_human_21122021.rds")))
sig_network <- readRDS(url(paste0(zenodo_path,
                                  "signaling_network_human_21122021.rds")))
gr_network <- readRDS(url(paste0(zenodo_path, "gr_network_human_21122021.rds")))

# Aggregate the individual data sources in a weighted manner to obtain
# a weighted integrated signaling network
weighted_networks <- construct_weighted_networks(
  lr_network = lr_network,
  sig_network = sig_network,
  gr_network = gr_network,
  source_weights_df = optimized_source_weights_df)

# Downweigh the importance of signaling and gene regulatory hubs
# Use the optimized parameters of this
weighted_networks <- apply_hub_corrections(
  weighted_networks = weighted_networks,
  lr_sig_hub = hyperparameter_list[
    hyperparameter_list$parameter == "lr_sig_hub",]$avg_weight,
  gr_hub = hyperparameter_list[
    hyperparameter_list$parameter == "gr_hub",]$avg_weight)

# In this example, we will calculate target gene regulatory potential scores for
# TNF and the combination TNF+IL6
# To compute it for all 1248 ligands (~1 min):
# ligands <- as.list(unique(lr_network$from))
ligands <- list("TNF", c("TNF","IL6"))
ligand_target_matrix = construct_ligand_target_matrix(
  weighted_networks = weighted_networks,
  ligands = ligands,
  algorithm = "PPR",
  damping_factor = hyperparameter_list[
    hyperparameter_list$parameter == "damping_factor",]$avg_weight,
  ltf_cutoff = hyperparameter_list[
    hyperparameter_list$parameter == "ltf_cutoff",]$avg_weight)

# Using LIANA
devtools::install_github("saezlab/liana")

liana::show_resources()

liana_db <- liana::decomplexify(liana::select_resource("CellPhoneDB")[[1]])
liana_db <- rename(liana_db, from = source_genesymbol, to = target_genesymbol)
liana_db$source <- "liana"
liana_db <- select(liana_db, from, to, source)

# Change source weights data frame (but in this case all source weights are 1)
source_weights <- add_row(source_weights_df, source = "liana",
                          weight = 1, .before = 1)

# Construct weighted network as before
weighted_networks <- construct_weighted_networks(
  lr_network = liana_db,
  sig_network = sig_network,
  gr_network = gr_network,
  source_weights_df = source_weights)
