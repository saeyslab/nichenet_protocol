
#### BOX 1 ####
## Model construction
# The complete networks can be downloaded from Zenodo

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

# aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df = source_weights_df)


# downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                          lr_sig_hub = hyperparameter_list[hyperparameter_list$parameter == "lr_sig_hub",]$avg_weight,
                                          gr_hub = hyperparameter_list[hyperparameter_list$parameter == "gr_hub",]$avg_weight)

ligands = list("TNF",c("TNF","IL6"))
ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                      damping_factor = hyperparameter_list[hyperparameter_list$parameter == "damping_factor",]$avg_weight,
                                                      ltf_cutoff = hyperparameter_list[hyperparameter_list$parameter == "ltf_cutoff",]$avg_weight)


liana_db <- liana::decomplexify(liana::select_resource("Consensus")[[1]])
liana_db <- rename(liana_db, from = source_genesymbol, to = target_genesymbol)
liana_db$source <- "liana"
liana_db <- select(liana_db, from, to, source)



# Change source weights dataframe (but in this case all source weights are 1)
source_weights <- add_row(source_weights_df, source = "liana", weight = 1, .before = 1)



# Aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network

weighted_networks <- construct_weighted_networks(lr_network = liana_db,
                                                 sig_network = sig_network,
                                                 gr_network = gr_network,
                                                 source_weights_df = source_weights)

# Plotting Ebi3 regulatory potential between genes in and out of the GSOI
# targets <- setdiff(rownames(ligand_target_matrix), geneset_oi)
# df <- data.frame(val = ligand_target_matrix[,"Ebi3"], gene = rownames(ligand_target_matrix)) %>% mutate(in_target = !(gene %in% targets))
# ggplot(df, aes(x=in_target, y = val)) + geom_violin()

