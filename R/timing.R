library(microbenchmark)

step1 <- function(seuratObj){
  Idents(seuratObj) <- seuratObj$celltype
  seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")
  sender_celltypes <- c("CD4 T","Treg", "Mono", "NK", "B", "DC")
  receiver <- "CD8 T"
  condition_oi <- "LCMV"
  condition_reference <- "SS"

  list_expressed_genes_sender <- lapply(sender_celltypes, function(celltype) {get_expressed_genes(celltype, seuratObj, pct = 0.10)})
  expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))

  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj,
                                                  pct = 0.10)

  background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]

  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  seurat_obj_receiver <- SetIdent(seurat_obj_receiver,
                                  value = seurat_obj_receiver[["aggregate"]])


  DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                   ident.1 = condition_oi, ident.2 = condition_reference,
                                   min.pct = 0.10)

  geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= 0.05 & abs(DE_table_receiver$avg_log2FC) >= 0.25, ]
  geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in% rownames(ligand_target_matrix)]

  all_ligands <- unique(lr_network$from)
  all_receptors <- unique(lr_network$to)

  expressed_ligands <- intersect(all_ligands, expressed_genes_sender)
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)


  potential_ligands <- lr_network[lr_network$from %in% expressed_ligands &
                                    lr_network$to %in% expressed_receptors, ]

  potential_ligands <- unique(potential_ligands$from)
}

step2 <- function(geneset_oi, background_expressed_genes, ligand_target_matrix, potential_ligands, expressed_receptors, lr_network, weighted_networks_lr_sig){
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands)

  ligand_activities <- ligand_activities[order(ligand_activities$aupr_corrected, decreasing = TRUE), ]

  best_upstream_ligands <- top_n(ligand_activities, 30, aupr_corrected)$test_ligand

  active_ligand_target_links_df <- lapply(best_upstream_ligands, get_weighted_ligand_target_links,
                                          geneset = geneset_oi,
                                          ligand_target_matrix = ligand_target_matrix,
                                          n = 200)

  active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df))
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks_lr_sig)
}


step3 <- function(active_ligand_target_links_df, ligand_target_matrix, ligand_receptor_links_df,
                  seuratObj, sender_celltypes, condition_oi, condition_reference, best_upstream_ligands, ligand_activities){

  active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)

  order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links)))
  order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links))

  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

  p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target,  "Prioritized ligands","Predicted target genes",
                                                 color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")

  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df, best_upstream_ligands, order_hclust = "receptors")

  p_ligand_receptor_network <- make_heatmap_ggplot(t(vis_ligand_receptor_network[, order_ligands]), "Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  celltype_order <- levels(Idents(seuratObj))

  DE_table_top_ligands <- lapply(celltype_order[celltype_order %in% sender_celltypes],
                                 get_lfc_celltype,
                                 seurat_obj = seuratObj,
                                 condition_colname = "aggregate", condition_oi = condition_oi, condition_reference = condition_reference,
                                 min.pct = 0,
                                 features = best_upstream_ligands, celltype_col = "celltype")
  DE_table_top_ligands <- reduce(DE_table_top_ligands, full_join)
  DE_table_top_ligands <- column_to_rownames(DE_table_top_ligands, "gene")

  vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ])

  p_ligand_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc, "Prioritized ligands","LFC in Sender",
                                                 low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",
                                                 x_axis_position = "top", legend_position = "bottom", legend_title = "LFC") +
    theme(axis.text.y = element_text(face = "italic"))

  p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
                       features = rev(best_upstream_ligands), cols = "RdYlBu") +
    coord_flip() +
    theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots


  # ligand activity heatmap
  ligand_aupr_matrix <- column_to_rownames(ligand_activities, "test_ligand")
  ligand_aupr_matrix <- ligand_aupr_matrix[rev(best_upstream_ligands), "aupr_corrected", drop=FALSE]
  colnames(ligand_aupr_matrix) <- "AUPR"

  vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)
  p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "bottom", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
  p_ligand_aupr

  cowplot::plot_grid(p_ligand_aupr + theme(legend.title = element_blank()),
                     p_ligand_lfc + theme(axis.title.y = element_blank()),
                     p_dotplot + scale_y_discrete(position = "right") + theme(legend.position = "bottom", legend.box = "vertical", axis.title = element_blank(), axis.text.y = element_text(size=9), axis.text.x = element_text(size=9, angle=90, hjust=0)),
                     ncol = 3, rel_widths = c(1, 2, 2), align = "hv")

  # Circos plot
  ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj, best_upstream_ligands[1:20], "celltype", func.agg = mean, slot = "data")
  active_ligand_target_links_df$target_type <- "LCMV-DE"
  circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                             active_ligand_target_links_df)
  ligand_colors <- c("General" = "lawngreen", "NK" = "royalblue", "B" = "darkgreen", "Mono" = "violet", "DC" = "steelblue2")
  target_colors <- c( "LCMV-DE" = "tomato")
  vis_circos_obj <- prepare_circos_visualization(circos_links, ligand_colors = NULL, target_colors = target_colors)
  vis_circos_obj <- prepare_circos_visualization(circos_links, ligand_colors = ligand_colors, target_colors = NULL)
  vis_circos_obj <- prepare_circos_visualization(circos_links)
  vis_circos_obj <- prepare_circos_visualization(circos_links, ligand_colors = ligand_colors, target_colors = target_colors)
  draw_circos_plot(vis_circos_obj, transparency = TRUE)

  # Circos plot receptors
  lr_network_top_df <- rename(ligand_receptor_links_df, ligand=from, target=to)
  lr_network_top_df$target_type = "LCMV_CD8T_receptor"
  lr_network_top_df <- inner_join(lr_network_top_df, ligand_type_indication_df)
  receptor_colors <- c("LCMV_CD8T_receptor" = "darkred")

  vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df, ligand_colors = ligand_colors, target_colors = receptor_colors)
  draw_circos_plot(vis_circos_receptor_obj, transparency = TRUE, link.visible = TRUE)
}

step4 <- function(lr_network, seuratObj, condition_oi, expressed_ligands, expressed_receptors, sender_celltypes, receiver){
  lr_network_renamed <- rename(lr_network, ligand=from, receptor=to)

  # By default, ligand_condition_specificty and receptor_condition_specificty are 0
  prioritizing_weights = c("de_ligand" = 1,
                           "de_receptor" = 1,
                           "activity_scaled" = 2,
                           "exprs_ligand" = 1,
                           "exprs_receptor" = 1,
                           "ligand_condition_specificity" = 0,
                           "receptor_condition_specificity" = 0)

  # Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
  celltypes <- unique(seuratObj$celltype)
  DE_table <- calculate_de(seuratObj, celltype_colname = "celltype",
                           condition_colname = "aggregate", condition_oi = condition_oi,
                           features = union(expressed_ligands, expressed_receptors))

  # Average expression information - only for LCMV condition
  expression_info <- get_exprs_avg(seuratObj, "celltype", condition_colname = "aggregate", condition_oi = condition_oi)

  # Calculate condition specificity - only for datasets with two conditions!
  condition_markers <- FindMarkers(object = seuratObj, ident.1 = condition_oi, ident.2 = condition_reference,
                                   group.by = "aggregate", min.pct = 0, logfc.threshold = 0,
                                   features = union(expressed_ligands, expressed_receptors))
  condition_markers <- rownames_to_column(condition_markers, "gene")

  # Combine DE of senders and receivers -> used for prioritization
  processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                            senders_oi = sender_celltypes, receivers_oi = receiver)

  processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)

  processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_renamed)

  ligand_activities$rank <- rank(rev(ligand_activities$aupr_corrected))

  prioritized_table <- generate_prioritization_tables(processed_expr_table,
                                                      processed_DE_table,
                                                      ligand_activities,
                                                      processed_condition_markers,
                                                      prioritizing_weights = prioritizing_weights)

  p_mushroom <- make_mushroom_plot(prioritized_table, top_n = 30)
}


steps_df <- microbenchmark(step1(seuratObj),
               step2(geneset_oi, background_expressed_genes, ligand_target_matrix, potential_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig),
               step3(active_ligand_target_links_df, ligand_target_matrix, ligand_receptor_links_df,
                     seuratObj, sender_celltypes, condition_oi, condition_reference, best_upstream_ligands, ligand_activities),
               step4(lr_network, seuratObj, condition_oi, expressed_ligands, expressed_receptors, sender_celltypes, receiver),
               times = 10, unit = "s")

# Separately
microbenchmark(step1(seuratObj), times = 10, unit = "s")
microbenchmark(step2(geneset_oi, background_expressed_genes, ligand_target_matrix, potential_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig), times = 10, unit = "s")
# microbenchmark(step3(active_ligand_target_links_df, ligand_target_matrix, ligand_receptor_links_df,
#                      seuratObj, sender_celltypes, condition_oi, condition_reference, best_upstream_ligands, ligand_activities), times = 10, unit = "s")
# microbenchmark(step4(lr_network, seuratObj, condition_oi, expressed_ligands, expressed_receptors, sender_celltypes, receiver), times = 10, unit = "s")


# Box 1
box1_df <- microbenchmark(nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj,
  receiver = "CD8 T",
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS",
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks), times = 10, unit = "s")

# Box 2
lr_network <- readRDS("~/Documents/nichenet/multinichenet_files/lr_network_mouse_21122021.rds")
sig_network <- readRDS("~/Documents/nichenet/multinichenet_files/signaling_network_human_21122021.rds")
gr_network <- readRDS("~/Documents/nichenet/multinichenet_files/gr_network_human_21122021.rds")

box2 <- function(lr_network, sig_network, gr_network, n_ligands = 10){
  # aggregate the individual data sources in a weighted manner to obtain a weighted integrated signaling network
  weighted_networks = construct_weighted_networks(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df = source_weights_df)


  # downweigh the importance of signaling and gene regulatory hubs - use the optimized parameters of this
  weighted_networks = apply_hub_corrections(weighted_networks = weighted_networks,
                                            lr_sig_hub = hyperparameter_list[hyperparameter_list$parameter == "lr_sig_hub",]$avg_weight,
                                            gr_hub = hyperparameter_list[hyperparameter_list$parameter == "gr_hub",]$avg_weight)

  ligands <- as.list(unique(lr_network$from)[1:n_ligands])
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks, ligands = ligands, algorithm = "PPR",
                                                        damping_factor = hyperparameter_list[hyperparameter_list$parameter == "damping_factor",]$avg_weight,
                                                        ltf_cutoff = hyperparameter_list[hyperparameter_list$parameter == "ltf_cutoff",]$avg_weight)


}

box2_df <- microbenchmark(box2(lr_network, sig_network, gr_network, n_ligands = 1287), times = 10, unit = "s")

timing_df <- rbind(steps_df, box1_df, box2_df)
saveRDS(timing_df, "timings_df.rds")
