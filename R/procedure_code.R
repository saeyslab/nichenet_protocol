seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))

zenodo_path <- "https://zenodo.org/record/7074291/files/"
ligand_target_matrix <- readRDS(url(paste0(zenodo_path, "ligand_target_matrix_nsga2r_final_mouse.rds")))
lr_network <- readRDS(url(paste0(zenodo_path, "lr_network_mouse_21122021.rds")))
weighted_networks <- readRDS(url(paste0(zenodo_path, "weighted_networks_nsga2r_final_mouse.rds")))

#### Procedure (steps are separated by an empty line) ####
library(nichenetr)
library(tidyverse)
library(Seurat)

## Feature extraction ##
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# For v5 objects
if (as.numeric(substr(packageVersion("Seurat"), 1, 1)) >= 5 &
    as.numeric(substr(seuratObj@version, 1, 1)) >= 5){
  count_matrix <- GetAssayData(seuratObj, assay = "RNA", layer = "count")
  data_matrix <- GetAssayData(seuratObj, assay = "RNA", layer = "data")
  
  print(paste("All genes from raw count and normalized data are the same:", all(rownames(count_matrix) == rownames(data_matrix))))
  
  # Convert gene names
  new_names <- convert_alias_to_symbols(rownames(count_matrix), "mouse", verbose=FALSE)
  
  # Check for duplicated gene names
  doubles <- names(which(table(new_names) > 1))
  # Set the duplicated names back to their old names
  genes_remove <- names(which(names(new_names[new_names %in% doubles]) != (new_names[new_names %in% doubles])))
  new_names[genes_remove] <- genes_remove
  
  # Change rownames
  rownames(count_matrix) <- new_names
  rownames(data_matrix) <- new_names
  
  seuratObj <- CreateSeuratObject(count_matrix,
                                  meta.data = seuratObj@meta.data)
  seuratObj[["RNA"]]$data <- data_matrix
  
  rm(doubles, genes_remove, new_names)
}

Idents(seuratObj) <- seuratObj$celltype

receiver <- "CD8 T"

expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network[lr_network$to %in% expressed_receptors, ]
potential_ligands <- unique(potential_ligands$from)

sender_celltypes <- c("CD4 T","Treg", "Mono", "NK", "B", "DC")
list_expressed_genes_sender <- lapply(sender_celltypes, function(celltype) {get_expressed_genes(celltype, seuratObj, pct = 0.05)})
expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

condition_oi <- "LCMV"
condition_reference <- "SS"

seurat_obj_receiver <- subset(seuratObj, idents = receiver)
DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate",
                                 min.pct = 0.05)

geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= 0.05 & abs(DE_table_receiver$avg_log2FC) >= 0.25, ]
geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]

## Ligand activity analysis and downstream prediction ##
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands)
ligand_activities <- ligand_activities[order(ligand_activities$aupr_corrected, decreasing = TRUE), ]

ligand_activities_all <- ligand_activities
ligand_activities <- ligand_activities[ligand_activities$test_ligand %in% potential_ligands_focused, ]

best_upstream_ligands <- top_n(ligand_activities, 30, aupr_corrected)$test_ligand

active_ligand_target_links_df <- lapply(best_upstream_ligands, get_weighted_ligand_target_links,
                                        geneset = geneset_oi,
                                        ligand_target_matrix = ligand_target_matrix,
                                        n = 200)
active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df))

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors,
                                                               lr_network, weighted_networks$lr_sig)

## Visualizations ##
# Ligand-target heatmap
active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links)))
order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
(make_heatmap_ggplot(vis_ligand_target,
                     y_name = "Prioritized ligands", x_name = "Predicted target genes",
                     color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple"))

# Ligand-receptor heatmap
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df,
                                                                     best_upstream_ligands, order_hclust = "receptors")
(make_heatmap_ggplot(t(vis_ligand_receptor_network[, order_ligands]),
                     "Ligands","Receptors",
                     color = "mediumvioletred",
                     legend_title = "Prior interaction potential"))

# Ligand activity heatmap
ligand_aupr_matrix <- column_to_rownames(ligand_activities, "test_ligand")
ligand_aupr_matrix <- ligand_aupr_matrix[rev(best_upstream_ligands), "aupr_corrected", drop=FALSE]
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)
(make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands", "Ligand activity",
                     legend_title = "AUPR", color = "darkorange") +
    theme(axis.text.x.top = element_blank()))

# LFC heatmap
celltype_order <- levels(Idents(seuratObj))
DE_table_top_ligands <- lapply(celltype_order[celltype_order %in% sender_celltypes],
                               get_lfc_celltype,
                               seurat_obj = seuratObj,
                               condition_colname = "aggregate",
                               condition_oi = condition_oi, condition_reference = condition_reference,
                               min.pct = 0, logfc.threshold = 0,
                               features = best_upstream_ligands, celltype_col = "celltype")
DE_table_top_ligands <- reduce(DE_table_top_ligands, full_join)
DE_table_top_ligands <- column_to_rownames(DE_table_top_ligands, "gene")

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop=FALSE])
(make_threecolor_heatmap_ggplot(vis_ligand_lfc, "Prioritized ligands","LFC in Sender",
                                low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",
                                legend_title = "LFC"))

# Dot plot
DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") +
  coord_flip() +
  scale_y_discrete(position = "right")

# Line plot to compare rankings between agnostic and focused approach
make_line_plot(ligand_activities = ligand_activities_all,
               potential_ligands = potential_ligands_focused)

# Chord diagram (ligand-target)
ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj, best_upstream_ligands[1:20],
                                                        celltype_col = "celltype")

active_ligand_target_links_df$target_type <- "LCMV-DE"
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df, cutoff = 0.40)

ligand_colors <- c("General" = "#377EB8", "NK" = "#4DAF4A", "B" = "#984EA3", "Mono" = "#FF7F00", "DC" = "#FFFF33", "Treg" = "#F781BF", "CD8 T"= "#E41A1C")
target_colors <- c( "LCMV-DE" = "#999999")
vis_circos_obj <- prepare_circos_visualization(circos_links, ligand_colors = ligand_colors, target_colors = target_colors)

make_circos_plot(vis_circos_obj, transparency = TRUE)

# Chord diagram (ligand-receptor)
lr_network_top_df <- rename(ligand_receptor_links_df, ligand=from, target=to)
lr_network_top_df$target_type = "LCMV_CD8T_receptor"
lr_network_top_df <- inner_join(lr_network_top_df, ligand_type_indication_df)
receptor_colors <- c("LCMV_CD8T_receptor" = "#E41A1C")
vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df, ligand_colors = ligand_colors, target_colors = receptor_colors)
make_circos_plot(vis_circos_receptor_obj, transparency = TRUE, link.visible = TRUE)

# Ligand to target signaling
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))
ligands_oi <- "Ebi3"
targets_oi <- c("Irf1", "Irf9")
active_signaling_network <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_oi,
                                                      targets_all = targets_oi, weighted_networks = weighted_networks,
                                                      top_n_regulators = 4, minmax_scaling = TRUE)

signaling_graph <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands_oi, targets_all = targets_oi,
                                                     sig_color = "indianred", gr_color = "steelblue")
DiagrammeR::render_graph(signaling_graph, layout = "tree")

## Prioritization ##
lr_network_filtered <- filter(lr_network, from %in% potential_ligands_focused &
                                to %in% expressed_receptors)[, c("from", "to")]

# A. Wrapper function
info_tables <- generate_info_tables(
  seuratObj,
  celltype_colname = "celltype",
  senders_oi = sender_celltypes,
  receivers_oi = receiver,
  lr_network = lr_network_filtered,
  condition_colname = "aggregate",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  scenario = "case_control")

processed_DE_table <- info_tables$sender_receiver_de
processed_expr_table <- info_tables$sender_receiver_info
processed_condition_markers <- info_tables$lr_condition_de

# B. Step-by-step
DE_table <- FindAllMarkers(subset(seuratObj, subset = aggregate == "LCMV"),
                           min.pct = 0, logfc.threshold = 0, return.thresh = 1,
                           features = unique(unlist(lr_network_filtered)))

expression_info <- get_exprs_avg(seuratObj, "celltype",
                                 condition_colname = "aggregate", condition_oi = condition_oi,
                                 features = unique(unlist(lr_network_filtered)))

condition_markers <- FindMarkers(object = seuratObj, ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate", min.pct = 0, logfc.threshold = 0,
                                 features = unique(unlist(lr_network_filtered)))
condition_markers <- rownames_to_column(condition_markers, "gene")

processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_filtered,
                                          senders_oi = sender_celltypes, receivers_oi = receiver)
processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_filtered)
processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_filtered)

# Optional custom weights
prioritizing_weights = c("de_ligand" = 1,
                         "de_receptor" = 1,
                         "activity_scaled" = 1,
                         "exprs_ligand" = 1,
                         "exprs_receptor" = 1,
                         "ligand_condition_specificity" = 1,
                         "receptor_condition_specificity" = 1)

prioritized_table <- generate_prioritization_tables(processed_expr_table,
                                                    processed_DE_table,
                                                    ligand_activities,
                                                    processed_condition_markers,
                                                    scenario = "case_control")

prioritized_table$sender <- factor(prioritized_table$sender, levels = celltype_order)
make_mushroom_plot(prioritized_table, top_n = 30,
                   show_all_datapoints = TRUE, show_rankings = TRUE)

#### Reproducing Figure 3 ####
top_n <- 10
legend_title_size <- 9
legend_text_size <- 7
legend_key_size <- unit(0.5, "cm")

p_ligand_target_network_subset <-
  make_heatmap_ggplot(vis_ligand_target %>% .[(nrow(.)-top_n+1):nrow(.),1:(top_n*2)],
                                                      "Prioritized ligands","Predicted target genes", color = "purple",
                                                      legend_title = "Regulatory potential",
                                                      legend_position = "bottom") +
  scale_fill_gradient(low = "whitesmoke",  high = "purple", n.breaks = 3)  +
  theme(legend.title = element_text(vjust = 0.75, size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.key.size = legend_key_size,
        axis.ticks = element_blank())

p_ligand_receptor_network_subset <- make_heatmap_ggplot(t(vis_ligand_receptor_network[, order_ligands]) %>% .[(nrow(.)-top_n+1):nrow(.),] %>% .[, colSums(.) > 0],
                                                        "Ligands", "Receptors", color = "mediumvioletred",
                                                        legend_title = "Prior interaction potential",
                                                        legend_position = "bottom") +
  theme(legend.title = element_text(vjust = 0.75, size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.key.size = legend_key_size,
        axis.ticks = element_blank())

p_ligand_aupr_subset <- make_heatmap_ggplot(vis_ligand_aupr %>% .[(nrow(.)-top_n+1):nrow(.),, drop=FALSE],
                                            "Prioritized ligands","Ligand activity",
                                            legend_position = "bottom", legend_title = "AUPR") +
  scale_x_discrete(expand = expansion(add=c(0,0)), position = "top") +
  scale_fill_gradient(n.breaks = 4, low="whitesmoke", high = "darkorange") +
  theme(legend.title = element_text(vjust = 0.75, size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.key.size = legend_key_size,
        legend.justification = c(0.7, 0),
        axis.title.x.top = element_text(margin=margin(0, 0, -30, 0)),
        axis.text.x.top = element_blank(),
        axis.ticks = element_blank())

p_ligand_lfc_subset <- make_threecolor_heatmap_ggplot(vis_ligand_lfc %>% .[(nrow(.)-top_n+1):nrow(.),], "Prioritized ligands","Sender cell types",
                                                      low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",
                                                      legend_position = "bottom", legend_title = "LFC") +
  theme(legend.title = element_text(vjust = 0.75, size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.key.size = legend_key_size,
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle=90))

p_dotplot_subset <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
                            features = rownames(vis_ligand_lfc %>% .[(nrow(.)-top_n+1):nrow(.),]),
                            cols = "RdYlBu") +
  coord_flip() +
  scale_y_discrete(position = "right") +
  theme(legend.text = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size, vjust = 0.8),
        legend.key.size = legend_key_size,
        legend.position = "bottom",
        legend.justification = c(0.5, 0),
        legend.box = "vertical",
        legend.box.margin = margin(10, 0, 0, 0),
        legend.margin = margin(-0.4,0,0,0, unit="cm"),
        axis.text.x = element_text(size=9, angle=90, hjust=0),
        axis.text.y = element_text(size=9),
        axis.ticks = element_blank(),
        axis.title = element_blank())

p_lineplot_subset <- make_line_plot(ligand_activities_all, potential_ligands_focused,
                                    ranking_range = c(1, 10),
                                    inset_scale = 0.75) +
  theme(legend.text = element_text(size = legend_text_size),
        legend.key.size = legend_key_size,
        plot.title = element_text(hjust = 0, size = legend_title_size, margin=margin(0,0,0,0)))

circos_links_subset <- get_ligand_target_links_oi(ligand_type_indication_df,
                                                  active_ligand_target_links_df,
                                                  cutoff = 0.95)
vis_circos_obj_subset <- prepare_circos_visualization(circos_links_subset, ligand_colors = ligand_colors[names(ligand_colors) %in% unique(circos_links_subset$ligand_type)],
                                                      target_colors = target_colors)

library(grid)
par(bg = "transparent")
# Create legend
circos_legend <- ComplexHeatmap::Legend(
  labels = c(unique(circos_links_subset$ligand_type), "Target"),
  background = c(ligand_colors[names(ligand_colors) %in% unique(circos_links_subset$ligand_type)], target_colors),
  type = "point",
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"),
  labels_gp = gpar(fontsize = 8)
  )


circos_legend_grob <- grid.grabExpr(ComplexHeatmap::draw(circos_legend))

make_circos_plot(vis_circos_obj_subset, transparency = TRUE, args.circos.text = list(cex = 0.7))
p_circos_no_legend <- recordPlot()
p_circos <- cowplot::plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.1))

# Signaling graph
active_signaling_network_subset <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_oi,
                                                             targets_all = targets_oi, weighted_networks = weighted_networks,
                                                             top_n_regulators = 1, minmax_scaling = TRUE)
signaling_graph_subset <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_subset,
                                                            ligands_all = ligands_oi, targets_all = targets_oi,
                                                            sig_color = "indianred", gr_color = "steelblue")
# It's going to be a woozy from now on
# Required packages: DiagrammeRsvg, rsvg, XML, and magick
# Export the graph to svg (won't work if you don't have the first two packages)
signaling_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(signaling_graph_subset, layout = "tree", output = "graph"))
# The width and height parameters of `render_graph` don't apply to SVGs, so change it manually with XML
svg_xml <- XML::xmlParse(signaling_svg)
root <- XML::xmlRoot(svg_xml)
svg_attrs <- XML::xmlAttrs(root)
svg_attrs[["width"]] <- paste0(as.numeric(str_extract(svg_attrs[["width"]], "[0-9]+"))*4, "pt") # Quadruple the width
svg_attrs[["height"]] <- paste0(as.numeric(str_extract(svg_attrs[["height"]], "[0-9]+"))*4, "pt") # Quadruple the height
XML::xmlAttrs(root) <- svg_attrs
# Now we read in the svg with magick, and turn it to a ggplot
signaling_magick <- magick::image_read_svg(XML::saveXML(svg_xml))
signaling_ggplot <- magick::image_ggplot(signaling_magick)

# Mushroom plot
p_mushroom_subset <- make_mushroom_plot(prioritized_table, top_n = 10,
                                        show_all_datapoints = TRUE, show_rankings = TRUE,
                                        legend.key.size = legend_key_size*0.5,
                                        #legend.justification = c(1, 0.5),
                                        legend.position = c(1.05, 0.85),
                                        #legend.key.height = unit(0.4, 'cm'),
                                        #legend.key.width = unit(0.3, 'cm'),
                                        #legend.box = "vertical",
                                        legend.spacing.x = unit(0, 'cm'),
                                        #legend.margin = margin(0,0,0,-10, unit="cm"),
                                        legend.title = element_text(size=legend_title_size),
                                        legend.text = element_text(size=legend_text_size),
                                        axis.text.x = element_text(angle=90, hjust=0),
                                        axis.title.x = element_text(hjust=0.2),
                                        axis.ticks = element_blank())

# Align left and right column first, before assembling the plots
library(cowplot)
left_column <- align_plots(p_ligand_target_network_subset,
                           p_ligand_lfc_subset, align = 'v', axis = 'l')

first_row <- plot_grid(left_column[[1]],
                       NULL,
                       p_ligand_receptor_network_subset,
                       NULL,
                       p_ligand_aupr_subset,
                       NULL,
                       ncol = 6, rel_widths = c(2.75, 0.5, 2, 0.25, 0.75, 0.25), align = "h", labels = c('a', '',  'b', '', 'c', ''))

second_row_1 <- align_plots(left_column[[2]], p_dotplot_subset, align = "h", axis = "tb")
second_row <- plot_grid(second_row_1[[1]],
                        NULL,
                        second_row_1[[2]],
                        NULL,
                        p_lineplot_subset,
                        labels = c('d','', 'e', '', 'f'),
                        ncol = 5, rel_widths = c(0.75, 0.2, 0.8, 0.15, 1.5), align = "h", axis = "t")


third_row <- cowplot::plot_grid(NULL,
                                ggdraw(p_circos),
                                NULL,
                                signaling_ggplot,
                                p_mushroom_subset,
                                labels = c('', 'g', '', 'h', 'i'),
                                nrow = 1, rel_widths = c(0.2, 1.2, 0.2, 1 , 2),
                                scale = c(1, 1.1, 1, 0.9, 1))

combined_plot <- cowplot::plot_grid(first_row,
                                    NULL,
                                    second_row,
                                    NULL,
                                    third_row,
                                    nrow = 5,
                                    rel_heights = c(3, 0.5, 3, 0.5, 2.5))

save_plot("~/Pictures/protocols_paper/combined_plot.pdf", combined_plot,
           base_width = 12, base_height = 12, dpi = 300)
