---
title: "Protocols code"
author: "Chananchida Sang-aram & Ruth Seurinck"
date: "2023-10-25"
output:
  html_document:
    keep_md: yes
---



# Materials
## Equipment setup

### Input data

The main function of NicheNet only needs a list of gene names for the 1) gene set of interest, 2) expressed ligands, and 3) expressed receptors. We provide code to obtain these genes, but as minimal input users need a gene expression matrix with cell type annotations. We do recommend users to provide a Seurat (at least v3) object for additional helper functions, but this is not strictly necessary.

### Hardware

Windows, Linux or Mac OS

### Software

R (RStudio desktop is recommended)

NicheNet can be installed by running the following command in R:


```r
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools") 

# devtools::install_github("saeyslab/nichenetr", ref ="devel") 
#Fix make sure it is on devel branch for some functionalities
```

### Example data

As example expression data of interacting cells, we will use mouse NICHE-seq data from Medaglia et al. (2017) to explore intercellular communication in the T cell area in the inguinal lymph node before and 72 hours after lymphocytic choriomeningitis virus (LCMV) infection. Specifically, we will prioritize which ligands can best explain the downstream changes after LCMV infection in CD8 T cells as the receiver population. This dataset contains 13,541 genes and 5,027 cells from 6 cell populations: CD4 T cells (including regulatory T cells), CD8 T cells, B cells, NK cells, dendritic cells (DCs) and inflammatory monocytes. The data is available on Zenodo (<https://zenodo.org/record/3531889>).

To download the file with the command line, commands like wget or curl can be used (Linux/macOS):

`wget https://zenodo.org/record/3531889/files/seuratObj.rds`

To download the file within the R session, entering the following commands:


```r
options(timeout = 3600) # increase time limit for downloading the data 
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds")) 
```

If the file is downloaded within the R session, it will have to be downloaded again once the session is restarted.

### Human/Mouse networks

Three networks are required to run the NicheNet analysis: the ligand-target prior model, the ligand-receptor network, and the weighted ligand-receptor network. These files are stored in Zenodo (<https://zenodo.org/record/7074291/>) and similar to above, can be downloaded locally or in the R session. We recommend users to have these files locally for convenience.

#### Downloading data locally


```r
ZENODO_PATH=https://zenodo.org/record/7074291/files

# Human ligand-target, ligand-receptor, and weighted networks
wget $ZENODO_PATH/ligand_target_matrix_nsga2r_final.rds
wget $ZENODO_PATH/lr_network_human_21122021.rds
wget $ZENODO_PATH/weighted_networks_nsga2r_final.rds

# Mouse ligand-target, ligand-receptor, and weighted networks
wget $ZENODO_PATH/ligand_target_matrix_nsga2r_final_mouse.rds
wget $ZENODO_PATH/lr_network_mouse_21122021.rds
wget $ZENODO_PATH/weighted_networks_nsga2r_final_mouse.rds
```

#### Download the file in the R session, mouse example


```r
zenodo_path <- "https://zenodo.org/record/7074291/files/" 

ligand_target_matrix <- readRDS(url(paste0(zenodo_path, "ligand_target_matrix_nsga2r_final_mouse.rds"))) 
lr_network <- readRDS(url(paste0(zenodo_path, "lr_network_mouse_21122021.rds"))) 
weighted_networks <- readRDS(url(paste0(zenodo_path, "weighted_networks_nsga2r_final_mouse.rds"))) 


ligand_target_matrix[1:5, 1:5] # target genes in rows, ligands in columns
```

```
##               2300002M23Rik 2610528A11Rik 9530003J23Rik            a
## 0610005C13Rik  0.000000e+00  0.000000e+00  1.311297e-05 0.000000e+00
## 0610009B22Rik  0.000000e+00  0.000000e+00  1.269301e-05 0.000000e+00
## 0610009L18Rik  8.872902e-05  4.977197e-05  2.581909e-04 7.570125e-05
## 0610010F05Rik  2.194046e-03  1.111556e-03  3.142374e-03 1.631658e-03
## 0610010K14Rik  2.271606e-03  9.360769e-04  3.546140e-03 1.697713e-03
##                        A2m
## 0610005C13Rik 1.390053e-05
## 0610009B22Rik 1.345536e-05
## 0610009L18Rik 9.802264e-05
## 0610010F05Rik 2.585820e-03
## 0610010K14Rik 2.632082e-03
```

```r
head(lr_network)
```

```
##            from    to database   source
## 1 2300002M23Rik  Ddr1 omnipath omnipath
## 2 2610528A11Rik Gpr15 omnipath omnipath
## 3 9530003J23Rik Itgal omnipath omnipath
## 4             a  Atrn omnipath omnipath
## 5             a  F11r omnipath omnipath
## 6             a  Mc1r omnipath omnipath
```

```r
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
```

```
##            from     to     weight
## 1 0610010F05Rik    App 0.10989552
## 2 0610010F05Rik    Cat 0.06732398
## 3 0610010F05Rik   H1f2 0.06601048
## 4 0610010F05Rik Lrrc49 0.08288421
## 5 0610010F05Rik  Nicn1 0.08639550
## 6 0610010F05Rik  Srpk1 0.12290087
```

The ligand-target prior model is a matrix describing the potential that a ligand may regulate a target gene, and it is used to run the ligand activity analysis. The ligand-receptor network contains information on potential ligand-receptor bindings, and it is used to identify potential ligands. Finally, the weighted ligand-receptor network contains weights representing the potential that a ligand will bind to a receptor, and it is used for visualization.

These networks were translated from human to mouse gene names using one-to-one orthologs when feasible, and one-to-many conversion was allowed when necessary (for instance, when one human gene symbol corresponded to two mouse gene symbols). Users that are interested in building prior models for other organisms can either create an organism-specific model using data sources relevant to that organism, or use the existing human NicheNet model to convert human gene symbols to their corresponding one-to-one orthologs in the organism of interest. However, this decision depends on one hand, the availability of data for the organism of interest and on the other, the homology between humans and the organism of interest. For instance, using the human model and converting gene symbols might work for primates, but creating a new model from species-specific data sources is better suited for organisms like Drosophila.

# Procedure

Here, we describe the procedure for both the sender-focused and sender-agnostic approach, as shown here:


![](figure2.svg){width=75%}

As two conditions are present in this example dataset, the gene set of interest is chosen as the DE genes between these conditions in the receiver cell type.

## Feature extraction

1. Load required libraries.


```r
library(nichenetr) 
# devtools::load_all("~/nichenetr/")
library(Seurat) 
library(tidyverse) 
```

2. *(Optional)* For older Seurat objects, update it to be compatible with the currently installed Seurat version. For expression data with older gene symbols, convert them to more recent gene symbols. 


```r
seuratObj <- UpdateSeuratObject((seuratObj))
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

seuratObj
```

```
## An object of class Seurat 
## 13541 features across 5027 samples within 1 assay 
## Active assay: RNA (13541 features, 1575 variable features)
##  3 layers present: counts, data, scale.data
##  4 dimensional reductions calculated: cca, cca.aligned, tsne, pca
```

3. Set the cell type annotation column as the identity of the Seurat object.


```r
Idents(seuratObj) <- seuratObj$celltype 
```

4. Define a "receiver" cell population. The receiver cell population can only consist of one cell type.


```r
receiver <- "CD8 T" 
```

5. Determine which genes are expressed in the receiver cell population. By default, get_expressed_genes considers genes to be expressed if they have non-zero counts in at least 10% of the cell population (pct argument). Users are also free to define expressed genes differently in a way that fits their data. Here, we have lowered the pct parameter to 5% as some of the ligands and receptors are very lowly expressed.


```r
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj,  pct = 0.05)

# Preview
length(expressed_genes_receiver)
```

```
## [1] 3903
```

```r
head(expressed_genes_receiver)
```

```
## [1] "Atraid"        "Tmem248"       "Tsr3"          "Sf3b6"        
## [5] "0610010K14Rik" "0610012G03Rik"
```


6. Get a list of all receptors available in the ligand-receptor network, and define expressed receptors as genes that are in the ligand-receptor network and expressed in the receiver.


```r
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver) 

# Preview
length(expressed_receptors)
```

```
## [1] 113
```

```r
head(expressed_receptors)
```

```
## [1] "Itgal"   "Notch1"  "Tspan14" "Itga4"   "Itgb1"   "Il6ra"
```

7. Define the potential ligands as all ligands whose cognate receptors are expressed.


```r
potential_ligands <- lr_network[lr_network$to %in% expressed_receptors, ] 
potential_ligands <- unique(potential_ligands$from) 

# Preview
length(potential_ligands)
```

```
## [1] 483
```

```r
head(potential_ligands)
```

```
## [1] "9530003J23Rik" "Adam10"        "Adam11"        "Adam12"       
## [5] "Adam15"        "Adam17"
```

8. *(Optional)* For the sender-focused approach, define sender cell types and expressed genes in all populations combined. Then, filter potential ligands to those that are expressed in sender cells.


```r
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC") 
list_expressed_genes_sender <- lapply(sender_celltypes, function(celltype) {
    get_expressed_genes(celltype, seuratObj, pct = 0.05)
  }) 
expressed_genes_sender <- unique(unlist(list_expressed_genes_sender)) 
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Preview
length(expressed_genes_sender)
```

```
## [1] 8492
```

```r
head(expressed_genes_sender)
```

```
## [1] "Atraid"        "Tmem248"       "Sf3b6"         "0610010K14Rik"
## [5] "0610012G03Rik" "Sptssa"
```

```r
length(potential_ligands_focused)
```

```
## [1] 127
```

```r
head(potential_ligands_focused)
```

```
## [1] "Adam10" "Adam15" "Adam17" "Adam9"  "Adgre5" "Alcam"
```

9. Define the reference condition and condition of interest. The condition of interest is the condition after the CCC event has taken place, or the 'case' group in case-control studies. Here, it represents the condition after LCMV infection.



```r
condition_oi <- "LCMV" 
condition_reference <- "SS" 
```

10. Define the gene set of interest that represents the cell-cell communication event to be studied. First, create a new Seurat object that only contains the receiver cell type. Then, perform DE analysis between the treatment conditions within the receiver cell type. Finally, define the gene set of interest as significantly DE genes, i.e., genes with adjusted p-value lower than or equal to 0.05 and absolute log-fold change greater than 0.25.

By default, both genes that are up and downregulated are considered. Users can choose to focus on only one direction (typically upregulation) by removing the `abs()` function and adjusting the equality term to either \>= 0.25 or \<= -0.25 for up and downregulation, respectively. We recommend the gene set of interest to contain between 20 and 2000 genes for optimal ligand activity prediction. Moreover, the number of background genes should be sufficiently greater than those of the gene set of interest.


```r
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,  
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate",
                                 min.pct = 0.05) 

geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= 0.05 & abs(DE_table_receiver$avg_log2FC) >= 0.25, ] 
geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in% rownames(ligand_target_matrix)] 

# Preview
length(geneset_oi)
```

```
## [1] 260
```

```r
head(geneset_oi)
```

```
## [1] "Ifi27l2b" "Irf7"     "Ly6a"     "Stat1"    "Ly6c2"    "Ifit3"
```

11. Determine background genes as all the genes expressed in the receiver cell type that are also in the ligand-target matrix.


```r
background_expressed_genes <- expressed_genes_receiver[ 
expressed_genes_receiver %in% rownames(ligand_target_matrix)] 

# Preview
length(background_expressed_genes)
```

```
## [1] 3476
```

```r
head(background_expressed_genes)
```

```
## [1] "Atraid"        "Tmem248"       "Tsr3"          "Sf3b6"        
## [5] "0610010K14Rik" "0610012G03Rik"
```

## Ligand activity analysis and downstream prediction

12. Perform the ligand activity analysis, then sort the ligands based on the area under the precision-recall curve (AUPR).


```r
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands) 

ligand_activities <- ligand_activities[order(ligand_activities$aupr_corrected, 	decreasing = TRUE), ] 

# Preview
dim(ligand_activities)
```

```
## [1] 483   5
```

```r
head(ligand_activities)
```

```
## # A tibble: 6 × 5
##   test_ligand auroc  aupr aupr_corrected pearson
##   <chr>       <dbl> <dbl>          <dbl>   <dbl>
## 1 Ifna1       0.714 0.433          0.358   0.498
## 2 Ifnb1       0.711 0.401          0.327   0.433
## 3 Ifnl3       0.683 0.392          0.317   0.433
## 4 Il27        0.682 0.391          0.316   0.445
## 5 Ifng        0.732 0.382          0.307   0.451
## 6 Ifnk        0.671 0.282          0.207   0.272
```

13. *(Optional)* If performing the sender-focused approach, subset the ligand activities to only contain expressed ligands.
**Note:** When using the sender-agnostic approach, simply replace ligand_activities with ligand_activities_all in Steps 14 and 20.


```r
ligand_activities_all <- ligand_activities 
ligand_activities <- ligand_activities[ligand_activities$test_ligand %in% potential_ligands_focused, ] 

# Preview
dim(ligand_activities)
```

```
## [1] 127   5
```

```r
head(ligand_activities)
```

```
## # A tibble: 6 × 5
##   test_ligand auroc  aupr aupr_corrected pearson
##   <chr>       <dbl> <dbl>          <dbl>   <dbl>
## 1 Il27        0.682 0.391          0.316   0.445
## 2 Ebi3        0.666 0.264          0.189   0.256
## 3 Tnf         0.671 0.205          0.131   0.249
## 4 Ptprc       0.660 0.198          0.124   0.168
## 5 H2-Eb1      0.656 0.195          0.120   0.182
## 6 Vsig10      0.649 0.194          0.119   0.170
```

14. Obtain the names of the top 30 ligands.


```r
best_upstream_ligands <- top_n(ligand_activities, 30, aupr_corrected)$test_ligand 
 
# Preview
length(best_upstream_ligands)
```

```
## [1] 30
```

```r
head(best_upstream_ligands)
```

```
## [1] "Il27"   "Ebi3"   "Tnf"    "Ptprc"  "H2-Eb1" "Vsig10"
```

15. Infer which genes in the gene set of interest have the highest regulatory potential for each top-ranked ligand. The function get_weighted_ligand_target_links will return genes that are in the gene set of interest and are the top `n` targets of a ligand (default: `n = 200`).


```r
active_ligand_target_links_df <- lapply(best_upstream_ligands,
                                        get_weighted_ligand_target_links, 
                                        geneset = geneset_oi, 
                                        ligand_target_matrix = ligand_target_matrix, 
                                        n = 200) 

active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df)) 

# Preview
dim(active_ligand_target_links_df)
```

```
## [1] 656   3
```

```r
head(active_ligand_target_links_df)
```

```
## # A tibble: 6 × 3
##   ligand target weight
##   <chr>  <chr>   <dbl>
## 1 Il27   Adar    0.163
## 2 Il27   B2m     0.170
## 3 Il27   Bst2    0.111
## 4 Il27   Calhm6  0.129
## 5 Il27   Cd274   0.111
## 6 Il27   Cxcl10  0.178
```

16. Similarly, identify which receptors have the highest interaction potential with the top-ranked ligands.


```r
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

# Preview
dim(ligand_receptor_links_df)
```

```
## [1] 54  3
```

```r
head(ligand_receptor_links_df)
```

```
## # A tibble: 6 × 3
##   from   to     weight
##   <chr>  <chr>   <dbl>
## 1 Adam17 Il6ra   0.447
## 2 Adam17 Itgb1   0.454
## 3 Adam17 Notch1  1.05 
## 4 App    Cd74    0.670
## 5 App    Sorl1   0.922
## 6 Ccl22  Ccr7    0.679
```

## Visualizations

Visualizations covered in this section include four heatmaps (Steps 17-22), a dot plot of cell type expression and percentage (Step 23), a line plot comparing ligand rankings between two approaches (Step 24), a chord diagram (Steps 25-29), and a signaling graph (Steps 30-31). ▲ CRITICAL Heatmaps depict the ligand-target regulatory potential (Steps 17-18), ligand-receptor interaction potential (Step 19), ligand activity (Step 20), and log-fold change of ligands between treatment conditions (Steps 21-22).

17. Prepare the weighted ligand-target data frame for visualization by transforming it into matrix. By default, regulatory potentials lower than the 25th percentile are set to zero for visualization clarity. This cutoff parameter can freely be tuned by the user.


```r
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25) 
```

18. Order the rows to follow the rankings of the ligands, and the columns alphabetically (Figure 3A).


```r
order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links))) 
order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links)) 

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

(make_heatmap_ggplot(vis_ligand_target, y_name = "Prioritized ligands", x_name = "Predicted target genes",
                     color = "purple", legend_title = "Regulatory potential") + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")) 
```

```
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

![](procedure_code_files/figure-html/visualizations-II-1.png)<!-- -->

19. Create a heatmap for ligand-receptor interactions (Figure 3B).


```r
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df, best_upstream_ligands,
  order_hclust = "receptors") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential")) 
```

![](procedure_code_files/figure-html/visualizations-III-1.png)<!-- -->

20. Create a heatmap of the ligand activity measure (Figure 3C).


```r
ligand_aupr_matrix <- column_to_rownames(ligand_activities, "test_ligand") 
ligand_aupr_matrix <- ligand_aupr_matrix[rev(best_upstream_ligands), "aupr_corrected", drop=FALSE] 
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  
```

![](procedure_code_files/figure-html/visualizations-IV-1.png)<!-- -->

21. For each cell type, compute the log-fold change of the top-ranked ligands between treatment conditions.


```r
celltype_order <- levels(Idents(seuratObj)) 

DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "aggregate",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- reduce(DE_table_top_ligands, full_join) 
```

```
## Joining with `by = join_by(gene)`
## Joining with `by = join_by(gene)`
## Joining with `by = join_by(gene)`
## Joining with `by = join_by(gene)`
## Joining with `by = join_by(gene)`
```

```r
DE_table_top_ligands <- column_to_rownames(DE_table_top_ligands, "gene") 
```

22. Create the heatmap (Figure 3D).


```r
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ]) 

(make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                "Prioritized ligands", "LFC in Sender",
                                low_color = "midnightblue", mid_color = "white",
                                mid = median(vis_ligand_lfc), high_color = "red",
                                legend_title = "LFC")) 
```

![](procedure_code_files/figure-html/visualizations-VI-1.png)<!-- -->

23. Create a dot plot showing the average expression of ligands per cell type, as well as the percentage of cells from the cell type expressing the ligands (Figure 3E).


```r
DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
        features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right") 
```

![](procedure_code_files/figure-html/visualizations-VII-1.png)<!-- -->

24. *(Optional)* Create a line plot comparing the rankings between the sender-agnostic and sender-focused approach (Figure 3F).


```r
(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
   theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))
```

![](procedure_code_files/figure-html/visualizations-VIII-1.png)<!-- -->

25. To create a ligand-target chord diagram, assign each ligand to a specific cell type. A ligand is only assigned to a cell type if that cell type is the only one to show an average expression of that ligand that is higher than the mean + one standard deviation across all cell types. Otherwise, it is assigned to "General".


```r
ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj, best_upstream_ligands[1:20], celltype_col = "celltype") 

# Preview
dim(ligand_type_indication_df)
```

```
## [1] 20  2
```

```r
head(ligand_type_indication_df)
```

```
##   ligand_type ligand
## 1       CD8 T  Clcf1
## 2        Treg    Tnf
## 3           B H2-Eb1
## 4           B  H2-M3
## 5           B H2-T24
## 6           B  H2-Oa
```


26. Using the weighted ligand-target data frame from Step 15, group target genes and filter out the lowest 40% of the regulatory potentials. In this case, there is only one grouping of target genes (DE genes after LCMV infection), but users can define multiple target gene groups if applicable. In case the resulting chord diagram is still overcrowded, users may adjust the `cutoff` parameter to filter out even more ligand-target links.


```r
active_ligand_target_links_df$target_type <- "LCMV-DE" 
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df, cutoff = 0.40) 
```

```
## Joining with `by = join_by(ligand)`
```

```r
# Preview
dim(circos_links)
```

```
## [1] 485   5
```

```r
head(circos_links)
```

```
## # A tibble: 6 × 5
##   ligand target weight target_type ligand_type
##   <chr>  <chr>   <dbl> <chr>       <chr>      
## 1 Il27   Adar    0.163 LCMV-DE     Mono       
## 2 Il27   B2m     0.170 LCMV-DE     Mono       
## 3 Il27   Bst2    0.111 LCMV-DE     Mono       
## 4 Il27   Calhm6  0.129 LCMV-DE     Mono       
## 5 Il27   Cd274   0.111 LCMV-DE     Mono       
## 6 Il27   Cxcl10  0.178 LCMV-DE     Mono
```

27. Assign colors to cell types and target gene groups. Then, prepare the data frame for visualization: the function assigns colors to ligands and targets and calculates gaps between sectors of the chord diagram.


```r
ligand_colors <- c("General" = "#377EB8", "NK" = "#4DAF4A", "B" = "#984EA3",
                   "Mono" = "#FF7F00", "DC" = "#FFFF33", "Treg" = "#F781BF",
                   "CD8 T"= "#E41A1C") 
target_colors <- c("LCMV-DE" = "#999999") 

vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors) 
```

```
## Joining with `by = join_by(ligand_type)`
## Joining with `by = join_by(target_type)`
```

28. Draw the chord diagram (Figure 3G).


```r
make_circos_plot(vis_circos_obj, transparency = FALSE,  args.circos.text = list(cex = 0.5)) 
```

![](procedure_code_files/figure-html/visualizations-XII-1.png)<!-- -->

29. To create a ligand-receptor chord diagram, perform Steps 26-28 using the weighted ligand-receptor data frame from Step 16. As `prepare_circos_visualization` accesses "target" and "target_type" columns, it is necessary to rename the columns accordingly even though the data frame contains receptor and not target gene information. When drawing the plot, the argument `link.visible` = TRUE is also necessary for making all links visible, since no cutoff is used to filter out ligand-receptor interactions.


```r
lr_network_top_df <- rename(ligand_receptor_links_df, ligand=from, target=to) 
lr_network_top_df$target_type = "LCMV_CD8T_receptor" 
lr_network_top_df <- inner_join(lr_network_top_df, ligand_type_indication_df) 
```

```
## Joining with `by = join_by(ligand)`
```

```r
receptor_colors <- c("LCMV_CD8T_receptor" = "#E41A1C") 

vis_circos_receptor_obj <- prepare_circos_visualization(lr_network_top_df,
                                                        ligand_colors = ligand_colors,
                                                        target_colors = receptor_colors) 
```

```
## Joining with `by = join_by(ligand_type)`
## Joining with `by = join_by(target_type)`
```

```r
make_circos_plot(vis_circos_receptor_obj, transparency = TRUE,
                 link.visible = TRUE,  args.circos.text = list(cex = 0.8)) 
```

![](procedure_code_files/figure-html/visualizations-XIII-1.png)<!-- -->

30. To create a signaling graph, first download the ligand-transcription factor matrix. Then, extract the most highly weighted paths from the ligand to the target genes of interest. The number of regulators that are extracted can be adjusted using `top_n_regulators`. By setting `minmax_scaling = TRUE`, we perform min-max scaling to make the weights between the signaling and gene regulatory network more comparable. Additionally, it is possible to check which data sources support the inferred pathway by using the function `infer_supporting_datasources`. This would require separate signaling and gene regulatory networks as input (see Box 3 for code to download these networks).


```r
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds")) 
ligands_oi <- "Ebi3" 
targets_oi <- c("Irf1", "Irf9") 

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 4, minmax_scaling = TRUE) 
```

31. Convert the data frames into a DiagrammeR object, and render the signaling graph (Figure 3H).


```r
signaling_graph <- diagrammer_format_signaling_graph(
  signaling_graph_list = active_signaling_network,
  ligands_all = ligands_oi, targets_all = targets_oi,
  sig_color = "indianred", gr_color = "steelblue") 

DiagrammeR::render_graph(signaling_graph, layout = "tree") 
```

```
## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
## `.name_repair` is omitted as of tibble 2.0.0.
## ℹ Using compatibility `.name_repair`.
## ℹ The deprecated feature was likely used in the DiagrammeR package.
##   Please report the issue at
##   <https://github.com/rich-iannone/DiagrammeR/issues>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```{=html}
<div class="grViz html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-0dd1ba32e0541469b5ed" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-0dd1ba32e0541469b5ed">{"x":{"diagram":"digraph {\n\ngraph [layout = \"neato\",\n       outputorder = \"edgesfirst\",\n       bgcolor = \"white\"]\n\nnode [fontname = \"Helvetica\",\n      fontsize = \"10\",\n      shape = \"circle\",\n      fixedsize = \"true\",\n      width = \"0.5\",\n      style = \"filled\",\n      fillcolor = \"aliceblue\",\n      color = \"gray70\",\n      fontcolor = \"gray50\"]\n\nedge [fontname = \"Helvetica\",\n     fontsize = \"8\",\n     len = \"1.5\",\n     color = \"gray80\",\n     arrowsize = \"0.5\"]\n\n  \"1\" [label = \"Ebi3\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#CD5C5C\", pos = \"7.5,7!\"] \n  \"2\" [label = \"Stat1\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"7.5,6!\"] \n  \"3\" [label = \"Stat2\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"3.5,3!\"] \n  \"4\" [label = \"Stat3\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"8,5!\"] \n  \"5\" [label = \"Stat4\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"8.5,4!\"] \n  \"6\" [label = \"Stat5a\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#7F7F7F\", pos = \"8,2!\"] \n  \"7\" [label = \"Irf1\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#4682B4\", pos = \"10.5,1!\"] \n  \"8\" [label = \"Irf9\", style = \"filled\", width = \"0.75\", fontcolor = \"white\", fillcolor = \"#4682B4\", pos = \"6,1!\"] \n\"1\"->\"2\" [penwidth = \"1.25977306006657\", color = \"indianred\"] \n\"1\"->\"3\" [penwidth = \"1.13167686673422\", color = \"indianred\"] \n\"1\"->\"4\" [penwidth = \"1.23806443285084\", color = \"indianred\"] \n\"1\"->\"5\" [penwidth = \"1.43281572060901\", color = \"indianred\"] \n\"1\"->\"6\" [penwidth = \"1.14009611891571\", color = \"indianred\"] \n\"2\"->\"3\" [penwidth = \"1.40533373915141\", color = \"indianred\"] \n\"2\"->\"4\" [penwidth = \"1.54144367250614\", color = \"indianred\"] \n\"2\"->\"5\" [penwidth = \"0.987020557480735\", color = \"indianred\"] \n\"2\"->\"6\" [penwidth = \"1.34717911930313\", color = \"indianred\"] \n\"3\"->\"2\" [penwidth = \"1.75\", color = \"indianred\"] \n\"3\"->\"6\" [penwidth = \"0.752397933486639\", color = \"indianred\"] \n\"4\"->\"2\" [penwidth = \"1.43118032085197\", color = \"indianred\"] \n\"4\"->\"3\" [penwidth = \"1.27328859647951\", color = \"indianred\"] \n\"4\"->\"5\" [penwidth = \"1.02801783517307\", color = \"indianred\"] \n\"4\"->\"6\" [penwidth = \"1.27932306138063\", color = \"indianred\"] \n\"5\"->\"2\" [penwidth = \"0.947717727085479\", color = \"indianred\"] \n\"5\"->\"3\" [penwidth = \"0.953340730930064\", color = \"indianred\"] \n\"5\"->\"4\" [penwidth = \"0.971892004543157\", color = \"indianred\"] \n\"5\"->\"6\" [penwidth = \"0.904230698245332\", color = \"indianred\"] \n\"6\"->\"3\" [penwidth = \"0.75\", color = \"indianred\"] \n\"6\"->\"4\" [penwidth = \"1.22412329170677\", color = \"indianred\"] \n\"6\"->\"5\" [penwidth = \"0.921665045908184\", color = \"indianred\"] \n\"2\"->\"7\" [penwidth = \"1.54631395817282\", color = \"steelblue\"] \n\"2\"->\"8\" [penwidth = \"1.34428158131059\", color = \"steelblue\"] \n\"3\"->\"7\" [penwidth = \"1.21853545540254\", color = \"steelblue\"] \n\"3\"->\"8\" [penwidth = \"1.00841001829415\", color = \"steelblue\"] \n\"4\"->\"7\" [penwidth = \"1.75\", color = \"steelblue\"] \n\"4\"->\"8\" [penwidth = \"1.32165980574267\", color = \"steelblue\"] \n\"5\"->\"7\" [penwidth = \"1.11251821968624\", color = \"steelblue\"] \n\"5\"->\"8\" [penwidth = \"0.75\", color = \"steelblue\"] \n\"6\"->\"7\" [penwidth = \"1.0743405567136\", color = \"steelblue\"] \n\"6\"->\"8\" [penwidth = \"0.909277567466469\", color = \"steelblue\"] \n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>
```

## Prioritization of ligand-receptor pairs

32. Filter the ligand-receptor network to only contain expressed interactions. 


```r
lr_network_filtered <- filter(lr_network, from %in% potential_ligands_focused &
                                to %in% expressed_receptors)[, c("from", "to")]
```

33. Calculate the values required for prioritization, including DE between cell types, average expression, and DE between conditions. To use the wrapper function, follow option A. To calculate these values step-by-step, follow option B.


A. **Using the wrapper function.**

(i) Run the wrapper function. 


```r
info_tables <- generate_info_tables( 
    seuratObj, 
    celltype_colname = "celltype", 
    senders_oi = sender_celltypes, 
    receivers_oi = receiver, 
    lr_network = lr_network_filtered, 
    condition_colname = "aggregate", 
    condition_oi = condition_oi, 
    condition_reference = condition_reference, 
    scenario = "case_control"
) 
```

```
## condition_* is given. Only cells from that condition will be considered in cell type specificity calculation.
```

```
## Calculating cluster CD8 T
```

```
## Calculating cluster CD4 T
```

```
## Calculating cluster B
```

```
## Calculating cluster Treg
```

```
## Calculating cluster NK
```

```
## Calculating cluster Mono
```

```
## Calculating cluster DC
```

(ii) Assign the output of the wrapper function to variables. 


```r
processed_DE_table <- info_tables$sender_receiver_de  
processed_expr_table <- info_tables$sender_receiver_info  
processed_condition_markers <- info_tables$lr_condition_de

# Preview
dim(processed_DE_table)
```

```
## [1] 1272   13
```

```r
head(processed_DE_table)
```

```
##   sender receiver ligand receptor lfc_ligand lfc_receptor
## 1     DC    CD8 T   Ccl5    Cxcr3   6.432043   0.16714791
## 2   Mono    CD8 T   Lyz2    Itgal   5.493265  -0.01687003
## 3     DC    CD8 T  H2-M2     Cd8a   3.416479   1.94059972
## 4     DC    CD8 T Cxcl16    Cxcr6   4.182085   0.54826454
## 5   Mono    CD8 T  Cxcl9    Cxcr3   4.328801   0.16714791
## 6   Mono    CD8 T  Cxcl9     Dpp4   4.328801   0.16416445
##   ligand_receptor_lfc_avg  p_val_ligand  p_adj_ligand p_val_receptor
## 1                3.299595  1.893317e-25  2.563740e-21   7.758812e-05
## 2                2.738198 1.728697e-160 2.340828e-156   4.973381e-02
## 3                2.678539 1.017174e-272 1.377355e-268  5.250531e-206
## 4                2.365175 1.138617e-243 1.541801e-239   5.987787e-21
## 5                2.247975 3.834954e-124 5.192911e-120   7.758812e-05
## 6                2.246483 3.834954e-124 5.192911e-120   6.628900e-04
##   p_adj_receptor pct_expressed_sender pct_expressed_receiver
## 1   1.000000e+00                1.000                  0.042
## 2   1.000000e+00                0.933                  0.188
## 3  7.109745e-202                0.429                  0.659
## 4   8.108063e-17                0.929                  0.089
## 5   1.000000e+00                0.547                  0.042
## 6   1.000000e+00                0.547                  0.148
```

```r
dim(processed_expr_table)
```

```
## [1] 10535     7
```

```r
head(processed_expr_table)
```

```
## # A tibble: 6 × 7
##   sender receiver ligand receptor avg_ligand avg_receptor ligand_receptor_prod
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                <dbl>
## 1 DC     Mono     B2m    Tap1           216.         8.59                1856.
## 2 DC     NK       B2m    Klrd1          216.         7.43                1607.
## 3 DC     B        B2m    Tap1           216.         7.35                1588.
## 4 DC     Treg     B2m    Tap1           216.         7.18                1552.
## 5 Mono   Mono     B2m    Tap1           158.         8.59                1353.
## 6 DC     DC       B2m    Tap1           216.         5.91                1277.
```

```r
dim(processed_condition_markers)
```

```
## [1] 215   9
```

```r
head(processed_condition_markers)
```

```
##   ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_avg  p_val_ligand
## 1 H2-Ab1      Cd4  2.4021254   0.11569357               1.2589095  4.424390e-06
## 2 Cxcl10     Dpp4  1.6066163   0.35175421               0.9791853  6.700636e-29
## 3    B2m     Tap1  0.7071427   1.13931050               0.9232266 6.936359e-174
## 4 H2-T22    Klrd1  1.5223370  -0.05659737               0.7328698 1.006291e-111
## 5 H2-T23    Klrd1  1.4651999  -0.05659737               0.7043013 1.789643e-114
## 6 Cxcl10    Cxcr3  1.6066163  -0.25400642               0.6763049  6.700636e-29
##    p_adj_ligand p_val_receptor p_adj_receptor
## 1  5.991066e-02   5.634068e-02   1.000000e+00
## 2  9.073332e-25   1.170731e-06   1.585287e-02
## 3 9.392524e-170   3.585450e-52   4.855057e-48
## 4 1.362618e-107   6.202530e-01   1.000000e+00
## 5 2.423356e-110   6.202530e-01   1.000000e+00
## 6  9.073332e-25   1.918372e-06   2.597667e-02
```
 

(B) **Running step-by-step calculations.**

(i) Calculate DE between cell types within the condition of interest.


```r
DE_table <- FindAllMarkers(subset(seuratObj, subset = aggregate == "LCMV"),
                           min.pct = 0, logfc.threshold = 0, return.thresh = 1,
                           features = unique(unlist(lr_network_filtered))) 
```

```
## Calculating cluster CD8 T
```

```
## Calculating cluster CD4 T
```

```
## Calculating cluster Treg
```

```
## Calculating cluster B
```

```
## Calculating cluster NK
```

```
## Calculating cluster Mono
```

```
## Calculating cluster DC
```

(ii) Calculate average expression of each gene per cell type.


```r
expression_info <- get_exprs_avg(seuratObj, "celltype",
                                 condition_colname = "aggregate",
                                 condition_oi = condition_oi,
                                 features = unique(unlist(lr_network_filtered)))
```

(iii) Calculate DE between conditions.


```r
condition_markers <- FindMarkers(object = seuratObj,
                                 ident.1 = condition_oi, ident.2 = condition_reference,
                                 group.by = "aggregate",
                                 min.pct = 0, logfc.threshold = 0,
                                 features = unique(unlist(lr_network_filtered)))

condition_markers <- rownames_to_column(condition_markers, "gene") 
```

(iv) Process the data frames from Steps (i)-(iii) to follow the same format.


```r
processed_DE_table <- process_table_to_ic(
  DE_table,
  table_type = "celltype_DE",
  lr_network_filtered,
  senders_oi = sender_celltypes,
  receivers_oi = receiver) 

processed_expr_table <- process_table_to_ic(
  expression_info,
  table_type = "expression",
  lr_network_filtered) 

processed_condition_markers <- process_table_to_ic( 
  condition_markers,
  table_type = "group_DE",
  lr_network_filtered) 
```


34. Generate the prioritization table containing rankings of cell-type-specific, ligand-receptor interactions. The "case_control" scenario sets all weights to one, while the "one_condition" scenario sets condition specificity to zero and the remaining weights to one.


```r
prioritized_table <- generate_prioritization_tables(
  processed_expr_table,
  processed_DE_table,
  ligand_activities,
  processed_condition_markers,
  scenario = "case_control") 
```

```
## Joining with `by = join_by(sender, receiver, ligand, receptor)`
## Joining with `by = join_by(sender, ligand, lfc_ligand, p_val_ligand)`
## Joining with `by = join_by(ligand)`
## Joining with `by = join_by(receiver, receptor, lfc_receptor, p_val_receptor)`
## Joining with `by = join_by(sender, ligand, avg_ligand)`
## Joining with `by = join_by(receiver, receptor, avg_receptor)`
## Joining with `by = join_by(ligand)`
## Joining with `by = join_by(receptor)`
```

```r
# Preview
dim(prioritized_table)
```

```
## [1] 1272   51
```

```r
head(prioritized_table)
```

```
## # A tibble: 6 × 51
##   sender receiver ligand receptor lfc_ligand lfc_receptor ligand_receptor_lfc_…¹
##   <chr>  <chr>    <chr>  <chr>         <dbl>        <dbl>                  <dbl>
## 1 NK     CD8 T    Ptprc  Dpp4          0.596        0.164                  0.380
## 2 Mono   CD8 T    Ptprc  Dpp4          0.438        0.164                  0.301
## 3 Mono   CD8 T    Cxcl10 Dpp4          4.27         0.164                  2.22 
## 4 Mono   CD8 T    Cxcl9  Dpp4          4.33         0.164                  2.25 
## 5 Treg   CD8 T    Ptprc  Dpp4          0.282        0.164                  0.223
## 6 Mono   CD8 T    Cxcl11 Dpp4          2.36         0.164                  1.26 
## # ℹ abbreviated name: ¹​ligand_receptor_lfc_avg
## # ℹ 44 more variables: p_val_ligand <dbl>, p_adj_ligand <dbl>,
## #   p_val_receptor <dbl>, p_adj_receptor <dbl>, pct_expressed_sender <dbl>,
## #   pct_expressed_receiver <dbl>, avg_ligand <dbl>, avg_receptor <dbl>,
## #   ligand_receptor_prod <dbl>, lfc_pval_ligand <dbl>,
## #   p_val_ligand_adapted <dbl>, scaled_lfc_ligand <dbl>,
## #   scaled_p_val_ligand <dbl>, scaled_lfc_pval_ligand <dbl>, …
```

Users may provide custom weights in a named vector to the argument `prioritizing_weights`. The vector must contain the following names, which correspond to the following criteria:

- `de_ligand`: upregulation of the ligand in a sender cell type compared to other cell types
- `de_receptor`: upregulation of the receptor in a receiver cell type
- `exprs_ligand`: average expression of the ligand in the sender cell type
- `exprs_receptor`: average expression of the receptor in the receiver cell type
- `ligand_condition_specificity`: condition-specificity of the ligand across all cell types
- `receptor_condition_specificity`: condition-specificity of the receptor across all cell types


```r
prioritizing_weights <- c("de_ligand" = 1,
                         "de_receptor" = 1,
                         "activity_scaled" = 1,
                         "exprs_ligand" = 1,
                         "exprs_receptor" = 1,
                         "ligand_condition_specificity" = 1,
                         "receptor_condition_specificity" = 1) 
```

35. Create a mushroom plot depicting ligand expression on one semicircle, and receptor expression on the other (Figure 3I).


```r
legend_adjust <- c(0.8, 0.7)
make_mushroom_plot(prioritized_table, top_n = 30,
                   show_all_datapoints = TRUE,
                   true_color_range = TRUE,
                   show_rankings = TRUE,
                   legend.title = element_text(size=8),
                   legend.text = element_text(size=8),
                   legend.key.size = unit(4, 'mm')) +
  theme(legend.justification = legend_adjust,
        axis.title.x = element_text(hjust = 0.25))
```

![](procedure_code_files/figure-html/prioritization-X-1.png)<!-- -->

## Box 1. Wrapper functions

To streamline the NicheNet analysis, we introduce three wrapper functions that automate Steps 5-23. The function `nichenet_seuratobj_aggregate` calculates the gene set of interest as the DE genes between two conditions within the receiver cell type. This function can be used to replicate the analysis in this paper as follows:


```r
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj,
  receiver = "CD8 T",
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
  condition_colname = "aggregate",
  condition_oi = "LCMV",
  condition_reference = "SS",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks) 
```

```
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"
```

```r
# Preview
names(nichenet_output)
```

```
##  [1] "ligand_activities"                     
##  [2] "top_ligands"                           
##  [3] "top_targets"                           
##  [4] "top_receptors"                         
##  [5] "ligand_target_matrix"                  
##  [6] "ligand_target_heatmap"                 
##  [7] "ligand_target_df"                      
##  [8] "ligand_expression_dotplot"             
##  [9] "ligand_differential_expression_heatmap"
## [10] "ligand_activity_target_heatmap"        
## [11] "ligand_receptor_matrix"                
## [12] "ligand_receptor_heatmap"               
## [13] "ligand_receptor_df"                    
## [14] "geneset_oi"                            
## [15] "background_expressed_genes"
```

Additionally, the sender-agnostic approach can be explicitly run by setting `sender = "undefined"`.

The resulting object is a list comprising various components, including the gene set of interest and background genes that were used for the analysis (`geneset_oi` and `background_expressed_genes`), output from the ligand activity analysis (`ligand_activities`), targets and receptors corresponding to the identified ligands (`top_targets` and `top_receptors`), and visualizations (`ligand_target_heatmap`, `ligand_expression_dotplot`, `ligand_receptor_heatmap`, etc.). If the sender-focused approach is used, the line plot from Step 24 will also be generated.

Another wrapper function, `nichenet_seuratobj_cluster_de`, calculates DE genes between two cell types as the gene set of interest. Additionally, when filtering for potential ligands, we only consider expressed ligands whose receptors are expressed by the "reference" receiver cell type. This function should only be used for specific biological scenarios, as shown in Figure 2. An example of using this function for the scenario where cell types differentiation occurs due to its niche is as follows:


```r
nichenet_seuratobj_cluster_de(
  seurat_obj = seurat_obj,
  receiver_affected = differentiated_celltype,
  receiver_reference = progenitor_cell,
  sender = niche_celltypes,
  ligand_target_matrix, lr_network,
  weighted_networks)
```

The final wrapper function, `nichenet_seuratobj_aggregate_cluster_de`, combines the aforementioned wrappers. The gene set of interest is calculated as the DE genes between the affected receiver cell type under the condition of interest and the reference receiver cell type in the reference condition.

## Box 2. Assessing the quality of predicted ligands

We can assess the quality of prioritized ligands by constructing a random forest model built using the top 30 predicted ligands. The aim is to evaluate its ability to predict if a particular gene belongs to the target gene set. Using kfold cross-validation, 1/k of the target gene set is isolated as “unseen” data. The performance of the model can then be evaluated using classification metrics like AUPR and AUROC, or using Fisher’s exact test on the confusion matrix.  


```r
# Define cross-validation folds and number of iterations
k <- 3
n <- 10

# Build random forest model and obtain prediction values
predictions_list <- lapply(1:n,
                           assess_rf_class_probabilities,
                           folds = k,
                           geneset = geneset_oi,
                           background_expressed_genes = background_expressed_genes,
                           ligands_oi = best_upstream_ligands,
                           ligand_target_matrix = ligand_target_matrix)

# Get classification metrics of the models, then calculate mean across all rounds
performances_cv <- bind_rows(lapply(predictions_list, classification_evaluation_continuous_pred_wrapper))
colMeans(performances_cv)
```

```
##                  auroc                   aupr         aupr_corrected 
##              0.8037847              0.4955904              0.4209850 
##        sensitivity_roc        specificity_roc mean_rank_GST_log_pval 
##              0.3759151              0.4278696            136.9676112 
##       pearson_log_pval      spearman_log_pval                pearson 
##            262.4489418             61.5030547              0.5399001 
##               spearman 
##              0.2765244
```

```r
# Calculate fraction of target genes and non-target genes
# that are among the top 5% predicted targets
fraction_cv <- bind_rows(lapply(predictions_list,
                                calculate_fraction_top_predicted,
                                quantile_cutoff = 0.95),
                         .id = "round")

# In this case, ~30% of target genes are in the top targets
# compared to ~1% of the non-target genes
mean(filter(fraction_cv, true_target)$fraction_positive_predicted)
```

```
## [1] 0.4488462
```

```r
mean(filter(fraction_cv, !true_target)$fraction_positive_predicted)
```

```
## [1] 0.01813953
```

```r
# Perform Fischer's exact test
lapply(predictions_list,
       calculate_fraction_top_predicted_fisher,
       quantile_cutoff = 0.95)
```

```
## [[1]]
## [1] 4.66469e-96
## 
## [[2]]
## [1] 1.20294e-102
## 
## [[3]]
## [1] 5.17928e-94
## 
## [[4]]
## [1] 1.20294e-102
## 
## [[5]]
## [1] 5.608708e-101
## 
## [[6]]
## [1] 7.513216e-104
## 
## [[7]]
## [1] 5.608708e-101
## 
## [[8]]
## [1] 1.20294e-102
## 
## [[9]]
## [1] 4.66469e-96
## 
## [[10]]
## [1] 3.881505e-88
```

```r
# Get which genes had the highest prediction values
top_predicted_genes <- lapply(1:n, get_top_predicted_genes, predictions_list)
top_predicted_genes <- reduce(top_predicted_genes, full_join, by = c("gene","true_target"))
top_predicted_genes
```

```
## # A tibble: 241 × 12
##    gene   true_target predicted_top_target_round1 predicted_top_target_round2
##    <chr>  <lgl>       <lgl>                       <lgl>                      
##  1 Gimap9 FALSE       TRUE                        TRUE                       
##  2 H2-Q4  FALSE       TRUE                        TRUE                       
##  3 Gbp4   TRUE        TRUE                        TRUE                       
##  4 Gbp9   TRUE        TRUE                        TRUE                       
##  5 Ifi203 TRUE        TRUE                        TRUE                       
##  6 Ifi209 TRUE        TRUE                        TRUE                       
##  7 Ifi213 TRUE        TRUE                        TRUE                       
##  8 Ifi208 TRUE        TRUE                        TRUE                       
##  9 Mndal  TRUE        TRUE                        TRUE                       
## 10 Ifi206 TRUE        TRUE                        TRUE                       
## # ℹ 231 more rows
## # ℹ 8 more variables: predicted_top_target_round3 <lgl>,
## #   predicted_top_target_round4 <lgl>, predicted_top_target_round5 <lgl>,
## #   predicted_top_target_round6 <lgl>, predicted_top_target_round7 <lgl>,
## #   predicted_top_target_round8 <lgl>, predicted_top_target_round9 <lgl>,
## #   predicted_top_target_round10 <lgl>
```


## Box 3. Constructing your own prior model

As the NicheNet prior model was constructed by integrating ligand-receptor, signaling, and gene regulatory databases, it is possible to replace each of these networks with external data sources. Constructing a customized prior model thus requires three directed networks, represented as data frames comprising three columns: `from`, `to`, and `source`. The `source` column contains the originating database of the interaction. As the reliability of each database can vary, we optimized the weights of each data source based on our validation procedure (as explained in Comparison with other methods). These optimized weights (`optimized_source_weights_df`), along with hyperparameters for model construction (`hyperparameter_list`), are provided in the NicheNet package. Key hyperparameters include correction factors for dominant hubs in the ligand-signaling and gene regulatory networks, as well as the central damping factor of the Personalized PageRank algorithm, the network propagation mechanism used to determine the ligand-target regulatory potential.


```r
zenodo_path <- "https://zenodo.org/record/7074291/files/"
lr_network <- readRDS(url(paste0(zenodo_path, "lr_network_human_21122021.rds")))
sig_network <- readRDS(url(paste0(zenodo_path, "signaling_network_human_21122021.rds")))
gr_network <- readRDS(url(paste0(zenodo_path, "gr_network_human_21122021.rds")))

# Aggregate the individual data sources in a weighted manner to obtain
# a weighted integrated signaling network
weighted_networks <- construct_weighted_networks(
  lr_network = lr_network,
  sig_network = sig_network,
  gr_network = gr_network,
  source_weights_df = rename(optimized_source_weights_df, weight = avg_weight))

# Downweigh the importance of signaling and gene regulatory hubs 
# Use the optimized parameters of this 
weighted_networks <- apply_hub_corrections(
  weighted_networks = weighted_networks,
  lr_sig_hub = hyperparameter_list[hyperparameter_list$parameter == "lr_sig_hub",]$avg_weight,
  gr_hub = hyperparameter_list[hyperparameter_list$parameter == "gr_hub",]$avg_weight
  ) 
```

In this example, we will calculate target gene regulatory potential scores for TNF and the combination TNF+IL6. 


```r
# To compute it for all 1248 ligands (~1 min):
# ligands <- as.list(unique(lr_network$from))

ligands <- list("TNF", c("TNF","IL6"))
ligand_target_matrix <- construct_ligand_target_matrix(
  weighted_networks = weighted_networks,
  ligands = ligands,
  algorithm = "PPR",
  damping_factor = hyperparameter_list[hyperparameter_list$parameter == "damping_factor",]$avg_weight,
  ltf_cutoff = hyperparameter_list[hyperparameter_list$parameter == "ltf_cutoff",]$avg_weight
  )
```


A frequent use case is one where users are interested in replacing the ligand-receptor network in NicheNet with those of expression permutation tools in order to make their results more comparable. To achieve this, we recommend employing LIANA (Dimitrov et al., 2022), a comprehensive CCC framework that integrates both resources and computational algorithms for ligand-receptor interaction inference. The `show_resources` function is used to check which resources are present in LIANA, and `select_resource` returns a data frame of the interactions in that resource. The `decomplexify` function of LIANA is necessary for this integration, as it separate receptors into their respective subunits. Note that unlike before, `source_weights_df` only represents unoptimized weights, where the weight of every data source is 1.


```r
if (!requireNamespace("liana", quietly=TRUE))
  devtools::install_github("saezlab/liana")

liana_db <- liana::decomplexify(liana::select_resource("CellPhoneDB")[[1]])
liana_db <- rename(liana_db, from = source_genesymbol, to = target_genesymbol)
liana_db$source <- "liana"
liana_db <- select(liana_db, from, to, source)

# Change source weights data frame (but in this case all source weights are 1)
source_weights <- add_row(source_weights_df,
                          source = "liana",
                          weight = 1, .before = 1)

# Construct weighted network as before
weighted_networks <- construct_weighted_networks(
  lr_network = liana_db,
  sig_network = sig_network,
  gr_network = gr_network,
  source_weights_df = source_weights)
```

What packages did I use?


```r
sessionInfo()
```

```
## R version 4.3.2 (2023-10-31)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS Stream 8
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Asia/Bangkok
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] forcats_1.0.0      stringr_1.5.0      dplyr_1.1.4        purrr_1.0.2       
##  [5] readr_2.1.2        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.4     
##  [9] tidyverse_1.3.1    SeuratObject_5.0.1 Seurat_4.4.0       nichenetr_2.1.0   
## 
## loaded via a namespace (and not attached):
##   [1] IRanges_2.34.1              progress_1.2.3             
##   [3] nnet_7.3-19                 goftest_1.2-3              
##   [5] vctrs_0.6.5                 spatstat.random_3.2-2      
##   [7] digest_0.6.33               png_0.1-8                  
##   [9] shape_1.4.6                 proxy_0.4-27               
##  [11] OmnipathR_3.9.6             ggrepel_0.9.4              
##  [13] deldir_2.0-2                parallelly_1.36.0          
##  [15] MASS_7.3-60                 reprex_2.0.1               
##  [17] reshape2_1.4.4              httpuv_1.6.13              
##  [19] foreach_1.5.2               BiocGenerics_0.46.0        
##  [21] withr_2.5.2                 xfun_0.41                  
##  [23] ggpubr_0.6.0                ellipsis_0.3.2             
##  [25] survival_3.5-7              memoise_2.0.1              
##  [27] zoo_1.8-12                  GlobalOptions_0.1.2        
##  [29] V8_4.3.3                    pbapply_1.7-2              
##  [31] Formula_1.2-5               prettyunits_1.2.0          
##  [33] promises_1.2.1              httr_1.4.7                 
##  [35] rstatix_0.7.2               globals_0.16.2             
##  [37] fitdistrplus_1.1-11         rstudioapi_0.15.0          
##  [39] miniUI_0.1.1.1              generics_0.1.3             
##  [41] base64enc_0.1-3             dir.expiry_1.8.0           
##  [43] curl_5.2.0                  S4Vectors_0.38.1           
##  [45] zlibbioc_1.46.0             ScaledMatrix_1.8.1         
##  [47] polyclip_1.10-6             randomForest_4.7-1.1       
##  [49] GenomeInfoDbData_1.2.10     xtable_1.8-4               
##  [51] doParallel_1.0.17           evaluate_0.23              
##  [53] S4Arrays_1.2.0              hms_1.1.3                  
##  [55] GenomicRanges_1.52.0        irlba_2.3.5.1              
##  [57] colorspace_2.1-0            filelock_1.0.2             
##  [59] visNetwork_2.1.2            ROCR_1.0-11                
##  [61] reticulate_1.34.0           readxl_1.4.3               
##  [63] spatstat.data_3.0-3         magrittr_2.0.3             
##  [65] lmtest_0.9-40               later_1.3.2                
##  [67] lattice_0.21-9              spatstat.geom_3.2-7        
##  [69] future.apply_1.11.0         scuttle_1.10.2             
##  [71] scattermore_1.2             shadowtext_0.1.2           
##  [73] cowplot_1.1.2               matrixStats_1.2.0          
##  [75] RcppAnnoy_0.0.21            class_7.3-22               
##  [77] Hmisc_5.1-0                 pillar_1.9.0               
##  [79] nlme_3.1-163                iterators_1.0.14           
##  [81] caTools_1.18.2              compiler_4.3.2             
##  [83] beachmat_2.16.0             stringi_1.7.6              
##  [85] gower_1.0.1                 tensor_1.5                 
##  [87] SummarizedExperiment_1.30.2 lubridate_1.9.3            
##  [89] devtools_2.4.3              plyr_1.8.9                 
##  [91] crayon_1.5.2                abind_1.4-5                
##  [93] locfit_1.5-9.8              haven_2.4.3                
##  [95] sp_2.1-2                    modelr_0.1.8               
##  [97] codetools_0.2-19            recipes_1.0.7              
##  [99] BiocSingular_1.16.0         bslib_0.6.1                
## [101] e1071_1.7-14                GetoptLong_1.0.5           
## [103] plotly_4.10.0               mime_0.12                  
## [105] splines_4.3.2               circlize_0.4.15            
## [107] DiagrammeRsvg_0.1           Rcpp_1.0.11                
## [109] basilisk_1.12.1             dbplyr_2.1.1               
## [111] sparseMatrixStats_1.12.2    cellranger_1.1.0           
## [113] knitr_1.45                  utf8_1.2.4                 
## [115] clue_0.3-64                 fs_1.6.3                   
## [117] listenv_0.9.0               checkmate_2.3.1            
## [119] DelayedMatrixStats_1.22.5   logger_0.2.2               
## [121] pkgbuild_1.4.3              ggsignif_0.6.4             
## [123] Matrix_1.6-4                statmod_1.5.0              
## [125] tzdb_0.4.0                  tweenr_2.0.2               
## [127] pkgconfig_2.0.3             tools_4.3.2                
## [129] cachem_1.0.8                viridisLite_0.4.2          
## [131] rvest_1.0.2                 DBI_1.1.3                  
## [133] fastmap_1.1.1               rmarkdown_2.11             
## [135] scales_1.3.0                grid_4.3.2                 
## [137] usethis_2.2.2               ica_1.0-3                  
## [139] liana_0.1.13                broom_0.7.12               
## [141] sass_0.4.8                  patchwork_1.1.3            
## [143] dotCall64_1.1-1             carData_3.0-5              
## [145] RANN_2.6.1                  rpart_4.1.21               
## [147] farver_2.1.1                yaml_2.3.8                 
## [149] MatrixGenerics_1.12.3       DiagrammeR_1.0.10          
## [151] foreign_0.8-85              cli_3.6.2                  
## [153] stats4_4.3.2                leiden_0.3.9               
## [155] lifecycle_1.0.4             caret_6.0-94               
## [157] uwot_0.1.16                 Biobase_2.60.0             
## [159] bluster_1.10.0              lava_1.7.3                 
## [161] sessioninfo_1.2.2           backports_1.4.1            
## [163] BiocParallel_1.34.2         timechange_0.2.0           
## [165] gtable_0.3.4                rjson_0.2.21               
## [167] ggridges_0.5.5              progressr_0.14.0           
## [169] parallel_4.3.2              pROC_1.18.5                
## [171] limma_3.56.2                edgeR_3.42.4               
## [173] jsonlite_1.8.8              bitops_1.0-7               
## [175] assertthat_0.2.1            Rtsne_0.17                 
## [177] spatstat.utils_3.0-4        BiocNeighbors_1.18.0       
## [179] metapod_1.8.0               jquerylib_0.1.4            
## [181] highr_0.10                  dqrng_0.3.2                
## [183] timeDate_4032.109           lazyeval_0.2.2             
## [185] shiny_1.7.1                 htmltools_0.5.7            
## [187] sctransform_0.4.0           rappdirs_0.3.3             
## [189] basilisk.utils_1.12.1       glue_1.6.2                 
## [191] spam_2.10-0                 XVector_0.40.0             
## [193] RCurl_1.98-1.12             scran_1.28.2               
## [195] gridExtra_2.3               igraph_1.2.11              
## [197] R6_2.5.1                    SingleCellExperiment_1.22.0
## [199] fdrtool_1.2.17              labeling_0.4.3             
## [201] cluster_2.1.4               pkgload_1.3.3              
## [203] GenomeInfoDb_1.36.1         ipred_0.9-14               
## [205] DelayedArray_0.26.7         tidyselect_1.2.0           
## [207] htmlTable_2.4.1             ggforce_0.4.1              
## [209] xml2_1.3.6                  car_3.1-2                  
## [211] future_1.33.0               ModelMetrics_1.2.2.2       
## [213] rsvd_1.0.5                  munsell_0.5.0              
## [215] KernSmooth_2.23-22          data.table_1.14.10         
## [217] htmlwidgets_1.6.2           ComplexHeatmap_2.16.0      
## [219] RColorBrewer_1.1-3          rlang_1.1.2                
## [221] spatstat.sparse_3.0-3       spatstat.explore_3.2-1     
## [223] remotes_2.4.2               ggnewscale_0.4.9           
## [225] fansi_1.0.6                 hardhat_1.3.0              
## [227] prodlim_2023.08.28
```
