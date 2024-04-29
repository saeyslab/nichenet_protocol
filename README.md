# NicheNet Protocol

This repository contains the R code found in https://arxiv.org/abs/2404.16358.

## Installation

To install NicheNet and its dependencies, enter the following code in R:

```
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

devtools::install_github("saeyslab/nichenetr")
```

Installation typically takes a few minutes, depending on the number of dependencies that has already been installed on your PC. nichenetr was tested on both Windows and Linux (most recently tested R version: R 4.3.2).

For the most up-to-date changes on NicheNet, please refer to our main GitHub repository: https://github.com/saeyslab/nichenetr.


## Content
### Input files

| File                    | Mouse                                                                                                   | Human                                                                                             |
|-------------------------|---------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| Example seurat object   | [Link](https://zenodo.org/record/3531889/files/seuratObj.rds)*                                           | N/A                                                                                               |
| Ligand-receptor network | [Link](https://zenodo.org/records/7074291/files/lr_network_mouse_21122021.rds?download=1)*               | [Link](https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds?download=1)‡         |
| Ligand-target matrix    | [Link](https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds?download=1)* | [Link](https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds?download=1) |
| Weighted networks       | [Link](https://zenodo.org/records/7074291/files/weighted_networks_nsga2r_final_mouse.rds?download=1)†    | [Link](https://zenodo.org/records/7074291/files/weighted_networks_nsga2r_final.rds?download=1)    |
| Ligand-TF matrix        | [Link](https://zenodo.org/records/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds?download=1)†     | [Link](https://zenodo.org/records/7074291/files/ligand_tf_matrix_nsga2r_final.rds?download=1)     |
| Signaling network       | [Link](https://zenodo.org/records/7074291/files/signaling_network_mouse_21122021.rds?download=1)        | [Link](https://zenodo.org/records/7074291/files/signaling_network_human_21122021.rds?download=1)‡  |
| Gene regulatory network | [Link](https://zenodo.org/records/7074291/files/gr_network_mouse_21122021.rds?download=1)               | [Link](https://zenodo.org/records/7074291/files/gr_network_human_21122021.rds?download=1)‡         |

\* Required to run the NicheNet ligand activity analysis and target gene prediction\
† Required to create certain visualizations (Steps 16, 30)\
‡ Required for Box 1 code

### Code files

* `procedure_code.*`: code from the _Procedure_ section provided in various file formats. The .html and .md files facilitate viewing on the web, while the .Rmd and .R files allow for easier modification. The code to reproduce Figure 3 is provided at the end of `procedure_code.R`.
* `box_code.R`: code from Boxes 1-3
* `timing.R`: code used to time the functions
