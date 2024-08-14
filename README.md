# ExPCom
This repository provides a tool to infer cell-cell communnication using a scRNAseq dataset. ExPCom allows to replicate the analysis found in Brunner et al. (2024). The functions calculate communication scores for any set of ligand-receptor pairs (e.g. CellChat Database) based on expression product. The L*R expression product method allows to obtain a continuous value by multiplying all possible combinations of the expression of both interacting molecules between two given cell types.

<p align="center">
  <img width="1100"  src="TanyCom.png">
  <center>Illustrative picture of the original workflow used in https://www.biorxiv.org/content/10.1101/2023.07.06.547914v1</center>
</p>


## Installation & System Requirements
ExPCom is supported in Windows, MAC OS and Linux and requires only a standard computer with enough RAM to support the in-memory operations.

ExPCom requires to have installed the package devtools, Seurat, stringr, reshape2 & multcomp. 
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("devtools","Seurat", "stringr", "reshape2", "multcomp")
library(devtools)

devtools::install_github("dalodriguez/ExPCom")

library(ExPCom)
```

## Example
<ul>
  <li><b>Calculating CCC in a dataset with a single condition</b></li>
</ul>

Combined scores can be calculated in a single condition using the ComScores() function. The following arguments are required. (1) A SeuratObject: a sample dataset can be loaded this github repository. The dataset is a downsample of the Brunner et al. (2024) paper. (2) The name of the metadata column containing the idents. (3) A database with ligand and receptor pairs. The cellchat database is available in the repository. Any custom database can be used as long as it is a dataframe with the following format: (4) a threshold, the minimum percentage of cells expressing the cells in source and target populations. Default is set to 0.

<table style="border:1px solid black;margin-left:auto;margin-right:auto;">
  <tr>
    <th>Ligand</th>
    <th>Receptor</th>
    <th>LR</th>
  </tr>
  <tr>
    <td>Kiss1</td>
    <td>Kiss1r</td>
    <td>Kiss1-Kiss1r</td>
  </tr>
    <tr>
    <td>...</td>
    <td>...</td>
    <td>...</td>
  </tr>
    <tr>
    <td>Lig</td>
    <td>Rec</td>
    <td>Lig-Rec</td>
  </tr>
</table>

```

data("LRPairs_CellChat_db")
data("TanSeurat")

ComScores <- ComScores(Seurat, idents = "cell_type", database = database, threshold = 0)

```


The ComScores() function returns a dataframe cotainning all possible cell-cell communication scores (e.g. tanycyte-neurons). 

<ul>
  <li><b>Calculating CCC in a dataset with several conditions</b></li>
</ul>

By using an integrated seurat object containing two or more conditions, differential communication analysis can be performed by calculating communication scores for each condition and infer the dynamics of cell-cell communication using the CondComScores() and DifComScores() functions.

The CondComScores() function simply calculates communication scores as the ComScores() function but taking into account conditions. The following arguments are required. (1) A SeuratObject (2) The name of the metadata column containing the idents. (3) The name of the metadata column containing the conditions. (4) A database with ligand and receptor pairs. (5) a threshold, the minimum percentage of cells expressing the cells in source and target populations. Default is set to 0.

```
database <- readRDS("CellChatDB.rds")
Seurat <- readRDS("TanSeurat.rds")

CondComScores <- CondComScores(Seurat, idents = "cell_type", condition = "orig.ident",database = database, threshold = 0)

```
  
The DifComScores() perfom a contrast between conditions using ANOVA and Tukey posthoc. The following arguments are required. (1) A SeuratObject (2) The name of the metadata column containing the idents. (3) The name of the metadata column containing the conditions. (4) The output of the CondComScores function, (5) a vector, order of the conditions to be compared for the contrast. 

```

DifComScores <- DifComScores(Seurat, idents = "cell_type", CondComScores = CondComScores, condition = "orig.ident",order = c("fed_int","fast12_int","fast24_int") )

```

Notes:
1)  Downstream analysis requires to correct p values for multiple comparisons. 
2)  The function requires to be adapted for parallelization for large datasets 






## Visualization
The results the DifComScores() function can be used to evaluate the differences of communication between conditions by calculating the delta of communication scores obtained or by looking at the dynamics of cell-cell communication. 

<p align="center">
  <img width="800"  src="CCC.png">
  <center>Example of the dynamics of cell-cell commuication</center>
</p>


## Suggestions and contributions 
Please use github issue tracker to report coding related issues or contact me directly, https://dlopez-rodriguez.ch/Contact/

## How to cite?
Brunner, M.&#185;, Lopez-Rodriguez, D.&#185;, Estrada-Meza, J. et al. Fasting induces metabolic switches and spatial redistributions of lipid processing and neuronal interactions in tanycytes. Nat Commun 15, 6604 (2024). https://doi.org/10.1038/s41467-024-50913-w

&#185; Co-First Author


## License
This project is covered under the Creative Commons Zero v1.0 Universal Licnse






   ![ExPCom](https://visitor-badge.laobi.icu/badge?page_id=dalodriguez.ExPCom)

