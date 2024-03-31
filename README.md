# ExPCom
This package provides a tool to infer cell-cell communnication using a scRNAseq dataset. The functions calculate communication scores for any set of ligand-receptor pairs (e.g. CellChat Database) based on expression product. The L*R expression product method allows to obtain a continuous value by multiplying all possible combinations of the expression of both interacting molecules between two given cell types.

<p align="center">
  <img width="1100"  src="TanyCom.png">
  <center>Illustrative picture of the original workflow used in https://www.biorxiv.org/content/10.1101/2023.07.06.547914v1</center>
</p>


## Installation 
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

Combined scores can be calculated in a single condition using the ComScores() function. The following arguments are required. (1) A SeuratObject: a sample dataset can be download from this github repository. The dataset is a downsample of the Brunner et al. (2024) paper. (2) The name of the metadata column containing the idents. (3) A database with ligand and receptor pairs. The database database.csv can also be found in the repository. Any custom database can be used as long as it is a dataframe with the following format: (4) a threshold, the minimum percentage of cells expressing the cells in source and target populations. Default is set to 0.

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

The ComScores() function returns a dataframe cotainning all possible cell-cell communication scores (e.g. tanycyte-neurons). 
```
database = path_to_database

CCC <- ComScores(seurat_object,  idents= levels(seurat_object), database)

```
<ul>
  <li><b>Calculating CCC in a dataset with several conditions</b></li>
</ul>

By using an integrated seurat object containing two or more conditions, we can calculate communication scores for each condition and infer the dynamics of cell-cell communication using the DifComScores() function.

The DifComScores() function requires to define the name of the seurat_object meta.data column containing the conditions and cell types. As bellow, the database and the name of the idents to be included can be reduced if required. 

```
conditions = "condition"
cell_type = "cell_type"
database = path_to_database

CCC <- DifComScores(seurat_object, conditions= conditions, cell_type = cell_type)
```




## Visualization
The results the DifComScores() function can be used to evaluate the differences of communication between conditions by calculating the delta of communication scores obtained or by looking at the dynamics of cell-cell communication. 

<p align="center">
  <img width="800"  src="CCC.png">
  <center>Example of the dynamics of cell-cell commuication</center>
</p>


## Suggestions and contributions 
Please use github issue tracker to report coding related issues or contact us directly, https://dalodriguez.github.io/Contact/

## How to cite?
M. Brunner&#185;, D. Lopez-Rodriguez&#185;, A. Messina, B. Thorens, F Santoni&#178;, F. Langlet&#178;. Pseudospatial transcriptional gradient analysis of hypothalamic ependymal cells: towards a new tanycyte classification. BioRxiv preprint. https://doi.org/10.1101/2023.07.06.547914

&#185; Co-First Author

&#178; Co-Last Author

https://www.biorxiv.org/content/10.1101/2023.07.06.547914v1



<p align="center">
  <a href="#">
     <img src="https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fsqjin%2FCellChat&labelColor=%233499cc&countColor=%2370c168](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fdalodriguez%2Ftest2&label=%20-%20&labelColor=%23000000&countColor=%23d9e3f0&style=flat-square)https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fdalodriguez%2Ftest2&label=%20-%20&labelColor=%23000000&countColor=%23d9e3f0&style=flat-square" />
   </a>
</p>

