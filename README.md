# The test2 package
This package provides blablabla for cell-cell communnication 

### Installation 
test2 can be installed with:
```
devtools::install_github("dalodriguez/test2")
```


### Examples

**Calculating CCC**  

In this example, we'll calculate CCC in a single condition
```
CCC <- ComScores(seurat)

```


In this example, we'll calculate CCC in several conditions
```

conditions = "condition"
cell_type= "cell_type"

CCC <- DifComScores(seurat, conditions= conditions, cell_type = cell_type)

```
