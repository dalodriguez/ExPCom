# The test2 package
This package provides blablabla for cell-cell communnication 

### Installation 
test2 can be installed with:
```
devtools::install_github("dalodriguez/test2")
library(Testing)
```


### Examples

**Calculating CCC**  

In this example, we'll calculate CCC in a single condition
```
database= path_to_database
CCC <- ComScores(seurat_object, database)

```


In this example, we'll calculate CCC in several conditions
```

conditions = "condition"
cell_type= "cell_type"
database= path_to_database

CCC <- DifComScores(seurat_object, conditions= conditions, cell_type = cell_type)

```
