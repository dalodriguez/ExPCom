#' ComScores: calculate communication scores
#' @param Seurat a Seurat object
#' @param idents a vector containing the list of idents to calculate communication scores, default is levels(Seurat)
#' @param database the cellchat database or any other CCC database containing a set of LR pairs
#' @return A list of dataframes including communication scores for each LR pair and interaction
#' @import Seurat stringr reshape2
#' @export

ComScores <- function(Seurat, idents = levels(Seurat), database) {
  # Prepare interactions:
    Interactions <- c()
    for (i in 1:length(idents)){
      Int1 <- rep(idents[i], (length(idents))-1) ; Int2 <- idents ; Int2 <-  Int2[-i]
      Int <- paste0(Int1,"-",Int2)
      Interactions <- c(Interactions, Int)}
  # Load database:
    database <- cbind(database, setNames(rep(list(""), length(Interactions)), c(Interactions)))
    for (i in 4:ncol(database)){database[i] <- colnames(database)[i]}
  # Calculate Cell-cell communication scores based on LR products:
    CCC <- list()
    for (i in 1:length(Interactions)){
      print(Interactions[i])
      CS <- database[,c(1,2,3,grep(Interactions[i], colnames(database)))]
      colnames(CS)[4] <- 'Interaction'
      CS$Cell1 <- str_split(CS$LR, "-", simplify = TRUE)[, 1]
      CS$Cell2 <- str_split(CS$LR, "-", simplify = TRUE)[, 2]
      CS$ct1  <- str_split(Interactions[i], "-", simplify = TRUE)[, 1]
      CS$ct2  <-str_split(Interactions[i], "-", simplify = TRUE)[, 2]
      CS$ComScore <- NA
      for (t in 1:20){ #nrow(CS) FOR TESTING 1:20
        #print(paste(Interactions[i], "LR_",t))
        tryCatch({
          cells1 = WhichCells(Seurat, idents = CS$ct1[t])
          cells2 = WhichCells(Seurat, idents = CS$ct2[t])
          expr1 = FetchData(Seurat, cells = cells1, CS$Ligand[t])
          expr2 = FetchData(Seurat, cells = cells2, CS$Receptor[t])
          LR <- outer(expr1[,1], expr2[,1], FUN = "*")
          CS$ComScore[t] <- mean(LR)
        }, error = function(e) {} )}
      CCC[[i]] <- CS
    }
    names(CCC) <- Interactions
    return(CCC)
}


#' DifComScores: calculate communication scores for each condition in a Seurat Object
#'
#' @param Seurat a Seurat object
#' @param idents a vector containing the list of idents to calculate communication scores, default is levels(Seurat)
#' @param conditions string, specify column in seurat meta.data containing the conditions
#' @param cell_type  string, specify column in seurat meta.data containing the idents
#' @param database the cellchat database or any other CCC database containing a set of LR pairs
#' @return A list of dataframes including communication scores for each LR pair and interaction
#' @import Seurat stringr reshape2
#' @export

DifComScores <- function(seurat, idents = levels(seurat), conditions, cell_type, database) {
    # Prepare interactions:
      Interactions <- c()
      for (a in 1:length(idents)){
        Int1 <- rep(idents[a], (length(idents))-1) ; Int2 <- idents ; Int2 <-  Int2[-a]
        Int <- paste0(Int1,"-",Int2)
        Interactions <- c(Interactions, Int)}
    # Load database:
      database <- cbind(database, setNames(rep(list(""), length(Interactions)), c(Interactions)))
      for (u in 4:ncol(database)){database[u] <- colnames(database)[u]}
    # Calculate Cell-cell communication scores based on LR products:
      CCC <- list()
      for (i in 1:length(Interactions)){
        print(Interactions[i])
        CS <- database[,c(1,2,3,grep(Interactions[i], colnames(database)))]
        colnames(CS)[4] <- 'Interaction'
        CS$Cell1 <- str_split(CS$LR, "-", simplify = TRUE)[, 1]
        CS$Cell2 <- str_split(CS$LR, "-", simplify = TRUE)[, 2]
        CS$ct1  <- str_split(Interactions[i], "-", simplify = TRUE)[, 1]
        CS$ct2  <-str_split(Interactions[i], "-", simplify = TRUE)[, 2]
        position = which(colnames(seurat@meta.data)==conditions)
        cond <- unique(seurat@meta.data[,position])
        for(c in 1:length(cond)){
          CS[, ncol(CS) + 1] <- NA
          names(CS)[ncol(CS)] <- paste0("ComScore_", cond[c])}
        for (t in 1:length(cond)){
          Idents(seurat) <- seurat@meta.data[,position]
          seurat_sub <- subset(seurat, idents = cond[t])
            for (k in 1:5){ #nrow(CS) FOR TESTING 1:20
              tryCatch({
                  ct_pos = which(colnames(seurat_sub@meta.data)==cell_type)
                  Idents(seurat_sub) <- seurat_sub@meta.data[,ct_pos]
                  cells1 = WhichCells(seurat_sub, idents = CS$ct1[k])
                  cells2 = WhichCells(seurat_sub, idents = CS$ct2[k])
                  expr1 = FetchData(seurat_sub, cells = cells1, CS$Ligand[k])
                  expr2 = FetchData(seurat_sub, cells = cells2, CS$Receptor[k])
                  LR <- outer(expr1[,1], expr2[,1], FUN = "*")
                  pos <- which(colnames(CS)== paste0("ComScore_", cond[t]))
                  CS[k,pos] <- mean(LR)
                    }, error = function(e) {} )}
                  }
                CCC[[i]] <- CS
                }

    names(CCC) <- Interactions
    return(CCC)
}














