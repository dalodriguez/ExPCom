#' ComScores: calculate communication scores
#' @param Seurat a Seurat object
#' @param idents the name of the Seurat column metadata that contains the idents for which the communication will be calculated.
#' @param database a dataframe with two columns containing ligand and receptors, default is the cellchat database.
#' @param threshold minimum percentage of cells expressing the cells in source and target populations. Default is set to 0

#' @return A dataframe including communication scores for each LR pair and interaction
#' @import Seurat stringr reshape2
#' @export
#'

ComScores <- function(Seurat, idents = "cell_type", database = database, threshold = 0){

    if (is.null(idents) == T){message("Add idents, the column metadata that contains the idents for which the communication will be calculated")}
    if (is.null(database) == T){message("Add a database of LR pairs")}

  # Extracting the count matrix from the Seurat object:
    Matrix = Seurat[["RNA"]]@data

  # Asign idents as default Seurat levels:
    idx = which(idents == colnames(Seurat@meta.data))
    Idents(Seurat) = Seurat@meta.data[,idx]
    message(paste("...The following idents will be used for computing cell-cell communication:"))
    message(cat(levels(Seurat),sep=", "))

  # Calculating the expression of ligand receptors in each population:
  preComm <- data.frame()

  for (ids in 1:length(levels(Seurat))){
    message(paste0('...Calculating LRs expression ',levels(Seurat)[ids], ": ", ids,' over ',length(levels(Seurat))))
    target <- which(Seurat@meta.data[,idx] ==levels(Seurat)[ids])
    others <- which(Seurat@meta.data[,idx] !=levels(Seurat)[ids])

    values <- apply(Matrix, 1, function(rows){
                    # Calculating the proportion of cells expressing the cells in the target population against the other populations:
                      count_average <- sum(rows[target])/length(others)
                    # Calculate the number of cells expressing the ligand/receptor in the target population
                      cell_number <- sum(rows[target] != 0)
                    # Calculate the proportion of cells expressing the ligand/receptor in the target population
                    cell_proportion <- (sum(rows[target] != 0)/length(target))*100 #The proportion of cells expressing the ligand/receptor in the cluster
                    c(count_average, cell_number, cell_proportion)})

    values <- as.data.frame(t(values))
    colnames(values) <- c( "count_average", "cell_number", "cell_proportion")
    values$cluster <- levels(Seurat)[ids]
    values$gene <- rownames(values) ; rownames(values) <- NULL
    values$L <- NA
    values$L[values$gene %in% database$Ligand] <- "L"
    values$R <- NA
    values$R[values$gene %in% database$Receptor] <- "R"
    values <- values[,c(4,5,6,7,1,2,3)]
    #values <- values[,c(9:12,1:8)]
    preComm <- rbind(preComm, values)
  }

  ##
  CommScores <- data.frame()
message("Starting Cell-Cell Communication:")
  for(ids in 1:length(levels(Seurat))){
    message(paste0("...Calculating CommScores for Source ", levels(Seurat)[ids]))
    ligands <- preComm[preComm$cluster == levels(Seurat)[ids],]
    ligands <- ligands[(ligands$cell_proportion > threshold),]
    ligands <- ligands[!is.na(ligands$L),]

    for (rec in 1:length(levels(Seurat))){
      tryCatch({
        message("     ...with Target ", levels(Seurat)[rec],": " ,rec, "/", length(levels(Seurat)))
        receptors <- preComm[preComm$cluster == levels(Seurat)[rec],]
        receptors <- receptors[(receptors$cell_proportion > threshold),]
        receptors <- receptors[!is.na(receptors$R),]

        pairs <- which((database$Ligand %in% ligands$gene) & (database$Receptor %in% receptors$gene))
        ligands <- ligands[ligands$gene %in% database$Ligand[pairs],]
        receptors <- receptors[receptors$gene %in% database$Receptor[pairs],]
        LRs <- data.frame()

        for (n in 1:nrow(ligands)){
          list <- database$Receptor[database$Ligand %in% ligands$gene[n]]
          list <- list[list %in% receptors$gene]

          for (r in list) {
            rows <- data.frame(ligands[n,c(1,2,5:7)], receptors[receptors$gene %in% r,c(1,2,5:7)])
            colnames(rows) <-c("cell.from","ligand", "lig.count_average", "lig.cell_number","lig.cell_prop", "cell.to", "receptor", "rec.count_average", "rec.cell_number", "rec.cell_prop")
            LRs <- rbind(LRs, rows)
          }
        }
        LRs <- LRs[,c(1,6,2,7,3,4,5,8,9,10)]
        LRs$CommScores <- NA

        for (t in 1:nrow(LRs)){
          LRs$CommScores[t]  <- mean(outer(Matrix[LRs$ligand[t], rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == LRs$cell.from[t]]],
                                           Matrix[LRs$receptor[t], rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == LRs$cell.to[t]]],
                                           FUN = "*"))
                              }


        CommScores <- rbind(CommScores, LRs)
      },
      error=function(e) {message('Error!')},
      warning=function(w) {message('Warning!')}
      )

    }
  }
  return(CommScores)
}



#' CondComScores: calculate communication scores for each condition in a Seurat Object
#'
#' @param Seurat a Seurat object
#' @param idents a vector containing the list of idents to calculate communication scores, default is levels(Seurat)
#' @param conditions string, specify column in seurat meta.data containing the conditions
#' @param database the cellchat database or any other CCC database containing a set of LR pairs
#' @return A dataframe including communication scores for each LR pair and interaction
#' @import Seurat stringr 
#' @export

CondComScores <- function(Seurat, idents = "cell_type", condition = "orig.ident",database = database, threshold = 0){

  if (is.null(idents) == T){message("ERROR: enter idents, the column metadata that contains the idents for which the communication will be calculated")}
  if (is.null(condition) == T){message("ERROR: enter condition, the column metadata that contains the conditions for which the communication will be calculated")}
  if (is.null(database) == T){message("ERROR: enter a database of LR pairs")}

  # Determine how many groups will be compared & subset the Matrix
  codx = which(condition  == colnames(Seurat@meta.data))
  if (length(codx) == 0){ message("ERROR: enter a valid metadata column") }
  if (length(unique(Seurat@meta.data[,codx])) == 1){ message("ERROR: there is only one condition, select the right metadata column.")}

        Matrix <- list()
        for (x in 1:length(unique(Seurat@meta.data[,codx]))){
          Matrix[[x]] <- Seurat[["RNA"]]@data[,rownames(Seurat@meta.data)[Seurat@meta.data[,codx] == unique(Seurat@meta.data[,codx])[x]]]
        }
        names(Matrix) <- unique(Seurat@meta.data[,codx])

  # Asign idents as default Seurat levels:
  idx = which(idents == colnames(Seurat@meta.data))
  Idents(Seurat) = Seurat@meta.data[,idx]
  message(paste("...The following idents will be used for computing cell-cell communication:"))
  message(cat(levels(Seurat),sep=", "))

  # Calculating the expression of ligand receptors in each population:
  preComm <- data.frame()
  
  for(t in 1:length(Matrix)){
      for (ids in 1:length(levels(Seurat))){
        message(paste0('...Calculating LRs expression in ',levels(Seurat)[ids], ": ", ids,' over ',length(levels(Seurat))))
        Mtx = Matrix[[t]]
        target <- which(Seurat@meta.data[colnames(Mtx),idx] ==levels(Seurat)[ids])  ###  !
        others <- which(Seurat@meta.data[colnames(Mtx),idx] !=levels(Seurat)[ids])  ###  !
        values <- apply(Mtx, 1, function(rows){
          # Calculating the proportion of cells expressing the cells in the target population against the other populations:
          count_average <- sum(rows[target])/length(others)
          # Calculate the number of cells expressing the ligand/receptor in the target population
          cell_number <- sum(rows[target] != 0)
          # Calculate the proportion of cells expressing the ligand/receptor in the target population
          cell_proportion <- (sum(rows[target] != 0)/length(target))*100 #The proportion of cells expressing the ligand/receptor in the cluster
          c(count_average, cell_number, cell_proportion)})

        values <- as.data.frame(t(values))
        colnames(values) <- c( "count_average", "cell_number", "cell_proportion")
        values$cluster <- levels(Seurat)[ids]
        values$gene <- rownames(values) ; rownames(values) <- NULL
        values$L <- NA
        values$L[values$gene %in% database$Ligand] <- "L"
        values$R <- NA
        values$R[values$gene %in% database$Receptor] <- "R"
        values$cond = names(Matrix)[t]

        values <- values[,c(4,5,6,7,1,2,3,8)]
        #values <- values[,c(9:12,1:8)]
        preComm <- rbind(preComm, values)
      }
  }
  
  message("Starting Cell-Cell Communication:")
  CondComScores <- list()
  for (x in 1:length(Matrix)){ 
      message(names(Matrix)[x])
      Mtx = Matrix[[x]]
      levels = unique(Seurat@meta.data$cell_type[Seurat@meta.data$orig.ident == names(Matrix)[x]])
      LR_commScores <- data.frame()
      
      for(ids in 1:length(levels)){
          message(paste0("...Calculating CommScores for Source ", levels[ids]))
          ligands <- preComm[preComm$cluster == levels[ids],]
          ligands <- ligands[ligands$cond == names(Matrix)[x],]
          ligands <- ligands[(ligands$cell_proportion > threshold),]
          ligands <- ligands[!is.na(ligands$L),]
    
          for (rec in 1:length(levels)){
            tryCatch({
              message("     ...with Target ", levels[rec],": " ,rec, "/", length(levels))
              receptors <- preComm[preComm$cluster == levels(Seurat)[rec],]
              receptors <- receptors[receptors$cond == names(Matrix)[x],]
              receptors <- receptors[(receptors$cell_proportion > threshold),]
              receptors <- receptors[!is.na(receptors$R),]
              
              pairs <- which((database$Ligand %in% ligands$gene) & (database$Receptor %in% receptors$gene))
              ligands <- ligands[ligands$gene %in% database$Ligand[pairs],]
              receptors <- receptors[receptors$gene %in% database$Receptor[pairs],]
              
              lists <- lapply(ligands$gene, function(ligand){ 
                list <- database$Receptor[database$Ligand %in% ligand]
                list <- list[list %in% receptors$gene] })
              names(lists) <- ligands$gene
              
              LRs <- lapply(lists, function(sublist){
                      recp <- lapply(sublist, function(element){
                                      row <- data.frame(receptors[receptors$gene %in% element,c(1,2,7)], row.names = NULL) })
                      recp <- do.call(rbind, recp)
                })
              
              LRs <- do.call(rbind, LRs)

              colnames(LRs) <-c( "cell.to", "receptor", "rec.cell_prop")
              
              ligands$lig_rep <- NA
              for (z in 1:length(lists)){ ligands$lig_rep[z]  <- length(lists[[z]])}
              lgnd <- as.data.frame(lapply(ligands[,c(1,2,7)], rep, ligands$lig_rep))
              
              LRs <- cbind(lgnd, LRs)
              colnames(LRs)[1:3] <-c( "cell.from", "ligand", "lig.cell_prop")
              rownames(LRs) <- NULL
              
              LRs <- LRs[,c(1,4,2,5,3,6)]
              LRs$CommScores <- NA
              LRs$CommScores <- apply(LRs, 1, function(rows){
                      CommScores <- mean(outer(Mtx[rows[3],rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == rows[1] &  Seurat@meta.data[,codx]  == names(Matrix)[x] ]], 
                                         Mtx[rows[4], rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == rows[2] &  Seurat@meta.data[,codx]  == names(Matrix)[x] ]], FUN = "*"))
                       c(CommScores )})
              
              colnames(LRs) <- paste0(colnames(LRs),"_",names(Matrix)[x])
              LR_commScores <- rbind(LR_commScores, LRs)
            },
            error=function(e) {message('Error!')},
            warning=function(w) {message('Warning!')}
            )
          }
      }
      CondComScores[[x]] <- LR_commScores
  }
  return(CondComScores)
  }

#' DifComScores: perform differential communication analysis between each condition 
#'
#' @param Seurat a Seurat object
#' @param idents a vector containing the list of idents to calculate communication scores, default is levels(Seurat)
#' @param conditions string, specify column in seurat meta.data containing the conditions
#' @param CondComScores the output of the CondComScores function
#' @param order a vector, order of the conditions to be compared for the contrast
#' @return A dataframe 
#' @import Seurat stringr reshape2 multcomp
#' @export

DifComScores <- function(Seurat, idents = "cell_type",  condition = "orig.ident",  CondComScores,  order = c("fed_int","fast12_int","fast24_int")){

  # Position of the Seurat metadata column containing the idents:
    idx = which(idents == colnames(Seurat@meta.data))
    if (length(idx) == 0){ message("ERROR: enter a valid metadata column") }
  
  # Position of the Seurat metadata column containing the condition:
    codx = which(condition  == colnames(Seurat@meta.data))
    if (length(codx) == 0){ message("ERROR: enter a valid metadata column") }
    if (length(unique(Seurat@meta.data[,codx])) == 1){ message("ERROR: there is only one condition, select the right metadata column.")}
  
  # Split matrix by conditions into a list:
    Matrix <- list()
    for (x in 1:length(unique(Seurat@meta.data[,codx]))){
      Matrix[[x]] <- Seurat[["RNA"]]@data[,rownames(Seurat@meta.data)[Seurat@meta.data[,codx] == unique(Seurat@meta.data[,codx])[x]]]}
    names(Matrix) <- unique(Seurat@meta.data[,codx])
  
  # Based on the results of the CondComScores function, obtain the unique LRpairs
    LRlist <- list()
    for (lrs in 1:length(CondComScores)){LRlist[[lrs]] <- unique(paste0(CondComScores[[lrs]][,3], "-",CondComScores[[lrs]][,4] ))}
    intersect <- as.data.frame(table(stack(setNames(LRlist, seq_along(LRlist)))))
    intersect <- intersect[!duplicated(intersect$values),]
  
  #  Differential communication analysis:
    message("Starting Differential Cell-Cell Communication Analysis:")
    message(paste0('... Total of LR pairs: ',length(intersect$values)))
    
    Results <- data.frame(data.frame(matrix(ncol = c(2+length(Matrix)+ncol(combn(names(Matrix),2))), nrow = 0)))
    colnames(Results)[1:2] <- c("LRpair","SourceTarget")
    for(k in 1:(length(Matrix))){ colnames(Results)[c(k+2)] <- paste0("CommScore_",names(Matrix)[k]) }
    for(z in 1:ncol(combn(names(Matrix),2))){ colnames(Results)[ncol(Results)-ncol(combn(names(Matrix),2))+z] = paste0(combn(names(Matrix),2)[1,z],"-",combn(names(Matrix),2)[2,z]) }
    
    for (y in 1:length(intersect$values)){ 
      message(paste0("    ...",y, " ",as.character(intersect$values[y])))
      
      for (w in 1:length(Matrix)){ 
        
            Cm <-  CondComScores[[w]][CondComScores[[w]][,3] == strsplit(as.character(intersect$values[y]), "-")[[1]][1] & CondComScores[[w]][,4] == strsplit(as.character(intersect$values[y]), "-")[[1]][2],]
           
             if (nrow(Cm) > 0){   
              
                for (g in unique(unique(paste0(Cm[,1],"-",Cm[,2]))) ){   
                    
                        Comm <- Cm[Cm[,1] == strsplit(g, "-")[[1]][1] & Cm[,2] == strsplit(g, "-")[[1]][2],]
                        ComRes <- data.frame(LRpair = as.character(intersect$values[y]), SourceTarget = paste0(Comm[1,1],"-",Comm[1,2]), data.frame(matrix(NA, nrow = 1, ncol = c(ncol(Results)-2)  )))
                        colnames(ComRes)[3:ncol(ComRes)]  <- colnames(Results)[3:ncol(Results)] 
                        
                        #ComRes[,which(paste0("CommScore_",names(Matrix)[w]) == colnames(ComRes))] = Comm[,7]
            
                        expr <-outer(Matrix[[w]][ strsplit(as.character(intersect$values[y]), "-")[[1]][1],
                                                  rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == strsplit(g, "-")[[1]][1] &  Seurat@meta.data[,codx]  == names(Matrix)[w] ]], 
                                     Matrix[[w]][ strsplit(as.character(intersect$values[y]), "-")[[1]][2], 
                                                  rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == strsplit(g, "-")[[1]][2] &  Seurat@meta.data[,codx]  == names(Matrix)[w] ]], FUN = "*")
                        
                        expr <- na.omit(melt(expr)) ;  expr$cond <- names(Matrix)[w]
                        
                        ComRes[,which(paste0("CommScore_",names(Matrix)[w]) == colnames(ComRes))] = mean(expr$value)
                        
                        
                        for (b in (1:length(Matrix))[-w]){
               
                            Comm2 <-  CondComScores[[b]][CondComScores[[b]][,3] == strsplit(as.character(intersect$values[y]), "-")[[1]][1] & CondComScores[[b]][,4] == strsplit(as.character(intersect$values[y]), "-")[[1]][2],]
                            
                              if (nrow(Comm2) > 0 ){ 
                                    Comm2 <- Comm2[Comm2[,1] == strsplit(g, "-")[[1]][1] & Comm2[,2] == strsplit(g, "-")[[1]][2],]
                                    #ComRes[,which(paste0("CommScore_",names(Matrix)[b]) == colnames(ComRes))] = Comm2[,7]  
                
                                    exprbis <-outer(Matrix[[b]][ strsplit(as.character(intersect$values[y]), "-")[[1]][1],
                                                              rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == strsplit(g, "-")[[1]][1] &  Seurat@meta.data[,codx]  == names(Matrix)[b] ]], 
                                                 Matrix[[b]][ strsplit(as.character(intersect$values[y]), "-")[[1]][2], 
                                                              rownames(Seurat@meta.data)[Seurat@meta.data[, idx] == strsplit(g, "-")[[1]][2] &  Seurat@meta.data[,codx]  == names(Matrix)[b] ]], FUN = "*")
                                    
                                    exprbis <- na.omit(melt(exprbis)) ;  exprbis$cond <- names(Matrix)[b]
                                    
                                    ComRes[,which(paste0("CommScore_",names(Matrix)[b]) == colnames(ComRes))] = mean(exprbis$value)
                                    
                                    expr <- rbind(expr,exprbis )
                                     } 
                                  }
                        if ( length(unique(expr$cond)) == length(Matrix)){ 
                          expr$cond <- factor(expr$cond, levels= order ) 
                          res.aov <- aov(value ~ cond, data = expr)  ; res2 <- summary(res.aov) 
                          if (res2[[1]][["Pr(>F)"]][1] > 0.05 | res2[[1]][["Pr(>F)"]][1] == "NaN"  ){ } 
                          else { 
                            post_test <-summary(glht(res.aov, linfct = mcp(cond = "Tukey"))) 
                            ComRes[,(ncol(ComRes)-ncol(combn(names(Matrix),2))+1):ncol(ComRes)] = post_test[["test"]][["pvalues"]]
                          } 
                        } 
                      }  
               Results <- rbind(Results,ComRes)
                    } 
              }
          }
  Results <- Results[!duplicated(Results),]
  return(Results)
    }
