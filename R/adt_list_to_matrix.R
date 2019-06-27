#' @export

adt_list_to_matrix <- function(adt_list){
  
  bad_rows <- c("bad_struct_adt",
                "no_match_adt",
                "total_reads_adt",
                "total_reads_adt_pf",
                "der_bad_struct_adt",
                "der_no_match_adt",
                "der_total_reads_adt",
                "der_total_reads_adt_pf",
                "der_gex_umi_sum",
                "der_unique_gene_count",
                "der_HTO_max1_v_HTO_min",
                "der_HTO_max2_v_HTO_min")
  
  for(iter in 1:length(adt_list)){
    
    #remove metadata rows
    adt_list[[iter]] <- adt_list[[iter]][!rownames(adt_list[[iter]])%in%bad_rows,]
    
    #remove ".IH" suffix
    rownames(adt_list[[iter]]) <- sub(x=rownames(adt_list[[iter]]),pattern=".IH",replacement="")
    
    #remove "adt_" prefix
    rownames(adt_list[[iter]]) <- sub(x=rownames(adt_list[[iter]]),pattern="adt_",replacement="")
    
    #remove all punctuation
    rownames(adt_list[[iter]]) <- gsub(x=rownames(adt_list[[iter]]),pattern="[[:punct:]]",replacement="")
    
    #remove "5" prefix from 5 prime samples
    if (substr(rownames(adt_list[[iter]]),1,1)==5){
      rownames(adt_list[[iter]]) <- substr(rownames(adt_list[[iter]]),2,nchar(rownames(adt_list[[iter]])))
    }
  }
  
  markers <- unique(unlist(lapply(adt_list,rownames)))
  marker_mat <- matrix(NA,nrow=length(markers),ncol=sum(unlist(lapply(adt_list,ncol))),
                       dimnames=list(markers,unlist(lapply(adt_list,colnames))))
  for(iter in 1:length(adt_list)){
    iter_markers <- intersect(markers,rownames(adt_list[[iter]]))
    marker_mat[iter_markers,colnames(adt_list[[iter]])] <- as.matrix(adt_list[[iter]][iter_markers,])
  }
  return(marker_mat)
  
}
