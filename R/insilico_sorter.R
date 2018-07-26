#' @export
# cell_to_batch - optional - allows reporting of #gated out cells per batch
insilico_sorter=function(umitab,insilico_gating,cell_to_batch=NULL){
  scores=list()
  gated_out_umitabs=list()
  if (!is.null(cell_to_batch)){
    tot_tab=table(cell_to_batch)
    counts=matrix(0,length(tot_tab),length(insilico_gating)+1,dimnames=list(names(tot_tab),c(names(insilico_gating),"total")))
    counts[,"total"]=tot_tab
  }
  if (!is.null(insilico_gating)){
    for (i in 1:length(insilico_gating)){
      score_i=colSums(umitab[intersect(rownames(umitab),insilico_gating[[i]]$genes),])/colSums(umitab)
      insilico_gating[[i]]$mask=names(which(score_i>=insilico_gating[[i]]$interval[1]&score_i<=insilico_gating[[i]]$interval[2]))
      if (is.null(cell_to_batch)){
        message("Gating out ",length(setdiff(names(score_i),insilico_gating[[i]]$mask))," / ",ncol(umitab)," ",names(insilico_gating)[i]," barcodes")
      }else{
        tab_gated_out=table(cell_to_batch[setdiff(names(score_i),insilico_gating[[i]]$mask)])
        counts[names(tab_gated_out),names(insilico_gating)[i]]=tab_gated_out
      }
      gated_out_umitabs[[i]]=umitab[,setdiff(colnames(umitab),insilico_gating[[i]]$mask),drop=F]
      umitab=umitab[,insilico_gating[[i]]$mask,drop=F]
      scores[[i]]=score_i
    }
    names(scores)=names(insilico_gating)
    names(gated_out_umitabs)=names(insilico_gating)
  }
  if (!is.null(cell_to_batch)){
    print(counts)
  }
  return(list(umitab=umitab,scores=scores,gated_out_umitabs=gated_out_umitabs))
  
}






