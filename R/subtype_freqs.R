
plot_subtype_freqs=function(freq_norm,celltype,plot_legend=T,cex.names=1,cex.axis=1,cex.legend=1,cluster_set_name=""){
  alt_cols<-rep(c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#fffac8","#800000","#aaffc3","#808000","#ffd8b1","#000080","#808080",
                   "#FFFFFF","#000000"),10)
  m=freq_norm[[celltype]]
  #cols=brewer.pal(length(ordc),"Set3")
  cols=alt_cols[1:ncol(m)]
  if (plot_legend){
    
    barplot(t(m[,ncol(m):1]),col=rev(cols),las=2,cex.names = cex.names,cex.axis = cex.axis)
    mtext(cluster_set_name,side = 2,cex = cex.names,line = 4)
    plot.new()
    legend("bottomleft",legend=colnames(m),col=cols,pch=15,cex=cex.legend,bty = "n")
  }
  else{
    barplot(t(m),col=cols,las=2,cex.names=cex.names,cex.axis=1)
  }
}

get_cell_counts=function(cell_to_cluster,cell_to_sample,selected_samples){
  scRNA_tab=table(cell_to_cluster,cell_to_sample)
  scRNA_tab=scRNA_tab[,as.character(selected_samples)]
  
  scRNA_tab=t(scRNA_tab)
  return(scRNA_tab)
}

get_freqs=function(cell_to_cluster,cell_to_sample,selected_samples){
  scRNA_tab=get_cell_counts(cell_to_cluster,cell_to_sample,selected_samples)
  freqs=(scRNA_tab)/rowSums(scRNA_tab)
  
  return(freqs)
}

normalize_by_clusterset_frequency=function(cell_to_cluster,cell_to_sample,samples,cluster_sets,pool_subtype=T,reg=0){
  freqs=get_freqs(cell_to_cluster,cell_to_sample,samples)
  
  pool_subtype_freqs=function(one_subtype){
    return(rowSums(freqs[,unlist(one_subtype),drop=F]))
  }
  
  norm_one_clusterset=function(one_clusterset){
    if (length(intersect(colnames(freqs),unlist(one_clusterset)))==0){
      return(NULL)
    }
    if (pool_subtype==F){
      tot=rowSums(reg+freqs[,unlist(one_clusterset),drop=F])
      norm_subtypes_freqs=(reg+freqs[,unlist(one_clusterset),drop=F])/(tot)
    } 
    else {
      if (length(unlist(one_clusterset))==1){
        norm_subtypes_freqs=reg+freqs[,unlist(one_clusterset),drop=F]
      }
      else {
        subtypes_freqs=reg+sapply(one_clusterset,pool_subtype_freqs)
        tot=rowSums(subtypes_freqs)
        norm_subtypes_freqs=subtypes_freqs/(tot)
      }
      colnames(norm_subtypes_freqs)=names(one_clusterset)
    }
    return(norm_subtypes_freqs)
  }
  return(sapply(cluster_sets,norm_one_clusterset,simplify = F))
  
}

