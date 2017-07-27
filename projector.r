library(Matrix)

split_sparse=function(sparse_umitab,cell_to_cluster){
  clusters=sort(as.character(unique(cell_to_cluster)))
  l=list()
  for (cl in clusters){
    l[[cl]]=sparse_umitab[,cell_to_cluster==cl,drop=F]
  }
  
  return(l)
}


get_one_likelihood=function(model_v,umitab,reg){
  return(colSums(umitab*log2(reg+model_v)))
}

getLikelihood=function(umitab,models,reg){
#  res=apply(models,2,get_one_likelihood,umitab,reg)/colSums(umitab)
#  colnames(res)=colnames(models)
#  rownames(res)=colnames(umitab)
     return(t(umitab)%*%log2(reg+models)/colSums(umitab))
  return(res)
}

MAP=function(likelihood){
  v=colnames(likelihood)[apply(likelihood,1,which.max)]
  names(v)=rownames(likelihood)
  return(v)
}

update_models=function(umis,cluster){
  counts=sapply(split_sparse(umis,cluster),rowSums)
  models=t(t(counts)/colSums(counts))
  return(models)
}


chisq_genes=function(umitab,cell_to_cluster){
  cluster_tot=sapply(split(colSums(umitab),cell_to_cluster[colnames(umitab)]),sum)
  counts=sapply(split_sparse(umitab,cell_to_cluster[colnames(umitab)]),rowSums)
  arrcont=array(c(counts,matrix(cluster_tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts),dim=c(dim(counts),2))
  res=t(apply(arrcont,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}))
  rownames(res)=rownames(counts)
  res=res[!is.na(res[,1]),]
  res=cbind(res,adjp=p.adjust(res[,1],method="BH"))
  return(res)
}
