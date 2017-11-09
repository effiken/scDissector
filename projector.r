library(Matrix)

split_sparse=function(sparse_umitab,cell_to_cluster){
  clusters=sort(as.character(unique(cell_to_cluster)))
  l=list()
  for (cl in clusters){
    l[[cl]]=sparse_umitab[,cell_to_cluster==cl,drop=F]
  }
  
  return(l)
}

get_total_likelihood=function(ll){
  return(mean(apply(ll,1,max)))
}


get_one_likelihood=function(model_v,umitab,reg){
  return(colSums(umitab*log2(reg+model_v)))
}

getLikelihood=function(umitab,models,reg){
     return(t(umitab)%*%log2(reg+models)/colSums(umitab))
  return(res)
}

update_beta=function(umitab,models,noise_model,avg_numis_per_model,reg,max_noise_fraction=.75){
  if (nrow(models)!=length(noise_model)){
    stop("noise_models and models have different number of genes")
  }
  get_ll_b=function(numis_noise_b,models,noise_model,umitab,reg){
    alpha=pmin(numis_noise_b/avg_numis_per_model,max_noise_fraction)
    
    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
    ll_b=getLikelihood(umitab,adjusted_models,reg=reg)
    return(ll_b)
  }
  func_to_opt=function(x){
    tot_ll=get_total_likelihood(get_ll_b(x,models,noise_model,umitab,reg))
    return(tot_ll)
  }
  
  #   message("Updating noise params")

    optim_val=optimize(func_to_opt,c(0,1000),maximum = T,tol=5)
    beta_noise=optim_val$maximum
  
  return(beta_noise)
}


getOneBatchCorrectedLikelihood=function(umitab,models,noise_model,beta_noise=NULL,  avg_numis_per_model,reg,max_noise_fraction=.75){
  
  ll=matrix(NA,ncol(umitab),ncol(models))
  rownames(ll)=colnames(umitab)
  colnames(ll)=colnames(models)
  
  beta_noise=update_beta(umitab,models,noise_model,avg_numis_per_model,reg=reg)
  alpha=pmin(beta_noise/avg_numis_per_model,max_noise_fraction)
  adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
  ll[,colnames(adjusted_models)]=as.matrix(getLikelihood(umitab,adjusted_models,reg=reg))
  
  return(ll)
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
