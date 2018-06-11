library(Matrix)
library(Matrix.utils)


get_total_likelihood=function(ll){
  return(mean(apply(ll,1,max)))
}


# get_one_likelihood
#
# returns the log likelihood of all the cells to a single model
# model_v - a probabiliy vector
# umitab
# reg - regularization parameter

get_one_likelihood=function(model_v,umitab,reg){
  return(colSums(umitab*log2(reg+model_v)))
}


getLikelihood=function(umitab,models,reg){
  return(as.matrix(t(umitab)%*%log2(reg+models)/colSums(umitab)))
}


update_alpha_single_batch=function(umitab,models,noise_model,reg,max_noise_fraction=.75,max_ncells=5000){
  if (ncol(umitab)>max_ncells){
    umitab=umitab[,sample(colnames(umitab),size = max_ncells)]   
  }
 
  if (nrow(models)!=length(noise_model)){
    stop("noise_models and models have different number of genes")
  }
  get_ll_b=function(alpha,models,noise_model,umitab,reg){
    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
    ll_b=getLikelihood(umitab,adjusted_models,reg=reg)
    return(ll_b)
  }
  func_to_opt=function(x,bi){
    #  message(bi," ",round(x,digits=2))
    
    tot_ll=get_total_likelihood(get_ll_a(x,models,noise_model,umitab,reg))
    return(tot_ll)
  }
  
  optim_val=optimize(func_to_opt,c(0,max_noise_fraction),maximum = T,tol=10)
  
  return(optim_val$maximum)
}


getOneBatchCorrectedLikelihood=function(umitab,models,noise_model,alpha_noise=NULL,reg,max_noise_fraction=.75){
  
  ll=matrix(NA,ncol(umitab),ncol(models))
  rownames(ll)=colnames(umitab)
  colnames(ll)=colnames(models)
  
  ll_noise=matrix(NA,ncol(umitab),1)
  rownames(ll_noise)=colnames(umitab)
  
 
  adjusted_models=t((1-alpha_noise)*t(models)+alpha_noise*matrix(noise_model,ncol(models),nrow(models),byrow=T))
  ll[colnames(umitab),colnames(adjusted_models)]=getLikelihood(umitab,adjusted_models,reg=reg)
  ll_noise[,1]=getLikelihood(umitab,noise_model,reg=reg)[,1]
  return(list(ll=ll,ll_noise=ll_noise))
}


noiseEMsingleBatch=function(umitab,models,noise_model,avg_numis_per_model,reg,max_noise_fraction=.75,trace=T){
  
    beta_noise=update_beta_single_batch(umitab,models,noise_model,avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
if (trace){
    message("Est 1: ~",round(beta_noise), " noise UMIs/cell")
  }
  cell_to_cluster=rep("",ncol(umitab))
  nmoved=Inf
  i=0
  min_n_moved=min(10,ncol(umitab)/100)
  while (nmoved>=min_n_moved&&i<6){
    res_boll=getOneBatchCorrectedLikelihood(umitab,models=models,noise_model,beta_noise=beta_noise,  avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
    prev_cell_to_cluster=cell_to_cluster
    cell_to_cluster=MAP(res_boll$ll)
    tmptab=sapply(split(colSums(umitab),cell_to_cluster[colnames(umitab)]),mean)
    avg_numis_per_model[names(tmptab)]=tmptab
   
    beta_noise=update_beta_single_batch(umitab,models,noise_model,avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)

    nmoved=sum(cell_to_cluster!=prev_cell_to_cluster)
    if (i>0&trace){
      message("Est 2: ~",round(beta_noise), " noise UMIs/cell")
      message("inner iter ",i, " ",nmoved,"/",length(cell_to_cluster)," cells moved")
    }
    i=i+1
  }

  res_boll=getOneBatchCorrectedLikelihood(umitab,models=models,noise_model,beta_noise=beta_noise,  avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
  
  return(list(beta_noise=beta_noise,avg_numis_per_model=avg_numis_per_model,ll=res_boll$ll))
}




get_expected_noise_UMI_counts=function(umis,cluster,batch,noise_models,beta_noise,clusters){
 ngenes=nrow(noise_models)
  nmodels=length(clusters)
  nsamps=ncol(noise_models)
  tmp_tab=table(batch,cluster)
  ncells=matrix(0,nsamps,nmodels,dimnames = list(colnames(noise_models),clusters))
  ncells[rownames(tmp_tab),colnames(tmp_tab)]=tmp_tab

  #beta_noise is the inferred number of noise molecules/cell.
  #tab contains the number of cells per (sample,cluster)
  tot_noise_umis=matrix(beta_noise[colnames(noise_models)],nsamps,nmodels,dimnames = list(colnames(noise_models),clusters))*ncells
  arr_tot_noise_umis=array(tot_noise_umis,dim=c(nsamps,nmodels,ngenes))
  arr_tot_noise_umis=aperm(arr_tot_noise_umis,c(1,3,2))
  arr_noise_models=array(noise_models[,rownames(ncells)],dim=c(ngenes,nsamps,nmodels))
  arr_noise_models=aperm(arr_noise_models,c(2,1,3))
  expected_noise_counts=arr_noise_models*arr_tot_noise_umis
  dimnames(expected_noise_counts)=list(colnames(noise_models),rownames(noise_models),clusters)
  return(expected_noise_counts)
}




MAP=function(likelihood){
  v=colnames(likelihood)[apply(likelihood,1,which.max)]
  names(v)=rownames(likelihood)
  return(v)
}

update_models=function(umis,cluster){
  counts=aggregate(t(umis),cluster)
  models=t(counts/rowSums(counts))
  return(as.matrix(models))
}

update_models_debatched=function(umis,cluster,batch,noise_models,alpha_noise){
  raw_counts=t(aggregate(t(umis),cluster))
  numis_per_batch=sapply(split(colSums(umis),batch),sum)[colnames(noise_models)]

  expected_noise_counts=noise_models%*%((numis_per_batch*alpha_noise))
  adj_counts=pmax(raw_counts-expected_noise_counts,0)
  
  models=t(t(adj_counts)/colSums(adj_counts))
  return(as.matrix(models))
}




