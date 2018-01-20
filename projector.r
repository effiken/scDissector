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


update_beta_single_batch=function(umitab,models,noise_model,avg_numis_per_model,reg,max_noise_fraction=.75,max_ncells=5000){
  if (ncol(umitab)>max_ncells){
    umitab=umitab[,sample(colnames(umitab),size = max_ncells)]   
  }
 
  if (nrow(models)!=length(noise_model)){
    stop("noise_models and models have different number of genes")
  }
  get_ll_b=function(numis_noise_b,models,noise_model,umitab,reg){
    alpha=pmin(numis_noise_b/avg_numis_per_model,max_noise_fraction)
    
    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
    ll_b=getLikelihood(umitab,adjusted_models,reg=reg)
    return(ll_b)
  }
  func_to_opt=function(x,bi){
    #  message(bi," ",round(x,digits=2))
    tot_ll=get_total_likelihood(get_ll_b(x,models,noise_model,umitab,reg))
    return(tot_ll)
  }
  
  optim_val=optimize(func_to_opt,c(0,1000),maximum = T,tol=10)
  
  return(optim_val$maximum)
}


getOneBatchCorrectedLikelihood=function(umitab,models,noise_model,beta_noise=NULL,  avg_numis_per_model,reg,max_noise_fraction=.75){
  
  ll=matrix(NA,ncol(umitab),ncol(models))
  rownames(ll)=colnames(umitab)
  colnames(ll)=colnames(models)
  
  ll_noise=matrix(NA,ncol(umitab),1)
  rownames(ll_noise)=colnames(umitab)
  
  alpha=pmin(beta_noise/avg_numis_per_model,max_noise_fraction)
  adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
  ll[colnames(umitab),colnames(adjusted_models)]=getLikelihood(umitab,adjusted_models,reg=reg)
  ll_noise[,1]=getLikelihood(umitab,noise_model,reg=reg)[,1]
  return(list(ll=ll,ll_noise=ll_noise))
}

get_expected_noise_UMI_counts=function(umis,cluster,batch,noise_models,beta_noise,avg_numis_per_model){
 ngenes=nrow(noise_models)
  nmodels=length(avg_numis_per_model)
  nsamps=ncol(noise_models)
  tmp_tab=table(batch,cluster)
  tab=matrix(0,nsamps,nmodels,dimnames = list(colnames(noise_models),names(avg_numis_per_model)))
  tab[,colnames(tmp_tab)]=tmp_tab

  #beta_noise is the inferred number of noise molecules/cell.
  #tab contains the number of cells per (sample,cluster)
  tot_noise_umis=matrix(beta_noise,nsamps,nmodels,dimnames = list(colnames(noise_models),names(avg_numis_per_model)))*tab
  arr_tot_noise_umis=array(tot_noise_umis,dim=c(nsamps,nmodels,ngenes))
  arr_tot_noise_umis=aperm(arr_tot_noise_umis,c(1,3,2))
  arr_noise_models=array(noise_models[,rownames(tab)],dim=c(ngenes,nsamps,nmodels))
  arr_noise_models=aperm(arr_noise_models,c(2,1,3))
  expected_noise_counts=arr_noise_models*arr_tot_noise_umis
  dimnames(expected_noise_counts)=list(colnames(noise_models),rownames(noise_models),names(avg_numis_per_model))
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




