
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




noiseEMsingleBatch=function(umitab,models,noise_model,avg_numis_per_model,reg,max_noise_fraction=.75,trace=T){
  getOneBatchCorrectedLikelihood_beta=function(umitab,models,noise_model,beta_noise=NULL,  avg_numis_per_model,reg,max_noise_fraction=.75){
  
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

  
    beta_noise=update_beta_single_batch(umitab,models,noise_model,avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
if (trace){
    message("Est 1: ~",round(beta_noise), " noise UMIs/cell")
  }
  cell_to_cluster=rep("",ncol(umitab))
  nmoved=Inf
  i=0
  min_n_moved=min(10,ncol(umitab)/100)
  while (nmoved>=min_n_moved&&i<6){
    res_boll=getOneBatchCorrectedLikelihood_beta(umitab,models=models,noise_model,beta_noise=beta_noise,  avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
    prev_cell_to_cluster=cell_to_cluster
    cell_to_cluster=MAP(res_boll$ll)
    tmptab=sapply(split((Matrix::colSums(umitab)),cell_to_cluster[colnames(umitab)]),mean)
    avg_numis_per_model[names(tmptab)]=tmptab
   
    beta_noise=update_beta_single_batch(umitab,models,noise_model,avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)

    nmoved=sum(cell_to_cluster!=prev_cell_to_cluster)
    if (i>0&trace){
      message("Est 2: ~",round(beta_noise), " noise UMIs/cell")
      message("inner iter ",i, " ",nmoved,"/",length(cell_to_cluster)," cells moved")
    }
    i=i+1
  }

  res_boll=getOneBatchCorrectedLikelihood_beta(umitab,models=models,noise_model,beta_noise=beta_noise,  avg_numis_per_model,reg=reg,max_noise_fraction=max_noise_fraction)
  
  return(list(beta_noise=beta_noise,avg_numis_per_model=avg_numis_per_model,ll=res_boll$ll))
}




