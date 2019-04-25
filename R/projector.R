library(Matrix)
library(Matrix.utils)
library(reshape2)

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
    return(as.matrix(Matrix::t(umitab)%*%log2(reg+models))/Matrix::colSums(umitab))
}


#update_alpha_single_batch=function(umitab,models,noise_model,reg,max_noise_fraction=.75,max_ncells=5000){
#  if (ncol(umitab)>max_ncells){
#    umitab=umitab[,sample(colnames(umitab),size = max_ncells)]
#  }
#
#  if (nrow(models)!=length(noise_model)){
#    stop("noise_models and models have different number of genes")
#  }
#  get_ll_a=function(alpha,models,noise_model,umitab,reg){
#    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
#    ll_b=getLikelihood(umitab,adjusted_models,reg=reg)
#    return(ll_b)
#  }
#  func_to_opt=function(x){
#    message(round(x,digits=2))
#    tot_ll=get_total_likelihood(get_ll_a(x,models,noise_model,umitab,reg))
#    return(tot_ll)
#  }
#
#  optim_val=optimize(func_to_opt,c(1e-3,max_noise_fraction),maximum = T,tol=1e-3)
#
#  return(optim_val$maximum)
#}

update_alpha_single_batch=function(umitab,models,noise_model,cell_to_cluster=NULL,reg,max_noise_fraction=.2,max_ncells=5000){
  if (ncol(umitab)>max_ncells){
    umitab=umitab[,sample(colnames(umitab),size = max_ncells)]
  }

  if (nrow(models)!=length(noise_model)){
    stop("noise_models and models have different number of genes")
  }


  get_tot_likelihood_fixed_clusters=function(alpha){
    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
    ll=getLikelihood(umitab,adjusted_models,reg=reg)
    return(mean(ll[cbind(1:nrow(ll),match(cell_to_cluster,colnames(ll)))]))
  }

  get_tot_likelihood=function(alpha){
    adjusted_models=t((1-alpha)*t(models)+alpha*matrix(noise_model,ncol(models),nrow(models),byrow=T))
    return(get_total_likelihood(getLikelihood(umitab,adjusted_models,reg=reg)))
  }

  if (!is.null(cell_to_cluster)){
    optim_val=optimize(get_tot_likelihood_fixed_clusters,c(0,max_noise_fraction),maximum = T,tol=1e-3)
  }
  else{
    optim_val=optimize(get_tot_likelihood,c(0,max_noise_fraction),maximum = T,tol=1e-3)
  }
  return(optim_val$maximum)
}





get_expected_noise_UMI_counts_beta=function(umis,cluster,batch,noise_models,beta_noise,clusters){
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

get_expected_noise_UMI_counts_alpha=function(umis,cluster,batch,noise_models,alpha_noise,clusters){
  ngenes=nrow(noise_models)
  nmodels=length(clusters)
  nsamps=ncol(noise_models)

  raw_counts=Matrix::t(Matrix.utils::aggregate.Matrix(Matrix::t(umis),cluster))
  ag=aggregate(Matrix::colSums(umis),by=list(batch,cluster),sum)
  numis_per_batch=sapply(split(Matrix::colSums(umis),batch),sum)[colnames(noise_models)]
  numis_per_batch_cluster=matrix(0,nsamps,nmodels,dimnames = list(colnames(noise_models),clusters))
  tmp_numis_per_batch_cluster=invisible(acast(data.frame(batch=batch,cluster=cluster,numis=Matrix::colSums(umis)),batch~cluster,fun.aggregate=sum,value.var = "numis")[colnames(noise_models),colnames(raw_counts)])
  numis_per_batch_cluster[rownames(tmp_numis_per_batch_cluster),colnames(tmp_numis_per_batch_cluster)]=tmp_numis_per_batch_cluster
  tot_noise_umis=matrix(numis_per_batch_cluster*alpha_noise,nsamps,nmodels,dimnames = list(colnames(noise_models),clusters))
  arr_tot_noise_umis=array(tot_noise_umis,dim=c(nsamps,nmodels,ngenes))
  arr_tot_noise_umis=aperm(arr_tot_noise_umis,c(1,3,2))
  arr_noise_models=array(noise_models,dim=c(ngenes,nsamps,nmodels))
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
    counts=Matrix.utils::aggregate.Matrix(Matrix::t(umis),cluster)
    models=Matrix::t(counts/Matrix::rowSums(counts))
  return(as.matrix(models))
}

update_models_debatched=function(umis,cell_to_cluster,batch,noise_models,alpha_noise,make_plots=F,figure_prefix=""){
  clusters=unique(cell_to_cluster)
  clusters=as.character(clusters[order(as.numeric(clusters))])
  raw_counts=Matrix(0,nrow(umis),length(clusters),dimnames = list(rownames(umis),clusters))
  tmp_raw_counts=Matrix::t(aggregate.Matrix(Matrix::t(umis),cell_to_cluster))
  raw_counts[match(rownames(tmp_raw_counts),rownames(raw_counts)),match(colnames(tmp_raw_counts),colnames(raw_counts))]=tmp_raw_counts
  rm(tmp_raw_counts)
  gc()

  numis_per_batch_cluster=matrix(0,ncol(noise_models),length(clusters),dimnames = list(colnames(noise_models),clusters))
  tmp_numis_per_batch_cluster=invisible(acast(data.frame(batch=batch,cluster=cell_to_cluster,numis=Matrix::colSums(umis)),batch~cell_to_cluster,fun.aggregate=sum,value.var = "numis")[colnames(noise_models),colnames(raw_counts)])
  numis_per_batch_cluster[rownames(tmp_numis_per_batch_cluster),colnames(tmp_numis_per_batch_cluster)]=tmp_numis_per_batch_cluster
  rm(tmp_numis_per_batch_cluster)
  gc()
  # Substracting the expected noise UMI counts from the raw UMI counts:
  expected_noise_counts=noise_models%*%(numis_per_batch_cluster*alpha_noise)
  adj_counts=pmax(as.matrix(raw_counts-expected_noise_counts),0)
  if (make_plots){
    for (i in 1:ncol(raw_counts)){
      if (!all(expected_noise_counts==0)){
        png(paste(figure_prefix,"_",i,".png",sep=""))
        plot(expected_noise_counts[,i],raw_counts[,i],log="xy",main=i,xlab="expected noise UMIs",ylab="raw UMIs");abline(0,1)
        dev.off()
      }
    }
  }
  #identify(raw_counts[,i],expected_noise_counts[,i],labels = rownames(raw_counts))
  #browser()
  rm(expected_noise_counts)
  gc()
  models=t(t(adj_counts)/colSums(adj_counts))
  return(models)
}




