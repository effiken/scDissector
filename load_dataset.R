library(Matrix)
library(Matrix.utils)


# cell_to_batch - optional - allows reporting of #gated out cells per batch
insilico_sorter=function(umitab,insilico_gating,cell_to_batch=NULL){
  scores=list()
  gated_out_umitabs=list()
  if (!is.null(insilico_gating)){
    for (i in 1:length(insilico_gating)){
      score_i=colSums(umitab[intersect(rownames(umitab),insilico_gating[[i]]$genes),])/colSums(umitab)
      insilico_gating[[i]]$mask=names(which(score_i>=insilico_gating[[i]]$interval[1]&score_i<=insilico_gating[[i]]$interval[2]))
      if (is.null(cell_to_batch)){
        message("Gating out ",length(setdiff(names(score_i),insilico_gating[[i]]$mask))," / ",ncol(umitab)," ",names(insilico_gating)[i]," barcodes")
      }else{
        tab_gated_out=table(cell_to_batch[setdiff(names(score_i),insilico_gating[[i]]$mask)])
        tab_total=table(cell_to_batch[colnames(umitab)])
        message("#Gated-out ",names(insilico_gating)[i]," barcodes:")
        print(cbind(tab_gated_out,tab_total[names(tab_gated_out)]))
      }
      gated_out_umitabs[[i]]=umitab[,setdiff(colnames(umitab),insilico_gating[[i]]$mask),drop=F]
      umitab=umitab[,insilico_gating[[i]]$mask,drop=F]
      scores[[i]]=score_i
    }
    names(scores)=names(insilico_gating)
    names(gated_out_umitabs)=names(insilico_gating)
  }
  return(list(umitab=umitab,scores=scores,gated_out_umitabs=gated_out_umitabs))
  
}



load_dataset_and_model=function(model_fn,sample_fns,min_umis=250,model_version_name="",max_umis=25000){

  if(is.null(names(sample_fns))){
    names(sample_fns)=sapply(strsplit(sapply(strsplit(sample_fns,"/"),tail,1),"\\_|\\."),function(x){paste(x[c(-1,-length(x))],collapse="_")})
  }
  
  model<-new.env()
  message("Loading model ",model_fn)
 
  load(file=model_fn,model)
  fn_prefix=strsplit(model_fn,"\\.")[[1]][1]
  output<-list()
  output$scDissector_params<-list()
  output$scDissector_params$excluded_clusters<-c()

  clusterset_fn=paste(fn_prefix,"_clustersets.txt",sep="")
  if (file.exists(clusterset_fn)){
    cluster_sets_tab=read.table(file=clusterset_fn,header = T,stringsAsFactors = F)
    if (nrow(cluster_sets_tab)>0){
      output$scDissector_params$cluster_sets<-strsplit(cluster_sets_tab[,2],",")
      output$scDissector_params$excluded_cluster_sets<-cluster_sets_tab[cluster_sets_tab[,3],1]
      names(output$scDissector_params$cluster_sets)<-cluster_sets_tab[,1]
      if (length(output$scDissector_params$excluded_cluster_sets)>0){
        output$scDissector_params$excluded_clusters<-unlist(output$scDissector_params$cluster_sets[[output$scDissector_params$excluded_cluster_sets]])
      }
    }
  }
  else{
    output$scDissector_params$cluster_sets<-list()
  }

  annnot_fn=paste(fn_prefix,"_annots.txt",sep="")
  if (file.exists(annnot_fn)){
    a=read.delim(annnot_fn,header=F,stringsAsFactors = F,row.names = 1)
    clustAnnots<-a[,1]
    clustAnnots[is.na(clustAnnots)]<-""
    names(clustAnnots)<-rownames(a)

  #  if (ncol(annot_tab)>1){
  #    output$cluster_sets<-split(rownames(annots_tab),annots_tab[,2])
  #  }

  }
  else{
    clustAnnots<-rep("",ncol(model$models))
    names(clustAnnots)<-colnames(model$models)
  }    
  output$clustAnnots=clustAnnots

  order_fn=paste(fn_prefix,"_order.txt",sep="")
  if (file.exists(order_fn)){
    cluster_order=as.character(read.table(file=order_fn,header = F)[,1])
  }
  else{
    cluster_order=colnames(model$models)
  }

  cluster_order=cluster_order[!cluster_order%in%output$scDissector_params$excluded_clusters]
  if (is.null(model$avg_numis_per_model)){
    model$avg_numis_per_model=rep(mean(colSums(model$umitab)),ncol(model$models))
    names(model$avg_numis_per_model)=colnames(model$models)
    tmptab=sapply(split(colSums(model$umitab),model$cell_to_cluster[colnames(model$umitab)]),mean)
    model$avg_numis_per_model[names(tmptab)]=tmptab
  }

  
  samples=names(sample_fns)

  dataset<-new.env()
  dataset$ds=list()
  
  dataset$ds_numis=NULL
  dataset$ll<-c()
  dataset$cell_to_cluster<-c()
  dataset$numis_before_filtering=list()
  dataset$cell_to_sample<-c()
  
  dataset$alpha_noise=rep(NA,length(samples))
  dataset$beta_noise=rep(NA,length(samples))
  dataset$avg_numis_per_sample_model<-matrix(NA,length(samples),ncol(model$models),dimnames = list(samples,colnames(model$models)))
  names(dataset$alpha_noise)=samples
  names(dataset$beta_noise)=samples
  genes=rownames(model$models)
  message("")
  dataset$umitab<-Matrix(,length(genes),,dimnames = list(genes,NA))
  dataset$gated_out_umitabs<-list()
  dataset$noise_models=matrix(0,length(genes),length(samples),dimnames = list(genes,samples))
  dataset$counts<-array(0,dim=c(length(samples),nrow(model$models),ncol(model$models)),dimnames = list(samples,rownames(model$models),colnames(model$models)))
  
  
  
  for (sampi in samples){
    message("Loading sample ",sampi)
    i=match(sampi,samples)
    
    
    tmp_env=new.env()
    load(sample_fns[i],envir = tmp_env)
    tmp_env$umitab=tmp_env$umitab[,setdiff(colnames(tmp_env$umitab),tmp_env$noise_barcodes)]
  
    colnames(tmp_env$umitab)=paste(sampi,colnames(tmp_env$umitab),sep="_")
    for (ds_i in 1:length(tmp_env$ds)){
      tmp_env$ds[[ds_i]]=tmp_env$ds[[ds_i]][,setdiff(colnames(tmp_env$ds[[ds_i]]),tmp_env$noise_barcodes)]
      colnames(tmp_env$ds[[ds_i]])=paste(sampi,colnames(tmp_env$ds[[ds_i]]),sep="_")
    }
    if (sampi==samples[1]){
      for (ds_i in 1:length(tmp_env$ds_numis)){
        dataset$ds[[ds_i]]<-Matrix(,length(genes),,dimnames = list(genes,NA))
      }
     
    }
    if (is.null(dataset$ds_numis)){
      dataset$ds_numis=tmp_env$ds_numis
    }
    else{
      if (!all(dataset$ds_numis==tmp_env$ds_numis)){
        message("Warning! Some of the samples don't share the same downsampling UMIs values")
        dataset$ds_numis=intersect(dataset$ds_numis,tmp_env$ds_numis)
      }
    }
    tmp_env$numis_before_filtering=colSums(tmp_env$umitab)    

    if (is.null(model$insilico_gating)){
      umitab=tmp_env$umitab
    }
    else{
      
      is_res=insilico_sorter(tmp_env$umitab,model$insilico_gating)
      umitab=is_res$umitab
     
      for (score_i in names(model$insilico_gating)){
        dataset$insilico_gating_scores[[score_i]]=c(dataset$insilico_gating_scores[[score_i]],is_res$scores[[score_i]])
        if (is.null(dataset$gated_out_umitabs[[score_i]])){
          dataset$gated_out_umitabs[[score_i]]=list()
        }
        dataset$gated_out_umitabs[[score_i]][[sampi]]=is_res$gated_out_umitabs[[score_i]]
      }
     
     
      
    }
    barcode_mask=tmp_env$numis_before_filtering[colnames(umitab)]>min_umis&tmp_env$numis_before_filtering[colnames(umitab)]<max_umis
    dataset$min_umis=min_umis
    dataset$max_umis=max_umis
    umitab=umitab[,barcode_mask]
    dataset$noise_models[genes,sampi]=tmp_env$noise_model[genes,1]
    genemask=intersect(rownames(umitab),rownames(model$models))
    projection_genemask=setdiff(genemask,model$params$genes_excluded)   
    noise_model=dataset$noise_models[projection_genemask,sampi]
    noise_model=noise_model/sum(noise_model)
    message("Projecting ",ncol(umitab)," cells")
    genes=intersect(rownames(umitab),genes)
    if (is.null(model$alpha_noise)&is.null(model$beta_noise)){
      ll=getLikelihood(umitab[projection_genemask,],models =model$models[projection_genemask,],reg = model$params$reg)
    }
    else {
      if (!is.null(model$beta_noise)){
        avg_numis_per_model=model$avg_numis_per_model
        gobclle_res=betaNoiseEMsingleBatch(umitab=umitab[projection_genemask,],models=model$models[projection_genemask,],noise_model=noise_model,avg_numis_per_model=avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
        ll=gobclle_res$ll
        dataset$beta_noise[sampi]=gobclle_res$beta_noise
        dataset$avg_numis_per_sample_model[sampi,names(gobclle_res$avg_numis_per_model)]=gobclle_res$avg_numis_per_model
      }
      else if (!is.null(model$alpha_noise)){
        alpha_b=update_alpha_single_batch(umitab[projection_genemask,],model$models[projection_genemask,],noise_model,reg=model$params$reg)
        message("%Noise = ",round(100*alpha_b,digits=2))
        res_l=getOneBatchCorrectedLikelihood(umitab=umitab[projection_genemask,],model$models[projection_genemask,],noise_model,alpha_noise=alpha_b,reg=model$params$reg)

         ll=res_l$ll
        dataset$alpha_noise[sampi]=alpha_b
      }
      else {
        error("Noise parameter does not exist!")
      }
    }
    
    cell_to_cluster=MAP(ll)
    cells_to_include=names(cell_to_cluster)[!cell_to_cluster%in%output$scDissector_params$excluded_clusters]
  
    cell_to_cluster=cell_to_cluster[cells_to_include]
    ll=ll[cells_to_include,]
    dataset$ll<-rbind(dataset$ll,ll)
    umitab=umitab[,cells_to_include]
    tmp_models=model$models[,setdiff(colnames(model$models),output$scDissector_params$excluded_clusters)]
  
    dataset$cell_to_cluster<-c(dataset$cell_to_cluster,cell_to_cluster)
    tmpmod=update_models(umitab,cell_to_cluster)
   
    tmp_counts=as.matrix(t(aggregate.Matrix(t(umitab[genemask,]),cell_to_cluster,fun="sum")))
    dataset$counts[sampi,rownames(tmp_counts),colnames(tmp_counts)]=tmp_counts
    dataset$numis_before_filtering[[sampi]]=tmp_env$numis_before_filtering
  
    for (ds_i in 1:length(tmp_env$ds_numis)){
      dataset$ds[[ds_i]]<-cBind(dataset$ds[[ds_i]][genes,],tmp_env$ds[[ds_i]][genes,intersect(colnames(tmp_env$ds[[ds_i]]),cells_to_include)])
     }
    dataset$umitab<-cBind(dataset$umitab[genes,],umitab[genes,])
    cellids=colnames(umitab)
    
    cell_to_sampi=rep(sampi,length(cellids))
    names(cell_to_sampi)=cellids
    dataset$cell_to_sample<-c(dataset$cell_to_sample,cell_to_sampi)
    message("")
    rm(list=c("ll","cell_to_cluster","tmp_env","umitab"))

  }

  message("...")


  dataset$umitab<-dataset$umitab[,-1]
  for (ds_i in 1:length(dataset$ds)){
    dataset$ds[[ds_i]]=dataset$ds[[ds_i]][,-1]
  }
  dataset$samples=samples
 

  if (!is.null(model$noise_models)){
    if (!is.null(dataset$beta_noise)){
      dataset$noise_counts=get_expected_noise_UMI_counts_beta(dataset$umitab,dataset$cell_to_cluster,dataset$cell_to_sample,dataset$noise_models,dataset$beta_noise,colnames(model$models))
    }
    if (!is.null(dataset$alpha_noise)){
      dataset$noise_counts=get_expected_noise_UMI_counts_alpha(dataset$umitab,dataset$cell_to_cluster,dataset$cell_to_sample,dataset$noise_models,dataset$alpha_noise,colnames(model$models))
    }
  }
  
  if (all(is.na(dataset$beta_noise))){
    dataset$beta_noise=NULL
  }
  if (all(is.na(dataset$alpha_noise))){
    dataset$alpha_noise=NULL
  }
  model$model_filename=model_fn
  output$dataset=dataset
  rm("dataset")
  ncells_per_cluster<-rep(0,dim(model$models)[2])
  names(ncells_per_cluster)<-colnames(model$models)
  temptab=table(model$cell_to_cluster)
  temptab=temptab[setdiff(names(temptab),output$scDissector_params$excluded_cluster_sets)]
  ncells_per_cluster[names(temptab)]<-temptab
  output$ncells_per_cluster=ncells_per_cluster
  output$model=model
  output$cluster_order<-cluster_order
  output$default_clusters<-cluster_order
  output$loaded_model_version<-model_version_name
  return(output)
  
}






