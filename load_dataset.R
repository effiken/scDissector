

insilico_sorter=function(umitab,insilico_gating){
  scores=list()
  if (!is.null(insilico_gating)){
    for (i in 1:length(insilico_gating)){
      score_i=colSums(umitab[intersect(rownames(umitab),insilico_gating[[i]]$genes),])/colSums(umitab)
      insilico_gating[[i]]$mask=names(which(score_i>=insilico_gating[[i]]$interval[1]&score_i<=insilico_gating[[i]]$interval[2]))
      message("Gating out ",length(setdiff(names(score_i),insilico_gating[[i]]$mask))," / ",ncol(umitab)," ",names(insilico_gating)[i]," barcodes")
      umitab=umitab[,insilico_gating[[i]]$mask]
      scores[[i]]=score_i
    }
  }
  return(list(umitab=umitab,scores=scores))
}



load_dataset_and_model=function(model_fn,sample_fns,min_umis=250,model_version_name="",max_umis=25000){

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
    tmptab=sapply(split(colSums(model$umitab),model$cell_to_cluster[colnames(model$umitab)]),mean)
    model$avg_numis_per_model[names(tmptab)]=tmptab
  }

  samples=names(sample_fns)

  tmp_dataset=new.env()
  dataset<-new.env()
  tmp_dataset$ds=list()
  dataset$ds=list()
  
  dataset$ds_numis=NULL
  dataset$ll<-c()
  dataset$cell_to_cluster<-c()
  dataset$numis_before_filtering=list()
  dataset$cell_to_sample<-c()
  tmp_dataset$counts=list()
  tmp_dataset$insilico_gating_scores=list()
  tmp_dataset$noise_models=list()
  tmp_dataset$beta_noise=list()
  genes=rownames(model$models)
  message("")
  dataset$umitab<-Matrix(,length(genes),,dimnames = list(genes,NA))
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
      tmp_dataset$insilico_gating_scores[[sampi]]=is_res$scores
    }
    barcode_mask=tmp_env$numis_before_filtering[colnames(umitab)]>min_umis&tmp_env$numis_before_filtering[colnames(umitab)]<max_umis
    dataset$min_umis=min_umis
    dataset$max_umis=max_umis
    umitab=umitab[,barcode_mask]
    tmp_dataset$noise_models[[sampi]]=tmp_env$noise_model

    message("Projecting ",ncol(umitab)," cells")
    genemask=intersect(rownames(umitab),rownames(model$models))
    projection_genemask=setdiff(genemask,model$params$genes_excluded)   
    genes=intersect(rownames(umitab),genes)
   
    if (is.null(model$beta_noise)){
      ll=getLikelihood(umitab[projection_genemask,],models =model$models[projection_genemask,],reg = model$params$reg)#params$reg)
    }
    else 
      {
        noise_model=tmp_dataset$noise_models[[sampi]]
        
        avg_numis_per_model=model$avg_numis_per_model
        
        beta_noise=update_beta_single_batch(umitab[genemask,],model$models[genemask,],noise_model[genemask,],avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
        cell_to_cluster=rep("",ncol(umitab))
        nmoved=Inf
        i=0
        while (nmoved>=10&&i<4){
          res_boll=getOneBatchCorrectedLikelihood(umitab[genemask,],models=model$models[genemask,],noise_model[genemask,],beta_noise=beta_noise,  avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
          
          prev_cell_to_cluster=cell_to_cluster
          cell_to_cluster=MAP(res_boll$ll)
          tmptab=sapply(split(colSums(umitab[genemask,]),cell_to_cluster[colnames(umitab)]),mean)
          avg_numis_per_model[names(tmptab)]=tmptab
          beta_noise=update_beta_single_batch(umitab[genemask,],model$models[genemask,],noise_model[genemask,],avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
          
          nmoved=sum(cell_to_cluster!=prev_cell_to_cluster)
          if (i>0){
            message("~",round(beta_noise), " noise UMIs/cell")
            message("iter ",i, " ",nmoved,"/",length(cell_to_cluster)," cells moved")
          }
          i=i+1
        }
        
        model$avg_numis_per_model=avg_numis_per_model
        tmp_dataset$beta_noise[[sampi]]=beta_noise
        res_boll=getOneBatchCorrectedLikelihood(umitab[genemask,],models=model$models[genemask,],noise_model[genemask,],beta_noise=beta_noise,  avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
        ll=res_boll$ll
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
    
    tmp_dataset$counts[[sampi]]=tmp_models*0
    tmp_dataset$counts[[sampi]][genemask,colnames(tmp_counts)]=tmp_counts
    
    dataset$numis_before_filtering[[sampi]]=tmp_env$numis_before_filtering
  
    if (length(tmp_dataset$ds)==0){
      for (ds_i in 1:length(tmp_env$ds_numis)){
        tmp_dataset$ds[[ds_i]]=list()
      }
    }
  
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



  dataset$samples=samples
  dataset$randomly_selected_cells<-list()
  dataset$bulk_avg=matrix(0,length(genes),length(samples))
  dataset$noise_models=matrix(0,length(genes),length(samples))
  colnames(dataset$noise_models)=names(tmp_dataset$noise_models)
  rownames(dataset$noise_models)=genes
  dataset$beta_noise=unlist(tmp_dataset$beta_noise)
  names(dataset$beta_noise)=names(tmp_dataset$beta_noise)
  if (!is.null(model$insilico_gating)){
    for (score_i in 1:length(model$insilico_gating)){
      dataset$insilico_gating_scores[[score_i]]=tmp_dataset$insilico_gating_scores[[1]][[score_i]]
    }
  }
  rownames(dataset$bulk_avg)=genes
  colnames(dataset$bulk_avg)=samples
  
  for (sampi in samples){
    dataset$noise_models[genes,sampi]=tmp_dataset$noise_models[[sampi]][genes,1]
    
    if (!is.null(model$insilico_gating)){
      for (score_i in 1:length(model$insilico_gating)){
        dataset$insilico_gating_scores[[score_i]]=c(dataset$insilico_gating_scores[[score_i]],tmp_dataset$insilico_gating_scores[[sampi]][[score_i]])
      }
    }
  }
  
  dataset$counts<-array(0,dim=c(length(tmp_dataset$counts),nrow(tmp_dataset$counts[[1]]),ncol(tmp_dataset$counts[[1]])),dimnames = list(names(tmp_dataset$counts),rownames(tmp_dataset$counts[[1]]),colnames(tmp_dataset$counts[[1]])))
  for (si in 1:length(tmp_dataset$counts)){
    dataset$counts[si,,]=tmp_dataset$counts[[si]]
  }

  if (!is.null(model$noise_models)){
    dataset$noise_counts=get_expected_noise_UMI_counts(dataset$umitab,dataset$cell_to_cluster,dataset$cell_to_sample,dataset$noise_models,dataset$beta_noise,model$avg_numis_per_model)
  }
  output$dataset=dataset
  rm("tmp_dataset","dataset")
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






