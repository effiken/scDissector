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



load_dataset_and_model=function(model_fn,sample_fns,min_umis=250){

  model<-new.env()
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


  samples=names(sample_fns)

  tmp_dataset=new.env()
  tmp_dataset$umitab=list()
  tmp_dataset$ll=list()
  tmp_dataset$cell_to_cluster
  tmp_dataset$avg=list()
  tmp_dataset$ds=list()
  tmp_dataset$ds_numis=NULL
  tmp_dataset$counts=list()
  tmp_dataset$bulksum=list()
  tmp_dataset$insilico_gating_scores=list()
  genes=rownames(model$models)
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
  
    if (is.null(tmp_dataset$ds_numis)){
      tmp_dataset$ds_numis=tmp_env$ds_numis
    }
    else{
      if (!all(tmp_dataset$ds_numis==tmp_env$ds_numis)){
        message("Warning! Some of the samples don't share the same downsampling UMIs values")
        tmp_dataset$ds_numis=intersect(tmp_dataset$ds_numis,tmp_env$ds_numis)
      }
    }
    tmp_env$numis_before_filtering=colSums(tmp_env$umitab)    
  
    if (is.null(model$insilico_gating)){
      tmp_dataset$umitab[[sampi]]=tmp_env$umitab
    }
    else{
      is_res=insilico_sorter(tmp_env$umitab,model$insilico_gating)
      tmp_dataset$umitab[[sampi]]=is_res$umitab
      tmp_dataset$insilico_gating_scores[[sampi]]=is_res$scores
    }
    tmp_dataset$umitab[[sampi]]=tmp_dataset$umitab[[sampi]][,tmp_env$numis_before_filtering[colnames(tmp_dataset$umitab[[sampi]])]>min_umis]
    
    message("Projecting ",ncol(tmp_dataset$umitab[[sampi]])," cells")
    genemask=intersect(rownames(tmp_dataset$umitab[[sampi]]),rownames(model$models))
  
    genes=intersect(genemask,genes)
    #  print(length(genes))
    
    if (is.null(model$beta_noise)){
      tmp_dataset$ll[[sampi]]=getLikelihood(tmp_dataset$umitab[[sampi]][genemask,],models =model$models[genemask,],reg = 1e-5)#params$reg)
    }
    else {
      
      beta_noise=update_beta(tmp_dataset$umitab[[sampi]][genemask,],model$models[genemask,],noise_model,avg_numis_per_model,reg=reg,max_noise_fraction=.75){
      tmp_dataset$ll=getOneBatchCorrectedLikelihood(tmp_dataset$umitab[[sampi]][genemask,],models=model$models[genemask,],noise_model,beta_noise=beta_noise,  avg_numis_per_model,reg=reg,max_noise_fraction=.75)
    }
    tmp_dataset$cell_to_cluster[[sampi]]=MAP(tmp_dataset$ll[[sampi]])
  
    cells_to_include=names(tmp_dataset$cell_to_cluster[[sampi]])[!tmp_dataset$cell_to_cluster[[sampi]]%in%output$scDissector_params$excluded_clusters]
  
    tmp_dataset$cell_to_cluster[[sampi]]=tmp_dataset$cell_to_cluster[[sampi]][cells_to_include]
    tmp_dataset$ll[[sampi]]=tmp_dataset$ll[[sampi]][cells_to_include,]
    tmp_dataset$umitab[[sampi]]=tmp_dataset$umitab[[sampi]][,cells_to_include]
    tmp_models=model$models[,setdiff(colnames(model$models),output$scDissector_params$excluded_clusters)]
  
    tmp_dataset$avg[[sampi]]=matrix(0, nrow(tmp_dataset$umitab[[sampi]]),ncol(tmp_models),dimnames = list(rownames(tmp_dataset$umitab[[sampi]]),colnames(tmp_models)))
    tmpmod=update_models(tmp_dataset$umitab[[sampi]],tmp_dataset$cell_to_cluster[[sampi]])
    tmp_dataset$avg[[sampi]][,colnames(tmpmod)]=tmpmod
    tmp_counts=sapply(split_sparse(tmp_dataset$umitab[[sampi]][genemask,],tmp_dataset$cell_to_cluster[[sampi]]),rowSums)
    tmp_dataset$counts[[sampi]]=tmp_models*0
    tmp_dataset$counts[[sampi]][genemask,colnames(tmp_counts)]=tmp_counts
    tmp_dataset$bulksum[[sampi]]=rowSums(tmp_dataset$umitab[[sampi]][genemask,])
    tmp_dataset$numis_before_filtering[[sampi]]=tmp_env$numis_before_filtering
  
    if (length(tmp_dataset$ds)==0){
      for (ds_i in 1:length(tmp_env$ds_numis)){
        tmp_dataset$ds[[ds_i]]=list()
      }
    }
    for (ds_i in 1:length(tmp_env$ds_numis)){
      tmp_dataset$ds[[ds_i]][[sampi]]=tmp_env$ds[[ds_i]][,intersect(colnames(tmp_env$ds[[ds_i]]),cells_to_include)]
    }
    message("")
    rm("tmp_env")
  }

  message("One more minute...")

  dataset<-new.env()
  dataset$ds_numis=tmp_dataset$ds_numis
  dataset$umitab<-tmp_dataset$umitab[[samples[1]]][genes,]
  dataset$ll<-tmp_dataset$ll[[samples[1]]]
  dataset$cell_to_cluster<-tmp_dataset$cell_to_cluster[[samples[1]]]
  names(dataset$cell_to_cluster)=colnames(tmp_dataset$umitab[[samples[1]]])
  dataset$cell_to_sample<-rep(samples[1],ncol(tmp_dataset$umitab[[samples[1]]]))
  names(dataset$cell_to_sample)=colnames(tmp_dataset$umitab[[samples[1]]])
  dataset$numis_before_filtering=tmp_dataset$numis_before_filtering[[1]]
  dataset$avg<-tmp_dataset$avg
  dataset$counts<-tmp_dataset$counts
  dataset$samples=samples
  dataset$ds<-list()
  dataset$randomly_selected_cells<-list()
  dataset$bulk_avg=matrix(0,length(genes),length(samples))
  if (!is.null(model$insilico_gating)){
    for (score_i in 1:length(model$insilico_gating)){
      dataset$insilico_gating_scores[[score_i]]=tmp_dataset$insilico_gating_scores[[1]][[score_i]]
    }
  }
  rownames(dataset$bulk_avg)=genes
  colnames(dataset$bulk_avg)=samples
  dataset$bulk_avg[genes,samples[1]]=tmp_dataset$bulksum[[samples[1]]][genes]/sum(tmp_dataset$bulksum[[samples[1]]][genes])
 
  for (ds_i in 1:length(dataset$ds_numis)){
    ds_sampi=tmp_dataset$ds[[ds_i]][[samples[1]]][genes,]
    dataset$ds[[ds_i]]=ds_sampi
    dataset$randomly_selected_cells[[ds_i]]<-list()
    for (randomi in 1:length(params$nrandom_cells_per_sample_choices)){
      nrandom_cells=params$nrandom_cells_per_sample_choices[randomi]
      if (nrandom_cells=="All"||as.numeric(nrandom_cells)>=ncol(ds_sampi)){
        dataset$randomly_selected_cells[[ds_i]][[randomi]]<-colnames(ds_sampi)
      }
      else{
        dataset$randomly_selected_cells[[ds_i]][[randomi]]<-sample(colnames(ds_sampi),size=as.numeric(nrandom_cells),replace=F)
      }
    }
  }

  if (length(samples)>1){
    for (sampi in samples[-1]){
    
      dataset$bulk_avg[genes,sampi]=tmp_dataset$bulksum[[sampi]][genes]/sum(tmp_dataset$bulksum[[sampi]][genes])
      cellids=colnames(tmp_dataset$umitab[[sampi]][genes,])
      dataset$umitab<-cBind(dataset$umitab,tmp_dataset$umitab[[sampi]][genes,])
      dataset$ll<-rbind(dataset$ll,tmp_dataset$ll[[sampi]])
      dataset$cell_to_cluster<-c(dataset$cell_to_cluster,tmp_dataset$cell_to_cluster[[sampi]])
      cell_to_sampi=rep(sampi,length(cellids))
      names(cell_to_sampi)=cellids
      dataset$cell_to_sample<-c(dataset$cell_to_sample,cell_to_sampi)
      if (!is.null(model$insilico_gating)){
        for (score_i in 1:length(model$insilico_gating)){
          dataset$insilico_gating_scores[[score_i]]=c(dataset$insilico_gating_scores[[score_i]],tmp_dataset$insilico_gating_scores[[score_i]])
      }
    }
      dataset$numis_before_filtering<-c(dataset$numis_before_filtering,tmp_dataset$numis_before_filtering[[sampi]])
      for (ds_i in 1:length(dataset$ds_numis)){
        ds_sampi=tmp_dataset$ds[[ds_i]][[sampi]][genes,]
        dataset$ds[[ds_i]]=cBind(dataset$ds[[ds_i]],ds_sampi)
      
      
        for (randomi in 1:length(params$nrandom_cells_per_sample_choices)){
          nrandom_cells=params$nrandom_cells_per_sample_choices[randomi]
          if (nrandom_cells=="All"||as.numeric(nrandom_cells)>=ncol(ds_sampi)){
            dataset$randomly_selected_cells[[ds_i]][[randomi]]<-c(dataset$randomly_selected_cells[[ds_i]][[randomi]],colnames(ds_sampi))
          }
          else{
            dataset$randomly_selected_cells[[ds_i]][[randomi]]<-c(dataset$randomly_selected_cells[[ds_i]][[randomi]],sample(colnames(ds_sampi),size=as.numeric(nrandom_cells),replace=F))
          }
        }
      }
    }
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
  return(output)
}