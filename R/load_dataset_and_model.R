#' @export


get_cluster_set_tree=function(mat,nodes_to_add=NULL){
  
  if (is.null(nodes_to_add)){
    nodes_to_add=setdiff(mat$parent,mat$node)
  }
  
  tr=list()  
  for (node in nodes_to_add){
    if (any(mat$parent==node)){
      tr[[node]]= get_cluster_set_tree(mat,mat$node[mat$parent==node])
      if (!is.null(tr[[node]])){
        names(tr)[length(tr)]=node
      }
    }
    else{
      tr[[length(tr)+1]]=node
      names(tr)[length(tr)]=node
    }
  }
  
  #  names(tr)=nodes_to_add
  return(tr)
}

#' Load dataset and project onto a predefined model
#'  
#' @param model_fn model file-name
#' @param sample_fns vector containing compiled sample file-names 
#' @param min_umis lower UMI-count threshold 
#' @param model_version_name model version name
#' @param max_umis upper UMI-count threshold 
#' @param excluded_clusters [optional] vector with clusters to exclude. Cells associated with these clusters are not loaded.
#' @param ds_numis [optional] vector containing UMI counts to down-sample the dataset to. 
#' @param genes [optional ] vector specifying genes which will be included in the analysis. If NULL (default) all genes are included. In cases of inconsistencies between the gene lists included in different samples this vector should provide the desired common gene list.  
#' @param max_ncells_per_sample [optional ] Integer. limits the number of cells randomly selected and loaded per sample. If NA (default) all cells are loaded
#' @param lightweight [optional ] Boolean. If true, unnecessary obectes are not loaded. (F is the default) 
#' @return LDM object
#' @export
load_dataset_and_model<-function(model_fn,sample_fns,min_umis=250,model_version_name="",max_umis=25000,excluded_clusters=NA,ds_numis=NA,genes=NULL,max_ncells_per_sample=NA,lightweight=F){
  if (all(is.na(excluded_clusters))){
    excluded_clusters=c()
  }
    require(Matrix)
    require(Matrix.utils)
    if(is.null(names(sample_fns))){
        names(sample_fns)=sapply(strsplit(sapply(strsplit(sample_fns,"/"),tail,1),"\\_|\\."),function(x){paste(x[c(-1,-length(x))],collapse="_")})
    }
    
    model<-new.env()
    message("Loading model ",model_fn)
    
    load(file=model_fn,model)
    fn_prefix=strsplit(model_fn,"\\.")[[1]][1]
    output<-list()
    output$scDissector_params<-list()

    
    
    
    get_cluster_set_list=function(clusts,i){
      if (i==0){
        return(clusts)
      }
      else{
        l=split(clusts,a[clusts,i])
        return(lapply(l,get_cluster_set_list,i-1))
      }
    }
    
    
    
    annnot_fn=paste(fn_prefix,"_annots.txt",sep="")
    if (file.exists(annnot_fn)){
        a=read.delim(annnot_fn,header=F,stringsAsFactors = F,row.names = 1)
        clustAnnots<-a[,1]
        clustAnnots[is.na(clustAnnots)]<-""
        names(clustAnnots)<-rownames(a)
        
    #     if (ncol(a)>1){
    #        output$cluster_sets<-get_cluster_set_list(rownames(a),ncol(a))
    #      }
        
    }
    else{
        clustAnnots<-rep("unannotated",ncol(model$models))
        names(clustAnnots)<-colnames(model$models)
    }
    output$clustAnnots=clustAnnots
    
    cluster_sets_fn=paste(fn_prefix,"_cluster_sets.txt",sep="")
    if (file.exists(cluster_sets_fn)){
      a=read.delim(cluster_sets_fn,header=T,stringsAsFactors = F)
      
      output$cluster_sets<-get_cluster_set_tree(a)
      
    }
    else{
      nodes=names(clustAnnots)
      parents=as.character(clustAnnots)
      nodes[nodes==""]="unannotated"
      parents[parents==""]="unannotated"
      output$cluster_sets<-get_cluster_set_tree(data.frame(node=nodes,parent=parents))
    }
    
    
    order_fn=paste(fn_prefix,"_order.txt",sep="")
    if (file.exists(order_fn)){
        cluster_order=as.character(read.table(file=order_fn,header = F)[,1])
    }
    else{
        cluster_order=colnames(model$models)
    }
    

    if (is.null(model$avg_numis_per_model)){
        model$avg_numis_per_model=rep(mean(Matrix::colSums(model$umitab)),ncol(model$models))
        names(model$avg_numis_per_model)=colnames(model$models)
        tmptab=sapply(split(Matrix::colSums(model$umitab),model$cell_to_cluster[colnames(model$umitab)]),mean)
        model$avg_numis_per_model[names(tmptab)]=tmptab
    }
    
    if (lightweight){
      model$umitab=NULL
    }
    
 
    samples=names(sample_fns)
    
    dataset<-new.env()
    dataset$ds=list()
    
    dataset$ds_numis=NULL
    dataset$ll<-c()
    dataset$cell_to_cluster<-c()
    dataset$numis_before_filtering=list()
    dataset$cell_to_sample<-c()
    
    dataset$adt_by_sample=list()
    dataset$hto_by_sample=list()
    dataset$alpha_noise=rep(NA,length(samples))

    dataset$avg_numis_per_sample_model<-matrix(NA,length(samples),ncol(model$models),dimnames = list(samples,colnames(model$models)))
    names(dataset$alpha_noise)=samples
    if (is.null(genes)){
      genes=rownames(model$models)
    }
    else{
      genes=intersect(genes,rownames(model$models))
    }
    message("")
    dataset$umitab<-Matrix(,length(genes),,dimnames = list(genes,NA))
    dataset$gated_out_umitabs<-list()
    dataset$noise_models=matrix(0,length(genes),length(samples),dimnames = list(genes,samples))
    dataset$counts<-array(0,dim=c(length(samples),nrow(model$models),ncol(model$models)),dimnames = list(samples,rownames(model$models),colnames(model$models)))
    
    if (!is.null(model$alpha_noise)){
      init_alpha=rep(NA,length(samples))
      names(init_alpha)=samples
      init_alpha[names(model$alpha_noise)]=model$alpha_noise
    
      if (any(is.na(init_alpha))){
        init_alpha[is.na(init_alpha)]=median(model$alpha_noise,na.rm=T)
      }
    }
    for (sampi in samples){
        message("Loading sample ",sampi)
        i=match(sampi,samples)
      
        
        tmp_env=new.env()
        load(sample_fns[i],envir = tmp_env)
        tmp_env$umitab=tmp_env$umitab[,setdiff(colnames(tmp_env$umitab),tmp_env$noise_barcodes)]
        
        colnames(tmp_env$umitab)=paste(sampi,colnames(tmp_env$umitab),sep="_")
        if ((!is.na(max_ncells_per_sample))&ncol(tmp_env$umitab)>max_ncells_per_sample){
          tmp_env$umitab=tmp_env$umitab[,sample(colnames(tmp_env$umitab),size = max_ncells_per_sample,replace = F)]
        }
        
        if (!is.na(ds_numis)){
          ds_numis_sampi=intersect(as.character(ds_numis),as.character(tmp_env$ds_numis))
        }
        else{
          ds_numis_sampi=tmp_env$ds_numis
        }
        for (ds_i in 1:length(ds_numis_sampi)){
            tmp_env$ds[[ds_i]]=tmp_env$ds[[ds_i]][,setdiff(colnames(tmp_env$ds[[ds_i]]),tmp_env$noise_barcodes)]
            colnames(tmp_env$ds[[ds_i]])=paste(sampi,colnames(tmp_env$ds[[ds_i]]),sep="_")
        }
        if (sampi==samples[1]){
            for (ds_i in 1:length(ds_numis_sampi)){
                dataset$ds[[ds_i]]<-Matrix(,length(genes),,dimnames = list(genes,NA))
            }
            
        }
        if (is.null(dataset$ds_numis)){
            dataset$ds_numis=ds_numis_sampi
        }
        else{
            if (!all(dataset$ds_numis==ds_numis_sampi)){
                message("Warning! Some of the samples don't share the same downsampling UMIs values")
                dataset$ds_numis=intersect(dataset$ds_numis,ds_numis_sampi)
            }
        }
        tmp_env$numis_before_filtering=Matrix::colSums(tmp_env$umitab)
        
        
        if (is.null(model$params$insilico_gating)){
            umitab=tmp_env$umitab
        }
        else{
            
            is_res=insilico_sorter(tmp_env$umitab,model$params$insilico_gating)
            umitab=is_res$umitab
            
            for (score_i in names(model$params$insilico_gating)){
                dataset$insilico_gating_scores[[score_i]]=c(dataset$insilico_gating_scores[[score_i]],is_res$scores[[score_i]])
                if (is.null(dataset$gated_out_umitabs[[score_i]])){
                    dataset$gated_out_umitabs[[score_i]]=list()
                }
                dataset$gated_out_umitabs[[score_i]][[sampi]]=is_res$gated_out_umitabs[[score_i]]
            }
            
            
            
        }
  
        ### AL edit 4/1/19
        if (!is.null(tmp_env$adttab)){
          colnames(tmp_env$adttab) <- paste(sampi,colnames(tmp_env$adttab),sep="_")
          dataset$adt_by_sample[[sampi]]=tmp_env$adttab
          #dataset$hto_by_sample[[sampi]]=tmp_env$htottab
        }
        if (!is.null(tmp_env$htotab)){
          colnames(tmp_env$htotab) <- paste(sampi,colnames(tmp_env$htotab),sep="_")
          dataset$hto_by_sample[[sampi]]=tmp_env$htottab
        }
        
        
        barcode_mask=tmp_env$numis_before_filtering[colnames(umitab)]>min_umis&tmp_env$numis_before_filtering[colnames(umitab)]<max_umis
        dataset$min_umis=min_umis
        dataset$max_umis=max_umis
        umitab=umitab[,barcode_mask]
        dataset$noise_models[genes,sampi]=tmp_env$noise_model[genes,1]
      #  genemask=intersect(rownames(umitab),rownames(model$models))
        projection_genemask=setdiff(genes,model$params$genes_excluded)
        noise_model=dataset$noise_models[projection_genemask,sampi]
        noise_model=noise_model/sum(noise_model)
        models=model$models[projection_genemask,]
        models=t(t(models)/colSums(models))
        message("Projecting ",ncol(umitab)," cells")
  #      genes=intersect(rownames(umitab),genes)
        if (is.null(model$alpha_noise)&is.null(model$avg_numis_per_model)){
            ll=getLikelihood(umitab[projection_genemask,],models =models,reg = model$params$reg)
        }
        else {
           if (!is.null(model$alpha_noise)){
            #        alpha_b=init_alpha[sampi]
            #          print(round(alpha_b,digits=6))
            #          res_l=getOneBatchCorrectedLikelihood(umitab=umitab[projection_genemask,],models,noise_model,alpha_noise=alpha_b,reg=model$params$reg)
            #          cell_to_cluster=MAP(res_l$ll)
            #          alpha_b=update_alpha_single_batch( umitab[projection_genemask,],models,noise_model,cell_to_cluster =cell_to_cluster,reg=model$params$reg )
            alpha_b=update_alpha_single_batch( umitab[projection_genemask,],models,noise_model,reg=model$params$reg )
            
            message("%Noise = ",round(100*alpha_b,digits=2))
            res_l=getOneBatchCorrectedLikelihood(umitab=umitab[projection_genemask,],cbind(models,noise_model),noise_model,alpha_noise=alpha_b,reg=model$params$reg)
            #     cells_to_exclude=apply(res_l$ll,1,which.max)==ncol(models)+1
            #      if (sum(cells_to_exclude)>0){
            #        message("Excluding ",sum(cells_to_exclude)," noisy barcodes")
            #        alpha_b=update_alpha_single_batch( umitab[projection_genemask,!cells_to_exclude],models,noise_model,reg=model$params$reg )
            #        message("%Noise = ",round(100*alpha_b,digits=2))
            #        res_l=getOneBatchCorrectedLikelihood(umitab=umitab[projection_genemask,!cells_to_exclude],models,noise_model,alpha_noise=alpha_b,reg=model$params$reg)
            
            #      print(median(colSums(umitab[projection_genemask,cells_to_exclude])))
            #      print(median(colSums(umitab[projection_genemask,!cells_to_exclude])))
            #    }
            
            #    message("%Noise = ",round(100*alpha_b,digits=2))
            
            dataset$alpha_noise[sampi]=alpha_b
            ll=res_l$ll[,1:ncol(models)]  
          }
          else {
            if (!is.null(model$avg_numis_per_model)){
              avg_numis_per_model=model$avg_numis_per_model
              gobclle_res=noiseEMsingleBatch(umitab=umitab[projection_genemask,],models=models,noise_model=noise_model,avg_numis_per_model=avg_numis_per_model,reg=model$params$reg,max_noise_fraction=.75)
              ll=gobclle_res$ll
              dataset$beta_noise[sampi]=gobclle_res$beta_noise
              dataset$avg_numis_per_sample_model[sampi,names(gobclle_res$avg_numis_per_model)]=gobclle_res$avg_numis_per_model
            }
          }
         
        cell_to_cluster=MAP(ll)
        cells_to_include=names(cell_to_cluster)[!cell_to_cluster%in%excluded_clusters]
        
        cell_to_cluster=cell_to_cluster[cells_to_include]
        ll=ll[cells_to_include,]
        dataset$ll<-rbind(dataset$ll,ll)
        umitab=umitab[,cells_to_include]
        
        dataset$cell_to_cluster<-c(dataset$cell_to_cluster,cell_to_cluster)
        
        tmp_counts=as.matrix(Matrix::t(aggregate.Matrix(Matrix::t(umitab[genes,]),cell_to_cluster,fun="sum")))
        dataset$counts[sampi,rownames(tmp_counts),colnames(tmp_counts)]=tmp_counts
        dataset$numis_before_filtering[[sampi]]=tmp_env$numis_before_filtering
        
        for (ds_i in 1:length(ds_numis_sampi)){
            dataset$ds[[ds_i]]<-cbind(dataset$ds[[ds_i]][genes,],tmp_env$ds[[ds_i]][genes,intersect(colnames(tmp_env$ds[[ds_i]]),cells_to_include)])
        }
        dataset$umitab<-cbind(dataset$umitab[genes,],umitab[genes,])
        cellids=colnames(umitab)
        
        cell_to_sampi=rep(sampi,length(cellids))
        names(cell_to_sampi)=cellids
        dataset$cell_to_sample<-c(dataset$cell_to_sample,cell_to_sampi)
        message("")
        rm(list=c("ll","cell_to_cluster","tmp_env","umitab"))
        }
    }
    
    message("...")
    
    
    dataset$umitab<-dataset$umitab[,-1]
    for (ds_i in 1:length(dataset$ds)){
        dataset$ds[[ds_i]]=dataset$ds[[ds_i]][,-1]
    }
    dataset$samples=samples
    
    
    if (!is.null(model$noise_models)){
          dataset$noise_counts=get_expected_noise_UMI_counts_alpha(dataset$umitab,dataset$cell_to_cluster,dataset$cell_to_sample,dataset$noise_models,dataset$alpha_noise,colnames(model$models))
        
    }
    
    if (all(is.na(dataset$alpha_noise))){
        dataset$alpha_noise=NULL
    }
    model$model_filename=model_fn
    output$dataset=dataset
    rm("dataset")
    included_clusters=setdiff(colnames(model$models),excluded_clusters)
    ncells_per_cluster<-rep(0,dim(model$models)[2])
    names(ncells_per_cluster)<-colnames(model$models)
    temptab=table(model$cell_to_cluster)
  #  temptab=temptab[setdiff(names(temptab),output$scDissector_params$excluded_cluster_sets)]
    ncells_per_cluster[names(temptab)]<-temptab
    output$ncells_per_cluster=ncells_per_cluster
    
    model$models=model$models[,included_clusters]
    output$clustAnnots=output$clustAnnots[included_clusters]
    output$ncells_per_cluster=output$ncells_per_cluster[included_clusters]
    output$dataset$ll=output$dataset$ll[,included_clusters]
    
    output$model=model
    output$cluster_order<-intersect(cluster_order,included_clusters)
    output$default_clusters<-intersect(cluster_order,included_clusters)
    output$loaded_model_version<-model_version_name
    return(output)
    
}


downsample=function(u,min_umis,chunk_size=100){
  non_zero_mask=Matrix::rowSums(u)!=0
  all_op=1:nrow(u[non_zero_mask,])
  base_tab=rep(0,sum(non_zero_mask))
  names(base_tab)=all_op
  
  downsamp_one=function(x,n){
    tab=base_tab
    tab2=table(sample(rep(all_op,x),size = n,replace=F))
    tab[names(tab2)]=tab2
    return(tab)
  }
  
  cell_mask=colnames(u)[Matrix::colSums(u,na.rm=T)>min_umis]  
  print(paste("Downsampling ", length(cell_mask), " cells to ",min_umis," UMIs",sep=""))
  
  
  breaks=unique(c(seq(from=1,to = length(cell_mask),by = chunk_size),length(cell_mask)+1))
  ds=Matrix(0,nrow =nrow(u),length(cell_mask),dimnames = list(rownames(u),cell_mask))
  
  
  for (i in 1:(length(breaks)-1)){
      ds[non_zero_mask,breaks[i]:(breaks[i+1]-1),drop=F]=Matrix(apply(u[non_zero_mask,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis))
 }
  
  
  return(ds)
}



import_dataset_and_model<-function(model_version_name,clustering_data_path,umitab,cell_to_cluster,cell_to_sample,min_umis=250,max_umis=25000,ds_numis=c(200,500,1000,2000),insilico_gating=NULL,clustAnnots=NA){
 
  require(Matrix)
  require(Matrix.utils)

  output<-list()
  output$scDissector_params<-list()
  
  
  
  
  get_cluster_set_list=function(clusts,i){
    if (i==0){
      return(clusts)
    }
    else{
      l=split(clusts,a[clusts,i])
      return(lapply(l,get_cluster_set_list,i-1))
    }
  }
  
  
 
 
  clusters=unique(cell_to_cluster)
  samples=unique(cell_to_sample)
  
  dataset<-new.env()
  dataset$ds=list()
  
  dataset$ds_numis=NULL
  dataset$ll<-c()



  
  dataset$alpha_noise=rep(NA,length(samples))
  
  names(dataset$alpha_noise)=samples
  
  genes=rownames(umitab)
  message("")
 
  dataset$gated_out_umitabs<-list()
  dataset$counts<-array(0,dim=c(length(samples),length(genes),length(clusters)),dimnames = list(samples,genes,clusters))

  
  dataset$numis_before_filtering=Matrix::colSums(umitab)
    
    barcode_mask=dataset$numis_before_filtering[colnames(umitab)]>min_umis&dataset$numis_before_filtering[colnames(umitab)]<max_umis
    dataset$min_umis=min_umis
    dataset$max_umis=max_umis
    umitab=umitab[,barcode_mask]
    cell_to_cluster=cell_to_cluster[barcode_mask]
    cell_to_sample=cell_to_sample[barcode_mask]
   
     for (ds_i in 1:length(ds_numis)){
      dataset$ds[[ds_i]]=downsample(umitab,min_umis=ds_numis[ds_i])
    }
    
   
    
    if (is.null(insilico_gating)){
      umitab=umitab
    }
    else{
      
      is_res=insilico_sorter(umitab,insilico_gating)
      umitab=is_res$umitab
      
      for (score_i in names(insilico_gating)){
        dataset$insilico_gating_scores[[score_i]]=c(dataset$insilico_gating_scores[[score_i]],is_res$scores[[score_i]])
        if (is.null(dataset$gated_out_umitabs[[score_i]])){
          dataset$gated_out_umitabs[[score_i]]=list()
        }
        dataset$gated_out_umitabs[[score_i]][[sampi]]=is_res$gated_out_umitabs[[score_i]]
      }
    }
    
    dataset$cell_to_cluster<-cell_to_cluster
    dataset$cell_to_sample<-cell_to_sample
    dataset$umitab<-umitab
    
    for (sampi in samples){
      maski=cell_to_sample==sampi
      tmp_counts=as.matrix(Matrix::t(aggregate.Matrix(Matrix::t(umitab[,maski]),cell_to_cluster[maski],fun="sum")))
      dataset$counts[sampi,rownames(tmp_counts),colnames(tmp_counts)]=tmp_counts
    }
    
    models=apply(dataset$counts,2:3,sum)
    models=t(t(models)/colSums(models))
    
    
    
    dataset$samples=samples
    dataset$ds_numis=ds_numis
    if (all(is.na(clustAnnots))){
      clustAnnots<-rep("unannotated",ncol(models))
      names(clustAnnots)<-colnames(models)
    }
    
    output$clustAnnots=clustAnnots
    
    
    nodes=names(clustAnnots)
    parents=as.character(clustAnnots)
    nodes[nodes==""]="unannotated"
    parents[parents==""]="unannotated"
    output$cluster_sets<-get_cluster_set_tree(data.frame(node=nodes,parent=parents))
    
    
  
   
   cluster_order=colnames(models)
    
  model=new.env()
  model$models=models
    
  output$dataset=dataset
  rm("dataset")
  included_clusters=unique(cell_to_cluster)
  ncells_per_cluster<-rep(0,length(included_clusters))
  names(ncells_per_cluster)<-included_clusters
  temptab=table(cell_to_cluster)
  #  temptab=temptab[setdiff(names(temptab),output$scDissector_params$excluded_cluster_sets)]
  ncells_per_cluster[names(temptab)]<-temptab
  output$ncells_per_cluster=ncells_per_cluster
  
  
  output$clustAnnots=output$clustAnnots[included_clusters]
  output$ncells_per_cluster=output$ncells_per_cluster[included_clusters]
 
  output$model=model
  output$cluster_order<-intersect(cluster_order,included_clusters)
  output$default_clusters<-intersect(cluster_order,included_clusters)
  output$loaded_model_version<-model_version_name
  fn_prefix=paste(clustering_data_path,"/imported_model_",model_version_name,sep="")
  output$model$model_filename<-paste(fn_prefix,".rd",sep="")
  
  order_fn=paste(fn_prefix,"_order.txt",sep="")
  if (file.exists(order_fn)){
    output$cluster_order=as.character(read.table(file=order_fn,header = F)[,1])
  }
  else{
    output$cluster_order=colnames(output$model$models)
  }
  output$default_clusters<-output$cluster_order
  
  return(output)
  
}



load_seurat_rds=function(rds_file,model_name="",clustering_data_path="",min_umis=200,max_umis=25000,ds_numis=c(200,500,1000,2000),sample_ID_converter=NULL){
    a=readRDS(rds_file)

    if (unlist(attributes(a)$version)[1]<3){
      umitab=attributes(a)[["raw.data"]]
      
      cells=attributes(a)[["cell.names"]]
      umitab=umitab[,cells]
      cluster_factor=as.factor(attributes(a)[["ident"]])
      annots=levels(cluster_factor)
      names(annots)=1:length(annots)
      cell_to_cluster=as.numeric(cluster_factor)
      names(cell_to_cluster)=cells
      colnames(umitab)=cells
      cell_to_sample=attributes(a)$meta.data$orig.ident
      if (!is.null(sample_ID_converter)){
        cell_to_sample=sample_ID_converter[cell_to_sample]
      }
      names(cell_to_sample)=cells
      l=import_dataset_and_model(model_name,clustering_data_path=clustering_data_path,umitab=umitab,cell_to_cluster=cell_to_cluster,cell_to_sample=cell_to_sample,min_umis=min_umis,max_umis=max_umis,ds_numis=ds_numis,insilico_gating=NULL,clustAnnots=annots)
    }
  else {
    umitab=attributes(a)$assay$RNA@counts
    cluster_factor=attributes(a)$active.ident
    cells=colnames(umitab)
    cell_to_cluster=as.numeric(cluster_factor)
    names(cell_to_cluster)=names(cluster_factor)
    cell_to_sample=attributes(a)$meta.data$orig.ident
    if (!is.null(sample_ID_converter)){
      cell_to_sample=sample_ID_converter[cell_to_sample]
    }
    names(cell_to_sample)=cells
    l=import_dataset_and_model(model_name,clustering_data_path=clustering_data_path,umitab=umitab,cell_to_cluster=cell_to_cluster,cell_to_sample=cell_to_sample,min_umis=min_umis,max_umis=max_umis,ds_numis=ds_numis,insilico_gating=NULL)
  }
  return(l)
}





load_metacell_clustering=function(mc_rda,mat_rda,clustering_data_path,name="",min_umis=200,max_umis=25000,ds_numis=c(200,500,1000,2000),metacell_to_metacellset=NULL,sample_ID_converter=NULL){
  mc=new.env()
  mat=new.env()
  load(mc_rda,envir = mc)
  load(mat_rda,envir=mat)
  umitab=attributes(mat$object)$mat
  cell_to_cluster=as.character(attributes(mc$object)$mc)
  names(cell_to_cluster)=names(attributes(mc$object)$mc)
  cell_to_sample=as.character(attributes(mat$object)$cell_metadata$amp_batch_id)
  names(cell_to_sample)=rownames(attributes(mat$object)$cell_metadata)
  cells=intersect(intersect(names(cell_to_cluster),names(cell_to_sample)),colnames(umitab))
  umitab=umitab[,cells]
  cell_to_cluster=cell_to_cluster[cells]
  if (!is.null(metacell_to_metacellset)){
    cell_to_cluster=as.character(metacell_to_metacellset[cell_to_cluster])
  }
  names(cell_to_cluster)=cells
  cell_to_sample=cell_to_sample[cells]
  if (!is.null(sample_ID_converter)){
    cell_to_sample=sample_ID_converter[cell_to_sample]
  }
  names(cell_to_sample)=cells
  ldm=import_dataset_and_model(name,clustering_data_path=clustering_data_path,umitab=umitab,cell_to_cluster=cell_to_cluster,cell_to_sample=cell_to_sample,min_umis=min_umis,max_umis=max_umis,ds_numis=ds_numis,insilico_gating=NULL,clustAnnots=NULL)
  return(ldm)
  
}

create_clustering_data_dir=function(path,samples=c(),model_names=c()){
  if(dir.exists(path)){
    stop("Error! Cannot create ",path,". It already exists!")
  }
  dir.create(path)
  dir.create(paste(path,"/metadata",sep=""))
  samples_tab=matrix("",length(samples),3,dimnames = list(NULL,c("index","path","title")))
  samples_tab[,1]=samples
  model_versions_tab=matrix("",length(model_names),2,dimnames = list(NULL,c("title","path")))
  model_versions_tab[,1]=model_names
  sample_annots_tab=matrix("",length(samples),15,dimnames=list(NULL,c("sample_ID","amp_batch_ID","old_lib_name","Disease","Origin","tissue","status","biotherapy_status","Inclusion date (M/D/Y)","drug","Patient ID","GRID ID","COMPASS ID","CHEMISTRY","details")))
  sample_annots_tab[,1]=samples
  write.csv(file=paste(path,"/samples.csv",sep=""),samples_tab,row.names = F)
  write.csv(file=paste(path,"/model_versions.csv",sep=""),model_versions_tab,row.names = F)
  write.csv(file=paste(path,"/metadata/sample_annots.csv",sep=""),sample_annots_tab,row.names = F)
}

