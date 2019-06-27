############################################

get_loaded_samples=function(session){
  return(session$userData$dataset$samples)
}

get_ds_options=function(session){
  return(session$userData$dataset$ds_numis)
}

get_max_umis=function(session){
  return(session$userData$dataset$max_umis)
}

get_min_umis=function(session){
  return(session$userData$dataset$min_umis)
}

dataset_exists=function(session){
  return (!is.null(session$userData$dataset))
}


get_all_cells=function(session){
  return(names(session$userData$dataset$cell_to_cluster))
}

get_all_genes=function(session){
  return(rownames(session$userData$dataset$umitab))
}

get_all_clusters=function(session){
  return(colnames(session$userData$model$models))
}

get_ds_cells=function(session,ds_version){
  return(colnames(session$userData$dataset$ds[[ds_version]]))
}

select_cells=function(session,cells=NULL,clusters=NULL,samples=NULL){
  if (is.null(cells)){
    cells=get_all_cells(session)
  }
  if (is.null(samples)){
    return(cells[session$userData$dataset$cell_to_cluster[cells]%in%clusters])
  }
  else if (is.null(clusters)){
    return(cells[session$userData$dataset$cell_to_sample[cells]%in%samples])
  }
  else{
    return(cells[session$userData$dataset$cell_to_cluster[cells]%in%clusters&session$userData$dataset$cell_to_sample[cells]%in%samples])
  }
}

get_cell_to_sample=function(session,cells=NULL){
  if (is.null(cells)){
    return(session$userData$dataset$cell_to_sample)
  }
  else {
    return(session$userData$dataset$cell_to_sample[cells])
  }
}

get_cell_to_cluster=function(session,cells=NULL){
  if (is.null(cells)){
    return(session$userData$dataset$cell_to_cluster)
  }
  else {
    return(session$userData$dataset$cell_to_cluster[cells])
  }
}


get_umitab=function(session,cells=NULL,genes=NULL){
  if (is.null(genes)){
    genes=rownames(session$userData$dataset$umitab)    
  }
  else {
    genes=intersect(genes,rownames(session$userData$dataset$umitab))
  }
  if (is.null(cells)){
    cells=colnames(session$userData$dataset$umitab)
  }
  else {
    cells=intersect(cells,colnames(session$userData$dataset$umitab))
  }
  return(session$userData$dataset$umitab[genes,cells,drop=F])
}

get_models=function(session,clusters=NULL,genes=NULL){
  if (is.null(genes)){
    genes=rownames(session$userData$model$models)    
  }
  else{
    genes=intersect(genes,rownames(session$userData$model$models))
  }
  if (is.null(clusters)){
    clusters=colnames(session$userData$model$models)
  }
  else{
    clusters=intersect(clusters,colnames(session$userData$model$models))
  }
  return(session$userData$model$models[genes,clusters,drop=F])
}

get_noise_models=function(session,samples=NULL,genes=NULL){
  if (is.null(session$userData$model$noise_models)){
    return()
  }
  if (is.null(genes)){
    genes=rownames(session$userData$model$noise_models)    
  }
  else{
    genes=intersect(genes,rownames(session$userData$model$models))
  }
  if (is.null(samples)){
    samples=colnames(session$userData$model$noise_models)
  }
  else{
    samples=intersect(samples,colnames(session$userData$model$noise_models))
  }
  return(session$userData$model$noise_models[genes,samples,drop=F])
}


get_dstab=function(session,cells=NULL,genes=NULL,ds_version){
  if (is.null(genes)){
    genes=rownames(session$userData$dataset$ds[[ds_version]])    
  }
  if (is.null(cells)){
    cells=colnames(session$userData$dataset$ds[[ds_version]])
  }
  genes=intersect(genes,rownames(session$userData$dataset$ds[[ds_version]]))
  cells=intersect(cells,colnames(session$userData$dataset$ds[[ds_version]]))
  return(session$userData$dataset$ds[[ds_version]][genes,cells,drop=F])
}

get_counts_array=function(session,samples=NULL,genes=NULL,clusters=NULL){
  if (is.null(samples)){
    samples=1:dim(session$userData$dataset$counts)[1]
  }
  else {
    samples=intersect(samples,dimnames(session$userData$dataset$counts)[[1]])
  }
  if (is.null(genes)){
    genes=1:dim(session$userData$dataset$counts)[2]
  }
  else {
    genes=intersect(genes,dimnames(session$userData$dataset$counts)[[2]])
  }
  if (is.null(clusters)){
    clusters=1:dim(session$userData$dataset$counts)[3]
  }else{
    clusters=intersect(clusters,dimnames(session$userData$dataset$counts)[[3]])
  }
  session$userData$dataset$counts[samples,genes,clusters,drop=F]
}

get_noise_counts_array=function(session,samples=NULL,genes=NULL,clusters=NULL){
  if (is.null(session$userData$dataset$noise_counts)){
    return(null)
  }
  if (is.null(samples)){
    samples=1:dim(session$userData$dataset$noise_counts)[1]
  }
  else {
    samples=intersect(samples,dimnames(session$userData$dataset$noise_counts)[[1]])
  }
  if (is.null(genes)){
    genes=1:dim(session$userData$dataset$noise_counts)[2]
  }
  else {
    genes=intersect(genes,dimnames(session$userData$dataset$noise_counts)[[2]])
  }
  if (is.null(clusters)){
    clusters=1:dim(session$userData$dataset$noise_counts)[3]
  }else{
    clusters=intersect(genes,dimnames(session$userData$dataset$noise_counts)[[3]])
  }
  session$userData$dataset$noise_noise_counts[samples,genes,clusters,drop=F]
}

get_insilico_gating_scores=function(session,cells=NULL){
  return(session$userData$dataset$insilico_gating_scores)
}

get_numis_before_filtering=function(session,samp){
  return(session$userData$dataset$numis_before_filtering[[samp]])
}

get_clustering_params=function(session){
  return(session$userData$model$params)
}

