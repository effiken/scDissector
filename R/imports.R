import_seurat=function(rds_fn,model_name,sample_id){
  a=readRDS(rds_fn)
  umitab=attributes(a)[["raw.data"]]
  umitab=umitab[,cells]
  cells=attributes(a)[["cell.names"]]
  cluster_factor=as.factor(attributes(a)[["ident"]])
  annots=levels(cluster_factor)
  names(annots)=1:length(annots)
  cell_to_cluster=as.numeric(cluster_factor)
  names(cell_to_cluster)=cells
  colnames(umitab)=cells
  cell_to_sample=rep(sample_id,length(cells))
  names(cell_to_sample)=cells
  ldm=import_dataset_and_model(model_name,umitab=umitab,cell_to_cluster=cell_to_cluster,cell_to_sample=cell_to_sample,min_umis=250,max_umis=25000,ds_numis=c(200,500,1000,2000),insilico_gating=NULL,clustAnnots=annots){
  return(ldm)
}
#example
#default_model_and_dataset=  import_seurat("pbmc3k.rds","seurat_pbmc3k","pbmc3k")
# run_scDissector()