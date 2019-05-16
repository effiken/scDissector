#' preload scRNA sample and project onto a model
#'  
#' @param clustering_data_path path to the folder that cntains the compiled samples, model versions and metadata
#' @param model_name model name as defined in clustering_data_path/model_versions.csv
#' @param sample_names vector containing sample names as defined in clustering_data_path/samples.csv
#' @return ldm object
#' @examples
#' # ldm=load_ldm(clustering_data_path="path_example/clustering_data_exaple",model_name="model_name_example",sample_names=c("sample1_example","sample2_example","sample3_example"))
#' # ldm can be stored on disk (optional):
#' # save(ldm,file="ldm_path_example")
#' # If saved, in the next time you can load the ldm file instead of regenrating it
#' # load(file="ldm_path_example")
#' # run scDissector
#' # run_scDissector(preloaded_data = ldm,clustering_data_path =clustering_data_path)
#' @export
#
load_scDissector_data=function(clustering_data_path,model_name,sample_names, min_umis=250,max_umis=25000,ds_numis=NA,max_ncells_per_sample=NA){
#sample_annots file can be used to select specific samples by their metadata
  annots=read.csv(paste(clustering_data_path,"/metadata/","sample_annots.csv",sep=""),stringsAsFactors = F)

  sample_to_fn=read.csv(paste(clustering_data_path,"/samples.csv",sep=""),stringsAsFactors = F,row.names = 1)

  
  model_fn=paste(clustering_data_path,"/",read.csv(paste(clustering_data_path,"/model_versions.csv",sep=""),stringsAsFactors = F,row.names = 1)[model_name,1],sep="")
  sample_fns=paste(clustering_data_path,sample_to_fn[as.character(sample_names),1],sep="/")
# Loading the samples and projecting them onto the model
  ldm=load_dataset_and_model(model_fn = model_fn,sample_fns = sample_fns,min_umis=min_umis,max_umis=max_umis,ds_numis=ds_numis,max_ncells_per_sample=max_ncells_per_sample)
  return(ldm)
}
