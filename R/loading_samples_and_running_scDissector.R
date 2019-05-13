# The clustering data dir contains all the compiled samples, model versions and metadata
clustering_data_dir="path_example/clustering_data_exaple"

# model_name is associated with a model file by the model_versions.csv table in the clustering_data dir
model_name="model_name_example"
# sample_names is a vector containing the sample names which are associated with sample files by the samples.csv file in the clustering_data dir
sample_names=c("sample1_example","sample2_example","sample3_example")

#sample_annots file can be used to select specific samples by their metadata
annots=read.csv(paste(clustering_data_dir,"/metadata/","sample_annots.csv",sep=""),stringsAsFactors = F)

sample_to_fn=read.csv(paste(clustering_data_dir,"/samples.csv",sep=""),stringsAsFactors = F,row.names = 1)

model_fn=paste(clustering_data_dir,"/",read.csv(paste(clustering_data_dir,"/model_versions.csv",sep=""),stringsAsFactors = F,row.names = 1)[model_name,1],sep="")
sample_fns=paste(clustering_data_dir,sample_to_fn[sample_names,1],sep="/")
# Loading the samples and projecting them onto the model
ldm=load_dataset_and_model(model_fn = model_fn,sample_fns = sample_fns)

# ldm can be stored on disk (optional):
save(ldm,file="ldm_path_example")

# If saved, in the next time you can load the ldm file instead of regenrating it
load(file="ldm_path_example")

# run scDissector
run_scDissector(preloaded_data = ldm,clustering_data_path =clustering_data_dir)