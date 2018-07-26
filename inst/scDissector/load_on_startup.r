load_on_startup <- function(){
    print("not loading")
}

#model_fn=paste(scDissector_datadir,"pbmc_model__30.rd",sep="/")

#samples_fn=paste(scDissector_datadir,c("data_Fresh_pbmc_68k.rd","data_pbmc_6k_v1.rd","data_regultatory_t.rd" ,"data_naive_t.rd","data_cd_14_monocytes.rd","data_b_cells.rd"),sep="/")
#names(samples_fn)=c("pbmc_68k_v1","pbmc_6k_v1","Zheng_regultatory_t","Zheng_naive_cd4_t","Zheng_cd_14_monocytes","Zheng_b_cells")
#default_model_dataset=load_dataset_and_model(model_fn,samples_fn)
