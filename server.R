library(Matrix)
library(ggvis)
library(dplyr)
library(gplots)
source("projector.r")
source("gene_symbol_convertors.r")
source("load_dataset.R")
source("scPlots.R")
library("heatmaply")
set.seed(3505)

non_data_tabs=c("Gating","Basics","Clusters","Truth","QC","Clustering QC","Gene Modules","Samples")



colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))
sample_cols<<-rep(paste("#",read.table("sample_colors.txt",stringsAsFactors = F)[,1],sep=""),10)


print(getwd())
genesetsfile="gene_sets.txt"
if (file.exists(genesetsfile)){
  geneList_tmp<-read.table(file=genesetsfile,header=T,stringsAsFactors = F,row.names =1)
  geneList<<-geneList_tmp[,1]
  names(geneList)<<-rownames(geneList_tmp)
}
source("load_on_startup.r")

gsc=get_gene_symbol_convetors()
gene_symbol_old2new<<-gsc$old2new
gene_symbol_new2old<<-gsc$new2old

hide_all_tabs=function(){
  for (tab in non_data_tabs){
    hideTab(inputId = "inMain", target = tab)
  }
}

show_all_tabs=function(){
  hideTab(inputId = "inMain", target = "Data")
  for (tab in non_data_tabs){
    showTab(inputId = "inMain", target = tab)
  }
}

update_clusters=function(session,new_clusts,delete_old=F){
  if(delete_old){
    session$userData$scDissector_params$previous_clusters<-new_clusts
  }
  updateTextInput(session,"inClusters",value=paste(new_clusts,collapse=",")) 
}

update_genes=function(session,new_genes,delete_old=F){
  if(delete_old){
    session$userData$scDissector_params$previous_genes<-new_genes
  }
  updateTextInput(session,"inGenes",value=paste(new_genes,collapse=", ")) 
}



update_all= function(session,ldm){
  for (item in names(ldm)){
    session$userData[[item]]=ldm[[item]]
  }
  
  ncells_choices=as.numeric(setdiff(params$nrandom_cells_per_sample_choices,"All"))
  ncells_per_sample=ncells_choices[which.min(abs((2000/length(session$userdata$samples))-ncells_choices))]
  clust_title=paste(session$userData$cluster_order," - ",session$userData$clustAnnots[session$userData$cluster_order],sep="") 
  session$userData$scDissector_params$previous_clusters=session$userData$default_clusters

  update_clusters(session,session$userData$default_clusters,T)
  updateSelectInput(session,"inAnnotateClusterNum",choices = clust_title)
  updateSelectInput(session,"inQCClust",choices =clust_title)
  updateSelectInput(session,"inClustForDiffGeneExprsProjVsRef",choices = clust_title)
  updateTextInput(session,"inSamplesToShow",value = paste(session$userData$dataset$samples,collapse=", "))
  
  updateSelectInput(session,"inGatingSample",choices = session$userData$dataset$samples)
  updateSelectInput(session,"inGatingShowClusters",choices = clust_title)
  updateSelectInput(session,"input$inTruthNcellsPerSample",choices=params$nrandom_cells_per_sample_choices,selected =ncells_per_sample )
  updateSelectInput(session,"inQCDownSamplingVersion",choices=session$userData$dataset$ds_numis,selected = ifelse("1000"%in%session$userData$dataset$ds_numis,"1000",max(session$userData$dataset$ds_numis)))
  updateSelectInput(session,"inTruthDownSamplingVersion",choices=session$userData$dataset$ds_numis,selected = ifelse("1000"%in%session$userData$dataset$ds_numis,"1000",max(session$userData$dataset$ds_numis)))
  updateSelectInput(session,"inModulesDownSamplingVersion",choices=session$userData$dataset$ds_numis,selected = max(session$userData$dataset$ds_numis))
  if (is.null(ldm$model$noise_models)){
     updateSelectInput(session,"inModelOrAverage",choices=c("Model","Average"))
  }
  
  
  updateTextInput(session,"inMaxUmis",value = session$userData$dataset$max_umis)
  updateTextInput(session,"inMinUmis",value = session$userData$dataset$min_umis)
  message("Successfully finished loading.")
  
  session$userData$loaded_flag<-T
}

randomly_select_cells=function(ldm,nrandom_cells_per_sample_choices){
  randomly_selected_cells=list()
  for (ds_i in 1:length(ldm$dataset$ds_numis)){
    randomly_selected_cells[[ds_i]]<-list()
    for (sampi in ldm$dataset$samples){
      for (nrandom_cells in nrandom_cells_per_sample_choices){
        randomly_selected_cells[[ds_i]][[nrandom_cells]]=c()
        if (nrandom_cells=="All"||pmax(0,as.numeric(nrandom_cells),na.rm=T)>=ncol(ldm$dataset$ds[[ds_i]])){
          randomly_selected_cells[[ds_i]][[nrandom_cells]]<-c(randomly_selected_cells[[ds_i]][[nrandom_cells]],colnames(ldm$dataset$ds[[ds_i]]))
        }
        else{
          randomly_selected_cells[[ds_i]][[nrandom_cells]]<-c(randomly_selected_cells[[ds_i]][[nrandom_cells]],sample(colnames(ldm$dataset$ds[[ds_i]]),size=as.numeric(nrandom_cells),replace=F))
        }
      }
    }
  }
  return(randomly_selected_cells)
}




tab3_left_margin=12


  function(input, output, session) {
    
    
    
    session$userData$prev_xy_range_G=c(0,0,0,0)
    session$userData$prev_xy_range_C=c(0,0,0,0)
    session$userData$update_ingene_text_input=T
    session$userData$update_inclusts_text_input=T
    hide_all_tabs()    
    session$userData$loaded_flag<-F
    session$userData$loaded_model_version=NA
    session$userData$prev_inModelColorScale_rel<-c(-4,4)
    session$userData$prev_inModelColorScale_abs<-c(-7,-1)
    session$userData$click_flag<-T
    
    
    
    
    observeEvent(input$inDatapath,{
      if (input$inDatapath==""){
        return()
      }
      verfile=paste(input$inDatapath,"model_versions.csv",sep="/")
      samples_file=paste(input$inDatapath,"samples.csv",sep="/")
      read_flag=T
      if (!file.exists(verfile)){
        message(verfile ," not found")
        read_flag=F
      }
      if  (!file.exists(samples_file)){
        message(samples_file ," not found")
        read_flag=F
      }
    
      if (!read_flag){
        updateSelectInput(session,"inModelVer", "Model Version:",choices = "")
        updateSelectInput(session,"inProjectedDataset","Projected Version:",choices="")
        return()
      }
      session$userData$vers_tab<-read.csv(file=verfile,header=T,stringsAsFactors = F)
      session$userData$samples_tab<-read.csv(file=samples_file,header=T,stringsAsFactors = F)
      
   
      myGeneListFile= paste(input$inDatapath,"gene_sets.txt",sep="/")
   
      if (file.exists(myGeneListFile)){
        added_gene_list_tmp=read.delim(file=myGeneListFile,header=T,stringsAsFactors = F)
    
        if (length(added_gene_list_tmp[,1])!=length(unique(added_gene_list_tmp[,1]))){
          message(myGeneListFile, " was not loaded since gene sets names are not unique.")
        }
        else{
          names2=c(added_gene_list_tmp[,1],names(session$userData$geneList))
          session$userData$geneList<-c(added_gene_list_tmp[,2],session$userData$geneList)
          names(session$userData$geneList)=names2
        }
      }  
    
    
    mySampleSetsFile= paste(input$inDatapath,"sample_sets.txt",sep="/")
    samples_to_show=session$userData$samples_tab$index
    if (file.exists(mySampleSetsFile)){
      sample_sets_tab<-read.table(file=mySampleSetsFile,header=T,stringsAsFactors = F,row.names =1)
      session$userData$sample_sets<-strsplit(sample_sets_tab[,"samples"],",| ,|, ")
      names(session$userData$sample_sets)<-rownames(sample_sets_tab)
      samples_to_show=c(names(session$userData$sample_sets),samples_to_show)
     }
    session$userData$scDissector_datadir<-input$inDatapath
    
    updateSelectInput(session,"inGeneSets",choices = names(session$userData$geneList))
    updateSelectInput(session,"inModelVer", "Model Version:",choices =  session$userData$vers_tab$title)
    updateSelectInput(session,"inSampleToAdd", "Samples:",choices =samples_to_show)
    updateSelectInput(session,"inProjectedDataset","Projected Version:",choices =  session$userData$vers_tab$title)
  })
  
  observeEvent(input$inLoad,{
    if (input$inModelVer=="")
    {
      return()
    }
    
    i=match(input$inModelVer, session$userData$vers_tab$title)
    model_fn=paste(input$inDatapath, session$userData$vers_tab$path[ifelse(is.na(i),1,i)],sep="/")
    
    if (!file.exists(model_fn)){
      message(model_fn, " not found")
      return()
    }

    session$userData$loaded_model_file<-model_fn
    
    samples=strsplit(input$inSamples,",| ,|, ")[[1]]
    tmp_samples=as.list(samples)
    
    if (!is.null(session$userData$sample_sets)){
      mask=samples%in%names(session$userData$sample_sets)
      tmp_samples[mask]=session$userData$sample_sets[samples[mask]]
    }
    samples=unique(unlist(tmp_samples))
    if (!all(samples%in%session$userData$samples_tab$index)){
      message("Error: These samples could not be recognized: ",paste(samples[!(samples%in%session$userData$samples_tab$index)],collapse=","),". Loading failed!")
      message("")
      return()
    }
    
    sample_paths=paste(input$inDatapath,session$userData$samples_tab$path[match(samples,session$userData$samples_tab$index)],sep="/")
    names(sample_paths)=samples
    
    isolate({
      min_umis=as.numeric(input$inMinUmis)
      max_umis=as.numeric(input$inMaxUmis)
      })

    ldm=load_dataset_and_model(model_fn,sample_paths,min_umis = min_umis,max_umis = max_umis)
    
    randomly_selected_cells=randomly_select_cells(ldm,params$nrandom_cells_per_sample_choices)
    ldm$dataset$randomly_selected_cells=randomly_selected_cells
    show_all_tabs()
    update_all(session,ldm)
    
  })
  
  observeEvent(input$inAddSample,{
    s1=paste(unique(c(strsplit(input$inSamples,",|, | ,")[[1]],input$inSampleToAdd)),collapse=", ")
    updateTextInput(session,"inSamples",value=s1) 
  })
  
  observeEvent(input$inResetSamples,{
    samples=colnames(session$userData$dataset$samples)
    updateTextAreaInput(session,"inSamplesToShow",value = paste(samples,collapse=", "))
  })
  
  observeEvent(input$inAbsOrRel,{
    if (input$inAbsOrRel=="Absolute"){
      
      session$userData$prev_inModelColorScale_rel<-input$inModelColorScale
      
      updateSliderInput(session,"inModelColorScale",label="Log10(expression)",min = -10,max=-.5,step = .5,value = session$userData$prev_inModelColorScale_abs)
    }
    else if (input$inAbsOrRel=="Relative"){
      session$userData$prev_inModelColorScale_abs<-input$inModelColorScale
      updateSliderInput(session,"inModelColorScale",label="Log2(expression/mean)",min = -8,max = 8,step = 1,value = session$userData$prev_inModelColorScale_rel)
      
    }
  })
  
  genes_reactive <-reactive({
    xy_range <- event_data("plotly_relayout")
    
    if (!session$userData$loaded_flag){
      return()
    }
 
    genes=init_genes(input$inGenes)
    if ((!is.null(xy_range))&(sum(c('xaxis.range[0]','xaxis.range[1]')%in%names(xy_range))==2)){
      if(!all(session$userData$prev_xy_range_G==unlist(xy_range))){
        if (min(unlist(xy_range),na.rm=T)<0){
          session$userData$prev_xy_range_G=unlist(xy_range) 
          update_genes(session,genes,T)
        }
        else{
          session$userData$prev_xy_range_G=unlist(xy_range)
          genes=genes[ceiling(as.numeric(xy_range$'xaxis.range[0]')):floor(as.numeric(xy_range$'xaxis.range[1]'))]
          update_genes(session,genes,F)
        }
      }
    }
    
    return(genes)
  })
  
  clusters_reactive <-reactive({
    xy_range <- event_data("plotly_relayout")
    if (!session$userData$loaded_flag){
      return()
    }
    
    clusts=strsplit(input$inClusters,",|, | ,")[[1]]
    clusts=intersect(clusts,setdiff(colnames(session$userData$model$models),session$userData$scDissector_params$excluded_clusters))
    if ((!is.null(xy_range))&(sum(c('yaxis.range[0]','yaxis.range[1]')%in%names(xy_range))==2)){
      if (!all(session$userData$prev_xy_range_C==unlist(xy_range))){
        if (min(unlist(xy_range),na.rm=T)<0){
          session$userData$prev_xy_range_C=unlist(xy_range) 
          clusts=session$userData$scDissector_params$previous_clusters
        }
        else{
          session$userData$prev_xy_range_C=unlist(xy_range)
          clusts=rev(rev(clusts)[ceiling(as.numeric(xy_range$'yaxis.range[0]')):floor(as.numeric(xy_range$'yaxis.range[1]'))])
        }
        session$userData$update_inclusts_text_input=F
        update_clusters(session,clusts)
        
      }
    }
    return(clusts)
  })
  
  
  clusters_genes_sampples_reactive <-reactive({
    
    return(list(clusters=clusters_reactive(),genes=genes_reactive(),samples=samples_reactive()))
  })
  
  
  output$event <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover on a point!" else d
  })
  
  samples_reactive <-reactive({
    samples=strsplit(input$inSamplesToShow,",|, | ,")[[1]]
     if (is.null(session$userData$dataset)){
      return(c())
     }
    samples=intersect(samples,session$userData$dataset$samples)
    return(samples)
  })
  

  
  modules_reactive <-reactive({
    modules=strsplit(input$inModules,",")[[1]]
    return(modules)
  })
  
  
  ds_QC_reactive <-reactive({
    clust=strsplit(input$inQCClust," - ")[[1]][1]
    samples=samples_reactive()
    
    if (is.null(session$userData$dataset)){
      return()
    }
    dataset=session$userData$dataset
    ds=dataset$ds[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]]
    sampling_mask=dataset$randomly_selected_cells[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]][[match(input$inQCNcellsPerSample,params$nrandom_cells_per_sample_choices)]]
    cluster_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster==clust]
    return(ds[,intersect(sampling_mask,cluster_mask),drop=F])
  })
  
  observeEvent(input$inGeneSets,{
    geneSets=input$inGeneSets
    updateTextInput(session,"inGenes",value=as.character(session$userData$geneList[geneSets]))
    })
  
  
  init_genes=function(genes_string){
    if (length(genes_string)==0){
      return()
    }
    data_genes=rownames(session$userData$dataset$umitab)
 
    gene_strings_adj=adjust_gene_names(genes_string, data_genes)
 
    session$userData$gcol<-ifelse(toupper(gene_strings_adj)%in%tfs ,4,1)
    names(session$userData$gcol)<-toupper(gene_strings_adj)
    return(gene_strings_adj)
    
  }
  
 
  
  observeEvent(input$inAnnotateCluster, {
    session$userData$clustAnnots[strsplit(input$inAnnotateClusterNum," - ")[[1]][1]]<-input$inClustAnnot
    updateTextInput(session,"inClustAnnot",value="")
  })
  
  observeEvent(input$inSaveAnnot, {
    f=paste(input$inDatapath, session$userData$vers_tab$path[ session$userData$vers_tab$title==session$userData$loaded_model_version],sep="/")
 
    annot_fn=paste(strsplit(f,"\\.")[[1]][1],"_annots.txt",sep="")
    write.table(file=annot_fn,session$userData$clustAnnots,row.names=T,col.names=F,quote=F,sep="\t")
    message("Cluster annotations saved to ",annot_fn,".")
    
  })
  
  
  observeEvent(input$inClusterGenes, {
    inclusts=clusters_reactive()
    ingenes=genes_reactive()
    mat<-session$userData$model$models[match(ingenes,rownames(session$userData$model$models)),inclusts]
    if (length(ingenes)==0){
      return()
    }
    if (input$inReorderingMethod=="Hierarchical clustering"){
      mat2=mat[!rowSums(is.na(mat))==ncol(mat),]
      cormat=cor(t(mat2),use="comp")
      mask=rowSums(!is.na(cormat))>1
      cormat=cormat[mask,mask]
    
      ord=hclust(dist(1-cormat))$order
      genes=c(rownames(cormat)[ord],setdiff(rownames(mat),rownames(cormat)))
      update_genes(session,genes,F)
    }
    else if (input$inReorderingMethod=="Diagonal"){
        clusters=clusters_reactive()
        ord=order(apply(mat[,clusters],1,which.max))
        genes=rownames(mat)[ord]
        update_genes(session,genes,F)
      }
  })
  
  observeEvent(input$inClusterModules, {
    inclusts=clusters_reactive()
    modulemat<-module_counts_reactive()
    if(is.null(modulemat)){
      return()
    }
    inmodules=modules_reactive()
    modulemat<-modulemat[inmodules,inclusts]
    
    if (input$inReorderingModulesMethod=="Hierarchical clustering"){
      cormat=cor(t(modulemat),use="comp")
      mask=rowSums(!is.na(cormat))>1
      cormat=cormat[mask,mask]
      
      ord=hclust(dist(1-cormat))$order
      updateTextInput(session,"inModules",value=paste(rownames(cormat)[ord],collapse=","))
    }
    else if (input$inReorderingModulesMethod=="Diagonal"){
      clusters=clusters_reactive()
      ord=order(apply(modulemat[,clusters],1,which.max))
      updateTextInput(session,"inModules",value=paste(rownames(modulemat)[ord],collapse=","))
    }
  })

  
  observeEvent(input$inBlindChisqSelectedClusters, {
    message("screening for variable ",input$inSelectGenesFrom)
    clusters=clusters_reactive()
    samples=samples_reactive()
    if (input$inSelectGenesFrom=="All genes"){
      pref="Var"
      genes=rownames(session$userData$model$models)
    } else if (input$inSelectGenesFrom=="TFs"){
      pref="Tfs"
      genes=intersect(tfs,rownames(session$userData$model$models))
    } else if(input$inSelectGenesFrom=="Surface markers"){
      pref="Surf"
      genes=intersect(surface_markers,rownames(session$userData$model$models))
      
    }
    cells=names(session$userData$dataset$cell_to_cluster[session$userData$dataset$cell_to_cluster%in%clusters])
    genes=intersect(genes,rownames(session$userData$dataset$umitab))
    chisq_res2=chisq_genes(session$userData$dataset$counts[samples,genes,clusters,drop=F])
    if (nrow(chisq_res2)==0){
      message("Chi sq detected nothing..")
      return()
    }
    counts=apply(session$userData$dataset$counts[samples,,,drop=F],2:3,sum)
    avg=t(t(counts)/colSums(counts))

    mask=rownames(chisq_res2)%in%genes&chisq_res2[,3]<0.05&
      rowSums(session$userData$dataset$umitab[rownames(chisq_res2),cells]>1)>=10&
      (apply(avg[rownames(chisq_res2),]/(rowSums(counts[rownames(chisq_res2),])/sum(counts)),1,max)>4|
      apply(avg[rownames(chisq_res2),]/rowMeans(avg[rownames(chisq_res2),]),1,max)>4)
    isolate({
      ngenes_to_show=as.numeric(input$inNgenes)
    })
  
    a=chisq_res2[mask,]
    genes_to_show=head(rownames(a)[order(a[,2],decreasing=T)],ngenes_to_show)
    genes_to_show_comma_delimited=paste(genes_to_show,collapse=",")
    new_set_name=paste("Chisq_",pref,"_",ngenes_to_show,"_",date(),sep="")
    
    message("chisq done")

    if (!is.null(session$userData$geneList)){
      geneList=session$userData$geneList
      geneList[length(geneList)+1]=genes_to_show_comma_delimited
    }
    else{
      geneList=c(genes_to_show_comma_delimited)
    }
    names(geneList)[length(geneList)]=new_set_name
    
    session$userData$geneList<-geneList
    updateSelectInput(session,"inGeneSets",choices=names(session$userData$geneList),selected = names(session$userData$geneList)[length(session$userData$geneList)])
     
  })
  
  chisq_genes=function(counts){
    counts=apply(counts,2:3,sum)
    gene_mask=apply(counts,1,max)>3
    counts=counts[gene_mask,]
    cluster_tot=colSums(counts)
    arrcont=array(c(counts,matrix(cluster_tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts),dim=c(dim(counts),2))
    suppressWarnings({ res=t(apply(arrcont,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}))})
    rownames(res)=rownames(counts)
    res=res[!is.na(res[,1]),]
    res=cbind(res,adjp=p.adjust(res[,1],method="BH"))
    return(res)
  }   
  
  observeEvent(input$inFC, {
    if (input$inSelectGenesFrom=="All genes"){
      pref="Var"
      genes=rownames(session$userData$model$models)
    } else if (input$inSelectGenesFrom=="TFs"){
      pref="Tfs"
      genes=tfs
    } else if(input$inSelectGenesFrom=="Surface markers"){
      pref="Surf"
      genes=surface_markers
    }
    mask=rownames(session$userData$model$umitab)%in%genes
    
    clusts_fg=strsplit(input$inFC_fgClusts,",")[[1]]
    clusts_bg=strsplit(input$inFC_bgClusts,",")[[1]]
    mask_fg=session$userData$model$cell_to_cluster%in%clusts_fg
    mask_bg=session$userData$model$cell_to_cluster%in%clusts_bg
    reg=20
    sfg=reg+rowSums(session$userData$model$umitab[mask,mask_fg])
    sbg=reg+rowSums(session$userData$model$umitab[mask,mask_bg])
    fc=(sfg/sum(sfg))/(sbg/sum(sbg))
    
    isolate({
      ngenes_to_show=as.numeric(input$inNgenes)
    })
    
    genes_to_show=head(names(sort(fc,decreasing=T)),floor(ngenes_to_show/2))
    genes_to_show=c(genes_to_show,tail(names(sort(fc,decreasing=T)),ceiling(ngenes_to_show/2)))
    genes_to_show_comma_delimited=paste(genes_to_show,collapse=",")

    if (!is.null(session$userData$geneList)){
      geneList=session$userData$geneList
      geneList[length(geneList)+1]=genes_to_show_comma_delimited
    }
    else{
      geneList=genes_to_show_comma_delimited
    }
    names(geneList)[length(geneList)]=paste("FC_",pref,"_",ngenes_to_show,"_",paste(clusts_fg,collapse="_"),"_VS_",paste(clusts_bg,collapse="_"),sep="")
    
    session$userData$geneList<-geneList
     updateSelectInput(session,"inGeneSets",choices=names(session$userData$geneList),selected = names(session$userData$geneList)[length(session$userData$geneList)])
  })
  
  
  observeEvent(input$inResetClusters, {
     clust_title=paste(session$userData$default_clusters," - ",session$userData$clustAnnots[session$userData$default_clusters],sep="") 
    update_clusters(session,session$userData$default_clusters)
    updateSelectInput(session,"inQCClust",choices=clust_title)
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title)
    
  })
  
  observeEvent(input$inSaveClusterOrder, {
    clusters=clusters_reactive()
   
    order_fn=paste(strsplit(session$userData$loaded_model_file,"\\.")[[1]][1],"_order.txt",sep="")
    if (all(colnames(session$userData$model$models)%in%clusters)){
      session$userData$default_clusters<-clusters
      write.table(clusters,file=order_fn,quote = F,sep="\t",row.names = F,col.names = F)
      message("Order has been successfully saved.")
    }else{
      message("Order could not be saved since 1 or more clusters are missing.")
    }
  })
  
  observeEvent(input$inUndoClusterChanges, {
    update_clusters(session,session$userData$scDissector_params$previous_clusters,T)
  })
  
  observeEvent(input$detectDiffExrs, {
    if (is.null(session$userData$model$models)){
      return()
    }
    clusters=clusters_reactive()
    cell_mask=session$userData$model$cell_to_cluster%in%clusters
   #tot (=#mols per cluster)
     tot=sapply(split(colSums(session$userData$model$umitab[rownames(session$userData$model$ds),cell_mask]),session$userData$model$cell_to_cluster[cell_mask][colnames(session$userData$model$umitab)[cell_mask]]),sum)
 
    fc=session$userData$model$models[rownames(session$userData$model$ds),]/(rowSums(session$userData$model$ds)/sum(session$userData$model$ds))
    fc_mask=apply(fc,1,max)>input$inputMinFC&rowMeans(session$userData$model$ds)>input$inMinAvgExprs
    print(sum(fc_mask))
    counts=sapply(split(as.data.frame(t(session$userData$model$umitab[rownames(session$userData$model$ds)[fc_mask],cell_mask])),session$userData$model$cell_to_cluster[cell_mask][colnames(session$userData$model$umitab)[cell_mask]]),colSums,na.rm=T)
  
    arrcont=array(c(counts,matrix(tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts),dim=c(dim(counts),2))
    res=t(apply(arrcont,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}))
    rownames(res)=rownames(counts)
    res=res[!is.na(res[,1]),]
    res=cbind(res,adjp=p.adjust(res[,1],method="BH"))
    
    chisq.res<-res[res[,3]<0.05,]
    session$userData$chisq.res<-chisq.res[order(chisq.res[,2],decreasing=T),]
  
    fdrmask=res[,3]<input$inputFDR
    names(fdrmask)=rownames(res)
    
    signiftab<-log2(2^-20+session$userData$model$models[names(fdrmask),])
    session$userData$signiftab<-round(signiftab[order(apply(signiftab,1,which.max)),],digits=1)
    print("done.")
  })
  
  output$downloadSignifTable <- downloadHandler(
    filename = function() { paste("DiffExprs", '.csv', sep='') },
    
    content = function(file) {
      write.csv(session$userData$signiftab, file)
    }
  )
  outputOptions(output, 'downloadSignifTable', suspendWhenHidden=FALSE)  
  
  output$downloadExprsTable <- downloadHandler(
    filename = function() { paste("Exprs", '.csv', sep='') },
    content = function(file) {
      write.csv(session$userData$model$models, file)
    }
  )
  outputOptions(output, 'downloadExprsTable', suspendWhenHidden=FALSE)   
  
  
  
  
  output$downloadCellCountsTable <- downloadHandler(
    filename = function() { paste("cell_counts", '.csv', sep='') },
    
    content = function(file) { 
      cellcounts_tab=table(unlist(dataset$cell_to_cluster),rep(names(dataset$cell_to_cluster),sapply(dataset$cell_to_cluster,length)))
      write.csv(cellcounts_tab, file)
    }
  
  )
  outputOptions(output, 'downloadCellCountsTable', suspendWhenHidden=FALSE)
  
  observeEvent(input$inOrderClusters, {
    inclusts=clusters_reactive()
    ingenes=genes_reactive()
    
    if (length(ingenes)==0){
      return()
    }
    mat<-session$userData$model$models[match(ingenes,rownames(session$userData$model$models)),inclusts]
    if (input$inReorderingClustersMethod=="Hierarchical clustering"){
    cormat=cor(mat,use="comp")
    d1=dist(1-cormat)
    d1[is.na(d1)]=100
    reorderv1=hclust(d1)$order
    
    }
    else if (input$inReorderingClustersMethod=="Diagonal"){
      reorderv1=order(apply(mat[ingenes,],2,which.max))
    }
    else if (input$inReorderingClustersMethod=="Formula"){
      s=input$inOrderClusetersByGenes
      s1=substr(s,1,1)
      genes=strsplit(s,"[+-]")[[1]]
      c1=strsplit(s,"")[[1]]
      sig=c1[c1=="+"|c1=="-"]
      if (s1!="+"&s1!="-"){
        sig=c("+",sig)
      }
      #paste(colnames(model$models)[order(-1e6*model$models["Mpo",]+1e6*model$models["Car2",]+1e12*model$models["Hlf",])],collapse = ",")
      score=rep(0,length(inclusts))
      for (i in 1:length(genes)){
        message(sig[i],genes[i])
  
        score=score+(10^(5*(length(genes)-i)))*ifelse(sig[i]=="+",1,-1)*(mat[genes[i],inclusts])
      }
      
      reorderv1=order(score)
    }
    clust_title=paste(inclusts," - ",session$userData$clustAnnots[inclusts],sep="") 
    
    update_clusters(session,inclusts[reorderv1])
    updateSelectInput(session,"inQCClust",choices=clust_title[reorderv1])
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title[reorderv1])
  })
  
  
  observeEvent(input$inBirdDefineClusterSet,{
    clusts_to_add=strsplit(input$inBirdAddClusters,",| ,|, ")[[1]]
    if (!all(clusts_to_add%in%colnames(session$userData$model$models))){
      updateTextInput(session,"inBirdAddClusters","Clusters:",value = paste(clusts_to_add[clusts_to_add%in%colnames(session$userData$model$models)],collapse = ","))  
      message("Warning! Cluster list contained unknown clusters. Cluster-set was not defined.")
      return()
    }
    intersect_with_defined=intersect(clusts_to_add,unlist(session$userData$scDissector_params$cluster_sets))
    if (length(intersect_with_defined)>0){
      message("warning! Cluster list contained clusters that have been already defined as cluster sets: ",paste(intersect_with_defined,collapse=","),". Cluster-set was not defined.")
      return()
    }
    
    if (input$inBirdAddClusterSetName%in%names(session$userData$scDissector_params$cluster_sets)){
      message("warning! Name has already been used. Cluster-set was not defined.")
      return()
    }
    
    session$userData$scDissector_params$cluster_sets[[input$inBirdAddClusterSetName]]<-clusts_to_add
    updateSelectInput(session,"inBirdRemoveClusterSetSelect","Cluster Set",choices =names(session$userData$scDissector_params$cluster_sets)) 
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =names(session$userData$scDissector_params$cluster_sets)) 
   })
  
  observeEvent(input$inBirdUndefineClusterSet,{
    session$userData$scDissector_params$is_cluster_excluded[session$userData$scDissector_params$cluster_sets[[input$inBirdUndefineClusterSetSelect]]]=F
    session$userData$scDissector_params$cluster_sets[[input$inBirdUndefineClusterSetSelect]]<-NULL
    session$userData$scDissector_params$excluded_cluster_sets<-setdiff(session$userData$scDissector_params$excluded_cluster_sets,input$inBirdUndefineClusterSetSelect)
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(session$userData$scDissector_params$cluster_sets),session$userData$scDissector_params$excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =session$userData$scDissector_params$excluded_cluster_sets)
    updateSelectInput(session,"inBirdUndefineClusterSetSelect","Cluster Set",choices =names(session$userData$scDissector_params$cluster_sets)) 
  })
  
  observeEvent(input$inBirdSaveClusterSets,{
    f=paste(input$inDatapath, session$userData$vers_tab$path[ session$userData$vers_tab$title==session$userData$loaded_model_version],sep="/")
  
    clusters_set_tab=data.frame(name=names(session$userData$scDissector_params$cluster_sets),clusters=sapply(session$userData$scDissector_params$cluster_sets,paste,collapse = ","),is_excluded=ifelse(names(session$userData$scDissector_params$cluster_sets)%in%session$userData$scDissector_params$excluded_cluster_sets,"T","F"))
    clusterset_fn=paste(strsplit(f,"\\.")[[1]][1],"_clustersets.txt",sep="")
    write.table(file=clusterset_fn,clusters_set_tab,row.names=F,col.names=T,quote=F,sep="\t")
   
  })
  
  observeEvent(input$inBirdClusterSetExclude,{
    
    session$userData$scDissector_params$is_cluster_excluded[session$userData$scDissector_params$cluster_sets[[input$inBirdClusterSetsExcludeSelect]]]<-T
    session$userData$scDissector_params$excluded_cluster_sets<-unique(c(session$userData$scDissector_params$excluded_cluster_sets,input$inBirdClusterSetsExcludeSelect))
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(scDissector_params$cluster_sets),session$userData$scDissector_params$excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =scDissector_params$excluded_cluster_sets)
    
  })
  
  observeEvent(input$inBirdClusterSetInclude,{
    
    session$userData$scDissector_params$is_cluster_excluded[session$userData$scDissector_params$cluster_sets[[input$inBirdClusterSetsIncludeSelect]]]<-F
    session$userData$scDissector_params$excluded_cluster_sets<-setdiff(session$userData$scDissector_params$excluded_cluster_sets,input$inBirdClusterSetsIncludeSelect)
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(scDissector_params$cluster_sets),scDissector_params$excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =scDissector_params$excluded_cluster_sets)
  })
  
  ###########################################################
  modules_varmean_reactive=reactive({
    dataset=session$userData$dataset
    ds_i=match(input$inModulesDownSamplingVersion,dataset$ds_numis)
    
    if (!session$userData$loaded_flag)
    {
      return(data.frame(m=0,v=0,gene=""))
    }
    
    inclusts=clusters_reactive()
    cluster_cells_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster%in%inclusts]
    ds=dataset$ds[[ds_i]][,intersect(cluster_cells_mask,dataset$randomly_selected_cells[[ds_i]][[match("All",params$nrandom_cells_per_sample_choices)]])]
    
    ds_mean<-rowMeans(ds)
    genemask=ds_mean>10^input$inVarMeanXlim[1]&ds_mean<10^input$inVarMeanXlim[2]
    ds=ds[genemask,]
    ds_mean=ds_mean[genemask]
    message("Estimating variance for ",nrow(ds)," genes")
    ds_var<-apply(ds,1,var)
    df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
    rownames(df)=names(ds_mean)
    return(df)
  })
  
  reactiveLoess<-reactive({
    varmean_df=modules_varmean_reactive()
    m=varmean_df$m
    v=varmean_df$v
    x=log10(m)
    breaks=seq(min(x),max(x),.2)
    lv=log2(v/m)
    z=sapply(split(lv,cut(x,breaks)),min,na.rm=T)
    maskinf=is.infinite(z)
    z=z[!maskinf]
    b=breaks[-length(breaks)]
    b=b[!maskinf]
    lo=loess(z~b)
    return(lo)
  })
  
  
  geneModuleMask_reactive<-reactive({
    if (!session$userData$loaded_flag)
    {
      return()
    }
    df=modules_varmean_reactive()
    lo=reactiveLoess()
    lline=predict(lo,newdata =log10(df$m))
    geneModuleMask<-log10(df$m)>as.numeric(input$inVarMean_MeanThresh)&log2(df$v/df$m)>lline+as.numeric(input$inVarMean_varmeanThresh)
    geneModuleMask[is.na(geneModuleMask)]<-F
    names(geneModuleMask)=rownames(df)
    return(geneModuleMask)
    
  })
  
  
  output$varMeanThreshPlot <- renderPlot({
    if (!session$userData$loaded_flag)
    {
      return()
    }
    df=modules_varmean_reactive()
    lo=reactiveLoess()
    plot(log10(df$m),log2(df$v/df$m),xlab="Log10(mean)",ylab="log2(var/mean)",panel.first=grid())
    x1=seq(input$inVarMeanXlim[1],input$inVarMeanXlim[2],l=100)
    lline=predict(lo,newdata =x1)
    lines(x1,lline+as.numeric(input$inVarMean_varmeanThresh),col=2)
    abline(v=input$inVarMean_MeanThresh,col=2)
    
    #abline(h=input$inVarMean_varmeanThresh,col=2)
    
    n=sum(geneModuleMask_reactive())
    legend("topright", paste(n,"genes"), bty="n",text.col=2) 
    
  })
  
  
  module_cor_reactive=reactive({
    if (!session$userData$loaded_flag)
    {
      return()
    }
    inclusts=clusters_reactive()
    dataset=session$userData$dataset
    cluster_cells_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster%in%inclusts]
    
    ds_i=match(input$inModulesDownSamplingVersion,dataset$ds_numis)
    #  ds=dataset$ds[[ds_i]][,dataset$randomly_selected_cells[[ds_i]][[match("All",params$nrandom_cells_per_sample_choices)]]]
    ds=dataset$ds[[ds_i]][,intersect(cluster_cells_mask,dataset$randomly_selected_cells[[ds_i]][[match("All",params$nrandom_cells_per_sample_choices)]])]
    
    message("calculating gene-to-gene correlations..")
    cormat=cor(as.matrix(t(log2(.1+ds[names(which(geneModuleMask_reactive())),]))),use="comp")
    return(cormat)
  })
  
  
  
  observeEvent(input$inGetModules,{
    if (!session$userData$loaded_flag)
    {
      return()
    }
    geneModuleMask=geneModuleMask_reactive()
    if (is.null(geneModuleMask)){
      return()
    }
    cormat=module_cor_reactive()
    if (is.null(cormat)){
      return()
    }
    c2c=cutree(hclust(as.dist(1-cormat)),k=as.numeric(input$inNUmberOfGeneModules))
 
    session$userData$modules<-(split(names(c2c),c2c))
    updateTextInput(session,"inModules",value=paste(names(session$userData$modules),collapse=","))
    updateSelectInput(session,"inModuleSelect",label="Show Module:",choices=names(session$userData$modules))
    print(session$userData$modules)
  })
  
  output$textModuleGenes<- renderText({
    input$inModules
    if (!is.null(session$userData$modules)){
      paste(session$userData$modules[[input$inModuleSelect]],collapse=",")
    }
  })
  
  ###########################################################
  
  output$ncells_barplot <-renderPlot({
    input$inModelVer
    clusters=clusters_reactive()
    cells=names(session$userData$model$cell_to_cluster)[session$userData$model$cell_to_cluster%in%clusters]
    ncells=table(session$userData$model$cell_to_cluster[cells])[clusters]
    par(mar=c(1,7,2,2))
    barplot(ncells,xaxs = "i",log="y",ylab="#cells",names.arg ="")
  })
  
  output$UMI_boxplot <- renderPlot({
    input$inModelVer
    par(xaxs="i")
    clusters=clusters_reactive()
    cells=names(session$userData$model$cell_to_cluster)[session$userData$model$cell_to_cluster%in%clusters]
    par(mar=c(2,7,1,2))
    cells=intersect(cells,colnames(session$userData$model$umitab))
    boxplot(split(log10(colSums(session$userData$model$umitab[,cells])),session$userData$model$cell_to_cluster[cells])[clusters],las=2,ylab="log10(#UMIs)")
  })
  
  output$BatchHeatmap <- renderPlot({
    clusters=clusters_reactive()
    cells=names(session$userData$model$cell_to_cluster)[session$userData$model$cell_to_cluster%in%clusters]
     if(length(unique(session$userData$model$cell_to_batch[cells]))==1){
      return()
     }
    if (is.null(unique(session$userData$model$cell_to_batch[cells]))){
      return(NULL)
    }
    obs=table(session$userData$model$cell_to_cluster[cells],session$userData$model$cell_to_batch[cells])[clusters,]
    cs=colSums(obs)
    rs=rowSums(obs)
    expected=(matrix(cs,nrow(obs),ncol(obs),byrow=T)/sum(cs))*rs
    reg=5
    logr=log2((reg+obs)/(reg+expected))
    par(mar=c(10,7,4,2))
    
   # sort(apply(obs,1,function(x){chisq.test(x,simulate.p.value = 1000)$p.value}))
    image(logr[,ncol(logr):1],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=c(-100,seq(-2,2,l=99),100),main="Batch enrichment")
    mtext(text = session$userData$clustAnnots[rownames(logr)],side = 1,at = seq(0,1,l=dim(logr)[1]),las= 2,cex=1)
    # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =colnames(logr), side=2,  at=seq(1,0,l=dim(logr)[2]),las=2,cex=1)
    #print(obs)
    box()
  })
  
 
  output$gaiting_plots_dynamic <- renderUI({
    if (!is.null(session$userData$dataset$insilico_gating_scores)){
      nplots=length(session$userData$dataset$insilico_gating_scores)
      he=max(c(500,500*nplots,na.rm=T))
      plotOutput("gating_plots", width = "100%", height = he)
    }
  })
  
  
  output$gating_plots <- renderPlot({
    nplots=length(session$userData$dataset$insilico_gating_scores)
    samp=input$inGatingSample
    numis=session$userData$dataset$numis_before_filtering[[samp]]
    clust=strsplit(input$inGatingShowClusters," ")[[1]][1]
    cell_to_cluster=session$userData$dataset$cell_to_cluster[session$userData$dataset$cell_to_sample==samp]
    layout(matrix(1:(nplots),nplots,1))
    for (i in 1:nplots){
      mask=intersect(names(numis),names(session$userData$dataset$insilico_gating_scores[[i]]))
      plot(numis[mask],session$userData$dataset$insilico_gating_scores[[i]][mask],log="x",ylab=paste("fraction",names(session$userData$model$insilico_gating)[i]),xlab="#UMIs",col=ifelse(is.na(cell_to_cluster[mask]),"gray",ifelse(cell_to_cluster[mask]==clust,2,1)))
      points(numis[mask],session$userData$dataset$insilico_gating_scores[[i]][mask],pch=ifelse(cell_to_cluster[mask]==clust,20,NA),col=2)
      rect(xleft = session$userData$dataset$min_umis,session$userData$model$insilico_gating[[i]]$interval[1],session$userData$dataset$max_umis,session$userData$model$insilico_gating[[i]]$interval[2],lty=3,lwd=3,border=2)
    }

   
  })
  
  
  output$samples_enrichment <- renderUI({
    inclusts=clusters_reactive()
    insamples=samples_reactive()
    he=max(c(500,12*length(inclusts)),na.rm=T)
    if (length(insamples)>1&(length(inclusts)>1)){
      plotOutput("samples_enrichment_plot",height=he)
    }
    
  })
  
  output$cluster_sizes<-renderUI({
    inclusts=clusters_reactive()
    he=max(c(500,12*length(inclusts)),na.rm=T)
    plotOutput("cluster_sizes_plot",height=he)
  })
  
  
  output$avg_profile_plot <- renderUI({
    inclusts=clusters_reactive()
    insamples=samples_reactive()
    he=max(c(500,12*length(inclusts)),na.rm=T)
    if (length(insamples)>1&(length(inclusts)>1)){
      plotlyOutput("avg_heatmap_interactive", width = "100%", height = he)
    # d3heatmapOutput("avg_heatmap_interactive", width = "150%", height = he)
    }
    else{
 #   d3heatmapOutput("avg_heatmap_interactive", width = "150%", height = he)
      plotlyOutput("avg_heatmap_interactive", width = "100%", height = he)
    }
   })
  
  
  output$avg_module_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("avg_heatmap_modules", width = "100%", height = he)
  })
  
  output$sample_avg_profile_plot <- renderUI({
    he=max(c(300,12*length(samples_reactive())),na.rm=T)
    plotOutput("avg_heatmap_samples", width = "100%", height = he)
  })
  
  
  output$projection_avg_heatmap_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("avg_heatmap_projection",width="100%",height=he)
  })
  
  output$projection_barplot_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("projection_barplot",width="100%",height=he)
  })
  
          
 
  
  output$samples_enrichment_plot <- renderPlot({
    inclusts=clusters_reactive()
    insamples=samples_reactive()
    par(mar=c(7,0,1.6,0))
    tab=matrix(0,ncol(session$userData$model$models),length(session$userData$dataset$samples))
    rownames(tab)=colnames(session$userData$model$models)
    colnames(tab)=session$userData$dataset$samples
    tmptab=table(session$userData$dataset$cell_to_cluster,session$userData$dataset$cell_to_sample)
    tab[rownames(tmptab),colnames(tmptab)]=tmptab
    tab=t(t(tab)/colSums(tab))
    tab=(tab/rowSums(tab))[inclusts,]
    barplot(t(tab[nrow(tab):1,]),col=sample_cols[1:ncol(tab)],horiz =T,yaxs = "i",names.arg=NULL,main="Samples",axes=F)
    
  })
  
  output$cluster_sizes_plot <- renderPlot({
    inclusts=clusters_reactive()
    insamples=samples_reactive()
    if (input$inModelOrAverage=="Model"){
      tab=session$userData$ncells_per_cluster
    }
    else{
      tab=session$userData$ncells_per_cluster*0
      tmp_tab=table(session$userData$dataset$cell_to_cluster[session$userData$dataset$cell_to_sample%in%insamples])
      tab[names(tmp_tab)]=tmp_tab
    }
    precentage=100*tab[inclusts]/sum(tab[setdiff(names(tab),session$userData$scDissector_params$excluded_clusters)])
   
    par(mar=c(5.9,.3,0.5,0))
    barplot(rev(precentage),horiz=T,border=F,col="gray",xaxs="i",xlab="Cells ( % )",cex.axis = .8)
  })
  
  
  # A reactive expression with the ggvis plot
  output$avg_heatmap <- renderPlot({
  ##### Don't delete!!
    input$inAnnotateCluster
    input$inModelVer
    input$inBirdClusterSetExclude
    input$inBirdClusterSetInclude
  #####
    cgs=clusters_genes_sampples_reactive()
    
    inclusts=cgs$clusters
    ingenes=cgs$genes
    insamples=cgs$samples
     if (!session$userData$loaded_flag){
      return()
    }
  
    # Lables for axes
    if (length(ingenes)==0){
        return()
    }
    if (length(inclusts)==0){
      return()
    }

    zlim=input$inModelColorScale
    if (length(insamples)>1&(length(inclusts)>1)){
      layout(matrix(1:2,1,2),widths=c(1,10))
      par(mar=c(7,1,1,.1))
      tab=matrix(0,ncol(session$userData$model$models),length(session$userData$dataset$samples))
      rownames(tab)=colnames(session$userData$model$models)
      colnames(tab)=session$userData$dataset$samples
      tmptab=table(session$userData$dataset$cell_to_cluster,session$userData$dataset$cell_to_sample)
      tab[rownames(tmptab),colnames(tmptab)]=tmptab
      tab=t(t(tab)/colSums(tab))
      tab=(tab/rowSums(tab))[inclusts,]
      barplot(t(tab[nrow(tab):1,]),col=sample_cols[1:ncol(tab)],horiz =T,yaxs = "i",names.arg=NULL,main="Samples",axes=F)
    }
    par(mar=c(7,tab3_left_margin,1,9))
    if (input$inModelOrAverage=="Model"){
      mat<-session$userData$model$models[match(ingenes,rownames(session$userData$model$models)),]
    }
    else {
      gene_match=match(ingenes,dimnames(session$userData$dataset$counts)[[2]])
      mat<-apply(session$userData$dataset$counts[insamples,gene_match,inclusts,drop=F],2:3,sum)

      if (input$inModelOrAverage=="Batch-corrected Average"){
        if (!is.null(session$userData$dataset$noise_counts)){
          

          
          mat_noise=apply(session$userData$dataset$noise_counts[insamples,gene_match,inclusts,drop=F]*arr_weights,2:3,sum)
          mat<-pmax(mat-mat_noise,0)
        }

      }
      rownames(mat)=ingenes
      mat<-t(t(mat)/colSums(mat,na.rm=T))
    }
    isolate({
    mat1=mat[,inclusts,drop=F]
    abs_or_rel=input$inAbsOrRel
    }) 
    main_title=paste(input$inModelOrAverage,":",session$userData$loaded_model_version)
    genes=rownames(mat1)
    gene.cols=session$userData$gcol[toupper(rownames(mat1))]
    clusters=colnames(mat1)
    clusters_text=paste(" (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster[setdiff(names(session$userData$ncells_per_cluster),session$userData$scDissector_params$excluded_clusters)]),digits=1),"% )",sep="")
    annots=session$userData$clustAnnots[inclusts]

    plot_avg_heatmap(mat1,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute=abs_or_rel)
  
  })
  
#  output$avg_heatmap_interactive <-renderD3heatmap({
  output$avg_heatmap_interactive <-renderPlotly({
    
      ##### Don't delete!!
      input$inAnnotateCluster
      input$inModelVer
      input$inBirdClusterSetExclude
      input$inBirdClusterSetInclude
      #####
      
      cgs=clusters_genes_sampples_reactive()
      
      inclusts=cgs$clusters
      ingenes=cgs$genes
      insamples=cgs$samples
  
      if (!session$userData$loaded_flag){
        return()
      }

      # Lables for axes
      if (length(ingenes)==0){
        return()
      }
      if (length(inclusts)==0){
        return()
      }
      
      zlim=input$inModelColorScale
     
      par(mar=c(7,tab3_left_margin,1,9))
      if (input$inModelOrAverage=="Model"){
        mat<-session$userData$model$models[match(ingenes,rownames(session$userData$model$models)),,drop=F]
      }
      else {
        ncells_per_sample=table(session$userData$dataset$cell_to_sample)
        gene_match=match(ingenes,dimnames(session$userData$dataset$counts)[[2]])
        weights=ncells_per_sample[insamples]/sum(ncells_per_sample[insamples])
        arr_weights=array(weights,dim=c(length(insamples),length(gene_match),length(inclusts)))
        mat<-apply(session$userData$dataset$counts[insamples,gene_match,inclusts,drop=F]*arr_weights,2:3,sum)
        numis_per_sample=apply(session$userData$dataset$counts,1,sum)
        
        if (input$inModelOrAverage=="Batch-corrected Average"){
          if (!is.null(session$userData$dataset$noise_counts)){

            mat_noise=apply(session$userData$dataset$noise_counts[insamples,gene_match,inclusts,drop=F]*arr_weights,2:3,sum)
            mat<-pmax(mat-mat_noise,0)
          }
        }
        rownames(mat)=ingenes
        mat<-t(t(mat)/colSums(mat,na.rm=T))
      }
      isolate({
        mat1=mat[,inclusts,drop=F]
        abs_or_rel=input$inAbsOrRel
      }) 
      main_title=paste(input$inModelOrAverage,":",session$userData$loaded_model_version)
      genes=rownames(mat1)
      gene.cols=session$userData$gcol[toupper(rownames(mat1))]
      clusters=colnames(mat1)
      clusters_text=paste(" (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster[setdiff(names(session$userData$ncells_per_cluster),session$userData$scDissector_params$excluded_clusters)]),digits=1),"% )",sep="")
      annots=session$userData$clustAnnots[inclusts]
      return(plot_avg_heatmap_interactive(mat1,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute=abs_or_rel))
  })
  
  module_counts_reactive<-reactive({
    input$inGetModules
    if (is.null(session$userData$modules)){
      return()
    }
    samps=samples_reactive()
    if (length(samps)==0){
      return()
    }
    
    counts=apply(session$userData$dataset$counts[samps,,,drop=F],2:3,sum)
    
    return(t(sapply(session$userData$modules,function(modi){colSums(counts[modi,,drop=F])})/colSums(counts)))
  
    
  })

  output$avg_heatmap_modules <- renderPlot({
    ##### Don't delete!!
    input$inModelVer
    samps=samples_reactive()
    #####
    
    inclusts=clusters_reactive()
    inmodules=modules_reactive()
    if (!session$userData$loaded_flag){
      return()
    }
    if (is.null(session$userData$modules)){
      return()
    }
    modulemat<-module_counts_reactive()
    if(is.null(modulemat)){
      return()
    }
    modulemat<-modulemat[inmodules,inclusts]
    
    zlim=input$inModelColorScale
    par(mar=c(7,tab3_left_margin,1,9))
   
   
    
    mat1=modulemat[,inclusts]
    
    if (input$inAbsOrRel=="Relative"){
      mat_to_show=log2(1e-5+mat1/pmax(1e-5,rowMeans(mat1,na.rm=T)))
      break1=-1e6
      break2=1e6
    }
    else if (input$inAbsOrRel=="Absolute"){
      mat_to_show=log10(mat1)
      break1=0
      break2=1e-1
    }
  
    image(mat_to_show[,ncol(mat1):1],col=colgrad,breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main=session$userData$loaded_model_version)
  
    box()
    
    mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1)
    mtext(text =paste(" ",colnames(mat1)," (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster),digits=1),"% )",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =paste(session$userData$clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
  })
  
  output$avg_heatmap_samples <- renderPlot({
    inclusts=clusters_reactive()
    ingenes=genes_reactive()
    insamples=samples_reactive()
    print(insamples)
    if (!session$userData$loaded_flag){
      return()
    }
    
    # Lables for axes
    if (length(ingenes)==0){
      return()
    }
    
    zlim=input$inSamplesColorScale
    
    par(mar=c(7,7,1,9))
    mat1=t(apply(session$userData$dataset$counts[insamples,ingenes,inclusts],1:2,sum))
    mat1=t(t(mat1)/colSums(mat1))
      if (ncol(mat1)>1){
        mat_to_show=log2(1e-5+mat1/pmax(1e-5,rowMeans(mat1,na.rm=T)))
        break1=-1e6
        break2=1e6
      }
      else{
        return()
      }

    isolate({
      if (length(inclusts)==ncol(session$userData$model$models)){
        clusters_string="All Cells"
      }
      else{
        clusters_string=paste("Clusters ",paste(inclusts,collapse=","))
      }
      image(mat_to_show[,ncol(mat1):1,drop=F],col=colgrad,breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main=paste("Pooled Sample Averages over",clusters_string))
    })
    box()
    
    mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1,col=session$userData$gcol[toupper(rownames(mat1))])
    mtext(text =colnames(mat1),side=2,at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
  })
  
  
  
 
  output$colorLegendModel<-renderPlot({
    par(mar=c(2,1,1,1))
    image(matrix(1:100,100,1),breaks=0:100,col=colgrad,axes=F)
    axs=c(input$inModelColorScale[1],round(mean(input$inModelColorScale),digits=2),input$inModelColorScale[2])
    axis(side =1,at = c(0,.5,1),labels = axs)
  })
  
  output$colorLegendTruth<-renderPlot({
    par(mar=c(2,1,1,1))
    axs=c(input$inTruthColorScale[1],input$inTruthColorScale[2])
    if (!is.null(axs)){
      image(matrix(1:100,100,1),breaks=0:100,col=colgrad,axes=F)
      axis(side =1,at = c(0,1),labels = axs)
    }
  })
  
  output$truthplot <- renderPlot({
   
    zlim=input$inTruthColorScale
   
    cgs=clusters_genes_sampples_reactive()
    inclusts=cgs$clusters
    ingenes=cgs$genes
    insamples=cgs$samples
    nclust=length(inclusts)
    if (!session$userData$loaded_flag){
      return()
    }
    if (length(ingenes)==0){
        
      return()
    }

    samps_cl=list() 
    ds=session$userData$dataset$ds[[match(input$inTruthDownSamplingVersion,session$userData$dataset$ds_numis)]]
    
    cells_per_sample=as.numeric(input$inTruthNcellsPerSample)
    genes=intersect(ingenes,rownames(ds))
    cells_selected=session$userData$dataset$randomly_selected_cells[[match(input$inTruthDownSamplingVersion,session$userData$dataset$ds_numis)]][[match(input$inTruthNcellsPerSample,params$nrandom_cells_per_sample_choices)]]
    cell_mask=session$userData$dataset$cell_to_cluster[colnames(ds)]%in%inclusts & 
    session$userData$dataset$cell_to_sample[colnames(ds)]%in%insamples &colnames(ds)%in%cells_selected
    ds=ds[genes,cell_mask]
    ds=ds[,order(match(session$userData$dataset$cell_to_cluster[colnames(ds)],inclusts))]
    samps=session$userData$dataset$cell_to_sample[colnames(ds)]
    ncells=rep(0,length(inclusts))
    names(ncells)=inclusts
    tmp_ncells=sapply(split(colnames(ds),session$userData$dataset$cell_to_cluster[colnames(ds)]),length)
    ncells[names(tmp_ncells)]=tmp_ncells
   
    #ncells=sapply(ds_cl,function(x){n=ncol(x);if(is.null(n)){return(0)};return(n)})
   # names(ncells)=names(ds_cl)
    
    
    pmat=as.matrix(ds)[,ncol(ds):1]
    spacer_size=ceiling(dim(pmat)[2]/200)
    pmat2=log2(1+pmat)

    layout(matrix(1:2,1,2),widths=c(40,1))   
    par(mar=c(10,3,1,1))
   # image(pmat2[,ncol(pmat2):1],col=c("gray",colgrad),axes=F,breaks=c(-3e6,-1e6,seq(zlim[1],zlim[2],l=99),1e6))
    image(pmat2,col=c("gray",colgrad),axes=F,breaks=c(-3e6,-1e6,seq(zlim[1],zlim[2],l=99),1e6))
    
    
    box()
    if (input$inTruthShowSeparatorBars){
      abline(h=1-cumsum(ncells)/sum(ncells),col="gray")
    }
    mtext(text =rownames(pmat), side=1, at=seq(0,1,l=dim(pmat)[1]),las=2,cex=1,col=session$userData$gcol[toupper(rownames(pmat))])
    a=cumsum(ncells)
    b=a-floor(ncells[inclusts[inclusts%in%names(ncells)]]/2)
    mtext(text =inclusts, side=2, at=1-(b/a[length(ncells)]),las=2,cex=1,adj=1)
    par(mar=c(10,0,1,1))
    image(t(as.matrix(match(rev(samps),insamples))),axes=F,breaks=0:length(insamples)+.5,col=sample_cols[1:length(insamples)])
  })
    

    varmean_reactive <- reactive({
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      samples=samples_reactive()
      if (is.null(session$userData$dataset)){
        return(data.frame(m=0,varmean=0,gene=""))
      }
      dataset=session$userData$dataset

      ds=dataset$ds[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]]
      cluster_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster==clust]
      ds=ds[,intersect(colnames(ds),cluster_mask),drop=F]

      if (is.null(ds))(return(data.frame(m=0,varmean=0,gene=0)))
      if (ncol(ds)<2)(return(data.frame(m=0,varmean=0,gene=0)))
      
      m=rowMeans(ds)
      logm=log10(m)
      mask1=logm>input$mean[1]&logm<input$mean[2]
      v=apply(ds[mask1,],1,var)
      df=data.frame(m=logm[mask1],varmean=log2(v/m[mask1]),gene=rownames(ds)[mask1])
      mask2=df$varmean<input$varmean[2]&df$varmean>input$varmean[1]
      df[mask2,]
      #  
      # Apply filters
    })
    
    
    output$modelSampleLegend <- renderPlot({
    
      insamples=samples_reactive()
      par(mar=c(0,0,0,0))
      plot.new()
      leg=session$userData$samples_tab[match(insamples,session$userData$samples_tab$index),"title"]
      if(length(leg)==0){
        return()
      }
      ncol=ceiling(length(insamples)/10)
      leg[is.na(leg)]=insamples[is.na(leg)]
      leg[leg==""]=insamples[leg==""]
      legend("topleft",pch=15,col=sample_cols[1:length(insamples)],legend=leg,cex=1,xpd=T,ncol=ncol)
    })
    
    
    output$truthSampleLegend <- renderPlot({
      insamples=samples_reactive()
      par(mar=c(0,0,0,0))
      plot.new()
      leg=session$userData$samples_tab[match(insamples,session$userData$samples_tab$index),"title"]
      if(length(leg)==0){
        return()
      }
      ncol=ceiling(length(insamples)/10)
      leg[is.na(leg)]=insamples[is.na(leg)]
      legend("topleft",pch=15,col=sample_cols[1:length(insamples)],legend=leg,cex=1,xpd=T,ncol=ncol)
    })
    
    # Function for generating tooltip text
    gene_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$gene)) return(NULL)
      
      # Pick out the movie with this ID
      
      paste0("<b>", x$gene, "</b>")
    }
    
    click_tooltip <- function(x) {
      
    }
    
    
    
    # A reactive expression with the ggvis plot
    vis <- reactive({
      # Lables for axes
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      # Normally we could do something like props(x = ~BoxOffice, y = ~Reviews),
      # but since the inputs are strings, we need to do a little more work.
      #   xvar <- prop("x", as.symbol(input$xvar))
      #    yvar <- prop("y", as.symbol(input$yvar))
      dataset=session$userData$dataset
      varmean_reactive %>%
        ggvis(x = ~m, y = ~varmean) %>%
        layer_points(size := 50, size.hover := 200,
                     fillOpacity := 0.2, fillOpacity.hover := 0.5, key := ~gene) %>%
        add_tooltip(gene_tooltip, "hover")%>% 
        add_tooltip(click_tooltip, "click")%>% 
        add_axis("x", title = "log10(mean)") %>%
        add_axis("y", title = "log2(var/mean)") %>%
        add_axis("x", orient = "top", ticks = 0, title = paste(session$userData$loaded_model_version,clust),
                 properties = axis_props(
                   axis = list(stroke = "white"),
                   labels = list(fontSize = 0))) %>%
        scale_numeric("x", domain = c(input$mean[1], input$mean[2])) %>%
        scale_numeric("y", domain = c(input$varmean[1],input$varmean[2])) %>%
        set_options(width = 500, height = 500)
      
      
    })
    
    vis %>% bind_shiny("varmeanplot")

    output$cellcor<-renderPlot({
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      insamples=samples_reactive()
      ds=ds_QC_reactive()
      dataset=session$userData$dataset
      if (is.null(ds)){
        return()
      }
      if (ncol(ds)<2){
        message("Warning: Cluster ",clust, " has less than 2 cells above the UMIs threshold.")
        return()
      }

      v=apply(ds,1,var)
      m=rowMeans(ds)
      genemask=rownames(ds)[m>1e-1]
      var_genes=head(genemask[order(v[genemask]/m[genemask],decreasing=T)],200)
      z=ds[var_genes,]
      z=z/mean(z)
      
      cormat=cor(as.matrix(z),use="comp")
 
      ord=hclust(dist(1-cormat))$order
      samps=dataset$cell_to_sample[colnames(cormat)[ord]]
      layout(matrix(c(1:6),nrow = 3),heights = c(8,1,1),widths=c(9,1))
      par(mar=c(.5,2,2,2))
      image(cormat[ord,ord],col=greenred(100),breaks=c(-1,seq(-1,1,l=99),1),axes=F,main=paste("Cluster",clust,": ",ncol(cormat),"/",sum(dataset$cell_to_cluster==clust),"cells"))
      lab=paste("Single-cells (cluster ",clust,")",sep="")

      
      par(mar=c(.5,2,.5,2))
      gene1=match(input$inGene1,rownames(ds))
      gene2=match(input$inGene2,rownames(ds))
     
      if (!is.na(gene1)){
        x=t(ds[gene1,colnames(cormat),drop=F])
        x=log2(1e-5+x/(pmax(1e-5,mean(x,na.rm=T))))
        image(as.matrix(x),breaks=c(-1e6,seq(-3,3,l=99),1e6),col=colgrad,axes=F)
        mtext(2,text = input$inGene1)
      }
      else{
        plot.new()
      }
      if (!is.na(gene2)){
        x=t(ds[gene2,colnames(cormat),drop=F])
        x=log2(1e-5+x/(pmax(1e-5,mean(x,na.rm=T))))
        image(as.matrix(x),breaks=c(-1e6,seq(-3,3,l=99),1e6),col=colgrad,axes=F)
        mtext(2,text = input$inGene2)
      }
      else{
        plot.new()
      }

      par(mar=c(.5,.5,2,.5))
      image(t(as.matrix(match(samps,insamples))),axes=F,breaks=0:length(insamples)+.5,col=sample_cols[1:length(insamples)])
      mtext(side=3,text = "Sample",cex=.8)
  
    })
    
    click_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$gene)) return(NULL)
      
      if (session$userData$click_flag){
        updateTextInput(session,"inGene1",,x$gene)
      }
      else{
        updateTextInput(session,"inGene2",,x$gene)
      }
      session$userData$click_flag<-!session$userData$click_flag
    } 
    
    output$n_movies <- renderText({ nrow(movies()) })
    
    output$gene2 <- renderText({ paste("_________",input$inGene2,sep="")})
    output$gene1 <- renderText({ input$inGene1})  
   
    output$table <- renderTable({
      if (!session$userData$loaded_flag){
        return()
      }
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      dataset=session$userData$dataset
      ds=ds_QC_reactive()
      if (ncol(ds)==0){
        return()
      }
      mask=dataset$cell_to_cluster[colnames(ds)]%in%clust
      flag=sum(c(input$inGene1,input$inGene2)%in%rownames(ds))==2
      if (flag){
        dat=floor(log2(1+data.frame(t(as.matrix(ds[c(input$inGene1,input$inGene2),mask])))))
        tab=table(dat)
        x=matrix(0,1+max(dat[,1],na.rm=T),1+max(dat[,2],na.rm=T))
        rownames(x)=0:(nrow(x)-1)
        colnames(x)=0:(ncol(x)-1)
        x[rownames(tab),colnames(tab)]=tab
      }
      else{
        x=matrix(NA,2,2)
        rownames(x)=0:1
        colnames(x)=0:1
      }
      x
    },rownames=T,colnames=T,width=18)
    

    
    
    # A reactive expression with the ggvis plot
    output$avg_heatmap_projection <- renderPlot({
      ##### Don't delete!!
      #####
      insamples=samples_reactive()
      inclusts=clusters_reactive()
      ingenes=genes_reactive()
      
      
      samples1=strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]]
      samples2=strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]]

      if (!((samples1%in%insamples|all(session$userData$sample_sets[[samples1]]%in%insamples))&&input$inProjectSampleGroup2%in%insamples|all(session$userData$sample_sets[[samples2]]%in%insamples))){
         return()
      }
      
      # Lables for axes
      if (length(ingenes)==0){
        
        updateSelectInput(session,"inGeneSets",selected = names(session$userData$geneList)[1])
        updateTextInput(session,"inGenes",value = init_genes(session$userData$geneList[1]))
        
        return()
      }
      
  
      zlim=input$inModelColorScale
      par(mar=c(7,7,1,9))
      
      mat1<-session$userData$model$models[match(ingenes,rownames(session$userData$model$models)),inclusts]
   
      n=ncol(mat1)
      ncells_per_projected_cluster=rep(0,length(inclusts))
      names(ncells_per_projected_cluster)=inclusts
      t1=table(projected$cell_to_cluster)
      ncells_per_projected_cluster[names(t1)]=t1
      
      percents_ref=round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster),digits=1)
      percents_projected=round(100*ncells_per_projected_cluster[inclusts]/sum(ncells_per_projected_cluster),digits=1)
      
      
      layout(matrix(1:2,1,2))
        if (input$inProjectPlotType=="Tile"){
        matc=cbind(mat1,mat1p)[,rep(1:n,each=2)+rep(c(0,n),n)]
        isolate({
          image(log2(1e-5+matc/pmax(1e-5,rowMeans(matc)))[,ncol(matc):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=paste(session$userData$loaded_model_version,input$inProjectedDataset,sep="/"))
        })
        box()
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
     
        sapply(yv,function(y){arrows(xusr[1],y,xusr[2],y,length=0,col="gray")})
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=session$userData$gcol[toupper(rownames(mat1))])
        mtext(text =paste(" ",colnames(mat1)," (",percents_ref,"% /",percents_projected,"% )",sep=""), side=4, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(session$userData$clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
      }
      else if (input$inProjectPlotType=="Side by Side"){
        avg=rowMeans(cbind(mat1,mat1p))
        image(log2(1e-5+mat1/pmax(1e-5,avg))[,ncol(mat1):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=session$userData$loaded_model_version)
        box()
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=session$userData$gcol[toupper(rownames(mat1))])
       # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(" (",percents_ref,"%)",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(session$userData$clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        image(log2(1e-5+mat1p/pmax(1e-5,avg))[,ncol(mat1p):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=input$inProjectedDataset)
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=session$userData$gcol[toupper(rownames(mat1))])
       # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(" (",percents_projected,"%)",sep=""), side=4,  at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(session$userData$clustAnnots[inclusts]," ",sep=""), side=2,  at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        
        box()
      }
    })
    
    output$projection_barplot=renderPlot({
      ##### Don't delete!!
      #####
      insamples=samples_reactive()
      inclusters=clusters_reactive()
      ingenes=genes_reactive()
      
      dataset=session$userData$dataset
     
      samples1=as.list(strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]])
      samples2=as.list(strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]])
      if (!is.null(session$userData$sample_sets)){
        samples1[!(samples1%in%insamples)]=session$userData$sample_sets[unlist(samples1[!samples1%in%insamples])]
        samples2[!(samples2%in%insamples)]=session$userData$sample_sets[unlist(samples2[!samples2%in%insamples])]
      }
      samples1=unlist(samples1)
      samples2=unlist(samples2)
      
      if (is.null(samples1)||is.null(samples2)){
        return()
      }
      if (!(all(samples1%in%insamples)&&all(samples2%in%insamples))){
        return()
      }
      if (is.null(dataset)){
        return()
      }
      full_table=function(cell_to_cluster,clusters){
        tab=rep(0,length(clusters)+1)
        names(tab)=c("-1",clusters)
        tab2=table(cell_to_cluster)
        mask=intersect(names(tab2),clusters)
        tab[mask]=tab2[mask]
        tab["-1"]=sum(tab2[setdiff(names(tab2),clusters)])
        return(tab)
      }
      tab1=t(sapply(dataset$cell_to_cluster[dataset$cell_to_sample%in%samples1],full_table,inclusters))
      tab2=t(sapply(dataset$cell_to_cluster[dataset$cell_to_sample%in%samples2],full_table,inclusters))
      
      tot1=colSums(tab1)   
      tot2=colSums(tab2)
      
      freq1=tot1[setdiff(names(tot1),"-1")]/sum(tot1)
      freq2=tot2[setdiff(names(tot2),"-1")]/sum(tot2)
      
      percents1=round(100*freq1,digits=1)
      percents2=round(100*freq2,digits=1)
   
     
      par(mar=c(7,7,1,9))
      
      
      layout(matrix(1:2,1,2),widths=c(2,1))
      
      m=matrix(0,2,length(inclusters),dimnames=list(c("1","2"),inclusters))
      m[1,names(percents1)]=percents1
      m[2,names(percents2)]=percents2
      
      par(mar=c(5,10,1,1))
      barplot(m[,dim(m)[2]:1],beside=T,names.arg = rev(paste(session$userData$clustAnnots[inclusters]," - ",inclusters,sep="")),horiz = T,las=2)
      legend("bottomright",pch=15,col= gray.colors(2),legend=c("1","2"),border=T)
    #  reg=100*10/((ncol(model$umitab)+ncol(projected$umitab))/2)
      reg=1e-3
      barplot(rev(log2((reg+m[2,])/(reg+m[1,]))),names.arg = rev(paste(session$userData$clustAnnots[inclusters]," - ",inclusters,sep="")),horiz = T,las=2,xlab="log2(2/1)")
    })
    
    click_tooltip_gene_proj_vs_ref <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$gene)) return(NULL)
      
     
      updateTextInput(session,"Gene1ForExprsTableRefVsProj",,x$gene)
      
    } 
    
    ge_proj_vs_ref <- reactive({
      insamples=samples_reactive()
      clust=strsplit(input$inClustForDiffGeneExprsProjVsRef," - ")[[1]][1]
      samples1=as.list(strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]])
      samples2=as.list(strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]])
      dataset=session$userData$dataset
      if (is.null(session$userData$dataset)){
        return(data.frame(x=0,y=0,gene=0))
      }
      if (!is.null(session$userData$sample_sets)){
        samples1[!(samples1%in%insamples)]=session$userData$sample_sets[unlist(samples1[!samples1%in%insamples])]
        samples2[!(samples2%in%insamples)]=session$userData$sample_sets[unlist(samples2[!samples2%in%insamples])]
      }
      samples1=unlist(samples1)
      samples2=unlist(samples2)
      
      if(length(samples1)==0||length(samples2)==0||(!(all(samples1%in%insamples)&&all(samples2%in%insamples)))||is.null(dataset)){
        return(data.frame(x=0,y=0,gene=0))
      }
      counts1=apply(session$userData$dataset$counts[samples1,,,drop=F],2:3,sum)
      counts2=apply(session$userData$dataset$counts[samples2,,,drop=F],2:3,sum)
      
      model1=t(t(counts1)/colSums(counts1))
      model2=t(t(counts2)/colSums(counts2))
      
      df=data.frame(x=log10(model1[,clust]),y=log2(model2[,clust]/model1[,clust]),gene=rownames(session$userData$model$models))
      
    })
    # A reactive expression with the ggvis plot
    vis2 <- reactive({
    
      clust=strsplit(input$inClustForDiffGeneExprsProjVsRef," - ")[[1]][1]
      
      ge_proj_vs_ref %>%
        ggvis(x = ~x, y = ~y) %>%
        layer_points(size := 50, size.hover := 200,
                     fillOpacity := 0.2, fillOpacity.hover := 0.5, key := ~gene) %>%
        add_tooltip(gene_tooltip, "hover")%>% 
        add_tooltip(click_tooltip_gene_proj_vs_ref, "click")%>% 
        add_axis("x", title = "log10(mean)") %>%
        add_axis("y", title = "log2(samples_gr2/samples_gr1)") %>%
        add_axis("x", orient = "top", ticks = 0, title = paste(session$userData$loaded_model_version," ",clust),
                 properties = axis_props(
                   axis = list(stroke = "white"),
                   labels = list(fontSize = 0))) %>%
        scale_numeric("x", domain = c(-6,-1)) %>%
        scale_numeric("y", domain = c(-6,6)) %>%
       set_options(width = 500, height = 500)
      
      
    })
    
    vis2 %>% bind_shiny("DiffGeneExprsProjVsRef")

    
    output$Gene1ExprsTableRefVsProj<- renderTable({
      if (!session$userData$loaded_flag){
        return()
      }
      if (!exists("projected")){
        return()
      }
   
      g=input$Gene1ForExprsTableRefVsProj
      clust=strsplit(input$inClustForDiffGeneExprsProjVsRef," - ")[[1]][1]
      
      if ((g%in%rownames(ds))&(g%in%rownames(projected$ds))){
        mask1=session$userData$model$cell_to_cluster[colnames(ds)]==clust
        mask2=projected$cell_to_cluster[colnames(projected$ds)]==clust
        x1=table(floor(log2(1+data.frame(ds[g,mask1]))))
        x2=table(floor(log2(1+data.frame(projected$ds[g,mask2]))))
        z=max(c(as.numeric(names(x1)),as.numeric(names(x2))))
        x=matrix(0,2,1+z,dimnames=list(c("Reference","Projected"),(0:z)))
        x[1,names(x1)]=x1
        x[2,names(x2)]=x2
        }
      else {
        return(NULL)
      }
      print(x)
      x
    },rownames=T,colnames=T,width=18)
    
    output$ClusteringsComparisonTable<- renderTable({
      if (!session$userData$loaded_flag){
        return()
      }
      if (!all(names(session$userData$model$cell_to_cluster)==names(projected$source_cell_to_cluster))){
        return() 
      }
      tab=table(session$userData$model$cell_to_cluster,projected$source_cell_to_cluster)
      x=matrix(tab,nrow(tab),ncol(tab),dimnames=list(rownames(tab),colnames(tab)))
      print(x)
      x
      
    })
    
    
    click_tooltip_ll <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$cell)) return(NULL)
    } 
    
    cell_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$cell)) return(NULL)
      
      # Pick out the movie with this ID
      
      paste0("<b>", x$cell, "</b>")
    }
    
    
    if (exists("scDissector_datadir")){
      updateTextInput(session,"inDatapath",,scDissector_datadir)
    }
    
    if (exists("default_model_dataset")){
      
      hideTab(inputId = "inMain", target = "Data")
      ldm=default_model_dataset
      randomly_selected_cells=randomly_select_cells(ldm,params$nrandom_cells_per_sample_choices)
      ldm$dataset$randomly_selected_cells=randomly_selected_cells
      update_all(session,ldm)
      show_all_tabs()
      updateTabsetPanel(session, "inMain", selected = "Model")
    }
    
  }





