library(Matrix)
library(ggvis)
library(dplyr)
library(gplots)
library("heatmaply")
set.seed(3505)
non_data_tabs=c("MetaData","Gating","Basics","Clusters","Cells","QC","Clustering QC","Gene Modules","Samples")

#write.table(file="~/Documents/GitHub/scDissector/viridis_colors.txt",viridis(100),quote=T,row.names=F,col.names=F)

print(getwd())


#load_on_startup()


############################################

as_list_recursive=function(l){
  if (is.null(l)){
    return(list())
  }
  if (length(l)==0){
    return(l)
  }
  for (i in 1:length(l)){
    if (is.list(l[[i]])){
      l[[i]]=as_list_recursive(l[[i]])
    }
    else{
      item=l[[i]]
      names(l)[i]=item
    }
  }
  return(l)
}

as_cluster_sets_recursive=function(l){
  if (is.null(l)){
    return(l)
  }
  l2=list()
    for (i in 1:length(l)){
      if (!is.list(l[[i]])){
        l2[[i]]=names(l)[i]
      }
      else{
        l2[[i]]=as_cluster_sets_recursive(l[[i]])
      }
    }
  names(l2)=names(l)
  return(l2)
}

get_nodes=function(l){
  if (is.list(l)){
    return(setdiff(c(names(l),unlist(sapply(l,get_nodes))),unlist(l)))
  }
  else {
    return(c())
  }
}

remove_node=function(l,node){
  if (is.list(l)){
    if (node%in%names(l)){
      l2=l
      l2[[node]]=NULL
      if (is.list(l[[node]])){
        for (n in names(l[[node]])){
          l2[[as.character(n)]]=l[[node]][[as.character(n)]]
        }
      }
      return(l2)
    }
    else {
      return(sapply(l,remove_node,node,simplify = F))
    }
  }
  else{
    return(l)
  }
}

get_edges=function(l,name=NA){
  if (is.null(l)){
    return()
  }
  else{
    node=c()
    parent=c()
    for (i in 1:length(names(l))){
      if (!is.na(name)&!is.null(names(l))){
        parent[i]=name
        node[i]=names(l)[i]
      }
    }
    m=data.frame(node=node,parent=parent)
    for (i in length(l):1){
      if (!is.null(names(l))){
        m=rbind(get_edges(l[[i]],names(l)[i]),m)
      }
    }
    return(m)
  }
}


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
  session$userData$reactiveVars<-reactiveValues()
  session$userData$reactiveVars$clusters_genes_samples_reactive=list(clusters=c(),genes=c(),samples=c())
  session$userData$reactiveVars$clusters_samples_reactive=list(clusters=c(),samples=c())
  session$userData$sample_colors=default_sample_colors[1:length(get_loaded_samples(session))]
  ncells_choices=as.numeric(setdiff(params$nrandom_cells_per_sample_choices,"All"))
  ncells_per_sample=ncells_choices[which.min(abs((2000/length(get_loaded_samples(session)))-ncells_choices))]

  
  if (is.null(session$userData$clustAnnots)&&(is.null(ldm$cluster_sets))){
    annots=NULL
  }
  else{
    if (is.null(ldm$cluster_sets)){
      annots=session$userData$clustAnnots
    }
    else{
      edges=get_edges(ldm$cluster_sets)
      edges=edges[match(session$userData$default_clusters,edges[,"node"]),]
      annots=edges[,2]
      names(annots)=edges[,1]
    }
  }

  clust_title=paste(session$userData$cluster_order," - ",annots[session$userData$cluster_order],sep="") 
  print(clust_title)
#  clust_title=paste(session$userData$cluster_order," - ",session$userData$clustAnnots[session$userData$cluster_order],sep="") 
  session$userData$cluster_sets=ldm$cluster_sets
  session$userData$scDissector_params$previous_clusters=session$userData$default_clusters
  update_clusters(session,session$userData$default_clusters,T)
  #print("*****")
  #print(as_list_recursive(session$userData$cluster_sets))
  #updateTree(session,"clusters_sets_shinytree",data = as_list_recursive(session$userData$cluster_sets))
  updateSelectInput(session,"inAnnotateClusterNum",choices = clust_title)
  updateSelectInput(session,"inQCClust",choices =clust_title)
  updateTextInput(session,"inSamplesToShow",value = paste(get_loaded_samples(session),collapse=", "))
  updateSelectInput(session,"inClustForDiffGeneExprsProjVsRef",choices = clust_title)
  cluster_sets=unique(unlist(get_nodes(as_cluster_sets_recursive(session$userData$cluster_sets))))
  updateSelectInput(session,"categorizeSamplesBy",choices=colnames(session$userData$sample_annots),selected = "sample_ID")
  updateSelectInput(session,"removeClusterSetSelectInput",choices=cluster_sets)
  updateSelectInput(session,"inReorderSamplesBy",choices=c("All",cluster_sets))
  updateSelectInput(session,"inGatingSample",choices = get_loaded_samples(session))
  updateSelectInput(session,"inGatingShowClusters",choices = clust_title)
  updateSelectInput(session,"input$inTruthNcellsPer",choices=params$nrandom_cells_per_sample_choices,selected =ncells_per_sample )
  updateSelectInput(session,"inQCDownSamplingVersion",choices=get_ds_options(session),selected = ifelse("1000"%in%get_ds_options(session),"1000",max(get_ds_options(session))))
  updateSelectInput(session,"inTruthDownSamplingVersion",choices=get_ds_options(session),selected = ifelse("1000"%in%get_ds_options(session),"1000",max(get_ds_options(session))))
  updateSelectInput(session,"inModulesDownSamplingVersion",choices=get_ds_options(session),selected = max(get_ds_options(session)))
  updateTree(session,"clusters_sets_shinytree",data = as_list_recursive(ldm$cluster_sets))
  if (is.null(get_noise_models(session))){
     updateSelectInput(session,"inModelOrAverage",choices=c("Model parameters","Average"))
  }
  
  
  
  
  
  updateTextInput(session,"inMaxUmis",value = get_max_umis(session))
  updateTextInput(session,"inMinUmis",value = get_min_umis(session))
  message("Successfully finished loading.")
  
  session$userData$loaded_flag<-T
}


load_metadata=function(session,datapath){
  verfile=paste(datapath,"model_versions.csv",sep="/")
  samples_file=paste(datapath,"samples.csv",sep="/")
  read_flag=T
  if (!file.exists(verfile)){
    message(verfile ," not found")
    read_flag=F
  }
  
  if  (!file.exists(samples_file)){
    message(samples_file ," not found")
    read_flag=F
  }
 
  if (read_flag){
    session$userData$vers_tab<-read.csv(file=verfile,header=T,stringsAsFactors = F)
    session$userData$samples_tab<-read.csv(file=samples_file,header=T,stringsAsFactors = F)
  }
  
  myGeneListFile= paste(datapath,"gene_sets.txt",sep="/")
  
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
  
  sample_annots_fn=paste(datapath,"metadata","sample_annots.csv",sep="/")
  if (dir.exists(paste(datapath,"metadata",sep="/"))){
    if (file.exists(sample_annots_fn)){
      session$userData$sample_annots=read.csv(sample_annots_fn,stringsAsFactors = F)
      rownames(session$userData$sample_annots)=session$userData$sample_annots$sample_ID
      session$userData$sample_annots=session$userData$sample_annots
    }
  }
  
  
  mySampleSetsFile= paste(datapath,"sample_sets.txt",sep="/")
  samples_to_show=session$userData$samples_tab$index
  if (file.exists(mySampleSetsFile)){
    sample_sets_tab<-read.table(file=mySampleSetsFile,header=T,stringsAsFactors = F,row.names =1)
    session$userData$sample_sets<-strsplit(sample_sets_tab[,"samples"],",| ,|, ")
    names(session$userData$sample_sets)<-rownames(sample_sets_tab)
    session$userData$samples_to_show=c(names(session$userData$sample_sets),samples_to_show)
  }
  else{
    session$userData$samples_to_show=session$userData$samples_tab$index
  }
  

  
}




tab3_left_margin=12


  function(input, output, session) {
    
    
    if (!exists("use_plotly")){
      use_plotly=T
    }
    session$userData$prev_xy_range_G=c(0,0,0,0)
    session$userData$prev_xy_range_C=c(0,0,0,0)
    session$userData$update_ingene_text_input=T
    session$userData$update_inclusts_text_input=T
    hide_all_tabs()    
    session$userData$loaded_flag<-F
    session$userData$loaded_model_version=NA
    session$userData$prev_inModelColorScale_rel<-c(-4,4)
    session$userData$prev_inModelColorScale_abs<-c(-7,-1)
    session$userData$prev_inModelAbsOrRel="Relative"
    session$userData$click_flag<-T
    
    
    observeEvent(input$inDatapath,{
      if (input$inDatapath==""){
        return()
      }
      load_metadata(session,input$inDatapath)
    session$userData$scDissector_datadir<-input$inDatapath
    
    updateSelectInput(session,"inGeneSets",choices = names(session$userData$geneList))
    updateSelectInput(session,"inModelVer", "Model Version:",choices =  session$userData$vers_tab$title)
    updateSelectInput(session,"inSampleToAdd", "Samples:",choices =session$userData$samples_to_show)
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
      browser()
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

    ldm=scDissector::load_dataset_and_model(model_fn,sample_paths,min_umis = min_umis,max_umis = max_umis,lightweight = T)
    

    show_all_tabs()

    update_all(session,ldm)
   
  })
  
  
  observeEvent(input$inSamplesAddClusterSet,{
    cluster_sets=vis_freqs_cluster_sets_reactive()
    if (input$inReorderSamplesBy=="All"){
      added=names(cluster_sets_reactive())
    }
    else{
      added=input$inReorderSamplesBy
    }
    print(names(added))
    updateTextAreaInput(session,"inSamplesClusterSets",value=paste(unique(c(cluster_sets,added)),collapse = ", "))
  })
  
  observeEvent(input$inReorderSamples,{
    insamples=samples_reactive()
    cluster_sets=cluster_sets_reactive()

    freq_norm=normalize_by_clusterset_frequency(get_cell_to_cluster(session),get_cell_to_sample(session),insamples,cluster_sets,pool_subtype = T,reg = 0)
    
    reorderby=vis_freqs_cluster_sets_reactive()
  
    x=do.call(cbind,freq_norm[reorderby])
    x=pmax(x,0,na.rm=T)
    if (nrow(x)>=2){
      ord=get_order(seriate(as.dist(1-cor(t(x))),method = "OLO_complete"))
      updateTextInput(session,"inSamplesToShow",value = paste(insamples[ord],collapse=", "))
    }
  })
  
  observeEvent(input$inAddSample,{
    s1=paste(unique(c(strsplit(input$inSamples,",|, | ,")[[1]],input$inSampleToAdd)),collapse=", ")
    updateTextInput(session,"inSamples",value=s1) 
  })
  
  observeEvent(input$selectAllSamples,{
    proxy <- DT::dataTableProxy("mytable")
    DT::selectRows(proxy,selected = input$mytable_rows_all)
  })
  
  observeEvent(input$inAbsOrRel,{
    if (input$inAbsOrRel=="Absolute"){
      if (session$userData$prev_inModelAbsOrRel=="Relative"){
        session$userData$prev_inModelColorScale_rel<-input$inModelColorScale
      }
      session$userData$prev_inModelAbsOrRel="Absolute"
      session$userData$modelColorGrad=colgrad_abs
      updateSliderInput(session,"inModelColorScale",label="Log10(expression)",min = -10,max=-.5,step = .5,value = session$userData$prev_inModelColorScale_abs)
    }
    else if (input$inAbsOrRel=="Relative"){
      if (session$userData$prev_inModelAbsOrRel=="Absolute"){
        session$userData$prev_inModelColorScale_abs<-input$inModelColorScale
      }
      session$userData$prev_inModelAbsOrRel="Relative"
      session$userData$modelColorGrad=colgrad_rel
      updateSliderInput(session,"inModelColorScale",label="Log2(expression/mean)",min = -8,max = 8,step = 1,value = session$userData$prev_inModelColorScale_rel)
      
    }
  })
  
  
  observeEvent(input$inSamples,{
    update_clusters_genes_samples_reactive()
  })
  
  observeEvent(input$inClusters,{
    update_clusters_genes_samples_reactive()
  })
  
  observeEvent(input$inGenes,{
    update_clusters_genes_samples_reactive()
  })
  
  
  observe({
    genes=init_genes(input$inGenes)
    genes=init_genes(input$inGenes)
    clusts=strsplit(input$inClusters,",|, | ,")[[1]]
    clusts=intersect(clusts,setdiff(get_all_clusters(session),session$userData$scDissector_params$excluded_clusters))
    
    xy_range <- event_data("plotly_relayout")
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
    session$userData$reactiveVars$clusters_genes_samples_reactive=list(clusters=clusts,genes=genes,samples=samples_reactive())
  })
    
  samples_reactive <-reactive({
    samples=strsplit(input$inSamplesToShow,",|, | ,")[[1]]
    if (!dataset_exists(session)){
      return(c())
    }
    samples=intersect(samples,get_loaded_samples(session))
    return(samples)
  })
  
  
  
  cluster_sets_reactive <-reactive ({
   # print(input$clusters_sets_shinytree)
    if (is.null(input$clusters_sets_shinytree)){
      return(session$userData$cluster_set)
    }
    as_cluster_sets_recursive(input$clusters_sets_shinytree)
  })
  
  
  vis_freqs_cluster_sets_reactive <-reactive ({
    return(cluster_sets=unique(strsplit(input$inSamplesClusterSets,",|, | ,")[[1]]))

  })
  
  
  cluster_annots_reactive<-reactive({
    cluster_sets=cluster_sets_reactive()
    if (is.null(cluster_sets)){
      cluster_sets=session$userData$cluster_sets
    }
    if (is.null(session$userData$clustAnnots)&&(is.null(cluster_sets))){
      return()
    }
    else{
      if (is.null(cluster_sets)){
        return(session$userData$clustAnnots)
      }
      else{
        edges=get_edges(cluster_sets)
        edges=edges[match(session$userData$default_clusters,edges[,"node"]),]
        annots=edges[,2]
        names(annots)=edges[,1]
        return(annots)
      }
    }
  })
  
 
  update_clusters_genes_samples_reactive =function(){
    if (!session$userData$loaded_flag){
      return()
    }
 
    genes=init_genes(input$inGenes)
    clusts=strsplit(input$inClusters,",|, | ,")[[1]]
    clusts=intersect(clusts,setdiff(get_all_clusters(session),session$userData$scDissector_params$excluded_clusters))
    
  
    
    session$userData$reactiveVars$clusters_genes_samples_reactive=list(clusters=clusts,genes=genes,samples=samples_reactive())
    session$userData$reactiveVars$clusters_samples_reactive=list(clusters=clusts,samples=samples_reactive())
  }
  
  sample_annots_reactive<-reactive({
    return(session$userData$sample_annots[intersect(get_loaded_samples(session),rownames(session$userData$sample_annots)),])
  })
  
  
  cells_reactive <-reactive({
    inclusts=session$userData$reactiveVars$clusters_samples_reactive$clusters
    insamples=session$userData$reactiveVars$clusters_samples_reactive$samples
   
    return(sort(select_cells(session,cells=get_all_cells(session),clusters=inclusts,samples=insamples)))
  })
  
  sample_cells_by_cluster_reactive <-reactive({
    if (is.null(session$userData$randomly_selected_cells_by_cluster)){
      session$userData$randomly_selected_cells_by_cluster=list()
      for (i in 1:length(get_ds_options(session))){
        session$userData$randomly_selected_cells_by_cluster[[i]]=list()
      }
    }
    
    ds_i=match(input$inTruthDownSamplingVersion, get_ds_options(session))
    nrandom_cells=input$inTruthNcellsPer
    if (is.na(match(nrandom_cells,session$userData$randomly_selected_cells_by_cluster[[ds_i]]))){
      session$userData$randomly_selected_cells_by_cluster[[ds_i]][[nrandom_cells]]=c()
    }
    for (clusteri in unlist(session$userData$cluster_sets)){
      maski=select_cells(session,cells=get_ds_cells(session,ds_i),clusters=clusteri)
      if (length(maski)==0){
       next 
      }
      #     print(paste(nrandom_cells,sampi))
      if (nrandom_cells=="All"||pmax(0,as.numeric(nrandom_cells),na.rm=T)>=length(maski)){
        session$userData$randomly_selected_cells_by_cluster[[ds_i]][[nrandom_cells]]<-c(session$userData$randomly_selected_cells_by_cluster[[ds_i]][[nrandom_cells]],maski)
      }
      else{
        session$userData$randomly_selected_cells_by_cluster[[ds_i]][[nrandom_cells]]<-c(session$userData$randomly_selected_cells_by_cluster[[ds_i]][[nrandom_cells]],sample(maski,size=as.numeric(nrandom_cells),replace=F))
      }
    }
    return(session$userData$randomly_selected_cells_by_cluster)
  })
  
  
  common_sample_cells_func=function(downsampling_version,ncells_per_sample){
    if (is.null(session$userData$randomly_selected_cells_by_sample)){
      session$userData$randomly_selected_cells_by_sample=list()
      for (i in 1:(length(get_ds_options(session)))+1){
        session$userData$randomly_selected_cells_by_sample[[i]]=list()
      }
    }
    
    ds_i=match(downsampling_version, get_ds_options(session))
    if (is.na(match(ncells_per_sample,session$userData$randomly_selected_cells_by_sample[[ds_i]]))){
      session$userData$randomly_selected_cells_by_sample[[ds_i]][[ncells_per_sample]]=c()
    }
    for (sampi in get_loaded_samples(session)){
      maski=select_cells(session,cells=get_ds_cells(session,ds_i),samples=sampi)
      if (ncells_per_sample=="All"||pmax(0,as.numeric(ncells_per_sample),na.rm=T)>=length(maski)){
        session$userData$randomly_selected_cells_by_sample[[ds_i]][[ncells_per_sample]]<-c(session$userData$randomly_selected_cells_by_sample[[ds_i]][[ncells_per_sample]],maski)
      }
      else{
        session$userData$randomly_selected_cells_by_sample[[ds_i]][[ncells_per_sample]]<-c(session$userData$randomly_selected_cells_by_sample[[ds_i]][[ncells_per_sample]],sample(maski,size=as.numeric(ncells_per_sample),replace=F))
      }
    }
    return(session$userData$randomly_selected_cells_by_sample)
  }
  
  sample_cells_by_sample_reactive <-reactive({
    common_sample_cells_func(input$inTruthDownSamplingVersion,input$inTruthNcellsPer)
  })
  
  sample_cells_by_sample_reactive_QC <-reactive({
    common_sample_cells_func(input$inQCDownSamplingVersion,input$inQCNcellsPerSample)
  })
  

  
  output$event <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover on a point!" else d
  })
  
 
  

  
  modules_reactive <-reactive({
    modules=strsplit(input$inModules,",")[[1]]
    return(modules)
  })
  
  
  ds_cells_QC_reactive <-reactive({
    clust=strsplit(input$inQCClust," - ")[[1]][1]
    samples=samples_reactive()
    
    if (!dataset_exists(session)){
      return()
  }
    sampling_mask=sample_cells_by_sample_reactive_QC()[[match(input$inQCDownSamplingVersion,get_ds_options(session))]][[input$inQCNcellsPerSample]]
    return(select_cells(session,cells = sampling_mask,clusters = clust,samples = samples))#

  })
  
  observeEvent(input$inGeneSets,{
    geneSets=input$inGeneSets
    updateTextInput(session,"inGenes",value=as.character(session$userData$geneList[geneSets]))
    })
  
  
  init_genes=function(genes_string){
    if (length(genes_string)==0){
      return()
    }
    data_genes=get_all_genes(session)
 
    gene_strings_adj=unique(adjust_gene_names(genes_string, data_genes))
 
    session$userData$gcol<-ifelse(toupper(gene_strings_adj)%in%tfs ,4,1)
    names(session$userData$gcol)<-toupper(gene_strings_adj)
    return(gene_strings_adj)
    
  }
  
 
  
#  observeEvent(input$inAnnotateCluster, {
#    session$userData$clustAnnots[strsplit(input$inAnnotateClusterNum," - ")[[1]][1]]<-input$inClustAnnot
#    updateTextInput(session,"inClustAnnot",value="")
#  })
  
#  observeEvent(input$inSaveAnnot, {
#    annot_fn=paste(strsplit(session$userData$loaded_model_file,"\\.")[[1]][1],"_annots.txt",sep="")
#    write.table(file=annot_fn,session$userData$clustAnnots,row.names=T,col.names=F,quote=F,sep="\t")
#    message("Cluster annotations saved to ",annot_fn,".")
#    
#  })

    

  
  observeEvent(input$inClusterGenes, {
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    ingenes=cgs$genes
    mat<-get_models(session,clusters = inclusts,genes = ingenes)
    if (length(ingenes)==0){
      return()
    }
    if (input$inReorderingMethod=="Hierarchical clustering"){
      mat2=mat[!rowSums(is.na(mat))==ncol(mat),]
      cormat=cor(t(mat2),use="comp")
      mask=rowSums(!is.na(cormat))>1
      cormat=cormat[mask,mask]
    
      d1=dist(1-cormat)
      d1[is.na(d1)]=100
      #  reorderv1=hclust(d1)$order
      ord=get_order(seriate(d1,method="GW_complete"))
      
      genes=c(rownames(cormat)[ord],setdiff(rownames(mat),rownames(cormat)))
      update_genes(session,genes,F)
    }
    else if (input$inReorderingMethod=="Diagonal"){
        mat2=mat[!rowSums(is.na(mat))==ncol(mat),]
        ord=order(apply(mat2[,inclusts],1,which.max))
        genes=rownames(mat2)[ord]
        update_genes(session,genes,F)
      }
  })
  
  observeEvent(input$inClusterModules, {
    inclusts=session$userData$reactiveVars$clusters_genes_samples_reactive$clusters
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
      ord=order(apply(modulemat[,inclusts],1,which.max))
      updateTextInput(session,"inModules",value=paste(rownames(modulemat)[ord],collapse=","))
    }
  })

  
  observeEvent(input$inModulesAddGeneList,{
    module_ids=unlist(modules_reactive())
    modsl=gene_to_module_reactive()
    genes_to_show_comma_delimited=paste(unlist(modsl[module_ids]),collapse = ", ")
    new_set_name=paste("Modules_",input$inNUmberOfGeneModules,"_",date(),sep="")
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
    
    }
  )
  
  
  observeEvent(input$inResetModules,{
    updateTextInput(session,"inModules",value=paste(names(gene_to_module_reactive()),collapse=","))
  })
  
  observeEvent(input$inVarMeanScreenSelectedClusters, {
    numis=1000
    ds_i=match(as.character(1000),get_ds_options(session))
    if (is.na(ds_i)){
      ds_i=which.max(as.numeric(get_ds_options(session)))
    }
    if (!session$userData$loaded_flag)
    {
      return(data.frame(m=0,v=0,gene=""))
    }
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    ggs_res=get_genes_to_screen()
    genes=ggs_res$genes
    pref=ggs_res$pref
    inclusts=cgs$clusters
    insamples=cgs$samples
    cell_mask=select_cells(session,clusters = inclusts,samples = insamples)
    ds=get_dstab(session,cells = cell_mask,genes = genes,ds_version = ds_i)
    s1=Matrix::rowSums(ds,na.rm=T) 
    s2=Matrix::rowSums(ds^2,na.rm=T)
    mask1=s1/(numis*ncol(ds))>10^as.numeric(input$inMinExprForScreen[1])&s1/(numis*ncol(ds))<10^as.numeric(input$inMinExprForScreen[2])
    m1=s1[mask1]/(ncol(ds)*numis)
    m2=s2[mask1]/(ncol(ds)*numis)
    v1=m2-m1^2
    x=log10(m1)
    breaks=seq(min(x,na.rm=T),max(x,na.rm=T),.2)
    lv=log2(v1/m1)
    llv=split(lv,cut(x,breaks))
    mask_llv=sapply(llv,length)>0   
    z=sapply(llv[mask_llv],min,na.rm=T)
    b=breaks[-length(breaks)]
    b=b[mask_llv]
    lo=loess(z~b)
    ngenes_to_show=as.numeric(input$inNgenes)
    var_score=lv-predict(lo,newdata =x)
    if (ngenes_to_show>length(inclusts)){
      positive_var_gens=names(var_score)[pmax(var_score,0,na.rm=T)>0.001]
      umisum=as.matrix(aggregate.Matrix(t(ds[positive_var_gens,]),get_cell_to_cluster(session,cells=colnames(ds))))
      cluster_avg=(t(umisum/rowSums(umisum)))
      cluster_avg_normed=log2((1e-6+cluster_avg)/(1e-6+rowMeans(cluster_avg)))
      high_var_genes=as.vector(unlist(sapply(split(rownames(cluster_avg_normed),apply(cluster_avg_normed,1,which.max)),FUN=function(x){names(head(sort(var_score[x],decreasing = T),floor(ngenes_to_show/ncol(cluster_avg_normed))))})))
    }
    else{
      high_var_genes=names(head(sort(var_score,decreasing=T),ngenes_to_show))  
    }
#   
    genes_to_show_comma_delimited=paste(high_var_genes,collapse = ", ")
    new_set_name=paste("Varmean_",pref,"_",ngenes_to_show,"_",date(),sep="")
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
  
  
  function(cluster_avg,var_score){
    clusts=sample(colnames(cluster_avg),size=ncol(cluster_avg),replace = F)
    for (clust in clusts){
      
    }
  }
  
  observeEvent(input$inBlindChisqSelectedClusters, {
    message("screening for variable ",input$inSelectGenesFrom)
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    clusters=cgs$clusters
    samples=cgs$samples
    ggs_res=get_genes_to_screen()
    genes=ggs_res$genes
    pref=ggs_res$pref
    cells=select_cells(session,clusters = clusters)
    genes=intersect(genes,get_all_genes(session))
  
    counts_to_chisq=get_counts_array(session,samples = samples,genes = genes,clusters = clusters)
    noise_counts=get_noise_counts_array(session,samples = samples,genes = genes,clusters = clusters)
    if (!is.null(noise_counts)){
      counts_to_chisq=pmax(counts_to_chisq-noise_counts,0)
    }
    
    chisq_res2=chisq_genes(counts_to_chisq)
    if (nrow(chisq_res2)==0){
      message("Chi sq detected nothing..")
      return()
    }
    counts=apply(get_counts_array(session,samples = samples,genes = genes,clusters = clusters),2:3,sum)
    avg=t(t(counts)/colSums(counts))

    mask=rownames(chisq_res2)%in%genes&chisq_res2[,3]<0.01&
      rowSums(get_umitab(session,cells = cells,genes =rownames(chisq_res2))>1)>=10
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
  
  observeEvent(input$inAddClusterSetButton, {
    if (is.null(input$clusters_sets_shinytree)){
      l=session$userData$cluster_sets
    }
    else{
      l=cluster_sets_reactive()
    }
    l[[input$inAddClusterSet]]=list()
    updateTree(session,"clusters_sets_shinytree",data = as_list_recursive(l))
    updateSelectInput(session,"removeClusterSetSelectInput",choices=unique(get_nodes(l)))
    updateSelectInput(session,"inReorderSamplesBy",choices=c("All",unique(get_nodes(l))))
  })
  
  observeEvent(input$inRemoveClusterSetButton, {
    if (is.null(input$clusters_sets_shinytree)){
      l=session$userData$cluster_sets
    }
    else {
      l=cluster_sets_reactive()
    }
    l2=remove_node(l,input$removeClusterSetSelectInput)
    updateTree(session,"clusters_sets_shinytree",data = as_list_recursive(l2))
    updateSelectInput(session,"removeClusterSetSelectInput",choices=unique(unlist(get_nodes(l2))))
    updateSelectInput(session,"inReorderSamplesBy",choices=c("All",unique(unlist(get_nodes(l2)))))
  })
  
  observeEvent(input$saveClusterSetButtion,{
   
    m=get_edges(cluster_sets_reactive())
    #Bug fix AL 08/09/19
    cluster_set_fn=paste(substr(session$userData$loaded_model_file,1,nchar(session$userData$loaded_model_file)-3),"_cluster_sets.txt",sep="")
    #cluster_set_fn=paste(strsplit(session$userData$loaded_model_file,"\\.")[[1]][1],"_cluster_sets.txt",sep="")
    write.table(file=cluster_set_fn,m,row.names=F,col.names=T,quote=F,sep="\t")
    message("Cluster-sets were saved to ",cluster_set_fn,".")
  })
  
  observeEvent(input$reloadClusterSetButtion,{
    #Bug fix AL 08/27/19
    cluster_sets_fn=paste(substr(session$userData$loaded_model_file,1,nchar(session$userData$loaded_model_file)-3),"_cluster_sets.txt",sep="")
    #cluster_sets_fn=paste(strsplit(session$userData$loaded_model_file,"\\.")[[1]][1],"_cluster_sets.txt",sep="")
    if (file.exists(cluster_sets_fn)){
      a=read.delim(cluster_sets_fn,header=T,stringsAsFactors = F)
      l=get_cluster_set_tree(a)
      
      updateTree(session,"clusters_sets_shinytree",data = as_list_recursive(l))
      updateSelectInput(session,"removeClusterSetSelectInput",choices=unique(get_nodes(l)))
      updateSelectInput(session,"inReorderSamplesBy",choices=c("All",unique(get_nodes(l))))
    }
    
  })
  
  observeEvent(input$selectSamples,{
    #rows_all is ordered as the table
    
    samples=intersect(input$mytable_rows_all,input$mytable_rows_selected)
    
    samples=as.character(sample_annots_reactive()[samples,"sample_ID"])

      updateTextAreaInput(session,inputId = "inSamplesToShow",value = paste(samples,collapse=", "))
    
  })
  
  
  
  chisq_genes=function(counts){
    counts=apply(counts,2:3,sum)
    gene_mask=apply(counts,1,max)>3
    counts=counts[gene_mask,,drop=F]
    cluster_tot=colSums(counts)
    arrcont=array(c(counts,matrix(cluster_tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts),dim=c(dim(counts),2))
    suppressWarnings({ res=t(apply(arrcont,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}))})
    rownames(res)=rownames(counts)
    res=res[!is.na(res[,1]),]
    res=cbind(res,adjp=p.adjust(res[,1],method="BH"))
    return(res)
  }   
  
  get_genes_to_screen=function(){
   
    if (input$inSelectGenesFrom=="All genes"){
      pref="All"
      genes=get_all_genes(session)
    } else if (input$inSelectGenesFrom=="Without RPs"){
      pref="NoRps"
      genes=get_all_genes(session)
      genes=setdiff(genes,grep("^RP|^Rp",genes,val=T))
    }else if (input$inSelectGenesFrom=="TFs"){
      pref="Tfs"
      genes=tfs
    } else if(input$inSelectGenesFrom=="Surface markers"){
      pref="Surf"
      genes=surface_markers
    }
    return(list(genes=genes,pref=pref))
  }
  
  observeEvent(input$inFC, {
    ggs_res=get_genes_to_screen()
    genes=ggs_res$genes
    samples=ggs_res$samples
    pref=ggs_res$pref
    mask=intersect(get_all_genes(session),genes)
    
    clusts_fg=strsplit(input$inFC_fgClusts,",")[[1]]
    clusts_bg=strsplit(input$inFC_bgClusts,",")[[1]]
    
    mask_fg=select_cells(session,samples = samples,clusters=clusts_fg)
    mask_bg=select_cells(session,samples = samples,clusters=clusts_bg)
    reg=20
    sfg=reg+rowSums(get_umitab(session,cells = mask_fg,genes = mask))
    sbg=reg+rowSums(get_umitab(session,cells = mask_bg,genes = mask))
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
     clust_title=paste(session$userData$default_clusters," - ",cluster_annots_reactive()[session$userData$default_clusters],sep="") 
    update_clusters(session,session$userData$default_clusters)
    updateSelectInput(session,"inQCClust",choices=clust_title)
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title)
    
  })
  
  observeEvent(input$inSaveClusterOrder, {
    clusters=session$userData$reactiveVars$clusters_genes_samples_reactive$clusters
    order_fn=paste(strsplit(session$userData$loaded_model_file,"\\.")[[1]][1],"_order.txt",sep="")
    if (all(get_all_clusters(session)%in%clusters)){
      session$userData$default_clusters<-clusters
      write.table(clusters,file=order_fn,quote = F,sep="\t",row.names = F,col.names = F)
      message("Order has been successfully saved.")
    }else{
      message("Order could not be saved since 1 or more clusters are missing.")
    }
  })
  
  

 # output$downloadExprsTable <- downloadHandler(
#    filename = function() { paste("Exprs", '.csv', sep='') },
#    content = function(file) {
#      write.csv(session$userData$model$models, file)
#    }
#  )
#  outputOptions(output, 'downloadExprsTable', suspendWhenHidden=FALSE)   
  
  
  
  observeEvent(input$inOrderClusters, {
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    ingenes=cgs$genes
    
    if (length(ingenes)==0){
      return()
    }
    mat<-get_models(session,clusters = inclusts,genes=ingenes)
    if (input$inReorderingClustersMethod=="Hierarchical clustering"){
    cormat=cor(mat,use="comp")
    d1=dist(1-cormat)
    d1[is.na(d1)]=100
  #  reorderv1=hclust(d1)$order
    reorderv1=get_order(seriate(d1,method="GW_complete"))
    }
    else if (input$inReorderingClustersMethod=="Diagonal"){
      reorderv1=order(apply(mat[ingenes,],2,which.max))
    }
    else if (input$inReorderingClustersMethod=="Cluster-sets"){
      cluster_sets=cluster_sets_reactive()
      if (is.null(cluster_sets)){
        reorderv1=1:length(inclusts)
      }
      else{
        x=function(l){if (!is.list(l)){return(l)}; for (i in 1:length(l)){return(sapply(l,x))}}
        y=function(l){if (!is.list(l)){return(intersect(inclusts,l))}else{sapply(l,y)}}
        reorderv1=order(match(inclusts,unlist(y(x(cluster_sets)))))
      }
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
      score=rep(0,length(inclusts))
      for (i in 1:length(genes)){
        message(sig[i],genes[i])
  
        score=score+(10^(5*(length(genes)-i)))*ifelse(sig[i]=="+",1,-1)*(mat[genes[i],inclusts])
      }
      
      reorderv1=order(score)
    }
    clust_title=paste(inclusts," - ",cluster_annots_reactive()[inclusts],sep="") 
    
    update_clusters(session,inclusts[reorderv1])
    updateSelectInput(session,"inQCClust",choices=clust_title[reorderv1])
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title[reorderv1])
  })
  
  

  ###########################################################
  
   
  
  
  
  modules_varmean_reactive=reactive({

    ds_i=match(input$inModulesDownSamplingVersion,get_ds_options(session))
    
    if (!session$userData$loaded_flag)
    {
      return(data.frame(m=0,v=0,gene=""))
    }
    
    cgs=session$userData$reactiveVars$clusters_samples_reactive
    inclusts=cgs$clusters
    insamples=cgs$samples
    cell_mask=select_cells(session,cells=get_all_cells(session),clusters = inclusts,samples=insamples)
    ds=get_dstab(session,cells = cell_mask,ds_version = ds_i)
 
    ds_mean<-Matrix::rowMeans(ds)
    genemask=ds_mean>10^input$inVarMeanXlim[1]&ds_mean<10^input$inVarMeanXlim[2]
    ds=ds[genemask,]
    ds_mean=ds_mean[genemask]
    message("Estimating variance for ",nrow(ds)," genes")
    s1=Matrix::rowSums(ds,na.rm=T) 
    s2=Matrix::rowSums(ds^2,na.rm=T)
    ds_mean=s1/ncol(ds)
    m2=s2/ncol(ds)
    ds_var=m2-ds_mean^2
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
  
  sample_to_category<-reactive({
    samps=samples_reactive()
    if (input$categorizeSamplesBy==""){
      v=samps
      names(v)=as.character(samps)
    }
    else{
      v=session$userData$sample_annots[samps,input$categorizeSamplesBy]
      names(v)=as.character(samps)
    }
    return(v)
  })
  
  
  getGeneModuleMask=function(){
    if (!session$userData$loaded_flag)
    {
      return()
    }
    df=modules_varmean_reactive()
    lo=reactiveLoess()
    lline=predict(lo,newdata =log10(df$m))
    isolate({
      geneModuleMask<-log10(df$m)>as.numeric(input$inVarMean_MeanThresh)&log2(df$v/df$m)>lline+as.numeric(input$inVarMean_varmeanThresh)
    })
    geneModuleMask[is.na(geneModuleMask)]<-F
    names(geneModuleMask)=rownames(df)
    return(geneModuleMask)
    
  }
  
 
  
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
    
    n=sum(getGeneModuleMask())
    legend("topright", paste(n,"genes"), bty="n",text.col=2) 
    
  })
  
  
  module_cor_reactive=reactive({
    if (!session$userData$loaded_flag)
    {
      return()
    }
    input$inModulesGetCormap
     
    isolate({
      cell_mask=cells_reactive()
      ds_i=match(input$inModulesDownSamplingVersion,get_ds_options(session))
    })
      ds= get_dstab(session,cells = cell_mask,ds_version = ds_i)
    
    message("calculating gene-to-gene correlations..")
    cormat=get_avg_gene_to_gene_cor(ds[names(which(getGeneModuleMask())),],get_cell_to_sample(session,cells=colnames(ds)))
    return(cormat)
  })
  
  
  gene_to_module_reactive=reactive({
    if (!session$userData$loaded_flag)
    {
      return()
    }
    if (input$inModulesGetCormap==0){
      return()
    }
    if (is.null(getGeneModuleMask())){
      return()
    }
    cormat=module_cor_reactive()
    if (is.null(cormat)){
      return()
    }
    c2c=cutree(hclust(as.dist(1-cormat)),k=as.numeric(input$inNUmberOfGeneModules))
    modsl=split(names(c2c),c2c)
    updateTextInput(session,"inModules",value=paste(names(modsl),collapse=","))
    updateSelectInput(session,"inModuleSelect",label="Show Module:",choices=names(modsl))
    print(modsl)
    return(modsl)
    
  })
 
 
  output$textModuleGenes<- renderText({
    modsl=gene_to_module_reactive()
    if (!is.null(modsl)){
      paste(modsl[[input$inModuleSelect]],collapse=",")
    }
  })
  
  output$inDownloadModuleList <- downloadHandler(
    filename = function() {
      paste("modules-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      
      isolate({modsl=gene_to_module_reactive()})
      tab=cbind(module_name=names(modsl),genes=sapply(modsl,paste,collapse=","))
      write.table(tab, file,row.names = F,col.names = T,quote =F)
    }
  )
  
  ###########################################################
  
  output$mytable = DT::renderDataTable({
    tab=sample_annots_reactive()
   # v=sprintf(
  #    '<input type="checkbox" name="%s" value="%s"/>',
  #    1:nrow(tab), T)
    if (is.null(tab)){
      removeUI("#selectAllSamples",immediate = T)
      removeUI("#selectSamples",immediate = T)
      removeUI("#mytable",immediate = T)
      return()
    }
    DT::datatable(tab,filter = "top",options = list(pageLength = 50),rownames = F,escape=F,selection=list(mode = 'multiple', selected =1:nrow(tab)))
    
  })
  
  
 
  
  
  output$ncells_barplot <-renderPlot({
    input$inModelVer
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    clusters=cgs$clusters
    samples=cgs$samples
    cells=select_cells(session,clusters = clusters,samples = samples)
    ncells=table(get_cell_to_cluster(session,cells =cells))[clusters]
    par(mar=c(1,7,2,2))
    ylim=c(1,10^ceiling(log10(max(ncells,na.rm=T))))
    barplot(ncells,xaxs = "i",log="y",ylab="#cells",names.arg ="",ylim=ylim)
  })
  
  output$UMI_boxplot <- renderPlot({
    input$inModelVer
    par(xaxs="i")
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    clusters=cgs$clusters
    samples=cgs$samples
    cells=select_cells(session,clusters = clusters,samples = samples)
    par(mar=c(2,7,1,2))
    u=get_umitab(session,cells = cells)
    boxplot(split(log10(colSums(u)),get_cell_to_cluster(session,cells=cells))[clusters],las=2,ylab="log10(#UMIs)")
  })
  
 
  output$gaiting_plots_dynamic <- renderUI({
    insilico_gating_scores=get_insilico_gating_scores(session)
    if (!is.null(insilico_gating_scores)){
      nplots=length(insilico_gating_scores)
      he=max(c(500,500*nplots,na.rm=T))
      plotOutput("gating_plots", width = "100%", height = he)
    }
  })
  
  output$subtype_freqs <- renderUI({
    cluster_sets=vis_freqs_cluster_sets_reactive()
    cluster_sets_length=sapply(cluster_sets_reactive(),length)[cluster_sets]
    cluster_sets=cluster_sets[pmax(cluster_sets_length,0,na.rm=T)>1]
    if(is.null(cluster_sets)){
      return()
    }

    if (length(cluster_sets)>0){
      he=300*length(cluster_sets)
      plotOutput("subtype_freqs_barplots", width = "100%", height = he)
    }
  })
  
  output$gating_plots <- renderPlot({
    clustering_params=get_clustering_params(session)
    insilico_gating_scores=get_insilico_gating_scores(session)
    nplots=length(insilico_gating_scores)
    samp=input$inGatingSample
    numis=get_numis_before_filtering(session,samp)
    clust=strsplit(input$inGatingShowClusters," ")[[1]][1]
    cell_to_cluster=get_cell_to_cluster(session,cells=select_cells(session,samples=samp))
    layout(matrix(1:(nplots),nplots,1))
    for (i in 1:nplots){
      mask=intersect(names(numis),names(insilico_gating_scores[[i]]))
      plot(numis[mask],insilico_gating_scores[[i]][mask],log="x",ylab=paste("fraction",names(clustering_params$insilico_gating)[i]),xlab="#UMIs",col=ifelse(is.na(cell_to_cluster[mask]),"gray",ifelse(cell_to_cluster[mask]==clust,2,1)))
      points(numis[mask],insilico_gating_scores[[i]][mask],pch=ifelse(cell_to_cluster[mask]==clust,20,NA),col=2)
      rect(xleft = get_min_umis(session),clustering_params$insilico_gating[[i]]$interval[1],get_max_umis(session),clustering_params$insilico_gating[[i]]$interval[2],lty=3,lwd=3,border=2)
    }

   
  })
  
  
  output$samples_enrichment <- renderUI({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    insamples=cgs$samples
    he=max(c(500,12*length(inclusts)),na.rm=T)
    if (length(insamples)>1&(length(inclusts)>1)){
      plotOutput("samples_enrichment_plot",height=he)
    }
    
  })
  
  output$cluster_sizes<-renderUI({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    he=max(c(500,12*length(inclusts)),na.rm=T)
    plotOutput("cluster_sizes_plot",height=he)
  })
  
  
  output$avg_profile_plot <- renderUI({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    insamples=cgs$samples
    he=max(c(500,12*length(inclusts)),na.rm=T)
    
    if (use_plotly){
     plotlyOutput("avg_heatmap_interactive", width = "100%", height = he)
    }
    else{
      plotOutput("avg_heatmap", width = "100%", height = he)
    }
  })
  
  
  output$avg_module_plot <- renderUI({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    he=max(c(500,12*length(inclusts)),na.rm=T)
    plotOutput("avg_heatmap_modules", width = "100%", height = he)
  })
  
  output$cor_gene_module_plot <- renderUI({
 
    if (input$inModulesGetCormap==0){
      return()
    }
    he=6*nrow(module_cor_reactive())
    plotOutput("genecor_heatmap_modules", width = "100%", height = he)
  })
  
  output$cor_module_plot <- renderUI({
    
    if (input$inModulesGetCormap==0){
      return()
    }
    he=14*length(gene_to_module_reactive())
    plotOutput("modulecor_heatmap_modules", width = "100%", height = he)
  })
  
  output$genecor_heatmap_modules <- renderPlot({
    cormat=module_cor_reactive()
    if (is.null(cormat)){
      return()
    }
    zbreaks=c(-1,seq(-1,1,l=99),1)
    cor_cols=colorRampPalette(c("blue","white","red"))(100)
    #  ord=hclust(as.dist(1-cormat))$order
    par(mar=c(3,3,3,3))
    ord=get_order(seriate(as.dist(1-cormat),method="OLO_complete"))
    image(cormat[ord,ord],col=cor_cols,breaks=zbreaks,axes=F)
    mtext(text = colnames(cormat)[ord],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    box()
  })

  output$modulecor_heatmap_modules <- renderPlot({
    modsl=gene_to_module_reactive()
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    if (is.null(modsl)){
     return()
    }
    ds_i=match(input$inModulesDownSamplingVersion,get_ds_options(session))

    mods=rep(1:length(modsl),sapply(modsl,length))
    names(mods)=unlist(modsl)
    modsums=aggregate(get_dstab(session,ds_version = ds_i,genes = names(mods),cells=select_cells(session = session,clusters = cgs$clusters,samples = cgs$samples)),groupings=mods,fun="sum")
    cormat=cor(t(as.matrix(modsums)))
    zbreaks=c(-1,seq(-1,1,l=99),1)
    cor_cols=colorRampPalette(c("blue","white","red"))(100)
    par(mar=c(3,3,3,3))
    ord=get_order(seriate(as.dist(1-cormat),method="OLO_complete"))
    image(cormat[ord,ord],col=cor_cols,breaks=zbreaks,axes=F)
    mtext(text = colnames(cormat)[ord],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    mtext(text = colnames(cormat)[ord],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
    box()
  
  })
  
  output$sample_avg_profile_plot <- renderUI({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    insamples=cgs$samples
    he=max(c(300,12*length(insamples)),na.rm=T)
    plotOutput("avg_heatmap_samples", width = "100%", height = he)
  })
  
  
 
  
  output$samples_enrichment_plot <- renderPlot({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
#    sample_cols=sample_colors_reactive()
    inclusts=cgs$clusters
    insamples=cgs$samples
    samp_to_cat=sample_to_category()[insamples]
    cats=as.character(unique(samp_to_cat))
    cat_cols= session$userData$sample_colors[1:length(cats)]
    names(cat_cols)=cats
    sample_cols=cat_cols[as.character(samp_to_cat)]
    par(mar=c(7,0,1.6,0))
    tab=matrix(0,length(get_all_clusters(session)),length(get_loaded_samples(session)))
    rownames(tab)=get_all_clusters(session)
    colnames(tab)=get_loaded_samples(session)
    cells=select_cells(session,samples = insamples)
    tmptab=table(get_cell_to_cluster(session,cells=cells),get_cell_to_sample(session,cells = cells))
    tab[rownames(tmptab),colnames(tmptab)]=tmptab
    tab=tab[,insamples,drop=F]
    tab=t(t(tab)/colSums(tab))
    tab=(tab/rowSums(tab))[inclusts,,drop=F]
    
 #   barplot(t(tab[nrow(tab):1,]),col=sample_cols[match(insamples,get_loaded_samples(session))],horiz =T,yaxs = "i",names.arg=NULL,main="Samples",axes=F)
    barplot(t(tab[nrow(tab):1,,drop=F]),col=sample_cols,horiz =T,yaxs = "i",names.arg=NULL,main="Samples",axes=F)
    
  })
  
  output$cluster_sizes_plot <- renderPlot({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    insamples=cgs$samples
    
    tab=session$userData$ncells_per_cluster*0
    tmp_tab=table(get_cell_to_cluster(session,select_cells(session,samples=insamples)))
    tab[names(tmp_tab)]=tmp_tab
    
    precentage=100*tab[inclusts]/sum(tab[setdiff(names(tab),session$userData$scDissector_params$excluded_clusters)])
    
    par(mar=c(6.9,.3,1.35,0))
    barplot(rev(precentage),horiz=T,border=F,col="gray",axes=F,yaxs="i",xlab="Cells ( % )",cex.axis = .8)
    axis(1)
  })
  


  
  output$avg_heatmap_interactive <-renderPlotly({

      ##### Don't delete!!
      input$inAnnotateCluster
      input$inModelVer
      input$inBirdClusterSetExclude
      input$inBirdClusterSetInclude
      #####
      cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
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
      if (input$inModelOrAverage=="Model parameters"){
        mat<-get_models(session,genes=ingenes,clusters = inclusts)
      }
      else {
        ncells_per_sample=table(get_cell_to_sample(session))[insamples]
        counts=get_counts_array(session,samples = insamples,genes = ingenes,clusters = inclusts)
        weights=ncells_per_sample/sum(ncells_per_sample)
        arr_weights=array(weights,dim=c(length(insamples),length(ingenes),length(inclusts)))
        numis_per_clust_mat=apply(counts[insamples,,inclusts,drop=F],c(1,3),sum,na.rm=T)
        numis_per_clust_arr=aperm(array(numis_per_clust_mat,dim=c(length(insamples),length(inclusts),length(ingenes))),c(1,3,2))
        
        if (input$inModelOrAverage=="Batch-corrected Average"){
          noise_counts=get_noise_counts_array(session,samples = insamples,genes = ingenes,clusters = inclusts)
          if (!is.null(noise_counts)){
            
            counts=counts-noise_counts
            counts<-pmax(counts,0)
            
            noise_numis_per_clust_mat=apply(noise_counts,c(1,3),sum,na.rm=T)
            noise_numis_per_clust_arr=aperm(array(noise_numis_per_clust_mat,dim=c(length(insamples),length(inclusts),length(ingenes))),c(1,3,2))
            numis_per_clust_arr=pmax(numis_per_clust_arr-noise_numis_per_clust_arr,0)
            
          }
        }
 
        mat=apply(counts,2:3,sum,na.rm=T)/apply(numis_per_clust_arr,2:3,sum,na.rm=T)
        rownames(mat)=ingenes
      }
      if (min(dim(mat))==0){
        return()
      }

      isolate({
        mat1=mat
        abs_or_rel=input$inAbsOrRel
      }) 
      
      main_title=paste(input$inModelOrAverage,":",session$userData$loaded_model_version)
      genes=rownames(mat1)
      gene.cols=session$userData$gcol[toupper(rownames(mat1))]
      clusters=colnames(mat1)
      clusters_text=paste(" (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster[setdiff(names(session$userData$ncells_per_cluster),session$userData$scDissector_params$excluded_clusters)]),digits=1),"% )",sep="")
      annots=cluster_annots_reactive()[inclusts]
   
      return(plot_avg_heatmap_interactive(mat1,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute=abs_or_rel,colgrad=session$userData$modelColorGrad))
      # plot_avg_heatmap(mat1,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute=abs_or_rel,reg=1e-6,colgrad =session$userData$modelColorGrad)
  })
  
  
  output$avg_heatmap <-renderPlot({
    
    ##### Don't delete!!
    input$inAnnotateCluster
    input$inModelVer
    input$inBirdClusterSetExclude
    input$inBirdClusterSetInclude
    #####
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
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
    if (input$inModelOrAverage=="Model parameters"){
      mat<-get_models(session,genes=ingenes,clusters = inclusts)
    }
    else {
      ncells_per_sample=table(get_cell_to_sample(session))[insamples]
      counts=get_counts_array(session,samples = insamples,genes = ingenes,clusters = inclusts)
      weights=ncells_per_sample/sum(ncells_per_sample)
      arr_weights=array(weights,dim=c(length(insamples),length(ingenes),length(inclusts)))
      numis_per_clust_mat=apply(counts[insamples,,inclusts,drop=F],c(1,3),sum,na.rm=T)
      numis_per_clust_arr=aperm(array(numis_per_clust_mat,dim=c(length(insamples),length(inclusts),length(ingenes))),c(1,3,2))
      
      if (input$inModelOrAverage=="Batch-corrected Average"){
        noise_counts=get_noise_counts_array(session,samples = insamples,genes = ingenes,clusters = inclusts)
        if (!is.null(noise_counts)){
          
          counts=counts-noise_counts
          counts<-pmax(counts,0)
          
          noise_numis_per_clust_mat=apply(noise_counts,c(1,3),sum,na.rm=T)
          noise_numis_per_clust_arr=aperm(array(noise_numis_per_clust_mat,dim=c(length(insamples),length(inclusts),length(ingenes))),c(1,3,2))
          numis_per_clust_arr=pmax(numis_per_clust_arr-noise_numis_per_clust_arr,0)
          
        }
      }
      mat=apply(counts,2:3,sum,na.rm=T)/apply(numis_per_clust_arr,2:3,sum,na.rm=T)
      rownames(mat)=ingenes
    }
    if (min(dim(mat))==0){
      return()
    }
    
    isolate({
      mat1=mat
      abs_or_rel=input$inAbsOrRel
    }) 
    
    main_title=paste(input$inModelOrAverage,":",session$userData$loaded_model_version)
    genes=rownames(mat1)
    gene.cols=session$userData$gcol[toupper(rownames(mat1))]
    clusters=colnames(mat1)
    clusters_text=paste(" (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster[setdiff(names(session$userData$ncells_per_cluster),session$userData$scDissector_params$excluded_clusters)]),digits=1),"% )",sep="")
    annots=cluster_annots_reactive()[inclusts]
    return(plot_avg_heatmap(mat1,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute=abs_or_rel,reg=1e-6,colgrad =session$userData$modelColorGrad))
  })
  
  
  module_counts_reactive<-reactive({
    modsl=gene_to_module_reactive()
    if (is.null(modsl)){
      return()
    }
    samps=samples_reactive()
    if (length(samps)==0){
      return()
    }
    
    counts=apply(get_counts_array(session,samples=samps),2:3,sum)
    totnumis=colSums(counts)
    return(t(sapply(modsl,function(modi){colSums(counts[modi,,drop=F])})/totnumis))
  
    
  })

  output$avg_heatmap_modules <- renderPlot({
    ##### Don't delete!!
    input$inModelVer
    
    #####
    
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    insamples=cgs$samples
    inmodules=modules_reactive()
    modsl=gene_to_module_reactive()
    if (!session$userData$loaded_flag){
      return()
    }
    if (is.null(modsl)){
      return()
    }
    if (length(inmodules)==0){
      return()
    }
    modulemat<-module_counts_reactive()
    if(is.null(modulemat)){
      return()
    }
    modulemat<-modulemat[inmodules,inclusts]
    
    zlim=input$inAvgModuleColorScale
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
  
    image(mat_to_show[,ncol(mat1):1],col=colgrad_rel,breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main=session$userData$loaded_model_version)
  
    box()
    
    mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1)
    mtext(text =paste(" ",colnames(mat1)," (n=",session$userData$ncells_per_cluster[inclusts]," ; ",round(100*session$userData$ncells_per_cluster[inclusts]/sum(session$userData$ncells_per_cluster),digits=1),"% )",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =paste(cluster_annots_reactive()[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
  })
  
  output$avg_heatmap_samples <- renderPlot({
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
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
    
    zlim=input$inSamplesColorScale
    
    par(mar=c(7,7,1,9))
    mat1=t(apply(get_counts_array(session,samples=insamples,genes = ingenes,clusters = inclusts),1:2,sum))
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
      if (length(inclusts)==length(get_all_clusters(session))){
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
    image(matrix(1:100,100,1),breaks=0:100,col=session$userData$modelColorGrad,axes=F)
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
    plot_truth_heatmap_wrapper()
  })
  
  output$downloadTruthPlot <- downloadHandler(
    filename = function() { paste(session$userData$loaded_model_version, '.png', sep='') },
    content = function(file) {
      png(file,1200,900)
      plot_truth_heatmap_wrapper()
      dev.off()
    }
  ,contentType = "image/png")
  
  
 
  
  plot_truth_heatmap_wrapper=function(){
    zlim=input$inTruthColorScale
    #sample_cols=sample_colors_reactive()
    cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    inclusts=cgs$clusters
    ingenes=cgs$genes
    insamples=cgs$samples
    score_genes=init_genes(input$inGenesScore)
    reverse_score_order=input$inScoreReverseOrder
    global_order=input$inGlobalOrder
    samp_to_cat=sample_to_category()[insamples]
    cats=as.character(unique(samp_to_cat))
    cat_cols= session$userData$sample_colors[1:length(cats)]
    names(cat_cols)=cats
    sample_cols=cat_cols[as.character(samp_to_cat)]
    nclust=length(inclusts)
    if (!session$userData$loaded_flag){
      return()
    }
    if (length(ingenes)==0){
        
      return()
    }

    samps_cl=list() 
    ds=get_dstab(session,ds_version = match(input$inTruthDownSamplingVersion,get_ds_options(session)))
    cells_per_sample=as.numeric(input$inTruthNcellsPer)
    genes=intersect(ingenes,rownames(ds))
    if (input$inTruthEqualRepresent=="equally_represent_samples"){
      cells_selected=sample_cells_by_sample_reactive()[[match(input$inTruthDownSamplingVersion,get_ds_options(session))]][[input$inTruthNcellsPer]]
    }
    else {
      cells_selected=sample_cells_by_cluster_reactive()[[match(input$inTruthDownSamplingVersion,get_ds_options(session))]][[input$inTruthNcellsPer]]
      
    }
    cell_mask=select_cells(session,cells = cells_selected,clusters = inclusts,samples = insamples)
    ds=ds[,cell_mask]
    plot_truth_heatmap(ds,get_cell_to_sample(session,cells=colnames(ds)),get_cell_to_cluster(session,cells=colnames(ds)),insamples,genes,inclusts,zlim,sample_cols=sample_cols,showSeparatorBars=input$inTruthShowSeparatorBars,score_genes=score_genes,reverse_score_order=reverse_score_order,global_order=global_order)

  }
  
  outputOptions(output, "downloadTruthPlot", suspendWhenHidden=FALSE)
 

    varmean_reactive <- reactive({
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      samples=samples_reactive()
      if (!dataset_exists(session)){
        return(data.frame(m=0,varmean=0,gene=""))
      }
    
      ds=get_dstab(session,cells = select_cells(session,cluster=clust,samples=samples),ds_version = match(input$inQCDownSamplingVersion,get_ds_options(session)))
      
      if (is.null(ds))(return(data.frame(m=0,varmean=0,gene=0)))
      if (ncol(ds)<2)(return(data.frame(m=0,varmean=0,gene=0)))
      
      m=Matrix::rowMeans(ds)
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
      samp_to_cat=sample_to_category()[insamples]
      cats=unique(samp_to_cat)
      cat_cols= session$userData$sample_colors[1:length(cats)]
      names(cat_cols)=cats
       par(mar=c(0,0,0,0))
      plot.new()
         leg=names(cat_cols)
       if(length(leg)==0){
          leg=rep("",length(cat_cols))
        }
        ncol=ceiling(length(cat_cols)/10)
        leg[is.na(leg)]=cat_cols[is.na(leg)]
        leg[leg==""]=cat_cols[leg==""]
        legend("topleft",pch=15,col=cat_cols,legend=leg,cex=1,xpd=T,ncol=ncol)
    })
    
    
    output$truthSampleLegend <- renderPlot({
      insamples=samples_reactive()
      samp_to_cat=sample_to_category()[insamples]
      cats=unique(samp_to_cat)
      cat_cols= session$userData$sample_colors[1:length(cats)]
      names(cat_cols)=cats
      par(mar=c(0,0,0,0))
      plot.new()
      leg=names(cat_cols)
      if(length(leg)==0){
        leg=rep("",length(cat_cols))
      }
      ncol=ceiling(length(cat_cols)/10)
      leg[is.na(leg)]=cat_cols[is.na(leg)]
      leg[leg==""]=cat_cols[leg==""]
      legend("topleft",pch=15,col=cat_cols,legend=leg,cex=1,xpd=T,ncol=ncol)
    })
    
    
    output$subtype_freqs_barplots <- renderPlot({
      insamples=samples_reactive()

      cluster_sets=cluster_sets_reactive()
      cluster_set_names=vis_freqs_cluster_sets_reactive()
      freq_norm=normalize_by_clusterset_frequency(get_cell_to_cluster(session),get_cell_to_sample(session),insamples,cluster_sets,pool_subtype = T,reg = 0.00)

      celltypes=intersect(cluster_set_names,names(freq_norm)[(!sapply(freq_norm,is.null))])
      celltypes=celltypes[sapply(freq_norm[celltypes],ncol)>1]
      layout(matrix(1:(length(celltypes)*2),length(celltypes),2,byrow = T))
      par(mar=c(8,6,2,1))
      for (cluster_set_name in celltypes){
        plot_subtype_freqs(freq_norm,cluster_set_name,plot_legend = T,cex.names=1.5,cex.axis=1.5,cex.legend=1.5,cluster_set_name=cluster_set_name)
      }
      
    })
    
   output$correlation_betwen_clusters <-renderPlot({

     cgs=session$userData$reactiveVars$clusters_genes_samples_reactive
    cluster_sets=cluster_sets_reactive()
     if (is.null(cluster_sets)){
       return()
     }
      inclusters=rev(unlist(cluster_sets))
    
     models=get_models(session,clusters=inclusters)
     gene_mask=apply(models,1,max)>5e-5
     m=log2(1e-5+models[gene_mask,])
     m=m-rowMeans(m)
     cormat=cor(m)
  #   d=as.dist(1-cormat)
  #   order=seriate(d,method = "GW")
  #   ord=get_order(order)
     
     image(cormat,col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
     box(lwd=2)
     mtext(rownames(cormat),1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.8,line =.5)
     mtext(rownames(cormat),2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.8,line =.5)
   })
    
    output$clusters_sets_shinytree <- renderTree({
      return(as_list_recursive(session$userData$cluster_sets))
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
    
    
    
    
    if (exists(".scDissector_preloaded_data")){

      hideTab(inputId = "inMain", target = "Data")
    }
    
    
    
    ########################################################################
    
    
    
    
    
    
    # A reactive expression with the ggvis plot
    vis <- reactive({
      # Lables for axes
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      # Normally we could do something like props(x = ~BoxOffice, y = ~Reviews),
      # but since the inputs are strings, we need to do a little more work.
      #   xvar <- prop("x", as.symbol(input$xvar))
      #    yvar <- prop("y", as.symbol(input$yvar))
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
      samp_to_cat=sample_to_category()[insamples]
      cats=unique(samp_to_cat)
      cat_cols= session$userData$sample_colors[1:length(cats)]
      names(cat_cols)=cats
      sample_cols=cat_cols[samp_to_cat]

      
      ds=get_dstab(session,cells=ds_cells_QC_reactive(),ds_version = match(input$inQCDownSamplingVersion,get_ds_options(session)))
      if (is.null(ds)){
        return()
      }
      if (ncol(ds)<2){
        message("Warning: Cluster ",clust, " has less than 2 cells above the UMIs threshold.")
        return()
      }

      s1=Matrix::rowSums(ds,na.rm=T) 
      s2=Matrix::rowSums(ds^2,na.rm=T)
      ds_mean=s1/ncol(ds)
      m2=s2/ncol(ds)
      ds_var=m2-ds_mean^2
      m=ds_mean
      v=ds_var
      genemask=rownames(ds)[m>1e-1]
      var_genes=head(genemask[order(v[genemask]/m[genemask],decreasing=T)],200)
      z=ds[var_genes,]
      z=z/mean(z)
      
      cormat=cor(as.matrix(z),use="comp")
 
      ord=hclust(dist(1-cormat))$order
      samps=get_cell_to_sample(session,cells=colnames(cormat)[ord])
      layout(matrix(c(1:6),nrow = 3),heights = c(8,1,1),widths=c(9,1))
      par(mar=c(.5,2,2,2))
      ntotcells=length(select_cells(session,clusters = clust,samples = insamples))
      image(cormat[ord,ord],col=greenred(100),breaks=c(-1,seq(-1,1,l=99),1),axes=F,main=paste("Cluster",clust,": ",ncol(cormat),"/",ntotcells,"cells"))
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
      insamples=samples_reactive()
      if (!session$userData$loaded_flag){
        return()
      }
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      ds=get_dstab(session,cells=select_cells(session,clusters = clust,samples = insamples),ds_version = match(input$inQCDownSamplingVersion,get_ds_options(session)))
      
      if (ncol(ds)==0){
        return()
      }

      flag=sum(c(input$inGene1,input$inGene2)%in%rownames(ds))==2
      if (flag){
        dat=floor(log2(1+data.frame(t(as.matrix(ds[c(input$inGene1,input$inGene2),])))))
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
    
    
    if (exists(".scDissector_clustering_data_path")){
      load_metadata(session,.scDissector_clustering_data_path)
      updateTextInput(session,"inDatapath",,.scDissector_clustering_data_path)
    }
    if (exists(".scDissector_preloaded_data")){
      
      hideTab(inputId = "inMain", target = "Data")
      ldm=.scDissector_preloaded_data
      if (!exists(".scDissector_clustering_data_path")){
        scDissector_datadir=""
      }
      session$userData$loaded_model_file<-ldm$model$model_filename

      update_all(session,ldm)
      show_all_tabs()
      updateTabsetPanel(session, "inMain", selected = "MetaData")
    }
    
    
    
  }





