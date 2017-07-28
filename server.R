library(ggvis)
library(dplyr)
library(gplots)
source("projector.r")
set.seed(3505)


clean_all=function(){
  if (exists("projected")){
    rm(projected,envir = .GlobalEnv)  
  }
  if (exists("ext_prof")){
    rm(ext_prof,envir= .GlobalEnv)
  }
  if (exists("modules")){
    rm(modules,envir=.GlobalEnv)
  }
  if (exists("model")){
    rm(model,envir =.GlobalEnv)
  }
  if (exists("dataset")){
    rm(dataset,envir =.GlobalEnv)
  }
}

clean_all()

cap <- function(x) {
  paste(toupper(substring(x, 1,1)), tolower(substring(x, 2)),sep="", collapse=" ")
}
colgrad=c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))
sample_cols<<-rep(paste("#",read.table("sample_colors.txt",stringsAsFactors = F)[,1],sep=""),10)

click_flag<<-T
loaded_flag<<-F
loaded_version=NA

prev_inModelColorScale_abs<<-c(-5,-2)
prev_inModelColorScale_rel<<-c(-2,4)
print(getwd())
genesetsfile="gene_sets.txt"
geneList<<-read.table(file=genesetsfile,header=T,stringsAsFactors = F,row.names =1)
hgnc<<-read.delim("hgnc_complete_set.txt",header = T,stringsAsFactors = F)
old_symbol=ifelse(hgnc[,"prev_symbol"]=="",hgnc[,"symbol"],hgnc[,"prev_symbol"])
l_old_symbol=strsplit(old_symbol,"\\|")
old_symbol2=unlist(l_old_symbol)

new_symbol=hgnc[,"symbol"]
new_symbol2=rep(new_symbol,sapply(l_old_symbol,length))

gene_symbol_old2new<<-new_symbol2
names(gene_symbol_old2new)<<-old_symbol2

gene_symbol_new2old<<-old_symbol
names(gene_symbol_new2old)<<-new_symbol

referenceSets<-read.csv(file="ReferenceProfiles/refereceSets.csv",header=T,stringsAsFactors = F)
tab3_left_margin=12
  function(input, output, session) {
   
    updateSelectInput(session,"inReferenceSets",label = "Reference Set:",choices = referenceSets$title)  
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
    
    if (read_flag){
      vers_tab<<-read.csv(file=verfile,header=T,stringsAsFactors = F)
      samples_tab<<-read.csv(file=samples_file,header=T,stringsAsFactors = F)
  
      scDissector_datadir<<-input$inDatapath
      updateSelectInput(session,"inModelVer", "Model Version:",choices = vers_tab$title)
      updateSelectInput(session,"inSampleToAdd", "Samples:",choices =samples_tab$index)
      updateSelectInput(session,"inProjectedDataset","Projected Version:",choices = vers_tab$title)
    }
    else{
      updateSelectInput(session,"inModelVer", "Model Version:",choices = "")
      updateSelectInput(session,"inProjectedDataset","Projected Version:",choices="")
      return()
    }

   
      myGeneListFile= paste(input$inDatapath,"gene_sets.txt",sep="/")
   
      if (file.exists(myGeneListFile)){
        geneList<<-rbind( read.table(file=myGeneListFile,header=T,stringsAsFactors = F,row.names =1),geneList)
      }
    
    updateSelectInput(session,"inGeneSets",choices = rownames(geneList))
    
    
    mySampleSetsFile= paste(input$inDatapath,"sample_sets.txt",sep="/")
    
    if (file.exists(mySampleSetsFile)){
      sample_sets_tab<-read.table(file=mySampleSetsFile,header=T,stringsAsFactors = F,row.names =1)
      sample_sets<<-strsplit(sample_sets_tab[,"samples"],",| ,|, ")
      names(sample_sets)<<-rownames(sample_sets_tab)
      
      updateSelectInput(session,"inSampleToAdd", "Samples:",choices = c(names(sample_sets),samples_tab$index))
    }
    
    tfs_file="tfs.csv"
    if (file.exists(tfs_file)){
      tfs<<-read.csv(file=tfs_file,header=F,stringsAsFactors = F)[,1]
      tfs2=paste(substring(tfs, 1,1), tolower(substring(tfs, 2)),sep="")
       tfs<<-c(tfs,tfs2)
    }
    else{
      message(tfs_file ," not found")
    }
    
    surface_markers_file="surface_markers.csv"
    if (file.exists(surface_markers_file)){
      surface_markers<<-read.csv(file=surface_markers_file,header=F,stringsAsFactors = F)[,1]
      surface_markers2=paste(substring(surface_markers, 1,1), tolower(substring(surface_markers, 2)),sep="")
      surface_markers<<-c(surface_markers,surface_markers2)
    }
    else{
      message(surface_markers_file ," not found")
    }
    
 
  })
  
  observeEvent(input$inLoad,{
    if (input$inModelVer=="")
    {
      return()
    }
    
    i=match(input$inModelVer,vers_tab$title)
    f=paste(input$inDatapath,vers_tab$path[ifelse(is.na(i),1,i)],sep="/")
    
    if (!file.exists(f)){
      message(f, " not found")
      return()
    }
  
    #rm(c())
    clean_all()
    message("Loading ",f)
    model<<-new.env(parent = globalenv())
    load(file=f,model)
    
    
    loaded_file<<-f
 
    loaded_version<<-input$inModelVer
  
    is_cluster_excluded<<-rep(F,ncol(model$models))
    names(is_cluster_excluded)<<-colnames(model$models)
    
    clusterset_fn=paste(strsplit(loaded_file,"\\.")[[1]][1],"_clustersets.txt",sep="")
    if (file.exists(clusterset_fn)){
      cluster_sets_tab=read.table(file=clusterset_fn,header = T,stringsAsFactors = F)
      if (nrow(cluster_sets_tab)>0){
        cluster_sets<<-strsplit(cluster_sets_tab[,2],",")
        names(cluster_sets)<<-cluster_sets_tab[,1]
      
        is_cluster_excluded[unlist(cluster_sets)]<<-rep(cluster_sets_tab[,3],sapply(cluster_sets,length))
        excluded_cluster_sets<<-cluster_sets_tab[cluster_sets_tab[,3],1]
        updateSelectInput(session,"inBirdUndefineClusterSetSelect","Cluster Set",choices =names(cluster_sets)) 
        updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(cluster_sets),excluded_cluster_sets)) 
        updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =excluded_cluster_sets) 
      }
    }
    else{
      cluster_sets<<-list()
    }
    
    
    
    
    annnot_fn=paste(strsplit(f,"\\.")[[1]][1],"_annots.txt",sep="")
    if (file.exists(annnot_fn)){
      a=read.delim(annnot_fn,header=F,stringsAsFactors = F,row.names = 1)
      clustAnnots<<-a[,1]
      clustAnnots[is.na(clustAnnots)]<<-""
      names(clustAnnots)<<-rownames(a)
    }
    else{
      clustAnnots<<-rep("",ncol(model$models))
      names(clustAnnots)<<-colnames(model$models)
    }    
   
    order_fn=paste(strsplit(loaded_file,"\\.")[[1]][1],"_order.txt",sep="")
    if (file.exists(order_fn)){
      cluster_order=as.character(read.table(file=order_fn,header = F)[,1])
    }
    else{
      cluster_order=colnames(model$models)
    }
    
   
    
    samples=strsplit(input$inSamples,",| ,|, ")[[1]]
    tmp_samples=as.list(samples)
  
    if (exists("sample_sets")){
      mask=samples%in%names(sample_sets)
      tmp_samples[mask]=sample_sets[samples[mask]]
    }
    samples=unique(unlist(tmp_samples))
    if (!all(samples%in%samples_tab$index)){
      message("Error: These samples could not be recognized: ",paste(samples[!(samples%in%samples_tab$index)],collapse=","),". Loading failed!")
      message("")
      return()
    }
    
    tmp_dataset=new.env()
    tmp_dataset$umitab=list()
    tmp_dataset$ll=list()
    tmp_dataset$cell_to_cluster
    tmp_dataset$avg=list()
    tmp_dataset$ds=list()
    tmp_dataset$ds_numis=NULL
    tmp_dataset$counts=list()
   
    genes=rownames(model$models)
    for (sampi in samples){
      message("Loading sample ",sampi)
      i=match(sampi,samples_tab$index)
      f=paste(input$inDatapath,samples_tab$path[i],sep="/")
      tmp_env=new.env()
      load(f,envir = tmp_env)
      
      colnames(tmp_env$umitab)=paste(samples[sampi],colnames(tmp_env$umitab),sep="_")
      for (ds_i in 1:length(tmp_env$ds)){
        colnames(tmp_env$ds[[ds_i]])=paste(samples[sampi],colnames(tmp_env$ds[[ds_i]]),sep="_")
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
      tmp_dataset$umitab[[sampi]]=tmp_env$umitab
      message("Projecting ",ncol(tmp_dataset$umitab[[sampi]])," cells")
      genemask=intersect(rownames(tmp_dataset$umitab[[sampi]]),rownames(model$models))
      
      genes=intersect(genemask,genes)
    #  print(length(genes))
      tmp_dataset$ll[[sampi]]=getLikelihood(tmp_dataset$umitab[[sampi]][genemask,],models =model$models[genemask,],reg = 1e-5)#params$reg)
      tmp_dataset$cell_to_cluster[[sampi]]=MAP(tmp_dataset$ll[[sampi]])
      tmp_dataset$avg[[sampi]]=matrix(0, nrow(tmp_dataset$umitab[[sampi]]),ncol(model$models),dimnames = list(rownames(tmp_dataset$umitab[[sampi]]),colnames(model$models)))
      tmpmod=update_models(tmp_dataset$umitab[[sampi]],tmp_dataset$cell_to_cluster[[sampi]])
      tmp_dataset$avg[[sampi]][,colnames(tmpmod)]=tmpmod
      tmp_counts=sapply(split_sparse(tmp_dataset$umitab[[sampi]][genemask,],tmp_dataset$cell_to_cluster[[sampi]]),rowSums)
      tmp_dataset$counts[[sampi]]=model$models*0
      tmp_dataset$counts[[sampi]][genemask,colnames(tmp_counts)]=tmp_counts
      
      if (length(tmp_dataset$ds)==0){
        for (ds_i in 1:length(tmp_env$ds_numis)){
          tmp_dataset$ds[[ds_i]]=list()
        }
      }
      for (ds_i in 1:length(tmp_env$ds_numis)){
        tmp_dataset$ds[[ds_i]][[sampi]]=tmp_env$ds[[ds_i]]
      }
      message("")
     
    }
    
    message("One more minute...")
    
    dataset<<-new.env()
    dataset$ds_numis=tmp_dataset$ds_numis
    dataset$umitab<<-tmp_dataset$umitab[[samples[1]]][genes,]
    dataset$ll<<-tmp_dataset$ll[[samples[1]]]
    dataset$cell_to_cluster<<-tmp_dataset$cell_to_cluster[[samples[1]]]
    names(dataset$cell_to_cluster)=colnames(tmp_dataset$umitab[[samples[1]]])
    dataset$cell_to_sample<<-rep(samples[1],ncol(tmp_dataset$umitab[[samples[1]]]))
    names(dataset$cell_to_sample)=colnames(tmp_dataset$umitab[[samples[1]]])
    dataset$avg<<-tmp_dataset$avg
    dataset$counts<<-tmp_dataset$counts
    dataset$samples=samples
    dataset$ds<<-list()
    dataset$randomly_selected_cells<<-list()
    
    for (ds_i in 1:length(tmp_env$ds_numis)){
      ds_sampi=tmp_dataset$ds[[ds_i]][[samples[1]]][genes,]
      dataset$ds[[ds_i]]=ds_sampi
      dataset$randomly_selected_cells[[ds_i]]<<-list()
      for (randomi in 1:length(params$nrandom_cells_per_sample_choices)){
        nrandom_cells=params$nrandom_cells_per_sample_choices[randomi]
        if (nrandom_cells=="All"||as.numeric(nrandom_cells)>=ncol(ds_sampi)){
          dataset$randomly_selected_cells[[ds_i]][[randomi]]<<-colnames(ds_sampi)
        }
        else{
          dataset$randomly_selected_cells[[ds_i]][[randomi]]<<-sample(colnames(ds_sampi),size=as.numeric(nrandom_cells),replace=F)
        }
      }
    }
   
    if (length(samples)>1){
      for (sampi in samples[-1]){
        cellids=colnames(tmp_dataset$umitab[[sampi]][genes,])
        dataset$umitab<-cBind(dataset$umitab,tmp_dataset$umitab[[sampi]][genes,])
        dataset$ll<-rbind(dataset$ll,tmp_dataset$ll[[sampi]])
        dataset$cell_to_cluster<<-c(dataset$cell_to_cluster,tmp_dataset$cell_to_cluster[[sampi]])
        cell_to_sampi=rep(sampi,length(cellids))
        names(cell_to_sampi)=cellids
        dataset$cell_to_sample<<-c(dataset$cell_to_sample,cell_to_sampi)
        
        for (ds_i in 1:length(dataset$ds_numis)){
          ds_sampi=tmp_dataset$ds[[ds_i]][[sampi]][genes,]
          dataset$ds[[ds_i]]=cBind(dataset$ds[[ds_i]],ds_sampi)
         
          
          for (randomi in 1:length(params$nrandom_cells_per_sample_choices)){
            nrandom_cells=params$nrandom_cells_per_sample_choices[randomi]
            if (nrandom_cells=="All"||as.numeric(nrandom_cells)>=ncol(ds_sampi)){
              dataset$randomly_selected_cells[[ds_i]][[randomi]]<<-c(dataset$randomly_selected_cells[[ds_i]][[randomi]],colnames(ds_sampi))
            }
            else{
              dataset$randomly_selected_cells[[ds_i]][[randomi]]<<-c(dataset$randomly_selected_cells[[ds_i]][[randomi]],sample(colnames(ds_sampi),size=as.numeric(nrandom_cells),replace=F))
            }
          }
        }
      }
    }
    
    ncells_per_cluster<<-rep(0,dim(model$models)[2])
    names(ncells_per_cluster)<<-colnames(model$models)
    temptab=table(model$cell_to_cluster)
    ncells_per_cluster[names(temptab)]<<-temptab
    
    clust_title=paste(cluster_order," - ",clustAnnots[cluster_order],sep="") 
    default_clusters<<-cluster_order
    
    updateTextInput(session,"inClusters",value=paste(cluster_order,collapse=",")) 
    updateSelectInput(session,"inAnnotateClusterNum",choices = clust_title)
#    updateSelectInput(session,"inTweezersFromCluster",choices = clust_title)
#    updateSelectInput(session,"inTweezersToCluster",choices = clust_title)
    updateSelectInput(session,"inQCClust",choices =clust_title)
#    updateSelectInput(session,"inTweezersLLX",choices = clust_title)
#    updateSelectInput(session,"inTweezersLLY",choices = clust_title)
    updateSelectInput(session,"inClustForDiffGeneExprsProjVsRef",choices = clust_title)
    updateTextInput(session,"inTruthSamples",value = paste(samples,collapse=","))
    updateSelectInput(session,"inTruthDownSamplingVersion",choices=dataset$ds_numis,selected = max(dataset$ds_numis))
    updateSelectInput(session,"inQCDownSamplingVersion",choices=dataset$ds_numis,selected = max(dataset$ds_numis))
    updateSelectInput(session,"inModulesDownSamplingVersion",choices=dataset$ds_numis,selected = max(dataset$ds_numis))
    message("Successfully finished loading.")
    loaded_flag<<-T
  })
  
  
  observeEvent(input$inAddSample,{
    s1=paste(unique(c(strsplit(input$inSamples,",|, | ,")[[1]],input$inSampleToAdd)),collapse=", ")
    updateTextInput(session,"inSamples",value=s1) 
  })
  
  
  observeEvent(input$inAbsOrRel,{
    if (input$inAbsOrRel=="Absolute"){
      
      prev_inModelColorScale_rel<<-input$inModelColorScale
      
      updateSliderInput(session,"inModelColorScale",label="Log10(expression)",min = -12,max=-.5,step = .5,value = prev_inModelColorScale_abs)
    }
    else if (input$inAbsOrRel=="Relative"){
      prev_inModelColorScale_abs<<-input$inModelColorScale
      updateSliderInput(session,"inModelColorScale",label="Log2(expression/mean)",min = -8,max = 8,step = 1,value = prev_inModelColorScale_rel)
      
    }
  })
  
  genes_reactive <-reactive({
    
    if (!loaded_flag){
      return()
    }
    
    return(init_genes(input$inGenes) )
  })
  
  clusters_reactive <-reactive({
    clusts=strsplit(input$inClusters,",")[[1]]
    clusts=setdiff(intersect(clusts,colnames(model$models)),names(which(is_cluster_excluded)))
    return(clusts)
  })
  
  truth_samples_reactive <-reactive({
    samples=strsplit(input$inTruthSamples,",")[[1]]
     if (!exists("dataset")){
      return(c())
    }
    samples=intersect(samples,dataset$samples)
    return(samples)
  })
  
  modules_reactive <-reactive({
    modules=strsplit(input$inModules,",")[[1]]
    return(modules)
  })
  
  
  ds_QC_reactive <-reactive({
    clust=strsplit(input$inQCClust," - ")[[1]][1]
    samples=truth_samples_reactive()
    if (!exists("dataset")){
      return()
    }
    
    ds=dataset$ds[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]]
    sampling_mask=dataset$randomly_selected_cells[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]][[match(input$inQCNcellsPerSample,params$nrandom_cells_per_sample_choices)]]
    cluster_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster==clust]
    return(ds[,intersect(sampling_mask,cluster_mask),drop=F])
  })
  
  observeEvent(input$inGeneSets,{
    geneSets=input$inGeneSets
    
 #   if (!loaded_flag){
#      return()
  #  }
    updateTextInput(session,"inGenes",value=geneList[geneSets,1])
    #  genes=strsplit(geneList[geneSets,1],",")[[1]]
  #  if (loaded_flag){
  #    init_genes(genes)
  #  }
    
    })
  
  
  init_genes=function(genes_string){
    genes=strsplit(genes_string,",|, | ,")[[1]]
    mask1=toupper(genes)%in%rownames(model$models)
    genes[mask1]=toupper(genes[mask1])
    mask2=unlist(sapply(genes,cap))%in%rownames(model$models)
    genes[mask2]=unlist(sapply(genes[mask2],cap))
    mask3=gene_symbol_old2new[genes]%in%rownames(model$models)
    genes[mask3]=gene_symbol_old2new[genes[mask3]]
    mask4=gene_symbol_new2old[genes]%in%rownames(model$models)
    genes[mask4]=gene_symbol_new2old[genes[mask4]]
    
   # genes=genes[(mask1|mask2|mask3)]
    
    gcol<<-ifelse(toupper(genes)%in%tfs ,4,1)
    names(gcol)<<-toupper(genes)
    return(genes)
  #  updateTextInput(session,"inGenes",value=genes)
    
  }
  
 
  
  observeEvent(input$inAnnotateCluster, {
    clustAnnots[strsplit(input$inAnnotateClusterNum," - ")[[1]][1]]<<-input$inClustAnnot
    updateTextInput(session,"inClustAnnot",value="")
  })
  
  observeEvent(input$inSaveAnnot, {
    f=paste(input$inDatapath,vers_tab$path[vers_tab$title==loaded_version],sep="/")
 
    annot_fn=paste(strsplit(f,"\\.")[[1]][1],"_annots.txt",sep="")
    write.table(file=annot_fn,clustAnnots,row.names=T,col.names=F,quote=F,sep="\t")
    message("Cluster annotations saved to ",annot_fn,".")
    
  })
  
  
  
  read_reference_set=function(){
    fn=paste("ReferenceProfiles/",referenceSets[referenceSets$title==input$inReferenceSets,"path"],sep="")
    message("Reading ",fn)
    ext_prof<<-read.delim(fn,stringsAsFactors = F,row.names = 1,header=T)
    message("Done.")
  }
  observeEvent(input$inSelectHighlyCorrelatedProfiles,{
    if (!exists("ext_prof")){
      read_reference_set()
    }
    clusts=clusters_reactive()
    genes=genes_reactive()
    gmask=intersect(genes,intersect(toupper(rownames(model$models)),rownames(ext_prof)))
  
    cormat=cor(model$models[match(gmask,toupper(rownames(model$models))),clusts],ext_prof[gmask,],method="spearman",use="comp")
    profs=unique(as.vector(apply(cormat,1,function(x){colnames(ext_prof)[order(x,decreasing=T)[1:2]]})))
  
    updateTextInput(session, "inRefProfiles",value=paste(profs,collapse=","))
  })
  
  observeEvent(input$inSelectAllReferenceProfiles,{
 
    read_reference_set()
    
    updateTextInput(session, "inRefProfiles",value=paste(colnames(ext_prof),collapse=","))
  })
    
  observeEvent(input$inClusterGenes, {
    if (input$inReorderingMethod=="Hierarchical clustering"){
      mat2=mat[!rowSums(is.na(mat))==ncol(mat),]
      cormat=cor(t(mat2),use="comp")
      mask=rowSums(!is.na(cormat))>1
      cormat=cormat[mask,mask]
    
      ord=hclust(dist(1-cormat))$order
      updateTextInput(session,"inGenes",value=paste(c(rownames(cormat)[ord],setdiff(rownames(mat),rownames(cormat))),collapse=","))
    }
    else if (input$inReorderingMethod=="Diagonal"){
        clusters=clusters_reactive()
        ord=order(apply(mat[,clusters],1,which.max))
        updateTextInput(session,"inGenes",value=paste(rownames(mat)[ord],collapse=","))
      }
  })
  
  observeEvent(input$inClusterModules, {
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
    print(modules[ord])
  })
  
  
#  observeEvent(input$inFindVarGenes, {
#    message("screening for variable genes")
#    clusters=clusters_reactive()
#    mask=cell_to_cluster[colnames(ds)]%in%clusters
#    l=split(as.data.frame(t(ds[,mask])),cell_to_cluster[colnames(ds)[mask]])
#    counts=t(sapply(l,colSums,na.rm=T))
    
#    vs=rowSums(sapply(l,function(x){rowSums((t(x)-colMeans(x))^2)}))/(sum(mask)-1)
#    gene_scores<<-1-vs/apply(ds[,mask],1,var)
     
#    m1=apply(counts[,names(gene_scores)],2,max)>=50&gene_scores>0.02
#    lgenes=split(names(m1)[m1],apply(model$models[names(m1),][m1,],1,which.max))
#    isolate({
#    ngenes_to_show=as.numeric(input$inNgenes)
#    })
  
  observeEvent(input$inBlindChisq, {
       message("screening for variable ",input$inSelectGenesFrom)
        clusters=clusters_reactive()
        chisq_res=model$chisq_res
        if (input$inSelectGenesFrom=="All genes"){
          pref="Var"
          genes=rownames(chisq_res)
        } else if (input$inSelectGenesFrom=="TFs"){
          pref="Tfs"
          genes=tfs
        } else if(input$inSelectGenesFrom=="Surface markers"){
          pref="Surf"
          genes=surface_markers
        }
        mask=rownames(chisq_res)%in%genes&chisq_res[,3]<0.05&(apply(model$models[rownames(chisq_res),]/(rowSums(model$umitab[rownames(chisq_res),])/sum(model$umitab)),1,max)>4|apply(model$models[rownames(chisq_res),]/rowMeans(model$models[rownames(chisq_res),]),1,max)>4)
         isolate({
          ngenes_to_show=as.numeric(input$inNgenes)
        })
        a=chisq_res[mask,]
          genes_to_show=head(rownames(a)[order(a[,2],decreasing=T)],ngenes_to_show)
    message("done")
  
    geneList2=rbind(geneList,paste(genes_to_show,collapse=","))
    rownames(geneList2)[1:nrow(geneList)]=rownames(geneList)
    rownames(geneList2)[nrow(geneList2)]=paste("Chisq_",pref,"_",ngenes_to_show,"_",date(),sep="")
    geneList<<-geneList2
    updateSelectInput(session,"inGeneSets",choices=rownames(geneList),selected = rownames(geneList)[nrow(geneList)])
    
  })
  
  
  
  
  observeEvent(input$inBlindChisqSelectedClusters, {
    message("screening for variable ",input$inSelectGenesFrom)
    clusters=clusters_reactive()
    if (input$inSelectGenesFrom=="All genes"){
      pref="Var"
      genes=rownames(model$models)
    } else if (input$inSelectGenesFrom=="TFs"){
      pref="Tfs"
      genes=intersect(tfs,rownames(model$models))
    } else if(input$inSelectGenesFrom=="Surface markers"){
      pref="Surf"
      genes=surface_markers
    }
    cells=names(model$cell_to_cluster[model$cell_to_cluster%in%clusters])
    chisq_res2=chisq_genes(model$umitab,model$cell_to_cluster,genes,cells)
    mask=rownames(chisq_res2)%in%genes&chisq_res2[,3]<0.05&(apply(model$models[rownames(chisq_res2),]/(rowSums(model$umitab[rownames(chisq_res2),])/sum(model$umitab)),1,max)>4|apply(model$models[rownames(chisq_res2),]/rowMeans(model$models[rownames(chisq_res2),]),1,max)>4)
    isolate({
      ngenes_to_show=as.numeric(input$inNgenes)
    })
    a=chisq_res2[mask,]
    genes_to_show=head(rownames(a)[order(a[,2],decreasing=T)],ngenes_to_show)
    
    message("done")
    
    geneList2=rbind(geneList,paste(genes_to_show,collapse=","))
    rownames(geneList2)[1:nrow(geneList)]=rownames(geneList)
    rownames(geneList2)[nrow(geneList2)]=paste("Chisq_",pref,"_",ngenes_to_show,"_",date(),sep="")
    geneList<<-geneList2
    updateSelectInput(session,"inGeneSets",choices=rownames(geneList),selected = rownames(geneList)[nrow(geneList)])
    # updateTextInput(session,"inGenes",value=paste(genes_to_show,collapse=","))
    
  })
  
  chisq_genes=function(umitab,cell_to_cluster,genes,cells){
    
    umitab2=umitab[genes,cells]
    cluster_tot=sapply(split(colSums(umitab2),cell_to_cluster[colnames(umitab2)]),sum)
    counts=sapply(split_sparse(umitab2,cell_to_cluster[colnames(umitab2)]),rowSums)
    gene_mask=apply(counts,1,max)>3
    counts=counts[gene_mask,]
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
      genes=rownames(model$models)
    } else if (input$inSelectGenesFrom=="TFs"){
      pref="Tfs"
      genes=tfs
    } else if(input$inSelectGenesFrom=="Surface markers"){
      pref="Surf"
      genes=surface_markers
    }
    mask=rownames(model$umitab)%in%genes
    
    clusts_fg=strsplit(input$inFC_fgClusts,",")[[1]]
    clusts_bg=strsplit(input$inFC_bgClusts,",")[[1]]
    mask_fg=model$cell_to_cluster%in%clusts_fg
    mask_bg=model$cell_to_cluster%in%clusts_bg
    reg=50
    fc=(rowSums(reg+model$umitab[mask,mask_fg])/sum(reg+model$umitab[,mask_fg]))/(rowSums(reg+model$umitab[mask,mask_bg])/sum(reg+model$umitab[,mask_bg]))
    
    isolate({
      ngenes_to_show=as.numeric(input$inNgenes)
    })
    
    genes_to_show=head(names(sort(fc,decreasing=T)),floor(ngenes_to_show/2))
    genes_to_show=c(genes_to_show,tail(names(sort(fc,decreasing=T)),ceiling(ngenes_to_show/2)))
    
    geneList2=rbind(geneList,paste(genes_to_show,collapse=","))
    rownames(geneList2)[1:nrow(geneList)]=rownames(geneList)
    rownames(geneList2)[nrow(geneList2)]=paste("FC_",pref,"_",ngenes_to_show,"_",paste(clusts_fg,collapse="_"),"_VS_",paste(clusts_bg,collapse="_"),sep="")
    geneList<<-geneList2
     updateSelectInput(session,"inGeneSets",choices=rownames(geneList),selected = rownames(geneList)[nrow(geneList)])
  })
  
  
  observeEvent(input$inResetClusters, {
     clust_title=paste(default_clusters," - ",clustAnnots[default_clusters],sep="") 
    updateTextInput(session,"inClusters",value=paste(default_clusters,collapse=","))
    updateSelectInput(session,"inQCClust",choices=clust_title)
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title)
    
  })
  
  observeEvent(input$inSaveClusterOrder, {
    clusters=clusters_reactive()
   
    order_fn=paste(strsplit(loaded_file,"\\.")[[1]][1],"_order.txt",sep="")
    if (all(colnames(model$models)%in%clusters)){
      default_clusters<<-clusters
      write.table(clusters,file=order_fn,quote = F,sep="\t",row.names = F,col.names = F)
      message("Order has been successfully saved.")
    }else{
      message("Order could not be saved since 1 or more clusters are missing.")
    }
  })
  
  observeEvent(input$detectDiffExrs, {
    if (!exists('model$models')){
      return()
    }
    clusters=clusters_reactive()
    cell_mask=model$cell_to_cluster%in%clusters
   #tot (=#mols per cluster)
     tot=sapply(split(colSums(model$umitab[rownames(model$ds),cell_mask]),model$cell_to_cluster[cell_mask][colnames(model$umitab)[cell_mask]]),sum)
 
    fc=model$models[rownames(model$ds),]/(rowSums(model$ds)/sum(model$ds))
    fc_mask=apply(fc,1,max)>input$inputMinFC&rowMeans(model$ds)>input$inMinAvgExprs
    print(sum(fc_mask))
    counts=sapply(split(as.data.frame(t(model$umitab[rownames(model$ds)[fc_mask],cell_mask])),model$cell_to_cluster[cell_mask][colnames(model$umitab)[cell_mask]]),colSums,na.rm=T)
  
    arrcont=array(c(counts,matrix(tot,dim(counts)[1],dim(counts)[2],byrow=T)-counts),dim=c(dim(counts),2))
    res=t(apply(arrcont,1,function(x){unlist(chisq.test(x)[c("p.value","statistic")])}))
    rownames(res)=rownames(counts)
    res=res[!is.na(res[,1]),]
    res=cbind(res,adjp=p.adjust(res[,1],method="BH"))
    
    chisq.res<<-res[res[,3]<0.05,]
    chisq.res<<-chisq.res[order(chisq.res[,2],decreasing=T),]
  
    fdrmask=res[,3]<input$inputFDR
    names(fdrmask)=rownames(res)
    
    signiftab<<-log2(2^-20+model$models[names(fdrmask),])
    signiftab<<-round(signiftab[order(apply(signiftab,1,which.max)),],digits=1)
    print("done.")
  })
  
  output$downloadSignifTable <- downloadHandler(
    filename = function() { paste("DiffExprs", '.csv', sep='') },
    
    content = function(file) {
      write.csv(signiftab, file)
    }
  )
  outputOptions(output, 'downloadSignifTable', suspendWhenHidden=FALSE)  
  
  output$downloadExprsTable <- downloadHandler(
    filename = function() { paste("Exprs", '.csv', sep='') },
    content = function(file) {
      write.csv(model$models, file)
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
    
    if (input$inReorderingClustersMethod=="Hierarchical clustering"){
  
    cormat=cor(mat[,inclusts],use="comp")
    d1=dist(1-cormat)
    d1[is.na(d1)]=100
    reorderv1=hclust(d1)$order
    
    }
    else if (input$inReorderingClustersMethod=="Diagonal"){
      genes=genes_reactive()
      reorderv1=order(apply(mat[genes,],2,which.max))
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
      inclusts=clusters_reactive()
      score=rep(0,length(inclusts))
      for (i in 1:length(genes)){
        message(sig[i],genes[i])
  
        score=score+(10^(5*(length(genes)-i)))*ifelse(sig[i]=="+",1,-1)*(mat[genes[i],inclusts])
      }
      
      reorderv1=order(score)
    }
    clust_title=paste(inclusts," - ",clustAnnots[inclusts],sep="") 
    
    updateTextInput(session,"inClusters",value=paste(inclusts[reorderv1],collapse=","))
    updateSelectInput(session,"inQCClust",choices=clust_title[reorderv1])
    updateSelectInput(session,"inAnnotateClusterNum",choices=clust_title[reorderv1])
  })
  
  

  
 # observeEvent(input$inTweezersMoveCells, {
#    input_cells=strsplit(input$inTweezers_cells_to_move,",")[[1]]
#    from_cluster=strsplit(input$inTweezersFromCluster," - ")[[1]][1]
#    to_cluster=strsplit(input$inTweezersToCluster," - ")[[1]][1]
#    cells=intersect(colnames(model$umitab)[model$cell_to_cluster==from_cluster],input_cells)
#    if (length(cells)==0){
#      message("None of the input cells ",input$inTweezers_cells_to_move, " were found in cluster ", from_cluster)
#      return()
#    }
#    else {
#      message("Moving ",length(cells)," cells from cluster ", from_cluster, " to cluster ", to_cluster)
#      model$cell_to_cluster[cells]<<-to_cluster
#      model$models<<-update_model(model$umitab,model$cell_to_cluster)
#      ll<<-getLikelihood(model$umitab,model$models,reg=params$reg)
#      
#      message("Done.")
#    }
#  })
  
#  observeEvent(input$inTweezersSaveVersion, {
#    if (input$inTweezers_new_version_name%in%vers_tab$title){
#      message("Version name already exists. Version was not saved!")
#      return()
#    }
    
#    chisq_res<<-chisq_genes(model$umitab,model$cell_to_cluster)
#    fn=paste(input$inDatapath,input$inTweezers_new_version_filename,sep="/")
#    message("Saving ",input$inTweezers_new_version_name, " to ", fn)
#    save(model$umitab,ds,model$models,model$cell_to_cluster,ll,l_cells_per_cluster,cell_to_batch,chisq_res,params,file=fn)
#    vers<<-rbind(vers,c(input$inTweezers_new_version_name,input$inTweezers_new_version_filename))
#    write.csv(vers,file=paste(input$inDatapath,"model_versions.csv",sep="/"),row.names = F,quote=F)
    
#    annot_fn=paste(input$inDatapath,strsplit(input$inTweezers_new_version_filename,"\\.")[[1]][1],"_annots.txt",sep="")
#    write.table(file=annot_fn,clustAnnots,row.names=T,col.names=F,quote=F,sep="\t")
#    updateSelectInput(session,"inModelVer", "Clustering Version:",choices = vers_tab$title)
#    updateSelectInput(session,"inProjectedDataset","Projected Version:",choices = vers_tab$title)
#    message("Done.")
#  })
  ###########################################################
  observeEvent(input$inBirdDefineClusterSet,{
    clusts_to_add=strsplit(input$inBirdAddClusters,",| ,|, ")[[1]]
    if (!all(clusts_to_add%in%colnames(model$models))){
      updateTextInput(session,"inBirdAddClusters","Clusters:",value = paste(clusts_to_add[clusts_to_add%in%colnames(model$models)],collapse = ","))  
      message("Warning! Cluster list contained unknown clusters. Cluster-set was not defined.")
      return()
    }
    intersect_with_defined=intersect(clusts_to_add,unlist(cluster_sets))
    if (length(intersect_with_defined)>0){
      message("warning! Cluster list contained clusters that have been already defined as cluster sets: ",paste(intersect_with_defined,collapse=","),". Cluster-set was not defined.")
      return()
    }
    
    if (input$inBirdAddClusterSetName%in%names(cluster_sets)){
      message("warning! Name has already been used. Cluster-set was not defined.")
      return()
    }
    
    cluster_sets[[input$inBirdAddClusterSetName]]<<-clusts_to_add
    updateSelectInput(session,"inBirdRemoveClusterSetSelect","Cluster Set",choices =names(cluster_sets)) 
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =names(cluster_sets)) 
   })
  
  observeEvent(input$inBirdUndefineClusterSet,{
    is_cluster_excluded[cluster_sets[[input$inBirdUndefineClusterSetSelect]]]=F
    cluster_sets[[input$inBirdUndefineClusterSetSelect]]<<-NULL
    excluded_cluster_sets<<-setdiff(excluded_cluster_sets,input$inBirdUndefineClusterSetSelect)
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(cluster_sets),excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =excluded_cluster_sets)
    updateSelectInput(session,"inBirdUndefineClusterSetSelect","Cluster Set",choices =names(cluster_sets)) 
  })
  
  observeEvent(input$inBirdSaveClusterSets,{
    f=paste(input$inDatapath,vers_tab$path[vers_tab$title==loaded_version],sep="/")
  
    clusters_set_tab=data.frame(name=names(cluster_sets),clusters=sapply(cluster_sets,paste,collapse = ","),is_excluded=ifelse(names(cluster_sets)%in%excluded_cluster_sets,"T","F"))
    clusterset_fn=paste(strsplit(f,"\\.")[[1]][1],"_clustersets.txt",sep="")
    write.table(file=clusterset_fn,clusters_set_tab,row.names=F,col.names=T,quote=F,sep="\t")
   
  })
  
  observeEvent(input$inBirdClusterSetExclude,{
    
    is_cluster_excluded[cluster_sets[[input$inBirdClusterSetsExcludeSelect]]]<<-T
    excluded_cluster_sets<<-unique(c(excluded_cluster_sets,input$inBirdClusterSetsExcludeSelect))
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(cluster_sets),excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =excluded_cluster_sets)
  })
  
  observeEvent(input$inBirdClusterSetInclude,{
    
    is_cluster_excluded[cluster_sets[[input$inBirdClusterSetsIncludeSelect]]]<<-F
    excluded_cluster_sets<<-setdiff(excluded_cluster_sets,input$inBirdClusterSetsIncludeSelect)
    updateSelectInput(session,"inBirdClusterSetsExcludeSelect","Cluster Set",choices =setdiff(names(cluster_sets),excluded_cluster_sets))
    updateSelectInput(session,"inBirdClusterSetsIncludeSelect","Cluster Set",choices =excluded_cluster_sets)
  })
  
  ###########################################################
  modules_varmean_reactive=reactive({
    ds_i=match(input$inModulesDownSamplingVersion,dataset$ds_numis)
    
     if (!loaded_flag)
    {
      return(data.frame(m=0,v=0,gene=""))
    }
    
    
    
     ds=dataset$ds[[ds_i]][,dataset$randomly_selected_cells[[ds_i]][[match("All",params$nrandom_cells_per_sample_choices)]]]
    
    ds_mean<-rowMeans(ds)
    message("Estimating variance for ",nrow(ds)," genes")
    ds_var<-apply(ds,1,var)
    return(data.frame(m=ds_mean,v=ds_var,gene=rownames(ds)))
  })
  
  geneModuleMask_reactive<-reactive({
    if (!loaded_flag)
    {
      return()
    }
    df=modules_varmean_reactive()
    
    geneModuleMask<-log10(df$m)>as.numeric(input$inVarMean_MeanThresh)&log2(df$v/df$m)>as.numeric(input$inVarMean_varmeanThresh)
    geneModuleMask[is.na(geneModuleMask)]<-F
    return(geneModuleMask)
    
  })
   
  
  output$varMeanThreshPlot <- renderPlot({
    if (!loaded_flag)
    {
      return()
    }
    df=modules_varmean_reactive()
    plot(log10(df$m),log2(df$v/df$m),xlab="Log10(mean)",ylab="log2(var/mean)",panel.first=grid())
    abline(v=input$inVarMean_MeanThresh,col=2)
    
    abline(h=input$inVarMean_varmeanThresh,col=2)
    
    n=sum(geneModuleMask_reactive())
    legend("topright", paste(n,"genes"), bty="n",text.col=2) 
    
  })
  
  
  module_cor_reactive=reactive({
    if (!loaded_flag)
    {
      return()
    }
    ds_i=match(input$inModulesDownSamplingVersion,dataset$ds_numis)
    ds=dataset$ds[[ds_i]][,dataset$randomly_selected_cells[[ds_i]][[match("All",params$nrandom_cells_per_sample_choices)]]]
    message("calculating gene-to-gene correlations..")
    cormat=cor(as.matrix(t(log2(.1+ds[geneModuleMask_reactive(),]))),use="comp")
    return(cormat)
  })
  
  
  observeEvent(input$inGetModules,{
    if (!loaded_flag)
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
 
    modules<<-(split(names(c2c),c2c))
    updateTextInput(session,"inModules",value=paste(names(modules),collapse=","))
    updateSelectInput(session,"inModuleSelect",label="Show Module:",choices=names(modules))
    print(modules)
  })
  
  output$textModuleGenes<- renderText({
    input$inModules
    if (exists("modules")){
      paste(modules[[input$inModuleSelect]],collapse=",")
    }
  })
  
  ###########################################################
  
  output$ncells_barplot <-renderPlot({
    input$inModelVer
    clusters=clusters_reactive()
    cells=names(model$cell_to_cluster)[model$cell_to_cluster%in%clusters]
    ncells=table(model$cell_to_cluster[cells])[clusters]
    par(mar=c(1,7,2,2))
    barplot(ncells,xaxs = "i",log="y",ylab="#cells",names.arg ="")
  })
  
  output$UMI_boxplot <- renderPlot({
    input$inModelVer
    par(xaxs="i")
    clusters=clusters_reactive()
    cells=names(model$cell_to_cluster)[model$cell_to_cluster%in%clusters]
    par(mar=c(2,7,1,2))
    boxplot(split(log10(colSums(model$umitab[,cells])),model$cell_to_cluster[cells])[clusters],las=2,ylab="log10(#UMIs)")
  })
  
  output$BatchHeatmap <- renderPlot({
    clusters=clusters_reactive()
    cells=names(model$cell_to_cluster)[model$cell_to_cluster%in%clusters]
    if(length(unique(model$cell_to_batch[cells]))==1){
      return()
    }
    obs=table(model$cell_to_cluster[cells],model$cell_to_batch[cells])[clusters,]
    cs=colSums(obs)
    rs=rowSums(obs)
    expected=(matrix(cs,nrow(obs),ncol(obs),byrow=T)/sum(cs))*rs
    reg=5
    logr=log2((reg+obs)/(reg+expected))
    par(mar=c(10,7,4,2))
    
   # sort(apply(obs,1,function(x){chisq.test(x,simulate.p.value = 1000)$p.value}))
    image(logr[,ncol(logr):1],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=c(-100,seq(-2,2,l=99),100),main="Batch enrichment")
    mtext(text = clustAnnots[rownames(logr)],side = 1,at = seq(0,1,l=dim(logr)[1]),las= 2,cex=1)
    # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =colnames(logr), side=2,  at=seq(1,0,l=dim(logr)[2]),las=2,cex=1)
    #print(obs)
    box()
  })
  
  output$avg_profile_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("kellisogram", width = "150%", height = he)
   })
  
  
  output$avg_module_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("kellisogram_modules", width = "150%", height = he)
  })
  
  output$external_profiles_plot <- renderUI({
    external_profiles<<-strsplit(input$inRefProfiles,",")[[1]]
     if (length(external_profiles)==0){
      return()
    }
    
    he=12*length(external_profiles)
    plotOutput("kellisogram_external_profiles", width = "150%", height = he)
  })
  
  output$projection_kellisogram_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("kellisogram_projection",width="200%",height=he)
  })
  
  output$projection_barplot_plot <- renderUI({
    he=max(c(500,12*length(clusters_reactive())),na.rm=T)
    plotOutput("projection_barplot",width="200%",height=he)
  })
  
                  
  
  # A reactive expression with the ggvis plot
  output$kellisogram <- renderPlot({
  ##### Don't delete!!
    input$inAnnotateCluster
    input$inModelVer
    input$inBirdClusterSetExclude
    input$inBirdClusterSetInclude
  #####
    
    inclusts=clusters_reactive()
    ingenes=genes_reactive()#genes_reactive()
    insamples=truth_samples_reactive()
     if (!loaded_flag){
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
      tab=table(dataset$cell_to_cluster,dataset$cell_to_sample)
      tab=t(t(tab)/colSums(tab))[,insamples]
      tab=(tab/rowSums(tab))[intersect(inclusts,rownames(tab)),,drop=F]
      barplot(t(tab[nrow(tab):1,insamples]),col=sample_cols[1:ncol(tab)],horiz =T,yaxs = "i",names.arg=NULL,main="Samples",axes=F)
    }
    par(mar=c(7,tab3_left_margin,1,9))
    
    mat<<-model$models[match(ingenes,rownames(model$models)),]
    mat1=mat[,inclusts,drop=F]
    if (input$inAbsOrRel=="Relative"){
      if (ncol(mat1)>1){
        mat_to_show=log2(1e-5+mat1/pmax(1e-5,rowMeans(mat1,na.rm=T)))
        break1=-1e6
        break2=1e6
      }
      else{
        return()
      }
    }
    else if (input$inAbsOrRel=="Absolute"){
      mat_to_show=log10(mat1)
      break1=0
      break2=1e-1
    }
    isolate({
      image(mat_to_show[,ncol(mat1):1],col=colgrad,breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main=paste("Model:",loaded_version))
    })
    box()
   
    mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1,col=gcol[toupper(rownames(mat1))])
    mtext(text =paste(" ",colnames(mat1)," (n=",ncells_per_cluster[inclusts]," ; ",round(100*ncells_per_cluster[inclusts]/sum(ncells_per_cluster[names(which(!is_cluster_excluded))]),digits=1),"% )",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =paste(clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
  })
  
  
  module_counts_reactive<-reactive({
    input$inGetModules
    if (!exists("modules")){
      return()
    }

    counts=0
    for (i in 1:length(dataset$counts)){
      counts=counts+dataset$counts[[i]]
    }
    return(t(sapply(modules,function(modi){colSums(counts[modi,])})/colSums(counts)))
  
    
  })

  output$kellisogram_modules <- renderPlot({
    ##### Don't delete!!
    input$inModelVer
    #####

    inclusts=clusters_reactive()
    inmodules=modules_reactive()
   
    if (!loaded_flag){
      return()
    }
    if (!exists("modules")){
      return()
    }
    modulemat<<-module_counts_reactive()
    if(is.null(modulemat)){
      return()
    }
    modulemat<<-modulemat[inmodules,inclusts]
    
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
  
    image(mat_to_show[,ncol(mat1):1],col=colgrad,breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main=loaded_version)
  
    box()
    
    mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1)
    mtext(text =paste(" ",colnames(mat1)," (n=",ncells_per_cluster[inclusts]," ; ",round(100*ncells_per_cluster[inclusts]/sum(ncells_per_cluster),digits=1),"% )",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =paste(clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
  })
  
  
  
  output$kellisogram_external_profiles<- renderPlot({
    input$inRefProfiles
    ingenes=genes_reactive()
    
    if (!exists("ext_prof")){
      read_reference_set()
    }
    if (!loaded_flag){
      return()
    }
    
    # Lables for axes
    if (length(ingenes)==0){
      
      return()
    }
   
    zlim=c(-5,5)
    par(mar=c(1,tab3_left_margin,1,9))
    ext_profmat<-ext_prof[match(toupper(ingenes),toupper(rownames(ext_prof))),]
    mat1=as.matrix(ext_profmat[,external_profiles,drop=F])
    break1=-1e3
    break2=1e3
    isolate({
      image(mat1[,ncol(mat1):1,drop=F],col=colorRampPalette(c("blue","white","red"))(100),breaks=c(break1,seq(zlim[1],zlim[2],l=99),break2),axes=F,main="Reference Profiles")
    })
    box()
    
 #   mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las=2,cex=1,col=gcol[toupper(rownames(mat1))])
  #  mtext(text =paste(" ",colnames(mat1)," (n=",ncells_per_cluster[inclusts]," ; ",round(100*ncells_per_cluster[inclusts]/sum(ncells_per_cluster),digits=1),"% )",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    mtext(text =external_profiles, side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
    
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
    inclusts=clusters_reactive()
    ingenes=genes_reactive()
    insamples=truth_samples_reactive()
    nclust=length(inclusts)
    if (!loaded_flag){
      return()
    }
    if (length(ingenes)==0){
      updateSelectInput(session,"inGeneSets",selected = rownames(geneList)[1])
      updateTextInput(session,"inGenes",value = init_genes(geneList[1,1]))
      
      return()
    }

    samps_cl=list() 
    ds=dataset$ds[[match(input$inTruthDownSamplingVersion,dataset$ds_numis)]]
    
    cells_per_sample=as.numeric(input$inTruthNcellsPerSample)
    genes=intersect(ingenes,rownames(ds))
    cells_selected=dataset$randomly_selected_cells[[match(input$inTruthDownSamplingVersion,dataset$ds_numis)]][[match(input$inTruthNcellsPerSample,params$nrandom_cells_per_sample_choices)]]
                                                                                                            
    cell_mask=dataset$cell_to_cluster[colnames(ds)]%in%inclusts & 
              dataset$cell_to_sample[colnames(ds)]%in%insamples &
              colnames(ds)%in%cells_selected
    ds=ds[genes,cell_mask]
    ds=ds[,order(match(dataset$cell_to_cluster[colnames(ds)],inclusts))]
    samps=dataset$cell_to_sample[colnames(ds)]
    ncells=rep(0,length(inclusts))
    names(ncells)=inclusts
    tmp_ncells=sapply(split(colnames(ds),dataset$cell_to_cluster[colnames(ds)]),length)
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
    mtext(text =rownames(pmat), side=1, at=seq(0,1,l=dim(pmat)[1]),las=2,cex=1,col=gcol[toupper(rownames(pmat))])
    a=cumsum(ncells)
    b=a-floor(ncells[inclusts[inclusts%in%names(ncells)]]/2)
    mtext(text =inclusts, side=2, at=1-(b/a[length(ncells)]),las=2,cex=1,adj=1)
    par(mar=c(10,0,1,1))
    image(t(as.matrix(match(rev(samps),insamples))),axes=F,breaks=0:length(insamples)+.5,col=sample_cols[1:length(insamples)])
  })
    

    varmean_reactive <- reactive({
      clust=strsplit(input$inQCClust," - ")[[1]][1]
      samples=truth_samples_reactive()
      if (!exists("dataset")){
        return(data.frame(m=0,varmean=0,gene=""))
      }
      ds=dataset$ds[[match(input$inQCDownSamplingVersion,dataset$ds_numis)]]
      cluster_mask=names(dataset$cell_to_cluster)[dataset$cell_to_cluster==clust]
      ds=ds[,intersect(colnames(ds),cluster_mask),drop=F]

      if (is.null(ds))(return(data.frame(m=0,varmean=0,gene=0)))
      if (ncol(ds)<2)(return(data.frame(m=0,varmean=0,gene=0)))
      
      m=rowMeans(ds)
      v=apply(ds,1,var)
      
      df=data.frame(m=log10(m),varmean=log2(v/m),gene=rownames(ds))
      mask2=df$m>input$mean[1]&df$m<input$mean[2]&df$varmean<input$varmean[2]&df$varmean>input$varmean[1]
      df[mask2,]
      #  
      # Apply filters
    })
    
    
    output$modelSampleLegend <- renderPlot({
    
      insamples=truth_samples_reactive()
      par(mar=c(0,0,0,0))
      plot.new()
      leg=samples_tab[match(insamples,samples_tab$index),"title"]
      if(length(leg)==0){
        return()
      }
      ncol=ceiling(length(insamples)/10)
      leg[is.na(leg)]=insamples[is.na(leg)]

      legend("topleft",pch=15,col=sample_cols[1:length(insamples)],legend=leg,cex=1,xpd=T,ncol=ncol)
    })
    
    
    output$truthSampleLegend <- renderPlot({
      insamples=truth_samples_reactive()
      par(mar=c(0,0,0,0))
      plot.new()
      leg=samples_tab[match(insamples,samples_tab$index),"title"]
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
      varmean_reactive %>%
        ggvis(x = ~m, y = ~varmean) %>%
        layer_points(size := 50, size.hover := 200,
                     fillOpacity := 0.2, fillOpacity.hover := 0.5, key := ~gene) %>%
        add_tooltip(gene_tooltip, "hover")%>% 
        add_tooltip(click_tooltip, "click")%>% 
        add_axis("x", title = "log10(mean)") %>%
        add_axis("y", title = "log2(var/mean)") %>%
        add_axis("x", orient = "top", ticks = 0, title = paste(loaded_version,clust),
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
      insamples=truth_samples_reactive()
      ds=ds_QC_reactive()
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
     # mtext(3,lab)
    #  mtext(4,lab)
      
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
      
      if (click_flag){
        updateTextInput(session,"inGene1",,x$gene)
      }
      else{
        updateTextInput(session,"inGene2",,x$gene)
      }
      click_flag<<-!click_flag
    } 
    
    output$n_movies <- renderText({ nrow(movies()) })
    
    output$gene2 <- renderText({ paste("_________",input$inGene2,sep="")})
    output$gene1 <- renderText({ input$inGene1})  
   
    output$table <- renderTable({
      if (!loaded_flag){
        return()
      }
      clust=strsplit(input$inQCClust," - ")[[1]][1]
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
    output$kellisogram_projection <- renderPlot({
      ##### Don't delete!!
      #####
      insamples=truth_samples_reactive()
      inclusts=clusters_reactive()
      ingenes=genes_reactive()
      
      
      samples1=strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]]
      samples2=strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]]

      if (!((samples1%in%insamples|all(sample_sets[[samples1]]%in%insamples))&&input$inProjectSampleGroup2%in%insamples|all(sample_sets[[samples2]]%in%insamples))){
         return()
      }
      
      # Lables for axes
      if (length(ingenes)==0){
        
        updateSelectInput(session,"inGeneSets",selected = rownames(geneList)[1])
        updateTextInput(session,"inGenes",value = init_genes(geneList[1,1]))
        
        return()
      }
      
  
      zlim=input$inModelColorScale
      par(mar=c(7,7,1,9))
      
      mat1<-model$models[match(ingenes,rownames(model$models)),inclusts]
   
      n=ncol(mat1)
      ncells_per_projected_cluster=rep(0,length(inclusts))
      names(ncells_per_projected_cluster)=inclusts
      t1=table(projected$cell_to_cluster)
      ncells_per_projected_cluster[names(t1)]=t1
      
      percents_ref=round(100*ncells_per_cluster[inclusts]/sum(ncells_per_cluster),digits=1)
      percents_projected=round(100*ncells_per_projected_cluster[inclusts]/sum(ncells_per_projected_cluster),digits=1)
      
      
      layout(matrix(1:2,1,2))
        if (input$inProjectPlotType=="Tile"){
        matc=cbind(mat1,mat1p)[,rep(1:n,each=2)+rep(c(0,n),n)]
        isolate({
          image(log2(1e-5+matc/pmax(1e-5,rowMeans(matc)))[,ncol(matc):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=paste(loaded_version,input$inProjectedDataset,sep="/"))
        })
        box()
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
     
        sapply(yv,function(y){arrows(xusr[1],y,xusr[2],y,length=0,col="gray")})
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=gcol[toupper(rownames(mat1))])
        mtext(text =paste(" ",colnames(mat1)," (",percents_ref,"% /",percents_projected,"% )",sep=""), side=4, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
      }
      else if (input$inProjectPlotType=="Side by Side"){
        avg=rowMeans(cbind(mat1,mat1p))
        image(log2(1e-5+mat1/pmax(1e-5,avg))[,ncol(mat1):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=loaded_version)
        box()
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=gcol[toupper(rownames(mat1))])
       # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(" (",percents_ref,"%)",sep=""), side=4, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(clustAnnots[inclusts]," ",sep=""), side=2, at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        image(log2(1e-5+mat1p/pmax(1e-5,avg))[,ncol(mat1p):1],col=colgrad,breaks=c(-1e6,seq(zlim[1],zlim[2],l=99),1e6),axes=F,main=input$inProjectedDataset)
        yusr=par("usr")[3:4]
        xusr=par("usr")[1:2]
        yv=seq(yusr[1],yusr[2],l=dim(mat1)[2]+1)
        mtext(text = rownames(mat1),side = 1,at = seq(0,1,l=dim(mat1)[1]),las= 2,cex=1,col=gcol[toupper(rownames(mat1))])
       # mtext(text =paste(" ",colnames(mat1),sep=""), side=2, at=seq(1-(yusr[2]-1),-1*yusr[1],l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(" (",percents_projected,"%)",sep=""), side=4,  at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        mtext(text =paste(clustAnnots[inclusts]," ",sep=""), side=2,  at=seq(1,0,l=dim(mat1)[2]),las=2,cex=1)
        
        box()
      }
    })
    
    output$projection_barplot=renderPlot({
      ##### Don't delete!!
      #####
      insamples=truth_samples_reactive()
      inclusters=clusters_reactive()
      ingenes=genes_reactive()
      
      
      samples1=as.list(strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]])
      samples2=as.list(strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]])
      if (exists("sample_sets")){
        samples1[!(samples1%in%insamples)]=sample_sets[unlist(samples1[!samples1%in%insamples])]
        samples2[!(samples2%in%insamples)]=sample_sets[unlist(samples2[!samples2%in%insamples])]
      }
      samples1=unlist(samples1)
      samples2=unlist(samples2)
      
      if (is.null(samples1)||is.null(samples2)){
        return()
      }
      if (!(all(samples1%in%insamples)&&all(samples2%in%insamples))){
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
      barplot(m[,dim(m)[2]:1],beside=T,names.arg = rev(paste(clustAnnots[inclusters]," - ",inclusters,sep="")),horiz = T,las=2)
      legend("bottomright",pch=15,col= gray.colors(2),legend=c("1","2"),border=T)
    #  reg=100*10/((ncol(model$umitab)+ncol(projected$umitab))/2)
      reg=1e-3
      barplot(rev(log2((reg+m[2,])/(reg+m[1,]))),names.arg = rev(paste(clustAnnots[inclusters]," - ",inclusters,sep="")),horiz = T,las=2,xlab="log2(1/2)")
    })
    
    click_tooltip_gene_proj_vs_ref <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$gene)) return(NULL)
      
     
      updateTextInput(session,"Gene1ForExprsTableRefVsProj",,x$gene)
      
    } 
    
    ge_proj_vs_ref <- reactive({
      insamples=truth_samples_reactive()
      clust=strsplit(input$inClustForDiffGeneExprsProjVsRef," - ")[[1]][1]
      samples1=as.list(strsplit(input$inProjectSampleGroup1,",| ,|, ")[[1]])
      samples2=as.list(strsplit(input$inProjectSampleGroup2,",| ,|, ")[[1]])
      if (!exists("dataset")){
        return(data.frame(x=0,y=0,gene=0))
      }
      if (exists("sample_sets")){
        samples1[!(samples1%in%insamples)]=sample_sets[unlist(samples1[!samples1%in%insamples])]
        samples2[!(samples2%in%insamples)]=sample_sets[unlist(samples2[!samples2%in%insamples])]
      }
      samples1=unlist(samples1)
      samples2=unlist(samples2)
      
      if(length(samples1)==0||length(samples2)==0||(!(all(samples1%in%insamples)&&all(samples2%in%insamples)))){
        return(data.frame(x=0,y=0,gene=0))
      }
      
      counts1=model$models*0
      counts2=model$models*0
      for (sampi in 1:length(samples1)){
          tmp_counts=dataset$counts[[samples1[sampi]]]
         counts1[,colnames(tmp_counts)]=counts1[,colnames(tmp_counts)]+tmp_counts  
      }
      
      for (sampi in 1:length(samples2)){
        tmp_counts=dataset$counts[[samples2[sampi]]]
        counts2[,colnames(tmp_counts)]=counts2[,colnames(tmp_counts)]+tmp_counts  
      }
      
      model1=t(t(counts1)/colSums(counts1))
      model2=t(t(counts2)/colSums(counts2))
      
      df=data.frame(x=log10(model1[,clust]),y=log2(model2[,clust]/model1[,clust]),gene=rownames(model$models))
      
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
        add_axis("x", orient = "top", ticks = 0, title = paste(loaded_version," ",clust),
                 properties = axis_props(
                   axis = list(stroke = "white"),
                   labels = list(fontSize = 0))) %>%
        scale_numeric("x", domain = c(-6,-1)) %>%
        scale_numeric("y", domain = c(-6,6)) %>%
       set_options(width = 500, height = 500)
      
      
    })
    
    vis2 %>% bind_shiny("DiffGeneExprsProjVsRef")

    
    output$Gene1ExprsTableRefVsProj<- renderTable({
      if (!loaded_flag){
        return()
      }
      if (!exists("projected")){
        return()
      }
   
      g=input$Gene1ForExprsTableRefVsProj
      clust=strsplit(input$inClustForDiffGeneExprsProjVsRef," - ")[[1]][1]
      
      if ((g%in%rownames(ds))&(g%in%rownames(projected$ds))){
        mask1=model$cell_to_cluster[colnames(ds)]==clust
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
      if (!loaded_flag){
        return()
      }
      if (!all(names(model$cell_to_cluster)==names(projected$source_cell_to_cluster))){
        return() 
      }
      tab=table(model$cell_to_cluster,projected$source_cell_to_cluster)
      x=matrix(tab,nrow(tab),ncol(tab),dimnames=list(rownames(tab),colnames(tab)))
      print(x)
      x
      
    })
    
    
    click_tooltip_ll <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$cell)) return(NULL)
      
      
   #   updateTextInput(session,"Gene1ForExprsTableRefVsProj",,x$gene)
      
    } 
    
    cell_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x$cell)) return(NULL)
      
      # Pick out the movie with this ID
      
      paste0("<b>", x$cell, "</b>")
    }
    
#    ll_plot <- reactive({
#      clust=strsplit(input$inTweezersFromCluster," - ")[[1]][1]
#      llx=strsplit(input$inTweezersLLX," - ")[[1]][1]
#      lly=strsplit(input$inTweezersLLY," - ")[[1]][1]
   
#      print(llx)
#      print(lly)
#      if (is.na(clust)|!loaded_flag){return(data.frame(x=0,y=0,gene=0))}
      
#      print(llx)
#      print(lly)
#      df=data.frame(x=ll[model$cell_to_cluster==cluster,llx],y=ll[model$cell_to_cluster==cluster,lly],cell=rownames(ll)[model$cell_to_cluster==cluster])
      
#    })
    # A reactive expression with the ggvis plot
#    vis3 <- reactive({
      
#      clust=strsplit(input$inTweezersFromCluster," - ")[[1]][1]
      
#      ll_plot %>%
#        ggvis(x = ~x, y = ~y) %>%
#        layer_points(size := 50, size.hover := 200,
#                     fillOpacity := 0.2, fillOpacity.hover := 0.5, key := ~cell) %>%
#        add_tooltip(cell_tooltip, "hover")%>% 
#        add_tooltip(click_tooltip_ll, "click")%>% 
#        add_axis("x", title = "llx") %>%
#        add_axis("y", title = "lly") %>%
#        add_axis("x", orient = "top", ticks = 0, title = paste(loaded_version," ",clust),
#                 properties = axis_props(
#                   axis = list(stroke = "white"),
#                   labels = list(fontSize = 0))) %>%
#        scale_numeric("x", domain = c(-6,-1)) %>%
#        scale_numeric("y", domain = c(-6,6)) %>%
#        set_options(width = 500, height = 500)
#      
#      vis3 %>% bind_shiny("TweezersLikelihoodPlot")
#      
#    })
    
    
    if (exists("scDissector_datadir")){
      updateTextInput(session,"inDatapath",,scDissector_datadir)
    }
    
  }
#table(projected$source_cell_to_cluster,projected$cell_to_cluster)

#import_alternative_clustering=function(cell_to_cluster,fn){
#  umitab=umitab[,names(cell_to_cluster)]
#  model$models=update_model$models(umitab,cell_to_cluster)
#  l_cells_per_cluster=NULL
#  cell_to_batch=NULL
#  ds=downsample(umitab,params$min_umis)
#  save(umitab,ds,model$models,cell_to_cluster,ll,l_cells_per_cluster,cell_to_batch,chisq_res,params,file=fn)
#}





