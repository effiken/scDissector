plot_avg_heatmap=function(m,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute="Relative"){
  if (Relative_or_Absolute=="Relative"){
    if (ncol(m)>1){
      m=log2(1e-6+m/pmax(1e-6,rowMeans(m,na.rm=T)))
    }
    else{
      return()
    }
  }
  else if (Relative_or_Absolute=="Absolute"){
    m=log10(1e-6+m)
  }
  else{
    plot.new()
    return()
  }
  break1=min(m,na.rm=T)
  break2=max(m,na.rm=T)

  breaks=sort(c(break1,seq(zlim[1],zlim[2],l=99),break2))
  
  image(m[,ncol(m):1],col=colgrad,breaks=breaks,axes=F,main=main_title)
  box()
  mtext(text = genes,side = 1,at = seq(0,1,l=length(genes)),las=2,cex=1,col=gene.cols)
  mtext(text =paste(" ",clusters,clusters_text), side=4, at=seq(1,0,l=length(clusters)),las=2,cex=1)
  mtext(text =paste(annots," ",sep=""), side=2, at=seq(1,0,l=length(annots)),las=2,cex=1)
}

plot_avg_heatmap_interactive=function(m,zlim,main_title,genes,gene.cols,clusters,clusters_text,annots,Relative_or_Absolute="Relative"){
  if (Relative_or_Absolute=="Relative"){
    if (ncol(m)>1){
      m=log2(1e-6+m/pmax(1e-6,rowMeans(m,na.rm=T)))
    }
    else{
      return()
    }
  }
  else if (Relative_or_Absolute=="Absolute"){
    m=log10(1e-6+m)
  }
  else{
    plot.new()
    return()
  }
  break1=min(m,na.rm=T)
  break2=max(m,na.rm=T)
  
 # breaks=sort(c(break1,seq(zlim[1],zlim[2],l=99),break2))
  m[m<zlim[1]]=zlim[1]
  m[m>zlim[2]]=zlim[2]
 # image(m[,ncol(m):1],col=colgrad,breaks=breaks,axes=F,main=main_title)
#  box()
#  mtext(text = genes,side = 1,at = seq(0,1,l=length(genes)),las=2,cex=1,col=gene.cols)
#  mtext(text =paste(" ",clusters,clusters_text), side=4, at=seq(1,0,l=length(clusters)),las=2,cex=1)
#  mtext(text =paste(annots," ",sep=""), side=2, at=seq(1,0,l=length(annots)),las=2,cex=1)
  
  heatmaply(t(m),Rowv=F,Colv=F, scale = "none", colors = colgrad,hide_colorbar=T,label_names=c("Gene","Cluster","Value"),labRow=rev(paste(annots," (",colnames(m),")",sep="")),labCol=rownames(m),main_title=main_title,column_text_angle=90,margins =c(100,180,20,0),fontsize_row = 8,fontsize_col = 8)
  
}


plot_truth_heatmap=function(ldm,downSamplingVersion,ingenes,inclusts_list,insamples,cells_per_sample=1000,zlim=c(0,2),ShowSeparatorBars=T,cex.genes=1){
 

  samps_cl=list() 
  ds=ldm$dataset$ds[[match(downSamplingVersion,ldm$dataset$ds_numis)]]
  
 
  genes=intersect(ingenes,rownames(ds))

  cells_selected=ldm$dataset$randomly_selected_cells[[match(downSamplingVersion,ldm$dataset$ds_numis)]][[as.character(cells_per_sample)]]
  
  layout(matrix(1:(2*length(inclusts_list)),length(inclusts_list),2,byrow=T),widths=c(40,1))   
  
  for (i in 1:length(inclusts_list)){
   message(i)
   
    left_margins=5
    if (i==length(inclusts_list)){
      top_margin=.1
      bottom_margin=5
     
      show_gene_names=T
    }
    else{
      top_margin=.1
      bottom_margin=.1
      show_gene_names=F
    }
    
    inclusts=inclusts_list[[i]] 
    cell_mask=ldm$dataset$cell_to_cluster[colnames(ds)]%in%inclusts & 
      ldm$dataset$cell_to_sample[colnames(ds)]%in%insamples &colnames(ds)%in%cells_selected
      ds_i=ds[genes,cell_mask]
  
    ds_i=ds_i[,order(match(ldm$dataset$cell_to_cluster[colnames(ds_i)],inclusts))]
    samps=ldm$dataset$cell_to_sample[colnames(ds_i)]
    ncells=rep(0,length(inclusts))
    names(ncells)=inclusts
    tmp_ncells=sapply(split(colnames(ds_i),ldm$dataset$cell_to_cluster[colnames(ds_i)]),length)
    ncells[names(tmp_ncells)]=tmp_ncells
  
  #ncells=sapply(ds_cl,function(x){n=ncol(x);if(is.null(n)){return(0)};return(n)})
  # names(ncells)=names(ds_cl)
  
  
    pmat=as.matrix(ds_i)[,ncol(ds_i):1]
    spacer_size=ceiling(dim(pmat)[2]/200)
    pmat2=log2(1+pmat)  
  
    par(mar=c(bottom_margin,left_margins,top_margin,1))
    # image(pmat2[,ncol(pmat2):1],col=c("gray",colgrad),axes=F,breaks=c(-3e6,-1e6,seq(zlim[1],zlim[2],l=99),1e6))
    message(names(inclusts_list)[i])
    image(pmat2,col=c("gray",colgrad),axes=F,breaks=c(-3e6,-1e6,seq(zlim[1],zlim[2],l=99),1e6),ylab=names(inclusts_list)[i])
  
  
    box()
    if (ShowSeparatorBars){
      abline(h=1-cumsum(ncells)/sum(ncells),col="gray")
    }
    if (show_gene_names){
      mtext(text =rownames(pmat), side=1, at=seq(0,1,l=dim(pmat)[1]),las=2,cex=cex.genes,col=ldm$gcol[toupper(rownames(pmat))])
    }
    a=cumsum(ncells)
    b=a-floor(ncells[as.character(inclusts[inclusts%in%names(ncells)])]/2)
    mtext(text =inclusts, side=2, at=1-(b/a[length(ncells)]),las=2,cex=.7,adj=1)
    par(mar=c(bottom_margin,0,top_margin,1))
    image(t(as.matrix(match(rev(samps),insamples))),axes=F,breaks=0:length(insamples)+.5,col=sample_cols[1:length(insamples)])
  
  }
  
  
}
