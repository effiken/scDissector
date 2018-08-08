library(seriation)
fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
fisher.z2r <- function (z) {(exp(2 * z) - 1)/(1 + exp(2 * z))}


get_avg_gene_to_gene_cor=function(ds,cell_to_sample){
  zmat=matrix(0,nrow(ds),nrow(ds),dimnames=list(rownames(ds),rownames(ds)))
  samples=unique(cell_to_sample)
  for (samp in samples){
    print(samp)
    dsi=ds[Matrix::rowSums(ds[,cell_to_sample==samp])>5,cell_to_sample==samp]
    z=fisher.r2z(.99*cor(as.matrix(Matrix::t(dsi))))
    #    z=fisher.r2z(.99*sparse.cor(Matrix::t(dsi)))
    rm(dsi)
    z[is.na(z)]=0
    #   print(range(z))
    zmat[rownames(z),colnames(z)]=z+ zmat[rownames(z),colnames(z)]
    rm(z)
    gc()
  }
  return(fisher.z2r(zmat/length(samples)))
}

get_gene_cormap=function(ldm,ds_version="2000",cells=NA){
  
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  if (!all(is.na(cells))){
    ds=ds[,intersect(colnames(ds),cells)]
  }
  cormat=cor_analysis(ds,ldm)
  
  return(cormat)
}


save_gene_cor_map=function(cormat,modules_version,zbreaks=c(-1,seq(-.5,.5,l=99),1),ser_method="OLD_ward",cor_cols=colorRampPalette(c("blue","white","red"))(100)){
  pdf(paste(modules_version,"gene_cor.png",sep="_"),width=ncol(cormat)/12,height=ncol(cormat)/12)
  #  ord=hclust(as.dist(1-cormat))$order
  #  browser()
  par(mar=c(3,3,3,3))
  ord=get_order(seriate(as.dist(1-cormat),ser_method))
  image(cormat[ord,ord],col=cor_cols,breaks=zbreaks,axes=F)
  mtext(text = colnames(cormat)[ord],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  box()
  dev.off()
}


gene_cor_analysis=function(ldm,ds_version,min_varmean_per_gene=0.15,min_number_of_UMIs=50){
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  s1=rowSums(ds,na.rm=T) 
  mask1=s1>=min_number_of_UMIs&(!rownames(ds)%in%c(grep("RP",rownames(ds),val=T),grep("MT-",rownames(ds),val=T)))
  message(sum(mask1)," genes passed expression threshold")
  v1=apply(ds[mask1,],1,var,na.rm=T)
  m1=rowMeans(ds[mask1,],na.rm=T)
  x=log10(m1)
  breaks=seq(min(x,na.rm=T),max(x,na.rm=T),.2)
  lv=log2(v1/m1)
  llv=split(lv,cut(x,breaks))
  mask_llv=sapply(llv,length)>0   
  z=sapply(llv[mask_llv],min,na.rm=T)
  b=breaks[-length(breaks)]
  b=b[mask_llv]
  lo=loess(z~b)
  high_var_genes=names(which((lv-predict(lo,newdata =x))>min_varmean_per_gene))
  message(length(high_var_genes)," High var genes")
  cormat=get_avg_gene_to_gene_cor(log2(1+ds[high_var_genes,]),ldm$dataset$cell_to_sample[colnames(ds)]) 
}


parse_modules=function(ldm,cormat,ds_version,modules_version="",nmods=50,reg=1e-6,mod_size=4,min_mod_cor=0.1,zlim=c(-.9,.9),ord_viz_method="OLO_complete"){
  
  gene_mask2=names(which(apply(cormat,1,quantile,1-mod_size/ncol(cormat),na.rm=T)>=min_mod_cor))
  #gene_cor_map(cormat[gene_mask2,gene_mask2],modules_version=modules_version,ser_method="OLO_complete",zbreaks=c(-1,seq(zlim[1],zlim[2],l=99),1))
  ds=ldm$dataset$ds[[match(ds_version,ldm$dataset$ds_numis)]]
  mods=cutree(hclust(as.dist(1-cormat[gene_mask2,gene_mask2])),nmods)
  modsl=split(names(mods),mods)
  modsums=aggregate(ds[gene_mask2,],groupings=mods,fun="sum")
#  modclustsum=aggregate(t(modsums),groupings=ldm$dataset$cell_to_cluster[colnames(modsums)],fun="sum")
#  mod_freqs=t(modclustsum/rowSums(modclustsum))
#  mod_freqs_normed=log2((reg+mod_freqs)/(reg+rowMeans(mod_freqs)))
 # ord=get_order(seriate(as.dist(1-cor(t(as.matrix(mod_freqs_normed)))),method = "OLO"))

  pdf(paste("module_cor_",modules_version,".pdf",sep=""),width=nmods/10,height=nmods/10)
  cormat=cor(t(as.matrix(modsums)))
  zbreaks=c(-1,seq(zlim[1],zlim[2],l=99),1)
  cor_cols=colorRampPalette(c("blue","white","red"))(100)
  par(mar=c(3,3,3,3))
  ord=get_order(seriate(as.dist(1-cormat),method="OLO_complete"))
  cormat=cormat[ord,ord]
  colnames(cormat)=1:nmods
  rownames(cormat)=1:nmods
  modsl=modsl[ord]
  names(modsl)=1:nmods
  ord_viz=get_order(seriate(as.dist(1-cormat),method=ord_viz_method))
  image(cormat[ord_viz,ord_viz],col=cor_cols,breaks=zbreaks,axes=F)
  mtext(text = colnames(cormat)[ord_viz],side = 1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 3,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  mtext(text = colnames(cormat)[ord_viz],side = 4,at = seq(0,1,l=ncol(cormat)),las=2,cex=.5)
  box()
  dev.off()
  
#  mod_freqs=mod_freqs[ord,]
#  rownames(mod_freqs)=1:nmods
 # mod_freqs_normed=mod_freqs_normed[ord,]
#  rownames(mod_freqs_normed)=1:nmods
  
  write.table(file=paste("tables/",modules_version,"_modules.txt",sep=""),sapply(modsl,paste,collapse=","),row.names = T,col.names = F,quote=F)
#  write.csv(file=paste("tables/",modules_version,"_modules_log10_freq_per_cluster.csv",sep=""),log10(1e-10+as.matrix(mod_freqs)),row.names = T,quote=F)
#  open_plot(fn=paste(modules_version,"_modules_clusters_heatmap",sep=""),plot_type = "pdf",width = 10,height = 5)
#  par(mar=c(5,7,1,1))
#  image(as.matrix(mod_freqs_normed[,intersect(ldm$cluster_order,colnames(mod_freqs_normed))]),col=greenred(100),axes=F,breaks=c(-1e8,seq(-5,5,l=99),1e8))
#  mtext(text = ldm$clustAnnots[intersect(ldm$cluster_order,colnames(mod_freqs_normed))],side = 2,at = seq(0,1,l=ncol(mod_freqs_normed)),las=2,cex=.5)
#  mtext(text = paste("Module",rownames(mod_freqs_normed)),side = 1,at = seq(0,1,l=nrow(mod_freqs_normed)),las=2,cex=.5)
#  close_plot()
  return(modsl)
}


example=function(){
  cormat=gene_cor_analysis(ldm,"2000",0.1,50)
  save_gene_cor_map(cormat,"all_cells2",ser_method = "OLO_complete")
  parse_modules(ldm,cormat,"2000","all_cells2",nmods =50)
}
