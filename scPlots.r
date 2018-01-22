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