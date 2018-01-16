get_gene_symbol_convetors=function(path=""){
  hgnc<-read.delim(paste(path,"hgnc_complete_set.txt",sep=""),header = T,stringsAsFactors = F)
  old_symbol=ifelse(hgnc[,"prev_symbol"]=="",hgnc[,"symbol"],hgnc[,"prev_symbol"])
  l_old_symbol=strsplit(old_symbol,"\\|")
  old_symbol2=unlist(l_old_symbol)

  new_symbol=hgnc[,"symbol"]
  new_symbol2=rep(new_symbol,sapply(l_old_symbol,length))

  gene_symbol_old2new<-new_symbol2
  names(gene_symbol_old2new)<-old_symbol2

  gene_symbol_new2old<-old_symbol
  names(gene_symbol_new2old)<-new_symbol
  return(list(old2new=gene_symbol_old2new,new2old=gene_symbol_new2old))
}

cap <- function(x) {
  paste(toupper(substring(x, 1,1)), tolower(substring(x, 2)),sep="", collapse=" ")
}

adjust_gene_names=function(genes_string,data_genes){
  genes=strsplit(genes_string,",|, | ,")[[1]]
  mask1=toupper(genes)%in%data_genes
  genes[mask1]=toupper(genes[mask1])
  mask2=unlist(sapply(genes,cap))%in%data_genes
  genes[mask2]=unlist(sapply(genes[mask2],cap))
  mask3=gene_symbol_old2new[genes]%in%data_genes
  genes[mask3]=gene_symbol_old2new[genes[mask3]]
  mask4=gene_symbol_new2old[genes]%in%data_genes
  genes[mask4]=gene_symbol_new2old[genes[mask4]]
 
  return(genes)

}