params=list()
params$nrandom_cells_per_sample_choices=c(100,250,500,1000,2000,4000,"All")

tfs_file=system.file("extdata", "tfs.csv", package="scDissector")
surface_markers_file=system.file("extdata", "surface_markers.csv", package="scDissector")

if (file.exists(tfs_file)){
  tfs<-read.csv(file=tfs_file,header=F,stringsAsFactors = F)[,1]
  tfs2=paste(substring(tfs, 1,1), tolower(substring(tfs, 2)),sep="")
  tfs<-c(tfs,tfs2)
  rm(tfs2)
  rm(tfs_file)
}else{
  message(tfs_file ," not found")
}

if (file.exists(surface_markers_file)){
  surface_markers<-read.csv(file=surface_markers_file,header=F,stringsAsFactors = F)[,1]
  surface_markers2=paste(substring(surface_markers, 1,1), tolower(substring(surface_markers, 2)),sep="")
  surface_markers<-c(surface_markers,surface_markers2)
  rm(surface_markers2)
  rm(surface_markers_file)
}else{
  message(surface_markers_file ," not found")
}

colgrad_abs<<-read.table(system.file("extdata", "colors_viridis.txt", package="scDissector"),stringsAsFactors=F)[,1]
colgrad_rel<<-read.table(system.file("extdata", "colors_brewer_RdBu.txt", package="scDissector"),stringsAsFactors=F)[,1]

colgrad<<-read.table(system.file("extdata", "colors_paul.txt", package="scDissector"),stringsAsFactors=F)[,1]
#colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))

default_sample_colors<<-rep(paste("#",read.table(system.file("extdata", "sample_colors.txt", package="scDissector"),stringsAsFactors = F)[,1],sep=""),10)

### gene symbol converters
load_gene_symbol_converters=function(){
    hgnc<-read.delim(system.file("extdata", "hgnc_complete_set.txt", package="scDissector"),header = T,stringsAsFactors = F)
    old_symbol=ifelse(hgnc[,"prev_symbol"]=="",hgnc[,"symbol"],hgnc[,"prev_symbol"])
    l_old_symbol=strsplit(old_symbol,"\\|")
    old_symbol2=unlist(l_old_symbol)

    new_symbol=hgnc[,"symbol"]
    new_symbol2=rep(new_symbol,sapply(l_old_symbol,length))

    gene_symbol_old2new<-(new_symbol2)
    names(gene_symbol_old2new)<-old_symbol2

    gene_symbol_new2old<-old_symbol
    names(gene_symbol_new2old)<-new_symbol
    return(list(old2new=gene_symbol_old2new,new2old=gene_symbol_new2old))
}

gsc=load_gene_symbol_converters()
gene_symbol_old2new<<-gsc$old2new
gene_symbol_new2old<<-gsc$new2old
rm(gsc)
genesetsfile<<-system.file("extdata", "gene_sets.txt", package="scDissector")
geneList_tmp<-read.table(file=genesetsfile,header=T,stringsAsFactors = F,row.names =1)
geneList<-geneList_tmp[,1]
names(geneList)<-rownames(geneList_tmp)
rm(geneList_tmp)

