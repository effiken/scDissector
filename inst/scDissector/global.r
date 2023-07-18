params=list()
params$nrandom_cells_per_sample_choices=c(100,250,500,1000,2000,4000,"All")
source("ldm_interface.R")
tfs_file=system.file("extdata", "tfs.csv", package="scDissector")
if (tfs_file==""){
    tfs_file="../extdata/tfs.csv"
}
#surface_markers_file=system.file("extdata", "surface_markers.csv", package="scDissector")

surface_markers_file=system.file("extdata", "Martinez_Martin_cell_2018.csv", package="scDissector")
if (surface_markers_file==""){
    surface_markers_file="../extdata/surface_markers.csv"
}

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
  surface_markers<-unique(read.csv(surface_markers_file,skip=2,header=F,stringsAsFactors = F)[,1])
  surface_markers2=paste(substring(surface_markers, 1,1), tolower(substring(surface_markers, 2)),sep="")
  surface_markers<-c(surface_markers,surface_markers2)
  rm(surface_markers2)
  rm(surface_markers_file)
}else{
  message(surface_markers_file ," not found")
}

colgrad_abs_file=system.file("extdata", "colors_viridis.txt", package="scDissector")
if (colgrad_abs_file==""){
    colgrad_abs_file="../extdata/colors_viridis.txt"
}
colgrad_abs<<-read.table(colgrad_abs_file,stringsAsFactors=F)[,1]

colgrad_rel_file=system.file("extdata", "colors_brewer_RdBu.txt", package="scDissector")
if (colgrad_rel_file==""){
    colgrad_rel_file="../extdata/colors_brewer_RdBu.txt"
}
colgrad_rel<<-read.table(colgrad_rel_file,stringsAsFactors=F)[,1]

colgrad_file=system.file("extdata", "colors_paul.txt", package="scDissector")
if (colgrad_file==""){
    colgrad_file="../extdata/colors_paul.txt"
}
colgrad<<-read.table(colgrad_file,stringsAsFactors=F)[,1]
#colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))

default_sample_colors_file=system.file("extdata", "sample_colors.txt", package="scDissector")
if (default_sample_colors_file==""){
    default_sample_colors_file="../extdata/sample_colors.txt"
}
default_sample_colors<<-rep(paste("#",read.table(default_sample_colors_file,stringsAsFactors = F)[,1],sep=""),10)


gsc=load_gene_symbol_converters()
gene_symbol_old2new<<-gsc$old2new
gene_symbol_new2old<<-gsc$new2old
rm(gsc)
genesetsfile<<-system.file("extdata", "gene_sets.txt", package="scDissector")
if (genesetsfile==""){
    genesetsfile="../extdata/gene_sets.txt"
}
geneList_tmp<-read.table(file=genesetsfile,header=T,stringsAsFactors = F,row.names =1)
geneList<-geneList_tmp[,1]
names(geneList)<-rownames(geneList_tmp)
rm(geneList_tmp)

aggregate.Matrix=function(mat,by_v){
  return(fac2sparse(by_v) %*% mat)
}

