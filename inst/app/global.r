params=list()
params$nrandom_cells_per_sample_choices=c(100,250,500,1000,2000,4000,"All")

tfs_file="tfs.csv"
surface_markers_file="surface_markers.csv"

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


