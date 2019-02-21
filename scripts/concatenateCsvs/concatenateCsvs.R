#!/usr/bin/env Rscript
concatenateCsvs = function(filenames){ 
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {rbind(x,y)}, datalist)
}

if(exists("snakemake")){
  write.csv(concatenateCsvs(snakemake@input[["csvs"]]),snakemake@output[["csv"]],row.names=FALSE)
}else{
  #concatenateCsvs()
}

