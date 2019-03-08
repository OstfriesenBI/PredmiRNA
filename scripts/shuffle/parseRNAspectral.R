suppressWarnings(library(seqinr))

parseRNAspectral <- function(filename,fastasource){
  i <- read.delim(filename)
  seqs <- read.fasta(fastasource)
  commentvals<-names(seqs)
  i$ID <- NULL
  i$MFE <- NULL
  i$Len <- NULL
  colnames(i) <- paste("spectral", colnames(i), sep = "_")
  i$comment <- commentvals
  return(i)
}

if(exists("snakemake")){
  write.csv(parseRNAspectral(snakemake@input[[1]],snakemake@input[[2]]),file = snakemake@output[[1]],row.names = FALSE)
}else{
  data <- parseRNAspectral("data/real_izmir/stanley/0000001.spectral","data/real_izmir/split/real_izmir.fasta_chunk_0000001")
  print(head(data))
}
