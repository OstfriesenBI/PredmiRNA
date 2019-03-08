parseRNARandom <- function(filename,prefix){
  i <- read.delim(filename)[,-1];
  commentvals <- i$comment
  i$comment <- NULL
  colnames(i) <- paste(prefix, colnames(i), sep = "_")
  i$comment <- commentvals
  return(i)
}

if(exists("snakemake")){
  prefix <- paste("permuted",snakemake@wildcards["method"],snakemake@wildcards["shuffles"],sep = "_")
  write.csv(parseRNARandom(snakemake@input[[1]],prefix),file = snakemake@output[[1]],row.names = FALSE)
}else{
  data <- parseRNARandom("data/real_izmir/stanley/0000001-m-10.stats","test")
  View(head(data))
}
