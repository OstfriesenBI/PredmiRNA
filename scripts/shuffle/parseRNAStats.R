parseRNAStats <- function(filename){
  inputdata <- read.delim(filename);
  data.frame(comment=inputdata$comment,Q=inputdata$Q,normQ=inputdata$NQ,D=inputdata$D,normD=inputdata$ND)
}

if(exists("snakemake")){
  write.csv(parseRNAStats(snakemake@input[[1]]),file = snakemake@output[[1]],row.names = FALSE)
}else{
  data <- parseRNAStats("data/real_izmir/stanley/0000001.stats")
  print(head(data))
}
