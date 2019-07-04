library(dplyr)
library(stringr)
library(seqinr)

features_derived <- function(infile){

  csvdf <- read.csv(infile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  

  # Apply the function GC on the matrix of the data, which only contains the sequences  
  length <- apply(as.matrix(csvdf[, "sequence_nodm"]),1,nchar)
  freqs <- data.frame()
  for(i in csvdf[, "sequence_nodm"]){
    freqs <- rbind(freqs,data.frame(t(as.matrix(count(tolower(s2c(i)), wordsize = 2, freq = T, alphabet = s2c("acgu"))))))
  }
  NEFE <- csvdf[,"efe"]/length
  dG <-  csvdf[,"mfe"]/length
  n_stems <- apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,'[^\\(\\)]{3,}')})
  tot_bases <-  apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,'\\.')})
  loopregex='\\.{4,}'
  n_loops <- apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,loopregex)})
  gccontent <- apply(as.matrix(csvdf[, "sequence_nodm"]),1,function(rnaseq){GC(s2c(rnaseq))})
  mfeI1 <- dG/gccontent
  mfeI2 <- dG/n_stems
  mfeI3 <- dG/n_loops
  mfeI4 <- csvdf[,"mfe"]/tot_bases
  dP <- tot_bases/length
  diff <- dG-NEFE
  
  #
  # Custom
  #
  nefe_div_gccon=NEFE/gccontent
  nefe_times_dP=NEFE*dP
  nefe_times_dP_div_gccon=NEFE*dP/gccontent
  ua_div_gccon=freqs[["ua"]]/gccontent
  
  total <- cbind(comment=csvdf[,"comment"], length=length, nefe=NEFE,mfei1=mfeI1,mfei2=mfeI2,mfei3=mfeI3,mfei4=mfeI4,tot_bases,n_stems,n_loops,dP,dG,diff,gccontent,freqs,nefe_div_gccon,nefe_times_dP,nefe_times_dP_div_gccon,ua_div_gccon)
  total[is.na(total)] <- 0
  return(total)
}

if(exists("snakemake")){
	total<-features_derived(snakemake@input[[1]])
	write.csv(total, file=snakemake@output[[1]], row.names=FALSE)
}else{
	total<-features_derived("test.csv")
	write.csv(total,file="testout.csv",row.names=FALSE)
}
