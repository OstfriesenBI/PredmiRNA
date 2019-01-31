library(stringr)

features_derived <- function(infile,out){

  csvdf <- read.csv(infile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  

  # Apply the function GC on the matrix of the data, which only contains the sequences  
  length <- apply(as.matrix(csvdf[, "sequence"]),1,nchar)
  NEFE <- csvdf[,"efe"]/length
  dG <-  csvdf[,"mfe"]/length
  n_stems <- apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,'(\\(){3,}')})
  tot_bases <-  apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,'(\\()')})
  n_loops <- apply(as.matrix(csvdf[,"mfesecstructure"]),1, function(seq){str_count(seq,'(\\.){3,}')})
  # GC calc here
  #mfeI1 <- dG/gccontent
  mfeI2 <- dG/n_stems
  mfeI3 <- dG/n_loops
  mfeI4 <- csvdf[,"mfe"]/tot_bases
  dP <- tot_bases/length
  diff <- dG-NEFE
  total <- cbind(comment=csvdf[,"comment"], length=length, nefe=NEFE,mfei2=mfeI2,mfei3=mfeI3,mfei4=mfeI4,tot_bases,n_stems,n_loops,dP,dG,diff)
  # Writes Data Frame to a .csv file
  write.csv(total, file=out, row.names=FALSE)
}

features_derived("test.csv","testout.csv")
