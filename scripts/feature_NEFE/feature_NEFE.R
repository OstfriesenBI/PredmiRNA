
feature_NEFE <- function(infile,out){

  csvdf <- read.csv(infile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  

  # Apply the function GC on the matrix of the data, which only contains the sequences  
  length <- apply(as.matrix(csvdf[, "sequence"]),1,nchar)
  NEFE <- csvdf[,"efe"]/length

  total <- cbind(comment=csvdf[,"comment"], length=length, nefe=NEFE)

  # Writes Data Frame to a .csv file
  write.csv(total, file=out, row.names=FALSE)
}

feature_NEFE("test.csv","testout.csv")
