
feature_NEFE <- function(csvpath){

  csvdf <- read.csv(csvpath, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  

  # Apply the function GC on the matrix of the data, which only contains the sequences  
  gccontent <- apply(as.matrix(csvdf[, "sequence"]),1,nchar)
  
  total <- cbind(comment=csvdf[,"comment"], gccontent=gccontent)

  # Writes Data Frame to a .csv file
  write.csv(total, file="feature_gccontent.csv", row.names=FALSE)
}

# Call of the function to test, has to be removed for Snakemake (!)
feature_gccontent("test.csv")
