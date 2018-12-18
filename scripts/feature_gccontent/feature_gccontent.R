library(seqinr)

feature_gccontent <- function(csvpath){

  # Reads .csv file into the Data Fram "csvdf"
  csvdf <- read.csv(csvpath, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # Gsub replaces all matches of a string; in this case "U" is replaced with "T" to "simulate" a DNA
  csvdf[, "sequence"] <- gsub("U","T",csvdf[, "sequence"])

  # Apply the function GC on the matrix of the data, which only contains the sequences  
  gccontent <- apply(as.matrix(csvdf[, "sequence"]),1,function(rnaseq){GC(s2c(rnaseq))})
  
  # Merged the the comment and gccontent column to a new Data Frame
  total <- cbind(comment=csvdf[,"comment"], gccontent=gccontent)

  # Writes Data Frame to a .csv file
  write.csv(total, file="feature_gccontent.csv", row.names=FALSE)
}

# Call of the function to test, has to be removed for Snakemake (!)
feature_gccontent("test.csv")
