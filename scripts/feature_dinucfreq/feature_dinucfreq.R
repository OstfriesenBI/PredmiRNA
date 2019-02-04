library(seqinr)

feature_dinucfrq <- function(inputcsv){
  
  # Reads .csv file
  csvdf <- read.csv(inputcsv,header = TRUE, sep = ",",stringsAsFactors = FALSE)
  
  # Get Sequence Column from DF
  mySequence <- csvdf[, "sequence"]
  
  print(mySequence)
  
  # For-Loop to iterate over the sequences and apply func count()
  for(i in mySequence){
    print(count(tolower(s2c(i)), wordsize = 2, alphabet = s2c("acgu")))
  }
}

feature_dinucfrq("test.csv")