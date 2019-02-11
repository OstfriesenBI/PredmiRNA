library(seqinr)

feature_dinucfrq <- function(inputcsv){
  
  # Reads .csv file
  csvdf <- read.csv(inputcsv,header = TRUE, sep = ",",stringsAsFactors = FALSE)
  
  # Get Sequence Column from DF
  mySequence <- csvdf[, "sequence"]
  print(mySequence)
  
  # Output Dataframe
  countDF <- data.frame(Dinucs = character(),
                         Count = integer())

  
  # For-Loop to iterate over the sequences and apply func count()
  for(i in mySequence){
    countDF <- count(tolower(s2c(i)), wordsize = 2, alphabet = s2c("acgu"))
  }
  # Merge DFs
  outputDF <- merge(mySequence, countDF)
  # outputDF <- t(outputDF)
  outputDF <- outputDF[-c(1)]
  View(outputDF)
}

feature_dinucfrq("test.csv")