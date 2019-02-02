library(seqinr)

feature_dinucfrq <- function(inputcsv){
  
  # Reads .csv file
  csvdf <- read.csv(inputcsv,header = TRUE, sep = ",",stringsAsFactors = FALSE)
  
  # Get Sequence Column from DF
  mySequence <- csvdf[, "sequence"]
  
  print(mySequence)
  
  # PROBLEM: Dunno what is going on here rn
  countDinuc <- list()
  
  for(i in mySequence){
    tempCount <- count(i, wordsize=2, alphabet = s2c("acgu"))
    countDinuc <- append(tempCount)
  }
  # END of confusion
  
  print(countDinucs)
  
}

feature_dinucfrq("test.csv")