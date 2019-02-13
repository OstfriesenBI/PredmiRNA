library(seqinr)

feature_dinucfrq <- function(inputcsv){
  
  # Reads .csv file
  csvdf <- read.csv(inputcsv,header = TRUE, sep = ",",stringsAsFactors = FALSE)
  
  # Get Sequence Column from DF
  mySequence <- csvdf[, "sequence"]

  # Initialize Count Dataframe
  countDF <- data.frame(Dinucs = character(),
                         Count = integer())

  # Initialize DataFrame for results
  resultDF <- data.frame(Sequence = character(),
                         aa = integer(),
                         ac = integer(),
                         ag = integer(),
                         au = integer(),
                         ca = integer(),
                         cc = integer(),
                         cg = integer(),
                         cu = integer(),
                         ga = integer(),
                         gc = integer(),
                         gg = integer(),
                         gu = integer(),
                         ua = integer(),
                         uc = integer(),
                         ug = integer(),
                         uu = integer())
  
  # DinucArray
  dinucArray <- c("aa", "ac", "ag", "au", "ca", "cc", "cg", "cu", "ga", "gc", "gg", "gu", "ua", "uc", "ug", "uu")
  
  # For-Loop to iterate over the sequences and apply func count()
  for(i in mySequence){
    countDF <- count(tolower(s2c(i)), wordsize = 2, alphabet = s2c("acgu"))
    
  }
}

feature_dinucfrq("test.csv")