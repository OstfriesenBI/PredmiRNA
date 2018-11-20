library(Biostrings)
transferFastaTo <- function(in_path, out_path, realmiRNA) {


    fastaFile <- readRNAStringSet(in_path)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)

    df <- data.frame(seq_name, sequence, realmiRNA)
    
    write.csv(df,out_path,row.names=FALSE)
}

transferFastaTo("Test.fasta", "Datenbank.csv", 1)
