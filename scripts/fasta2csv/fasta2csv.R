library("Biostrings")
transferFastaTo <- function(in_path, out_path, real) {


    fastaFile <- readDNAStringSet(in_path)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)

    df <- data.frame(seq_name, sequence, real)
    write.csv(df,out_path,row.names=FALSE)
}

transferFastaTo("Test.fasta", "Datenbank.csv", 0)