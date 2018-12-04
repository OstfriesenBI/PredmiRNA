library(Biostrings)
# Loads the Biosting library into the script
# Note that the library has to be installed befor running the script

# Defines a function named "fasta2csv"
fasta2csv <- function(in_path, out_path, realmiRNA) {

    # Uses the function readRNAStringSet from Biostrings on given fasta file (in_path)
    fastaFile <- readRNAStringSet(in_path)

    # Assings the name and the sequence from the fasta file to a variable
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)

    # Creates a data frame
    df <- data.frame(seq_name, sequence, realmiRNA)
    
    # Writes the data frame to csv file to given path without extra row names
    write.csv(df,out_path,row.names=FALSE)
}

# Calls the function to test it, has to be swapped out late to be used with snakemake
fasta2csv("Test.fasta", "Datenbank.csv", 1)
