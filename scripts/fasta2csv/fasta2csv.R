library(Biostrings)
# Loads the Biosting library into the script
# Note that the library has to be installed befor running the script

# Defines a function named "fasta2csv"
fasta2csv <- function(in_path, out_path, realmiRNA) {

    # Uses the function readRNAStringSet from Biostrings on given fasta file (in_path)
    fastaFile <- readBStringSet(in_path)

    # Assings the name and the sequence from the fasta file to a variable
    comment = names(fastaFile)
    sequence = paste(fastaFile)
    
    lowercase = sapply(regmatches(sequence, gregexpr("[a-z]", sequence, perl=TRUE)), length)
    lengths = sapply(sequence,nchar)
    lowcomplexity = (lowercase/lengths)
    
    # Length until we a hit a stop codon
    stopLen = function(seq,codon){
     stoplens=regexpr(codon,seq)-1
     selector=stoplens==-2
     stoplens[selector]=nchar(seq[selector])
     stoplens
    }

    # Creates a data frame
    df <- data.frame(comment, sequence, realmiRNA,lowcomplexity,lenWoUGA=stopLen(sequence,"UGA"),lenWoUAG=stopLen(sequence,"UAG"),lenWoUAA=stopLen(sequence,"UAA"))
    
    # Writes the data frame to csv file to given path without extra row names
    write.csv(df,out_path,row.names=FALSE)
}

if(exists("snakemake")){
	real=0
	if(grepl(snakemake@params[["realmarker"]],snakemake@input[[1]])){
		real=1
	}else if(grepl(snakemake@params[["classmarker"]],snakemake@input[[1]])){
		real=-1
	}
	fasta2csv(snakemake@input[[1]], snakemake@output[[1]], real)
}else{
	# Calls the function to test it, has to be swapped out late to be used with snakemake
	fasta2csv("Test.fasta", "Datenbank.csv", 1)
}
