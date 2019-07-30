#!/usr/bin/env Rscript

#
# Converts a .csv dataframe into a Weka usable .arff
#

library(foreign)
set.seed(1)
infile="test.csv"
outfile="file.arff"
selected_feat = list("dP,mfei1")
if(exists("snakemake")){
	infile=snakemake@input[[1]]
	outfile=snakemake@output[[1]]
        selected_feat = snakemake@params[["sel"]]
}
data = read.csv(infile,stringsAsFactors = FALSE)
data <- data[data$realmiRNA==-1,]
comments <- data$comment
if(selected_feat!=list("all")){
        data <- data[,c(unlist(selected_feat),"realmiRNA")]
}
is.na(data) <- do.call(cbind,lapply(data, is.infinite))
data[sapply(data, is.character)] <- list(NULL)
data$realmiRNA <- factor(data$realmiRNA)
data <- data[,c(setdiff(colnames(data),c("realmiRNA")),"realmiRNA")]
#data$comment <- comments
write.arff(x=data, file=outfile, eol = "\n", relation = deparse(substitute(x)))
