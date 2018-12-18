library(foreign)
data = read.csv("test.csv",stringsAsFactors = FALSE)
write.arff(x=data, file= "file.arff", eol = "\n", relation = deparse(substitute(x)))


