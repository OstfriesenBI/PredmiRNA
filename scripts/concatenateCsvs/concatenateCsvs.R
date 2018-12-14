
concatenateCsvs = function(mypath){ 
  filenames=list.files(path=mypath, full.names=TRUE, pattern= "")
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {rbind(x,y)}, datalist)}

concatenateCsvs("mypath = Pfad zum Ordner mit den zu kombinierenden Dateien")

