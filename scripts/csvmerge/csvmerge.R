#setwd("C:\\Users\\filip\\OneDrive\\Documents\\GitHub\\PredmiRNA\\scripts\\csvmerge\\")
multmerge = function(filenames){ #Pfad im ordener wo die csv sind muss angegeben werden.
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  #lapply ruft eine funktion für jedes element aus der liste filenames, Das Ergebnis ist eine Liste die genauso groß ist mit dem output der funktion für jedes Element
  #file=x er liest aus den dateipfad x
  Reduce(function(x,y) {merge(x,y)}, datalist)} #reduce funktion produziert einen output der selber wieder als Input benutz werden kann. 
  #merge kombiniert nur 2, und mit reduce können quasi die folgenden csv Dateien zur Summe der bereits existierenden hinzugefügt werden.
#besonderheit bei R: letze funktion in der Gesamtfunktion ist der Rückgabewert (statt return)

if(exists("snakemake")){
	merged=multmerge(snakemake@input)
	write.csv(merged,file=snakemake@output[[1]], row.names=FALSE)
}else{
	filenames=list.files(path="example_files", full.names=TRUE, pattern= "") #.= beliebiges zeichen, + = whdGruppe, \csv. für dateinamen, es sollen nur csv einglesen werden
	mymergeddata <- multmerge(filenames)
	View(datalist)
}
