
multmerge = function(mypath){ #Pfad im ordener wo die csv sind muss angegeben werden.
  filenames=list.files(path=mypath, full.names=TRUE, pattern= "") #.= beliebiges zeichen, + = whdGruppe, \csv. für dateinamen, es sollen nur csv einglesen werden
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  #lapply ruft eine funktion für jedes element aus der liste filenames, Das Ergebnis ist eine Liste die genauso groß ist mit dem output der funktion für jedes Element
  #file=x er liest aus den dateipfad x
  Reduce(function(x,y) {merge(x,y)}, datalist)} #reduce funktion produziert einen output der selber wieder als Input benutz werden kann. 
  #merge kombiniert nur 2, und mit reduce können quasi die folgenden csv Dateien zur Summe der bereits existierenden hinzugefügt werden.
#besonderheit bei R: letze funktion in der Gesamtfunktion ist der Rückgabewert (statt return)
multmerge("example_files")


