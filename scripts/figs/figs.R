#!/usr/bin/env Rscript

library(foreign)
library(ggplot2)
library(reshape2)
library(ggfortify)

trainarff <- "data/all.arff"
outdir <- "outputfig"
classname  <- "realmiRNA"

if(exists("snakemake")){
  trainarff <- snakemake@input[["data"]]
  outdir <- snakemake@output[["outdir"]]
}

dir.create(outdir)

data <- read.arff(trainarff)

if(!exists("snakemake")){
  data<-data[,!grepl("permuted_.+",colnames(data))]
}

theme_set(theme_bw())

genPlots <- function(featurename,featurenames,bins=32){
  histo <- ggplot(data) + aes_string(x=featurename,fill=classname) + geom_histogram(bins = bins, position = "identity",alpha=0.45)
  otherfeaturenames <-  setdiff(featurenames,c(featurename))
  genPointPlot <- function(featurename2){
    ggplot(data) + aes_string(x=featurename,y=featurename2,color=classname) + geom_point()
  }
  pointplots=lapply(otherfeaturenames,genPointPlot)
  names(pointplots) <- otherfeaturenames
  res <- list(histo=histo,pointplots=pointplots)
  res
}
featurenames = setdiff(colnames(data),c(classname))
k<-lapply(featurenames,genPlots,featurenames=featurenames)
names(k) <- featurenames

# Create dirs
s<-lapply(names(k),function(x)dir.create(file.path(outdir,x)))
# Save histos
s<-lapply(names(k),function(y)ggsave(file.path(outdir,y,paste0(y,"_histo.png")),k[[y]]$histo,"png",width = 8,height = 8,units = "in"))
# Save point plots
s<-lapply(names(k),function(x)lapply(names(k[[x]]$pointplots),function(y)ggsave(file.path(outdir,x,paste0(x,"_2_",y,".png")),k[[x]]$pointplots[[y]],"png",width = 8,height = 8,units = "in")))


# Copied from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

nums <- unlist(lapply(data, is.numeric))
datan <- data[ , nums]

plotCorrHeatMap<-function(method){
  pear <- cor(datan, method = method)
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  pear<-reorder_cormat(pear)
  pear[lower.tri(pear)]<- NA
  ggplot(data = melt(pear), aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste0(method,"\nCorrelation")) +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90))+
    coord_fixed()
}
methods <- c("pearson","kendall","spearman")
corplots<-lapply(methods,plotCorrHeatMap)
names(corplots) <- methods
s<-lapply(methods, function(x)ggsave(file.path(outdir,paste0(x,"_cor.png")),corplots[[x]],"png",width = 10,height = 10,units = "in"))

pca <- prcomp(datan, center = TRUE,scale. = TRUE)
autoplot(pca, data = data, colour = classname)

ggsave(file.path(outdir,"pca_1_2_points.png"), autoplot(pca, data = data, colour = classname,x=1,y=2),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdir,"pca_1_2_eigen.png"), autoplot(pca, data = data, colour = classname,x=1,y=2,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdir,"pca_2_3_points.png"), autoplot(pca, data = data, colour = classname,x=2,y=3),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdir,"pca_2_3_eigen.png"), autoplot(pca, data = data, colour = classname,x=2,y=3,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdir,"pca_1_3_points.png"), autoplot(pca, data = data, colour = classname,x=1,y=3),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdir,"pca_1_3_eigen.png"), autoplot(pca, data = data, colour = classname,x=1,y=3,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5),"png",width = 10,height = 10,units = "in")

