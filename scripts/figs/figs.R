#!/usr/bin/env Rscript

library(foreign)
library(ggplot2)
library(reshape2)
library(Boruta)
library(ggfortify)

set.seed(1)

trainarff <- "data/models/Weka_Confirmed/Weka_Confirmed_train.arff"
outdirpca <- "outputfig"
outdirfeat <- "outputfig"
classname  <- "realmiRNA"
burotadata <- "boruta.csv"
burotaplot <- "boruta.png"

if(exists("snakemake")){
  trainarff <- snakemake@input[["data"]]
  outdirpca <-  snakemake@output[["outdir_pca"]]
  outdirfeat <-  snakemake@output[["outdir_feat"]]
  borutalog <-  snakemake@output[["borutalog"]]
  burotadata <-  snakemake@output[["borutadata"]]
  burotaplot <-  snakemake@output[["burotaplot"]]
}

theme_set(theme_bw())

dir.create(outdirpca)
dir.create(outdirfeat)

data <- read.arff(trainarff)

print("Running Boruta")

sink(borutalog)
data[is.na(data)] <- 0 
data.boruta <- Boruta(realmiRNA~., data = data, doTrace = 2,maxRuns=100)
print(data.boruta)
data.boruta <- TentativeRoughFix(data.boruta)
print(data.boruta)
formatres<-function(group)paste0(lapply(names(data.boruta$finalDecisio[data.boruta$finalDecision==group]),function(x)paste0("\"",x,"\"")),collapse = ",")
print("Confirmed")
cat(formatres("Confirmed"))
print("\nRejected")
cat(formatres("Rejected"))
sink()
data.boruta.df <- attStats(data.boruta)
data.boruta.df$feature <- rownames(data.boruta.df)
rownames(data.boruta.df) <- NULL
write.csv(data.boruta.df,burotadata)

png(burotaplot,width = 1280,height = 720)
plot(data.boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(data.boruta$ImpHistory),function(i)
  data.boruta$ImpHistory[is.finite(data.boruta$ImpHistory[,i]),i])
names(lz) <- colnames(data.boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(data.boruta$ImpHistory), cex.axis = 0.7)
dev.off()




genPlots <- function(featurename,featurenames,bins=32){
  histo <- ggplot(data) + aes_string(x=featurename,fill=classname) + geom_histogram(bins = bins, position = "identity",alpha=0.45) + theme(legend.position = "bottom")
  boxplot <- ggplot(data) + aes_string(x=classname,y=featurename) + geom_boxplot() + theme(legend.position = "bottom")
  # otherfeaturenames <-  setdiff(featurenames,c(featurename))
  # genPointPlot <- function(featurename2){
  #   ggplot(data) + aes_string(x=featurename,y=featurename2,color=classname) + geom_point() + theme(legend.position = "bottom")
  # }
  # pointplots=lapply(otherfeaturenames,genPointPlot)
  # names(pointplots) <- otherfeaturenames
  res <- list(histo=histo,boxplot=boxplot)#,pointplots=pointplots)
  res
}
featurenames = setdiff(colnames(data),c(classname))
print("Generating feat plots ...")
k<-lapply(featurenames,genPlots,featurenames=featurenames[!grepl("permuted_.+",featurenames)])
names(k) <- featurenames

print("Saving feat plots ...")
# Save histos
s<-lapply(names(k),function(y)ggsave(file.path(outdirfeat,paste0(y,"_histo.png")),k[[y]]$histo,"png",width = 8,height = 8,units = "in"))
s<-lapply(names(k),function(y)ggsave(file.path(outdirfeat,paste0(y,"_boxplot.png")),k[[y]]$boxplot,"png",width = 8,height = 8,units = "in"))
# Save point plots
#s<-lapply(names(k),function(x)lapply(names(k[[x]]$pointplots),function(y)ggsave(file.path(outdirfeat,x,paste0(x,"_2_",y,".png")),k[[x]]$pointplots[[y]],"png",width = 8,height = 8,units = "in")))


# Copied from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

nums <- unlist(lapply(data, is.numeric))
datan <- data[ , nums]
datan <- datan[ , 0!=unlist(lapply(datan,sd))]

plotCorrHeatMap<-function(method){
  pear <- cor(datan, method = method)
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  dd <- as.dist((1-pear)/2)
  hc <- hclust(dd)
  pear <-pear[hc$order, hc$order]
  pear_all <- pear
  pear[lower.tri(pear)]<- NA
  pl<-ggplot(data = melt(pear), aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste0(method,"\nCorrelation")) +
    theme_minimal(base_size = 6)+ 
    theme(axis.text.x = element_text(angle = 90))+
    coord_fixed()
  list(cormat=pear_all,dmat=dd,clstr=hc,plot=pl)
}
methods <- c("pearson","spearman")
print("Generating corr plots ...")
corplots<-lapply(methods,plotCorrHeatMap)
names(corplots) <- methods
print("Saving corr plots ...")
s<-lapply(methods, function(x)ggsave(file.path(outdirpca,paste0(x,"_cor.png")),corplots[[x]][["plot"]],"png",width = 10,height = 10,units = "in"))
for(m in methods){
  png(file.path(outdirpca,paste0(m,"_cor_vis_clust.png")),width = 1280*2,height = 720*2)
  par(cex=0.8, mar=c(5, 8, 4, 1))
  plot(corplots[[m]]$clstr,main=NULL,xlab="Feature",ylab="(1-Correlation)/2")
  dev.off()
  clust <- hclust(as.dist(1-abs(corplots[[m]]$cormat)))
  cutheight <- 0.33333
  grouplabels <- cutree(clust,h=cutheight)
  grouplabels <- as.data.frame(grouplabels)
  grouplabels$feature <- rownames(grouplabels)
  rownames(grouplabels) <- NULL
  groups <-merge(data.boruta.df, grouplabels,by="feature")
  groups <- groups[order(groups$grouplabels, -groups$medianImp),]
  selectedfeatures <- aggregate(.~grouplabels,groups,FUN=head,1)
  sink(file.path(outdirpca,paste0(m,"_cor_abs_clust_log.txt")))
  print("Height")
  print(cutheight)
  print("Best feature according to boruta median importance")
  cat(paste0(lapply(selectedfeatures$feature,function(x)paste0("\"",x,"\"")),collapse = ","))
  sink()
  png(file.path(outdirpca,paste0(m,"_cor_abs_clust.png")),width = 1280*2,height = 720*2)
  par(cex=0.8, mar=c(5, 8, 4, 1))
  plot(clust,main=NULL,xlab="Feature",ylab="(1-abs(Correlation))")
  dev.off()
}


print("Generating pca plots ...")
pca <- prcomp(datan, center = TRUE,scale. = TRUE)
#autoplot(pca, data = data, colour = classname)

print("Saving pca plots ...")
ggsave(file.path(outdirpca,"pca_1_2_points.png"), autoplot(pca, data = data, colour = classname,x=1,y=2) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdirpca,"pca_1_2_eigen.png"), autoplot(pca, data = data, colour = classname,x=1,y=2,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdirpca,"pca_2_3_points.png"), autoplot(pca, data = data, colour = classname,x=2,y=3) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdirpca,"pca_2_3_eigen.png"), autoplot(pca, data = data, colour = classname,x=2,y=3,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdirpca,"pca_1_3_points.png"), autoplot(pca, data = data, colour = classname,x=1,y=3) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")
ggsave(file.path(outdirpca,"pca_1_3_eigen.png"), autoplot(pca, data = data, colour = classname,x=1,y=3,loadings=T, loadings.label = TRUE, loadings.colour = 'blue',loadings.label.size = 5) + theme(legend.position = "bottom"),"png",width = 10,height = 10,units = "in")

