#!/usr/bin/env Rscript

library(foreign)
library(ggplot2)

featuresets <- c("All_Literature","Bor_Confirmed","Weka_Confirmed")
featuresets.labels <- c("Literatur","Boruta","Weka SVMeval")
datafiles <- as.list(paste0(paste0("data/figs/",featuresets),"/comparison.csv"))
featurefiles <- as.list(paste0(paste0(paste0("data/models/",featuresets),"/",featuresets),"_train.arff"))
cputime <- "time.png"
mem <-"mem.png"
fmeasures <- "fmeasure.png"
roc <- "roc.png"
selections <- "selections.csv"
featurecount <- "featurecount.png"
combinedoutputfile <- "combinedoutputfile.csv"

if(exists("snakemake")){
  featuresets <- snakemake@params[["featuresets"]]
  datafiles <- snakemake@input[["datafiles"]]
  featurefiles <- snakemake@input[["featurefiles"]]
  combinedoutputfile <- snakemake@output[["combinedoutputfile"]]
  cputime <- snakemake@output[["cputime"]]
  mem <- snakemake@output[["mem"]]
  fmeasures <- snakemake@output[["fmeasures"]]
  roc <- snakemake@output[["roc"]]
  selections <- snakemake@output[["selections"]]
  featurecount <- snakemake@output[["featurecount"]]
  featuresets.labels <- snakemake@params[["setlabels"]]
}


names(datafiles) <- featuresets
names(featurefiles) <- featuresets

fselection <- lapply(featuresets,function(x)colnames(read.arff(featurefiles[[x]])))
features <- unique(unlist(fselection))
featureset_sel <- cbind(Feature=features,do.call(cbind,lapply(fselection,function(x)as.data.frame(as.numeric(features %in% x)))))
colnames(featureset_sel)[2:(1+length(featuresets))]<-featuresets
write.csv(featureset_sel,selections)

data <- do.call(rbind,lapply(featuresets,function(x)cbind(`Feature Set`=x,read.csv(datafiles[[x]]))))
data$X <- NULL
data <- merge(data,data.frame(`FeatureSet`=featuresets,`feature_count`=sapply(featuresets,function(set)sum(featureset_sel[[set]]))),by.x="Feature Set",by.y="FeatureSet")
data$`Feature Set`=factor(data$`Feature Set`,levels = featuresets,labels = featuresets.labels)


p1<-ggplot(data) + aes(x=`Feature Set`,y=`F.Measure`,color=Model,group=Model) +geom_point()+geom_line() + theme(legend.position = "bottom")

ggsave(roc,p1+aes(y=ROC.Area)+ylab("ROC Area"),"png",width = 21,height = 15.3,units = "cm")
ggsave(fmeasures,p1+aes(y=F.Measure)+ylab("F Measure"),"png",width = 21,height = 15.3,units = "cm")
ggsave(mem,p1+aes(y=median_max_rss)+ylab("Used memory in MB"),"png",width = 14,height = 14,units = "cm")
ggsave(cputime,p1+aes(y=median_s)+ylab("CPU time in s"),"png",width = 14,height = 14,units = "cm")
ggsave(featurecount,ggplot(data[data$Model==data$Model[[1]],],aes(`Feature Set`,`feature_count`,group=1))+geom_point()+geom_line()+xlab("Feature Set")+ylab("Number of Features") + theme(legend.position = "bottom"),"png",width = 21,height = 15.3,units = "cm")

write.csv(data,combinedoutputfile)
