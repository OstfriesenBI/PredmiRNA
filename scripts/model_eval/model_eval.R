#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(ggfortify)

thfile = list("data/models/All_Literature/threshold/LibSVM.csv","data/models/All_Literature/threshold/NaiveBayes.csv")
log = list("data/models/All_Literature/LibSVM.log","data/models/All_Literature/NaiveBayes.log")
bench = list("data/models/All_Literature/LibSVM.benchmark.txt","data/models/All_Literature/NaiveBayes.benchmark.txt")

modelnames=list("libsvm","bayes")

roc="roc.png"
time="time.png"
mem="memory.png"
data="comparison.csv"
measure="measures.png"

if(exists("snakemake")){
  thfile = snakemake@input[["thfile"]]
  log = snakemake@input[["log"]]
  bench = snakemake@input[["bench"]]
  roc = snakemake@output[["roc"]]
  time = snakemake@output[["time"]]
  mem = snakemake@output[["mem"]]
  data = snakemake@output[["data"]]
  measure = snakemake@output[["measure"]]
  modelnames = snakemake@params[["modelnames"]]
}
theme_set(theme_bw())
modelnames <- unlist(modelnames)

threshvalues <- lapply(thfile,function(x)read.csv(x,stringsAsFactors = F))
threshvalues <- do.call(rbind,lapply(seq_along(modelnames),function(x)cbind(threshvalues[[x]],Model=modelnames[[x]])))
pl <- ggplot(threshvalues) + aes(X.False.Positive.Rate.,X.True.Positive.Rate.,color=Model) + geom_line() + geom_abline(slope=1,intercept=0) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + xlab("1-Specificity") + ylab("Sensitivity") + theme(legend.position = "bottom")
ggsave(roc,pl,"png",height = 16, width = 16, units = "cm")
rm(pl)

benchvals <- lapply(bench,function(x)read.delim(x))
benchvals <- do.call(rbind,lapply(seq_along(modelnames),function(x)cbind(benchvals[[x]],Model=modelnames[[x]])))
# RSS: RSS is the Resident Set Size and is used to show how much memory is allocated to that process and is in RAM. 
# It does not include memory that is swapped out. It does include memory from shared libraries as long as the pages from those libraries are actually in memory. 
# It does include all stack and heap memory.

# VSZ is the Virtual Memory Size. It includes all memory that the process can access, including memory that is swapped out, memory that is allocated, but not used, and memory that is from shared libraries.

# USS: Unique Set size, physical unshared memory

# PSS: Proportional Set Size, USS + process proportion of shared memory

pl <- ggplot(benchvals) + aes(Model,s) + geom_boxplot()+ xlab("Model") + ylab("CPU clock time in sec") + theme(legend.position = "bottom")
ggsave(time,pl,"png",height = 16, width = 16, units = "cm")
rm(pl)
                      
pl <- ggplot(benchvals) + aes(Model,max_rss) + geom_boxplot()+ xlab("Model") + ylab("Memory Usage (Resident Set Size) in MB") + theme(legend.position = "bottom")
ggsave(mem,pl,"png",height = 16, width = 16, units = "cm")
rm(pl)

benchvals <- benchvals %>% group_by(Model) %>% select_if(is.numeric) %>% summarise_all(median) %>% ungroup() %>% mutate(Model=as.character(Model))
colnames(benchvals) <- c("Model",paste0("median_",setdiff(colnames(benchvals),"Model")))

parseLogs <- function(x){
  cols <- c("TPR","FPR","Precision","Recall","F-Measure","MCC","ROC Area","RPC Area")
  pat<- paste0("Weighted Avg.",paste0(rep("\\s+([-+]?[0-9]*\\.?[0-9]+)",length(cols)),collapse = ""))
  vals<- t(as.matrix(as.numeric(str_match(x,pat)[,-1])))
  colnames(vals) <- cols
  vals
}

logtxt <- lapply(log, function(filename)readChar(filename, file.info(filename)$size))
logtxt <- cbind(Model=modelnames,as.data.frame(do.call(rbind,lapply(logtxt,parseLogs)))) %>% as_tibble()

alldata <- inner_join(logtxt%>% mutate_if(is.factor,as.character),benchvals%>% mutate_if(is.factor,as.character),by="Model")

molten_alldata<- alldata%>% mutate_if(is.factor,as.character)%>% rename(`Median RAM Usage in MB`=median_max_rss,`Median CPU time in s`=median_s) %>% melt(id.vars = c("Model"),measure.vars =c("Recall","Precision","F-Measure","ROC Area","Median RAM Usage in MB","Median CPU time in s"),value.name = "value", variable.name="Measure")

pl <- molten_alldata %>% group_by(Measure) %>%  mutate(multi=case_when(Measure == "Median RAM Usage in MB" |Measure == "Median CPU time in s" ~ 1, TRUE ~-1))%>% 
  mutate(value=value*multi) %>% mutate(Rank=rank(value)) %>% mutate(value=value*multi)  %>%
  ggplot() + aes(x=Model,y=Measure,fill=Rank) + geom_tile() +geom_text(aes(label=sprintf("%0.2f", round(value, digits = 2)))) +  scale_fill_gradient(low="green",high="red")
ggsave(measure,pl,"png",height = 14, width = 22, units = "cm")


write.csv(alldata,data)
