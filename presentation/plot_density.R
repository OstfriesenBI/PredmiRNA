library(ggplot2)
library(foreign)

df <- read.arff("2018-2-11-1-20-56.arff")

feature_names <- names(df)

pdf("plot_density.pdf")

for (name in feature_names) {
  plot <-ggplot(df, aes_string(x=name, color="class"))+geom_density()+ggtitle(name)
  print(plot)
}

dev.off()