library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(dplyr)


comb_tab <- read.table("analysis/genome_combine_summary.tab", sep="\t", header=T)             
comb_tab.phage <- comb_tab[comb_tab$IsPhage,]
ggplot(comb_tab.phage, aes(x=Flye.length, color=Vibrant.type, fill=Vibrant.type)) +
  geom_histogram()
