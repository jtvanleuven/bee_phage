library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)


##compare to other bee phage datasets
matthijnssen <- read.csv("raw_data/blastdb/matthijnssen_Supptable_S18.csv", stringsAsFactors = F)
acc <- matthijnssen[!matthijnssen$Genbank.accession == 'nan',]$Genbank.accession
#slow steps
#matthijnssen_seqs <- read.GenBank(acc,species.names=T)
#write.dna(matthijnssen_seqs, file="raw_data/blastdb/matthijnssen.fasta", format="fasta", colsep = "")
#blastn -db matthijnssen.fasta -query ../flye_meta.fasta -evalue 1E-10 -outfmt 6 > flye_meta_vs_matt.blastn
#downloaded entire engel dataset from github
#also downloaded high-quality viral contigs from genbank. start here
#blastn -db JAAOBB01.1.engel.fasta -query ../flye_meta.fasta -evalue 1E-10 -outfmt 6 > flye_meta_vs_engel.blastn
#make three-way blast

matt_blast <- read.csv("raw_data/blastdb/matthijnssen_vs_engel.blastn", sep="\t", stringsAsFactors = F)
engel_blast <- read.csv("raw_data/blastdb/flye_meta_vs_engel.blastn", sep="\t", stringsAsFactors = F)
names <- c("qseqid",
           "sseqid",
           "pident",
           "length",
           "mismatch",
           "gapopen",
           "qstart",
           "qend",
           "sstart",
           "send",
           "evalue",
           "bitscore")
names(matt_blast) <- names
names(engel_blast) <- names
