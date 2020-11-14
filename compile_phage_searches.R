library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)


##merge Tayor's PHASTER analysis
#Tayor copied contigs into PHASTER server by hand and saved output files
phaster_dirs <- list.files("raw_data/phaster/")
phaster_summaries <- data.frame(matrix(ncol=18, nrow=1))
headerline <-  c("CONTIG", "REGION" ,"REGION_LENGTH" ,"COMPLETENESS(score)","SPECIFIC_KEYWORD","REGION_POSITION","TRNA_NUM","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM", "HYPOTHETICAL_PROTEIN_NUM","PHAGE+HYPO_PROTEIN_PERCENTAGE","BACTERIAL_PROTEIN_NUM","ATT_SITE_SHOWUP","PHAGE_SPECIES_NUM","MOST_COMMON_PHAGE_NAME(hit_genes_count)","FIRST_MOST_COMMON_PHAGE_NUM","FIRST_MOST_COMMON_PHAGE_PERCENTAGE","GC_PERCENTAGE")
names(phaster_summaries) <- headerline
phaster_summaries$CONTIG <- str_replace(phaster_summaries$CONTIG,">","")
for(i in 1:length(phaster_dirs)){
  contig <- readLines(file(paste("raw_data/phaster/",phaster_dirs[i], "/detail.txt", sep=""),"r"),n=1)
  sum <- read.table(paste("raw_data/phaster/",phaster_dirs[i], "/summary.txt", sep=""),
                    fill = F, quote="", stringsAsFactors = F, skip = 33)
  sum <- cbind(rep(contig,nrow(sum)),sum)
  names(sum) <- headerline
  phaster_summaries <- rbind(phaster_summaries,sum)
  #genes <- read.table(paste("raw_data/phaster/",phaster_dirs[i], "/detail.txt", sep=""),
  #                    fill = R, stringsAsFactors = F, quote="")
  #genes <- genes[-3,]
}
phaster_summaries <- phaster_summaries[-1,]
#write.table(phaster_summaries,"analysis/phaster_summary.csv", quote = F, row.names = F, sep = "\t")
phaster_summaries_short <- phaster_summaries[,c("CONTIG","COMPLETENESS(score)", "SPECIFIC_KEYWORD", "TOTAL_PROTEIN_NUM", "PHAGE_HIT_PROTEIN_NUM", "ATT_SITE_SHOWUP", "MOST_COMMON_PHAGE_NAME(hit_genes_count)")]
phaster_smash <- aggregate(.~CONTIG, phaster_summaries_short, paste, collapse=",")

flye_info <- read.csv("raw_data/assembly_info_flye.txt", sep = "\t", stringsAsFactors = F)
names(flye_info) <- c("scaffold","length","cov","circ","repeat","mult","alt_grp","graph_path")
vibrant_pred <- read.csv("raw_data/VIBRANT_machine_assembly_flye_1000_filtered.tsv", sep="\t", stringsAsFactors = F)
row.names(vibrant_pred) <- vibrant_pred$scaffold
vibrant_pred$scaffold <- str_replace(vibrant_pred$scaffold,"_fragment_\\d","")
vibrant_pred_smash <- aggregate(prediction~scaffold, vibrant_pred, paste, collapse=",")
vibrant_qual <- read.csv("raw_data/VIBRANT_genome_quality_assembly_flye_1000_filtered.tsv", sep="\t", stringsAsFactors = F)
vibrant_qual$scaffold <- str_replace(vibrant_qual$scaffold,"_fragment_\\d","")
vibrant_qual_smash <- aggregate(.~scaffold, vibrant_qual, paste, collapse=",")
)

comb_tab <- merge(flye_info, vibrant_pred_smash, by="scaffold", all=T)
comb_tab <- merge(comb_tab, vibrant_qual_smash, by="scaffold", all=T)
comb_tab <- merge(comb_tab, phaster_smash, by=1, all=T)

#add blast seraches
matt_blast <- read.csv("raw_data/blastdb/flye_meta_vs_matthijnssen.blastn", sep="\t", stringsAsFactors = F)
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
matt_blast[,3:12] <- sapply(matt_blast[,3:12], as.numeric)
engel_blast[,3:12] <- sapply(engel_blast[,3:12], as.numeric)
comb_tab$matt_match <- NA
comb_tab$matt_match_len <- NA
comb_tab$matt_match_per <- NA
comb_tab$engel_match <- NA
comb_tab$engel_match_per <- NA
comb_tab$engel_match_len <- NA

for(i in 1:nrow(comb_tab)){
  if(comb_tab[i,]$scaffold %in% matt_blast$qseqid){
    mmatch <- matt_blast[matt_blast$qseqid==comb_tab[i,]$scaffold,] 
    mmatch_smash <- mmatch[,c("sseqid","pident","length","evalue", "bitscore")]
    mmatch_smash <- aggregate(.~sseqid, mmatch_smash, FUN=mean)
    best <- which(mmatch_smash$bitscore==max(mmatch_smash$bitscore))
    if(length(best) > 1){
      mmatch_smash[best[1],]$sseqid <- paste(mmatch_smash[best,]$sseqid,collapse=",")
      best <- best[1]
    }
    comb_tab[i,]$matt_match <- mmatch_smash[best,]$sseqid
    comb_tab[i,]$matt_match_per <- mmatch_smash[best,]$pident
    comb_tab[i,]$matt_match_len <- mmatch_smash[best,]$length
  }
  if(comb_tab[i,]$scaffold %in% engel_blast$qseqid){
    ematch <- engel_blast[engel_blast$qseqid==comb_tab[i,]$scaffold,] 
    ematch_smash <- ematch[,c("sseqid","pident","length","evalue", "bitscore")]
    ematch_smash <- aggregate(.~sseqid, ematch_smash, FUN=mean)
    best <- which(ematch_smash$bitscore==max(ematch_smash$bitscore))
    if(length(best) > 1){
      ematch_smash[best[1],]$sseqid <- paste(ematch_smash[best,]$sseqid,collapse=",")
      best <- best[1]
    }
    comb_tab[i,]$engel_match <- ematch_smash[best,]$sseqid
    comb_tab[i,]$engel_match_per <- ematch_smash[best,]$pident
    comb_tab[i,]$engel_match_len <- ematch_smash[best,]$length
  }
}

comb_tab$match <- "none"
comb_tab[which(!is.na(comb_tab$engel_match)),]$match <- "engel"
comb_tab[which(!is.na(comb_tab$matt_match)),]$match <- "matthijnssen"
comb_tab[which(!is.na(comb_tab$matt_match) & !is.na(comb_tab$engel_match)),]$match <- "both"

#write.table(comb_tab,"analysis/phaster_summary2.csv", quote = F, row.names = F, sep = "\t")
comb_tab <- read.csv("analysis/phaster_summary2.csv", sep="\t")





