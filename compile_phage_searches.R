library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)


##merge Tayor's PHASTER analysis
#Tayor copied contigs into PHASTER server by hand and saved output files
phaster_dirs <- list.files("raw_data/phaster/")
phaster_summaries <- data.frame(matrix(ncol=18, nrow=1))
headerline <-  c("CONTIG", "REGION" ,"REGION_LENGTH" ,"COMPLETENESS(score)","SPECIFIC_KEYWORD","REGION_POSITION","TRNA_NUM","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM", "HYPOTHETICAL_PROTEIN_NUM","PHAGE+HYPO_PROTEIN_PERCENTAGE","BACTERIAL_PROTEIN_NUM","ATT_SITE_SHOWUP","PHAGE_SPECIES_NUM","MOST_COMMON_PHAGE_NAME(hit_genes_count)","FIRST_MOST_COMMON_PHAGE_NUM","FIRST_MOST_COMMON_PHAGE_PERCENTAGE","GC_PERCENTAGE")
names(phaster_summaries) <- headerline
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
write.table(phaster_summaries,"analysis/phaster_summary.csv", quote = F, row.names = F, sep = "\t")
