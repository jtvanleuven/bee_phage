library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(readxl)
library(dplyr)
library(cowplot)


###fragment honey bee genome to run through PPRMeta (reviewer 2 of course :)
if(FALSE){
  assemblyStats <- read.csv('raw_data/assembly_info_flye.txt', fill = T, sep ='\t')
  assemblyStats <- assemblyStats[assemblyStats$length > 1000, ]   ###same filter as contigs
  assemblyStats$randbeeseq <- NA
  
  #GCA_003254395.2.fasta downloaded on Feb 10, 2022 from ENA
  bee <- read.fasta('../PPRMeta_beegenome/GCA_003254395.2.fasta')
  contiglens <- getLength(bee)
  for(i in 1:nrow(assemblyStats)){
    length <- assemblyStats[i,]$length
    picks <- which(contiglens > length*3)  ##make sure we pick from large bee contigs
    contig <- sample(picks, size = 1)
    start <- sample(1:getLength(bee[contig]), size = 1)
    if(start + length < getLength(bee[contig])){ ##make sure contig chromosome is long enough
      assemblyStats[i,]$randbeeseq <- paste(
        as.character(
          getFrag(bee[contig], begin = start, end = start+length)[[1]]),
        collapse = '')
    }else{
      assemblyStats[i,]$randbeeseq <- paste(
        as.character(
          getFrag(bee[contig], begin = start, end = start-length)[[1]]),  ##go backwards if not long enough
        collapse = '')
    }  
      
  }
  
  for(i in 1:nrow(assemblyStats)){
    write.fasta(assemblyStats[i,]$randbeeseq, i, '../PPRMeta_beegenome/randseqs.fasta', open = "a", nbchar = 60, as.string = FALSE)
  }
}



###take SI table 2 and modify it according to a few reviewer comments
### 1) add vHulk info
### 2) add checkV scores
### 3) add checkV clustering results. Note that I'm not removing contigs, just indicating when they belong to a cluster.

###simple lengths distribution plots 
comb_tab <- read.table("analysis/genome_combine_summary.tab", sep="\t", header=T)  

siTab2 <- read.table("revision/SI_table_2.tsv", sep="\t", fill = T)
names(siTab2) <- siTab2[1,]
siTab2 <- siTab2[-1,]

comb_tab <- siTab2

comb_tab$IsPhage <- as.logical(comb_tab$IsPhage)
comb_tab.phage <- comb_tab[comb_tab$IsPhage,]
comb_tab.phage$Vibrant.type <- str_replace(comb_tab.phage$Vibrant.type, 'lytic,lytic', 'lytic')
comb_tab.phage$Flye.length <- as.numeric(as.character(comb_tab.phage$Flye.length))
comb_tab.phage$Flye.length <- comb_tab.phage$Flye.length/1000
comb_tab.phage[which(is.na(comb_tab.phage$Vibrant.type)),]$Vibrant.type <- 'unidentified'
comb_tab.phage$Vibrant.type <- factor(comb_tab.phage$Vibrant.type, , levels=c('lytic','lysogenic','unidentified'))

###add checkV comparison
##compare completeness scores
drop <- c('\\d+','\\(','\\)')
comb_tab.phage$Phaster.COMPLETENESS.score <- str_remove_all(comb_tab.phage$Phaster.COMPLETENESS.score, paste(drop, collapse = '|'))
table(comb_tab.phage$Phaster.COMPLETENESS.score)
comb_tab.phage$Phaster.COMPLETENESS.score <- str_replace(comb_tab.phage$Phaster.COMPLETENESS.score, ',incomplete', '')
table(comb_tab.phage$Phaster.COMPLETENESS.score)
comb_tab.phage$Vibrant.Quality <- str_remove(comb_tab.phage$Vibrant.Quality, ',complete circular')

##checkV scores
#siTab2 <- read.table("revision/SI_table_2.tsv", sep="\t", fill = T)
#names(siTab2) <- siTab2[1,]
#siTab2 <- siTab2[-1,]

checkv <- read.csv('raw_data/checkV/out/quality_summary.tsv', sep="\t")
checkv_short <- checkv %>%
  select(contig_id, checkv_quality, completeness, warnings)
names(checkv_short) <- c('contig_id', 'checkV_quality', 'checkV_completeness', 'checkV_warning') 

comb_tab.phage <- merge(comb_tab.phage, checkv_short, by.x = 'Genome', by.y='contig_id')
#comb_tab.phage <- merge(comb_tab.phage, checkv_short, by.x = 'Genome', by.y='contig_id')


join1 <- paste(comb_tab.phage$Phaster.COMPLETENESS.score, comb_tab.phage$Vibrant.Quality, sep=',')
join2 <- paste(comb_tab.phage$Phaster.COMPLETENESS.score, comb_tab.phage$checkV_quality, sep=',')
join3 <- paste(comb_tab.phage$Vibrant.Quality, comb_tab.phage$checkV_quality, sep=',')

complete.tab <- data.frame(phaster=NA, vibrant=NA, table(join1))
complete.tab$phaster <- str_split(complete.tab$join1, ',', simplify = T)[,1]
complete.tab$vibrant <- str_split(complete.tab$join1, ',', simplify = T)[,2]
complete.tab$Freq <- as.numeric(complete.tab$Freq) 
complete.tab$phaster <- factor(complete.tab$phaster, 
                               levels = c('NA',  'questionable','incomplete', 'intact'))
complete.tab$vibrant <- factor(complete.tab$vibrant, 
                               levels = c('NA',  'low quality draft','medium quality draft', 'high quality draft'))

complete.tab.2 <- data.frame(phaster=NA, checkv=NA, table(join2))
complete.tab.2$phaster <- str_split(complete.tab.2$join2, ',', simplify = T)[,1]
complete.tab.2$checkv <- str_split(complete.tab.2$join2, ',', simplify = T)[,2]
complete.tab.2$Freq <- as.numeric(complete.tab.2$Freq) 
complete.tab.2$phaster <- factor(complete.tab.2$phaster, 
                                 levels = c('NA',  'questionable','incomplete', 'intact'))
complete.tab.2$checkv <- factor(complete.tab.2$checkv, 
                                levels = c('Not-determined',  'Low-quality','Medium-quality', 'High-quality', 'Complete'))


complete.tab.3 <- data.frame(vibrant=NA, checkv=NA, table(join3))
complete.tab.3$vibrant <- str_split(complete.tab.3$join3, ',', simplify = T)[,1]
complete.tab.3$checkv <- str_split(complete.tab.3$join3, ',', simplify = T)[,2]
complete.tab.3$Freq <- as.numeric(complete.tab.3$Freq) 
complete.tab.3$vibrant <- factor(complete.tab.3$vibrant, 
                                 levels = c('NA',  'low quality draft','medium quality draft', 'high quality draft'))
complete.tab.3$checkv <- factor(complete.tab.3$checkv, 
                                levels = c('Not-determined',  'Low-quality','Medium-quality', 'High-quality', 'Complete'))





p1 <- ggplot(complete.tab, aes(x=vibrant, y=phaster, fill=Freq)) +
  geom_tile() +
  theme_cowplot() +
  geom_text(aes(label = round(Freq, 1)), size=5) +
  scale_fill_gradient2(low='#dadaeb', high='#54278f') +
  labs(x='Vibrant', y='Phaster') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.1), legend.position = 'none', plot.margin = unit(c(0.1,0.35,0.1,0.1), 'in'))


p2 <- ggplot(complete.tab.2, aes(x=phaster, y=checkv, fill=Freq)) +
  geom_tile() +
  theme_cowplot() +
  geom_text(aes(label = round(Freq, 1)), size=5) +
  scale_fill_gradient2(low='#dadaeb', high='#54278f') +
  labs(x='Phaster', y='CheckV') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.1), legend.position = 'none', plot.margin = unit(c(0.1,0.35,0.1,0.1), 'in'))

p3 <- ggplot(complete.tab.3, aes(x=vibrant, y=checkv, fill=Freq)) +
  geom_tile() +
  theme_cowplot() +
  geom_text(aes(label = round(Freq, 1)), size=5) +
  scale_fill_gradient2(low='#dadaeb', high='#54278f') +
  labs(x='Vibrant', y='CheckV') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.1), legend.position = 'none', plot.margin = unit(c(0.1,0.35,0.1,0.1), 'in'))


plot_grid(p1,p2,p3, nrow=1)
#ggsave(filename = "plots/quality_revision.pdf", width = 15, height= 4.5, units = "in")
#ggsave(filename = "plots/quality_revision.png", width = 15, height= 4.5, units = "in")


#checkV clustering
checkv_contigs99 <- read.csv("raw_data/checkV/cluster/my_clusters_99.tsv", sep='\t', header = F)
checkv_contigs95 <- read.csv("raw_data/checkV/cluster/my_clusters_95.tsv", sep='\t', header = F)
checkv_contigs95$V3 <- str_replace_all(checkv_contigs95$V2, ',', '_')
checkv_contigs95$V3 <- paste(checkv_contigs95$V3, '_', sep="")
checkv_contigs99$V3 <- str_replace_all(checkv_contigs99$V2, ',', '_')
checkv_contigs99$V3 <- paste(checkv_contigs99$V3, '_', sep="")
comb_tab.phage$checkv_clust95 <- NA
comb_tab.phage$checkv_clust99 <- NA
for(i in 1:nrow(comb_tab.phage)){
  find <- paste(comb_tab.phage[i,]$Genome,'_', sep="")
  c95 <- checkv_contigs95[str_detect(checkv_contigs95$V3, find),]$V2
  c99 <- checkv_contigs99[str_detect(checkv_contigs99$V3, find),]$V2
  if(str_detect(c95, ",")){
    comb_tab.phage[i,]$checkv_clust95 <- c95
  }else{
    comb_tab.phage[i,]$checkv_clust95 <- NA
  }
  if(str_detect(c99, ",")){
    comb_tab.phage[i,]$checkv_clust99 <- c99
  }else{
    comb_tab.phage[i,]$checkv_clust99 <- NA
  }
}

vhulk <- read.csv('raw_data/vHULK_results.csv')  ###only runs on contigs >5kb
vhulkShort <- vhulk %>%
  select(BIN.genome, entropy)
names(vhulkShort) <- c('BIN.genome', 'vHulkEntropy')

comb_tab.phage2 <- merge(comb_tab.phage, vhulkShort, by.x='Genome', by.y='BIN.genome', all.x = T)
#write.table(comb_tab.phage2, file='revision/SI_table2_revision.tsv', sep="\t", quote = F, row.names = F)
add <- names(comb_tab.phage2)[which(!names(comb_tab.phage2) %in% names(comb_tab))]
comb_tab2 <- comb_tab[which(!comb_tab$Genome %in% comb_tab.phage2$Genome),]
comb_tab2$checkV_quality <- NA
comb_tab2$checkV_completeness <- NA
comb_tab2$checkV_warning <- NA
comb_tab2$checkv_clust95 <- NA
comb_tab2$checkv_clust99 <- NA
comb_tab2$vHulkEntropy <- NA
table(names(comb_tab.phage2) %in% names(comb_tab2))
comb_tab2 <- comb_tab2[,names(comb_tab.phage2)]
comb_tab2 <- rbind(comb_tab.phage2, comb_tab2)
#write.table(comb_tab2, file='revision/SI_table2_revision.tsv', sep="\t", quote = F, row.names = F)

###remove potentially duplicate contigs from list of phages
if(FALSE){
  get <- read.csv('revision/compressed_phage_contigs.txt', header = F, stringsAsFactors = F)
  for(i in 1:nrow(get)){
    get[i,] <- as.character(str_replace(get[i,], ">", ""))
  }
  
  contigsBefore <- read.fasta('raw_data/checkV/contigs_plasmids_phage.fa')
  find <- which(names(contigsBefore) %in% as.character(get$V1))
  contigsAfter <- contigsBefore[find]
  for(i in 1:length(contigsAfter)){
    write.fasta(paste(getSequence(contigsAfter[i])[[1]], collapse=''), names = names(contigsAfter)[i], file.out = 'revision/compressed_phage_contigs.fna', open = "a", nbchar = 60, as.string = FALSE)
  }
}


###trying to fix clustering figures and having trouble. redo cluster membership identfication and look for mess up
##this is the old run of vContact2!!!
#commmented out everyting
# eurocomp <- read.csv("analysis/vContact2_kbase_euro.csv")
# eurocomp$source <- NA
# eurocomp[str_detect(eurocomp$Genome, "JAAOBB"),]$source <- "Swiss"
# eurocomp[str_detect(eurocomp$Genome, "MN85"),]$source <- "Belgian"
# eurocomp[str_detect(eurocomp$Genome, "contig_|circular|scaffold_"),]$source <- "Texas"
# eurocomp[is.na(eurocomp$source),]$source <- "RefSeq"
# 
# comb_tab.phage2$VC_dbcheck <- NA
# comb_tab.phage2$VC_dbcheck_members <- NA
# for(i in 1:nrow(comb_tab.phage2)){
#   VC <- eurocomp[eurocomp$Genome==comb_tab.phage2[i,]$Genome,]$VC 
#   if(length(VC)==0){
#     VC <- "NaN"
#   }
#   comb_tab.phage2[i,]$VC_dbcheck <- VC
#   if(VC=="NaN"){
#     comb_tab.phage2[i,]$VC_dbcheck_members <- "NaN"
#   }else{
#     matches <- paste(eurocomp[eurocomp$VC==VC,]$Genome, collapse = ",")
#     comb_tab.phage2[i,]$VC_dbcheck_members <- matches
#   }
#   
# }
#write.table(comb_tab.phage2, file='revision/SI_table2_revision.tsv', sep="\t", quote = F, row.names = F)

  
  
  
