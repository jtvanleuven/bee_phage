library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(readxl)
library(dplyr)


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

#write.table(comb_tab,"analysis/phaster_summary2.csv", quote = F, row.names = F, sep = "\t")
comb_tab <- read.csv("analysis/phaster_summary2.csv", sep="\t")



#ran plasmids though what the phage and phaster. many had predicted phage so better add them. 
#only 4 found by phaster
phaster_p <- c('circular_pair_124', 'circular_pair_3549', 'circular_pair_1029', 'circular_pair_240')

deepvir <- read.csv("analysis/wtp_plasmids/deepvirfinder.txt", header = F)
metaphinder <- read.csv("analysis/wtp_plasmids/metaphinder.txt", header = F)
PPRmeta <- read.csv("analysis/wtp_plasmids/PPRmeta.txt", header = F)
vibrant <- read.csv("analysis/wtp_plasmids/vibrant.txt", header = F)
virfinder <- read.csv("analysis/wtp_plasmids/virfinder.txt", header = F)
virsorter_vir <- read.csv("analysis/wtp_plasmids/virsorter-virome.txt", header = F)
virsort <- read.csv("analysis/wtp_plasmids/virsorter.txt", header = F)
vibrant_vir <- read.csv("analysis/wtp_plasmids/vibrant-virome.txt", header = F)

find_phage <- list(deepvir=deepvir$V1, 
                   PPRmeta=PPRmeta$V1, 
                   vibrant=vibrant$V1, 
                   virfinder=virfinder$V1,
                   virsorter_virome=virsorter_vir$V1,
                   virsorter=virsort$V1,
                   vibrant_virome=vibrant_vir$V1,
                   phaster=phaster_p)
library(UpSetR)
upset(fromList(find_phage), order.by="freq", nsets = 9)

plasmid_tab <- data.frame(matrix(ncol=9, nrow=length(unique(unname(unlist(find_phage))))))
names(plasmid_tab) <- c('plasmid', 'deepvir', 'PPRmeta', 'vibrant', 'virfinder', 'virsorter_virome', 'virsorter', 'vibrant_virome', 'phaster')
plasmid_tab$plasmid <- unique(unname(unlist(find_phage)))
plasmid_tab$deepvir <- plasmid_tab$plasmid %in% deepvir$V1  
plasmid_tab$PPRmeta <- plasmid_tab$plasmid %in% PPRmeta$V1 
plasmid_tab$vibrant <- plasmid_tab$plasmid %in% vibrant$V1 
plasmid_tab$virfinder <- plasmid_tab$plasmid %in% virfinder$V1 
plasmid_tab$virsorter_virome <- plasmid_tab$plasmid %in% virsorter_vir$V1 
plasmid_tab$virsorter <- plasmid_tab$plasmid %in% virsort$V1 
plasmid_tab$vibrant_virome <- plasmid_tab$plasmid %in% vibrant_vir$V1
plasmid_tab$phaster <- plasmid_tab$plasmid %in% phaster_p

methods <- c('deepvir', 'PPRmeta', 'vibrant', 'virfinder', 'virsorter_virome', 'virsorter', 'vibrant_virome', 'phaster')
get <- which(apply(plasmid_tab[,methods], 1, any)) 

plasmids <- read.fasta("analysis/polished_1.fasta")
plasmids_phages <- plasmids[names(plasmids) %in% plasmid_tab$plasmid]
#write.fasta(plasmids_phages, names=names(plasmids_phages), file.out = 'analysis/plasmid_phage.fa')  

run <- F
all_phage <- read.fasta("analysis/contigs_plasmids_phage.fa")
if(run){
  for(i in 1:length(all_phage)) {
    if(length(all_phage[[i]]) > 5000){
      #write.fasta(all_phage[i], names = names(all_phage)[i], file.out = paste("analysis/vHULK/", names(all_phage)[i], ".fasta", sep=""))
    }
  }
}


##copied code from above to add plasmids to overview table.
##had to change the code a bit. When Taylor ran phaster, i think that she provided individual contigs.
###when I ran it on plasmids, I provided the whole multi-fasta file.
#phaster_dirs <- list.files("raw_data/phaster/")
headerline <-  c("CONTIG", "REGION" ,"REGION_LENGTH" ,"COMPLETENESS(score)","SPECIFIC_KEYWORD","REGION_POSITION","TRNA_NUM","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM", "HYPOTHETICAL_PROTEIN_NUM","PHAGE+HYPO_PROTEIN_PERCENTAGE","BACTERIAL_PROTEIN_NUM","ATT_SITE_SHOWUP","PHAGE_SPECIES_NUM","MOST_COMMON_PHAGE_NAME(hit_genes_count)","FIRST_MOST_COMMON_PHAGE_NUM","FIRST_MOST_COMMON_PHAGE_PERCENTAGE","GC_PERCENTAGE")
names(phaster_summaries) <- headerline
#contig <- readLines(file(paste("raw_data/phaster/",phaster_dirs[i], "/detail.txt", sep=""),"r"),n=1)
sum <- read.table("raw_data/phaster_plasmids/summary.txt", sep="", fill = F, quote="", stringsAsFactors = F, skip = 34)
#sum <- cbind(rep(contig,nrow(sum)),sum)
sum$V1 <- str_split(sum$V5, ":", simplify = T)[,1]
sum <- cbind(sum[,1], sum)
sum$CONTIG <- paste('>', sum$CONTIG, sep="")
names(sum) <- headerline
sum$REGION_POSITION <- str_split(sum$REGION_POSITION, ":", simplify = T)[,2]
sum$REGION <- 1
#write.table(sum,"analysis/phaster_summary_plasmids.csv", quote = F, row.names = F, sep = "\t")
phaster_summaries_short <- sum[,c("CONTIG","COMPLETENESS(score)", "SPECIFIC_KEYWORD", "TOTAL_PROTEIN_NUM", "PHAGE_HIT_PROTEIN_NUM", "ATT_SITE_SHOWUP", "MOST_COMMON_PHAGE_NAME(hit_genes_count)")]
phaster_summaries_short$CONTIG <- str_remove(phaster_summaries_short$CONTIG, ">")

flye_info <- read.csv("raw_data/assembly_info_flye_plasmids.txt", sep = "\t", stringsAsFactors = F)
names(flye_info) <- c("scaffold","length","cov")
vibrant_pred <- read.csv("raw_data/wtp_plasmids/VIBRANT_machine_polished_1_filtered.tsv", sep="\t", stringsAsFactors = F)
row.names(vibrant_pred) <- vibrant_pred$scaffold
vibrant_qual <- read.csv("raw_data/wtp_plasmids/VIBRANT_genome_quality_polished_1_filtered.tsv", sep="\t", stringsAsFactors = F)
vibrant_qual$scaffold <- str_replace(vibrant_qual$scaffold,"_fragment_\\d","")
vibrant_qual_smash <- aggregate(.~scaffold, vibrant_qual, paste, collapse=",")

comb_tab <- merge(flye_info, vibrant_pred, by="scaffold", all=T)
comb_tab <- merge(comb_tab, vibrant_qual_smash, by="scaffold", all=T)
comb_tab_plasmids <- merge(comb_tab, phaster_summaries_short, by=1, all=T)
#write.table(comb_tab_plasmids, "analysis/phaster_summary_plasmids.csv", row.names = F, quote = F, sep="\t")


#join plasmids and contigs and add vContact2
comb_tab <- read.csv("analysis/phaster_summary2.csv", sep="\t")
comb_tab <- comb_tab[,c('scaffold','length','cov','circ','prediction','type',"Quality","COMPLETENESS.score.","SPECIFIC_KEYWORD","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM","ATT_SITE_SHOWUP", "MOST_COMMON_PHAGE_NAME.hit_genes_count.")]

comb_tab_plasmids <- read.csv("analysis/phaster_summary_plasmids.csv", sep="\t")
comb_tab_plasmids$circ <- 'Y'
comb_tab_plasmids <- comb_tab_plasmids[,c('scaffold','length','cov','circ','prediction','type',"Quality","COMPLETENESS.score.","SPECIFIC_KEYWORD","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM","ATT_SITE_SHOWUP", "MOST_COMMON_PHAGE_NAME.hit_genes_count.")]

comb_tab_all <- rbind(comb_tab, comb_tab_plasmids)
names(comb_tab_all) <- c('Genome','length','cov','circ','prediction','type',"Quality","COMPLETENESS.score.","SPECIFIC_KEYWORD","TOTAL_PROTEIN_NUM","PHAGE_HIT_PROTEIN_NUM","ATT_SITE_SHOWUP", "MOST_COMMON_PHAGE_NAME.hit_genes_count.")

vhulk <- read.csv('raw_data/vHULK_results.csv')
vhulk.s <- vhulk[,c('BIN.genome','final_prediction')]
names(vhulk.s) <- c('Genome', 'vHulk.pred')
comb_tab_all <- merge(comb_tab_all, vhulk.s, by = 1, all=T)
genome2genome <- read.csv('raw_data/vContact2/genome_by_genome_overview.csv')
genome2genome.swiss <- genome2genome[str_detect(genome2genome$Genome, 'JAA'),]
table(genome2genome.swiss$Order)
genome2genome.belg <- genome2genome[str_detect(genome2genome$Genome, 'MN85'),]
table(genome2genome.belg$Order)
genome2genome.short <- genome2genome[genome2genome$Genome %in% comb_tab_all$Genome,c('Genome','VC', 'VC.Status', 'Size', 'VC.Subcluster', 'VC.Subcluster.Size', 'Quality')]
table(genome2genome.short$Order)
table(genome2genome.short$Genus)

clusters <- read.csv('raw_data/vContact2/viral_cluster_overview.csv')
clusters[str_detect(clusters$Members, 'contig_1897'),]  ##clustered
clusters[str_detect(clusters$Members, 'contig_2882'),]  ##overlapping
genome2genome.short[genome2genome.short$Genome=='contig_2882',]
clusters[clusters$VC == 'VC_198_0',]
clusters[clusters$VC == 'VC_222_0',]

tmp <- merge(comb_tab_all, genome2genome.short, by='Genome', all = T)
tmp <- tmp[,!names(tmp)=='ATT_SITE_SHOWUP']
names(tmp) <- c("Genome", "Flye.length", "Flye.cov", "Flye.circ" , "Vibrant.prediction", 
                "Vibrant.type", "Vibrant.Quality", "Phaster.COMPLETENESS.score", 
                "Phaster.SPECIFIC_KEYWORD", "Phaster.TOTAL_PROTEIN_NUM", "Phaster.PHAGE_HIT_PROTEIN_NUM", 
                "Phaster.MOST_COMMON_PHAGE_NAME.hit_genes_count", "vHulk.pred", 
                "VC", "VC.Status", "VC.Size", "VC.Subcluster", "VC.Subcluster.Size", "VC.Quality")
tmp$VC.members <- NA
tmp[which(tmp$VC.Subcluster == ""),]$VC.Subcluster <- NA
for(i in 1:nrow(tmp)){
  if(!is.na(tmp[i,]$VC.Subcluster == "")){
    tmp[i,]$VC.members <- clusters[which(clusters$VC == tmp[i,]$VC.Subcluster),]$Members
  }
}

phage <- read.fasta('analysis/all_phage.fasta')
tmp$IsPhage <- tmp$Genome %in% names(phage)
###not sure what this contig is. can't find it in assembly. remove it. contig_12251
tmp <- tmp[!tmp$Genome=='contig_12251',]

##add crispr blast
hits <- read.csv('analysis/crisprBlast.csv')
tmp$crisprhit <- 'NA'
for(i in 1:nrow(tmp)){
  found <- hits[hits$V1 == tmp[i,]$Genome,]
  species <- unique(str_extract(found$host, '\\w+\\s\\w+'))
  if(length(species)==1){
    tmp[i,]$crisprhit <- species
  }else if(length(species) > 1){
    genus <- unique(str_extract(species, '\\w+'))
    if(length(genus) == 1){
      tmp[i,]$crisprhit <- genus
    }else if(length(genus) == 1){
      print(paste(i, 'more than one genus', sep=' '))
    }
  }
}

##add info from belgian group
suptab <- read.csv('raw_data/blastdb/matthijnssen_Supptable_S18.csv', header = T)
tmp$mattInfo <- 'NA'
mattHits <- which(str_detect(tmp$VC.members, 'MN85'))
for(i in mattHits){
  get <- str_extract_all(tmp[i,]$VC.members,'MN85\\d+', simplify = T)
  fetch <- suptab[suptab$Genbank.accession %in% as.character(get), c('Cluster.Name', 'Host_crisprdb', 'Length')]
  fix <- paste(fetch, collapse = ',')
  fix <- str_remove_all(fix, 'c\\(')
  fix <- str_remove_all(fix, '\\)')
  fix <- str_remove_all(fix, '\\"')
  tmp[i,]$mattInfo <- fix
}                
           

##add info from the swiss group
swissSeqs <- read.fasta('raw_data/blastdb/JAAOBB01.1.engel.fasta')

swissInfo <- read_xlsx('raw_data/blastdb/HoneyBee-Virome-2020/pnas.2000228117.sd02.xlsx', skip = 2)
swissInfo <- swissInfo[-1,1:6]
names(swissInfo) <- c("VC", "Predicted Host", "Predicted LifeStyle", "Taxonomic Classification", "type","Total Length")

swissNames <- read.table('raw_data/blastdb/JAAOBB01.1.engel.fasta', fill=T, header=F, sep=';')
swissNames <- data.frame(swissNames[grep("^>", swissNames$V1),])
names(swissNames) <- 'header'
swissNames$contig <- str_extract(swissNames$header, 'JAAOBB\\d+')
swissNames$VC <- str_extract(swissNames$header, 'VC_\\d+')
swissNames[is.na(swissNames$VC),]$VC <- str_extract(swissNames[is.na(swissNames$VC),]$header, '\\w\\w7_\\d+')
#one cluter with funky naming convention
swissNames[swissNames$contig=='JAAOBB010000140',]$VC <- 'VC_6_1'

tmp$bonillaInfo <- 'NA'
bonillaHits <- which(str_detect(tmp$VC.members, 'JAAOBB'))
for(i in bonillaHits){
  get <- str_extract_all(tmp[i,]$VC.members,'JAAOBB\\d+', simplify = T)
  fetch <- unique(swissNames[swissNames$contig %in% as.character(get),]$VC)
  fix <- paste(fetch, collapse = ',')
  info <- swissInfo[swissInfo$VC %in% fix, c('Predicted Host', 'Predicted LifeStyle', 'Taxonomic Classification', 'Total Length')]
  info2 <- paste(info, collapse = ",")
  tmp[i,]$bonillaInfo <- info2
}    

#add DASH crispr info
#error in running
#(CrisprOpenDB_env) james@james-Precision-5530:~/Software/CrisprOpenDB/$ python CL_Interface.py -i contigs_plasmids_phage.fasta -m 0 -t -n 8 > contig_plasmids_phage_predictions.txt
#have to fix. got 451 contigs run through, but would not compile summary
#try to run on these 451
allcontigs <- read.fasta('~/Software/CrisprOpenDB/old/contigs_plasmids_phage.fasta', as.string = T)
run <- read.csv('~/Software/CrisprOpenDB/try1/done_list.txt', header = F)
run$V1 <- str_replace(run$V1, '.csv', '')
trimcontigs <- allcontigs[which(names(allcontigs) %in% run$V1)]
#write.fasta(sequences = trimcontigs, names = names(trimcontigs), file.out = '~/Software/CrisprOpenDB/contigs_run.fasta', as.string = T)

ntrun <- allcontigs[which(!names(allcontigs) %in% run$V1)]
#write.fasta(sequences = ntrun, names = names(ntrun), file.out = '~/Software/CrisprOpenDB/contigs_notrun.fasta', as.string = T)
#for(i in 1:length(ntrun)){
#  write.fasta(sequences = ntrun[i], 
#              names = names(ntrun[i]), 
#              file.out = paste('~/Software/CrisprOpenDB/', names(ntrun[i]), '.fasta', sep=""), 
#              as.string = T)
#}

genome_tab <- read.csv('analysis/genome_combine_summary.tab', sep='\t')
genome_tab <- genome_tab[!str_detect(names(genome_tab), 'CrisprDB')]
crisprDB <- read.csv('raw_data/crispr/results_CrisprOpenDB.txt', header = F)
crisprDB$V1 <- str_remove_all(crisprDB$V1, '\\(')
crisprDB$V1 <- str_remove_all(crisprDB$V1, '\'')
crisprDB$V2 <- str_remove_all(crisprDB$V2, '\'')
rownames(crisprDB) <- crisprDB$V1
crisprDBhits <- crisprDB[!str_detect(crisprDB$V2,'Sorry'),1:2]
names(crisprDBhits) <- c('Genome', "CrisprDB")
genome_tab <- merge(genome_tab, crisprDBhits, by='Genome', all.x = T)

write.table(genome_tab, file='analysis/genome_combine_summary.tab', quote = F, row.names = F, sep = "\t")



##############work on annotations using multiphate2
####real pain getting to work.
####key was to use -parse_seqids flag for makeblastdb
####otherwise ran into blast error
####still need to format contigs for running
###must be split and make this name list
##Genome 1 
##genome_file='Eb_P2.fasta'
##genome_type='phage'

contigs <- read.fasta('analysis/contigs_plasmids_phage.fa')
text <- vector(mode="character", length=length(contigs)*3)
for(i in 1:length(contigs)){
  header <- names(contigs)[i]
  #=(i-1)*3+1
  text[(i-1)*3+1] <- paste('Genome ', i, sep='')
  text[(i-1)*3+1+1] <- paste("genome_file=\'", header, ".fasta\'", sep='')
  text[(i-1)*3+1+2] <- "genome_type=\'phage\'"
  #write.fasta(contigs[i], names = names(contigs)[i], file.out = paste('analysis/multiphate2/', header, '.fasta', sep=''))
}
#write.table(text, 'analysis/multiphate2/genomes.txt', quote = F, col.names = F, row.names = F)

genome_tab <- read.table('analysis/genome_combine_summary.tab', header = T, sep = "\t")

contig_list <- list.dirs('analysis/multiphate2/', full.names = F)
contig_list <- unique(str_split(contig_list, '\\/', simplify = T)[,2])
contig_list <- contig_list[!contig_list=='']
multiphate <- vector(mode = "list", length = length(contig_list))
names(multiphate) <- contig_list

for(i in contig_list){
  annot <- read.gff(paste('analysis/multiphate2/PipelineOutput/', i, '/phate_sequenceAnnotation_main.gff', sep=''))
  annot <- annot %>%
    filter(type=='CDS')
  annot$prot <- NA
  annot$len <- abs(annot$end-annot$start)
  for(j in 1:nrow(annot)){
    genelist <- t(str_split(annot[j,]$attributes, ';', simplify = T))
    if(nrow(genelist) > 1){
      genelist <- str_split(genelist, 'protein', simplify = T)[,1]
      genelist <- str_split(genelist, '\\)', simplify=T)[,2]
      genelist <- str_remove_all(genelist, '\\[')
      genelist <- str_remove_all(genelist, '\\]')
      genelist <- str_squish(genelist)
      genelist <- genelist[-1] ##drop first one
      genelist <- paste(genelist, ' protein', sep='')
      genelist <- unique(genelist)
      annot_num <- length(genelist)
      genelist <- paste(genelist, collapse = ', ')
      annot[j,]$prot <- genelist
    }else{
      annot[j,]$prot <- NA
    }
  }
  multiphate[[i]] <- annot
}

multiphate.long <- matrix(nrow = 1, ncol=3)
for(i in 1:length(multiphate)){
  tmp <- cbind(rep(names(multiphate)[i], nrow(multiphate[[i]])),
               multiphate[[i]]$len, 
               multiphate[[i]]$prot)
  multiphate.long <- rbind(multiphate.long, tmp)
}
multiphate.long <- data.frame(multiphate.long[-1,])
names(multiphate.long) <- c('contig', 'protLen', 'protHits')
write.table(multiphate.long, 'analysis/multiphate_summary.tab', quote = F, row.names = F, sep = '\t')

multiphate.long$protLen <- as.numeric(multiphate.long$protLen)
hist(multiphate.long$protLen, breaks=50)
max(multiphate.long$protLen)
min(multiphate.long$protLen)
median(multiphate.long$protLen)
multiphate.long %>%
  count(protLen > 300)

multiphate.long$hypoflag <- str_detect(multiphate.long$protHits, 'hypothetical protein|Phage protein')
table(multiphate.long$hypoflag)
print('FALSE=annotated, TRUE=hypothetical')

annot.genes <- na.omit(multiphate.long[!multiphate.long$hypoflag,])

##29 phage have no functional annotation. 
length(unique(annot.genes$contig)) - length(unique(multiphate.long$contig))
##characterize these a little
darkgenomes <- unique(multiphate.long[!multiphate.long$contig %in% annot.genes$contig,]$contig)
darkgenomesAnnots <- multiphate.long[multiphate.long$contig %in% darkgenomes, ]
darkgenomesAnnots %>%
  count(contig)
#num genes from 3-31

darkgenomesAnnots %>%
  count(protLen > 300)

##common phage genes from 10.1038/s41467-020-18236-8
comm <- c('TAIL',
'TAIL FIBER',
'TAIL FIBRE',
'SCAFFOLDING PROTEIN',
'MAJOR CAPSID',
'MINOR CAPSID',
'\\w+ SPIKE',
'CAPSID',
'PORTAL',
'DNA-BINDING',
'DNA BINDING',
'TAPEMEASURE',
'TAPE MEASURE',
'REPLICATION',
'REPLICATION INITIATION',
'MORPHOGENESIS',
'PLASMID PROTEIN',
'\\w+YSIN',
'BASEPLATE',
'TAIL ASSEMBLY',
'HOLIN',
'RIBONUCLEOSIDE REDUCTASE',
'STRUCTURAL PROTEIN',
'RECOMBINATION',
'GDSL LIPASE', 
'TRASCRIPTIONAL REGULATOR',
'REPRESSOR',
'\\w+REPRESSOR',
'INTEGRON GENE CASSETTE',
'\\w+ESTRICTION',
'\\w+ PROTIEN',
'\\w+ASE',
'DNA \\w+ASE',
'\\w+ SYNTHASE',
'\\w+ ASSOCIATED PROTEIN',
'\\w+ REGULATOR',
'\\w+ TAIL \\w+',
'\\w+ TAIL',
'TAIL \\w+',
'\\w+ CAPSID')

functions <- data.frame(str_split(annot.genes$protHits, ',', simplify = T))
#annot.genes$concat <- apply(functions[,1:ncol(functions)], 1, paste, collapse=" ")

annot.genes$protHits <- str_replace(annot.genes$protHits, 'DNA-BINDING', 'DNA BINDING')
annot.genes$protHits <- str_replace(annot.genes$protHits, 'PBSX', 'PORTAL')

annot.genes$parseName <- NA
for(i in 1:nrow(annot.genes)){
  found <- unique(c(str_extract_all(toupper(annot.genes[i,]$protHits), comm, simplify = T)))
  found <- str_squish(found[!found==""])
  annot.genes[i,]$parseName <- paste(found, collapse = ', ')
}
annot.genes$parseName <- str_replace_all(annot.genes$parseName, ', ,', ',')
annot.genes$parseName <- str_replace_all(annot.genes$parseName, '\\|', ' ')
annot.genes$parseName <- str_replace_all(annot.genes$parseName, '_', ' ')
annot.genes$parseName <- str_squish(annot.genes$parseName)
write.table(annot.genes, 'analysis/multiphate2/PipelineOutput/simple_prot_list.tab', quote = F, row.names = F, sep='\t')

annot.genes <- read.csv('analysis/multiphate2/PipelineOutput/simple_prot_list.tab', sep='\t')
##must simplify
annot.genes$parseName <- str_squish(str_remove_all(annot.genes$parseName, '\\d+'))
annot.genes$simpName <- NA
for(i in 1:nrow(annot.genes)){
  hits <- str_squish(str_split(annot.genes[i,]$parseName, ',', simplify = T))
  if(length(hits) > 1){
    matrix <- matrix(nrow = length(hits), ncol=length(hits))
    for(j in 1:length(hits)){
      matrix[j,] <- str_detect(hits, hits[j])
    }
    matrix <- matrix*1
    get <- rowSums(matrix)-1 == 0
    annot.genes[i,]$simpName <- paste(hits[get], collapse = ', ')
  }else{
    annot.genes[i,]$simpName <- annot.genes[i,]$parseName
  }
}
annot.genes$parseName <- annot.genes$simpName
annot.genes <- annot.genes[,!names(annot.genes)=='simpName']

#lots of fixing
annot.genes$parseName <- str_replace(annot.genes$parseName, "RECOMBINASE, PHAGE INTEGRASE", "INTEGRASE")
annot.genes$parseName <- str_replace(annot.genes$parseName, "PHAGE INTEGRASE, RECOMBINASE", "INTEGRASE")
annot.genes$parseName <- str_replace(annot.genes$parseName, "PHAGE INTEGRASE", "INTEGRASE")
annot.genes$parseName <- str_replace(annot.genes$parseName, "RECOMBINASE", "INTEGRASE")
annot.genes$parseName <- str_replace(annot.genes$parseName, "TAPE MEASURE, PHAGE TAIL LENGTH", "TAPE MEASURE")
annot.genes$parseName <- str_replace(annot.genes$parseName, "SCAFFOLDING PROTEIN, PHAGE CAPSID", "SCAFFOLDING PROTEIN")
annot.genes$parseName <- str_replace(annot.genes$parseName, "AMIDASE, ENDOLYSIN", "ENDOLYSIN, AMIDASE") 
annot.genes$parseName <- str_replace(annot.genes$parseName, "23 TAIL FIBER, 25 TAIL FIBER, PHAGE TAIL FIBER", "TAIL FIBER")
annot.genes$parseName <- str_replace(annot.genes$parseName, "PHAGE TAIL FIBER, 21 TAIL FIBER", "TAIL FIBER")
annot.genes$parseName <- str_replace(annot.genes$parseName, "43 TAIL ASSEMBLY, PHAGE TAIL FIBER", "TAIL FIBER")
annot.genes[annot.genes$parseName == "PHAGE TAIL FIBER",]$parseName <- "TAIL FIBER"

sumtab <- table(annot.genes$parseName)
nrow(sumtab)
View(sumtab)

View(annot.genes[str_detect(toupper(annot.genes$protHits), 'BIOSYNTHESIS'),])



##########################analyze metaQUAST output for supplemental table 4
##read output file from metaQUAST
quast <- read.table("raw_data/metaQUAST/alignments_flye.tsv", fill = TRUE)
quast$V1 <- paste(str_split(quast$V1, "_", simplify = T)[,1], str_split(quast$V1, "_", simplify = T)[,2], sep="_")
quast.l <- melt(quast, id.vars = c("V1"))
quast.l <- quast.l[!quast.l$value=="", c("V1","value")]
table(quast.l$V1)  ##only 3 types of bacteria. Several species of each
##drop E. coli
quast.l <- quast.l[!quast.l$V1=="Escherichia_coli",]
names(quast.l) <- c('bact_hit', 'Genome')
#Lactobacillus, Bifidobacterium, and Gilliamella

###get contig coverage info
comb_tab <- read.table("analysis/genome_combine_summary.tab", sep="\t", header=T)             
comb_tab.phage <- comb_tab[comb_tab$IsPhage,]
comb_tab.phage$Vibrant.type <- str_replace(comb_tab.phage$Vibrant.type, 'lytic,lytic', 'lytic')
comb_tab.phage$Flye.length <- comb_tab.phage$Flye.length/1000
comb_tab.phage[which(is.na(comb_tab.phage$Vibrant.type)),]$Vibrant.type <- 'unidentified'
comb_tab.phage$Vibrant.type <- factor(comb_tab.phage$Vibrant.type, , levels=c('lytic','lysogenic','unidentified'))

tmp <- merge(quast.l, comb_tab, by="Genome")

##be careful here. quite a few temperate phages showing up here
tmp <- tmp %>%
  filter(IsPhage == FALSE) %>%
  filter(Flye.circ == 'N') %>%
  filter(Flye.length > 1000)
tmp <- tmp[is.na(tmp$Vibrant.type),]
tmp <- tmp[is.na(tmp$Phaster.COMPLETENESS.score),]

tmp$Flye.cov <- as.numeric(as.character(tmp$Flye.cov))
ggplot(tmp, aes(x=bact_hit, y=Flye.cov, colour=bact_hit)) +
  geom_violin() +
  geom_jitter(width = 0.2) +
  theme_classic()
  
ggplot(tmp, aes(x=Flye.cov, y=Flye.length, colour=bact_hit)) +
  geom_point() +
  theme_classic()


###this is a god damn mess. 
###tried aligning contigs to reference genomes using bwa mem
##parse sam file, filter by mapping quality and see if I get anywhere
#bwa mem refs/all.fna assembly_flye_combined.fasta > refs/assembly_flye_align_refs.sam
#samtools view -b -S assembly_flye_align_refs.sam > assembly_flye_align_refs.bam
#samtools sort -o assembly_flye_align_refs.sorted.bam assembly_flye_align_refs.bam

library(Rsamtools)
refs <- data.frame(acc=c('CP048272.1', 'NZ_JABZEC010000001.1', 'NZ_CP007445.1', 'CP009531.1', 'NZ_KQ034026.1'), 
                   name=c('Bifidobacterium asteroides strain ESL0447', 'Bombilactobacillus apium strain DCY120', 'Gilliamella apicola strain wkB1', 'Lactobacillus sp. wkB8', 'Lactobacillus melliventris strain Hma8'))       

aln <- scanBam("raw_data/refs/assembly_flye_align_refs.sorted.bam")
aln <- aln[[1]]
hist(aln$mapq, breaks=25)
