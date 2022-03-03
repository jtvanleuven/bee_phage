library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(readxl)


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


###prepare vContact2 genome-association file
prots <- read.fasta("raw_data/assembly_flye_1000_filtered.phages_combined.faa", as.string = T)
contigs <- paste(str_split(names(prots),"_",simplify = T)[,1],
                  str_split(names(prots),"_",simplify = T)[,2],sep="_")
contigs <- str_replace_all(contigs,",","")
header <- str_split(names(prots),"\t",simplify = T)[,1]
header <- str_replace(header, ",", "")
name <- str_split(names(prots),"\t",simplify = T)[,5]
name <- str_replace_all(name,",","")
name <- str_replace_all(name,";","")
name <- str_replace_all(name,"\\|","_")
name <- str_replace(name,"\"","")
name <- str_replace(name,"\\.","_")
tab <- cbind(header, contigs, name)
#write.table(tab, "analysis/vibrant_genes4vcontact2.csv", quote = F, row.names = F, col.names = F, sep=",")

prots2 <- prots
names(prots2) <- header
#write.fasta(prots2, names = header, file.out="analysis/assembly_flye_1000_filtered.phages_renamed.faa", as.string = T)

#pull out only contigs identified by at least 3 programs
upset_tab <- read.table("analysis/upset_table.tab", sep = "\t", header = T)
upset_tab[,25:33] <- upset_tab[,25:33] * 1
upset_tab$num_progs <- rowSums(upset_tab[,25:33])
viral_contigs <- upset_tab[upset_tab$num_progs > 3,]$scaffold

tab2 <- data.frame(tab, stringsAsFactors = F)
tab2 <- tab2[tab2$contigs %in% viral_contigs,]
#write.table(tab2, "analysis/vibrant_genes4vcontact2.csv", quote = F, row.names = F, col.names = F, sep=",")
#write.table(tab2, "../vContact2_try3/vibrant_genes4vcontact2.csv", quote = F, row.names = F, col.names = F, sep=",")

get <- which(names(prots2) %in% tab2$header)
prots.viral <- prots2[get]
#write.fasta(prots.viral, names = names(prots.viral), file.out="analysis/viral_prots.faa", as.string = T)
#vcontact2 --raw-proteins viral_prots.faa --rel-mode 'Diamond' --proteins-fp vibrant_genes4vcontact2.csv --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/james/miniconda3/envs/vContact2/bin/

contigs <- read.fasta("../assembly_flye_1000.fasta", as.string = T)
get <- which(names(contigs) %in% viral_contigs)
contigs.viral <- contigs[get]
#write.fasta(contigs.viral, names = names(contigs.viral), file.out = "analysis/viral_contigs.fna", as.string = T)




###combined viral contigs identified by >3 programs (includes plasmids) with european phage contigs.
##ran all though vContact2 using kbase
#import output file and compare
eurocomp <- read.csv("analysis/vContact2_kbase_euro.csv")
eurocomp_bee <- eurocomp[str_detect(eurocomp$Genome, "contig|JAAOBB|MN85|circular|scaffold"),]
eurocomp_bee$source <- NA
eurocomp_bee[str_detect(eurocomp_bee$Genome, "JAAOBB"),]$source <- "Swiss"
eurocomp_bee[str_detect(eurocomp_bee$Genome, "MN85"),]$source <- "Belgian"
eurocomp_bee[str_detect(eurocomp_bee$Genome, "contig|circular|scaffold"),]$source <- "Texas"


venn_list <- list(Swiss = eurocomp_bee[eurocomp_bee$source == "Swiss",]$VC,
                  Belgian = eurocomp_bee[eurocomp_bee$source == "Belgian",]$VC,
                  Texan = eurocomp_bee[eurocomp_bee$source == "Texas",]$VC)
library(VennDiagram)
venn.diagram(venn_list, filename = "plots/venn.tiff")

##need to deal with overlap
#this is actually a little tricky. Would like to have the numbers in the venn actually
#represent the number of contigs, but this gets complicated when clusters are linked by overlapping clusters
#going to try a heatmap (ggtile) that links each contig together
#clusters are links. x and y axes are contigs
all_phage <- read.fasta("analysis/all_phage.fasta")
contigs <- unique(eurocomp_bee$Genome)
##missing 1 contigs
names(all_phage)[!names(all_phage) %in% contigs]
length(all_phage$circular_read_641)
##it is only 2267bp in length. probably not a phage


##will be painful, but loop over everything and find connections
eurocomp_bee$overlap <- str_split(eurocomp_bee$VC.Status, " ", simplify = T)[,2]
eurocomp_bee$VC.Status <- str_split(eurocomp_bee$VC.Status, " ", simplify = T)[,1]
eurocomp_bee$VC.Status <- str_replace(eurocomp_bee$VC.Status,"Clustered/Singleton", "Singleton")
eurocomp_bee$overlap <- str_remove_all(eurocomp_bee$overlap, "\\(")
eurocomp_bee$overlap <- str_remove_all(eurocomp_bee$overlap, "\\)")
eurocomp_bee$overlap <- str_replace_all(eurocomp_bee$overlap, "/", ",")
table(eurocomp_bee$VC.Status)
plot.bar.tab <- eurocomp_bee %>%
  select(Genome,source,VC.Status) %>%
  group_by(source) %>%
  count(VC.Status)

ggplot(plot.bar.tab, aes(x=source, y=n, fill=VC.Status)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c("Belgian", "Texas", "Swiss")) +
  theme_classic() +
  scale_fill_npg() +
  labs(x="", y="# viral contigs")
#ggsave("plots/vc_stat_bar.pdf", width = 4, height = 3.5)
#ggsave("plots/vc_stat_bar.png", width = 4, height = 3.5)

bee_clusters <- eurocomp_bee[eurocomp_bee$VC.Status %in% c("Clustered","Overlap"),]
contigs <- unique(bee_clusters$Genome)
heat.tab <- data.frame(matrix(nrow=length(contigs), ncol = length(contigs)))
names(heat.tab) <- contigs
row.names(heat.tab) <- contigs
bee_clusters$VC <- str_split(bee_clusters$VC, "_", simplify = T)[,1]
bee_clusters$VC <- paste('VC_', bee_clusters$VC, sep='')
bee_clusters[bee_clusters$VC=='VC_NaN', ]$VC <- 'NaN'

for(i in 1:nrow(heat.tab)){
  contig <- names(heat.tab)[i]
  VC <- bee_clusters[bee_clusters$Genome == contig,]$VC
  overlap <- bee_clusters[bee_clusters$Genome == contig,]$overlap
  if(!VC=='NaN'){
    find <- bee_clusters[bee_clusters$VC==VC,]$Genome
    heat.tab[i,find] <- 3
    find.overlap <- bee_clusters[str_detect(bee_clusters$overlap, VC),]$Genome
    heat.tab[i,find.overlap] <- 2
  }
  if(!overlap==''){
    overlap <- str_split(overlap,",")[[1]]
    find <- bee_clusters[bee_clusters$VC %in% overlap,]$Genome
    heat.tab[i,find] <- 2
    for(j in overlap){
      find <- bee_clusters[str_detect(bee_clusters$overlap, j),]$Genome
      heat.tab[i,find] <- 1
    }
  }
}

heat.tab$genome <- row.names(heat.tab)
heat.long <- melt(heat.tab, id.vars = c('genome'))
heat.long <- heat.long[!heat.long$genome == heat.long$variable,]
heat.long$value <- as.factor(heat.long$value)
##do some clustering
heat.wide <- dcast(heat.long, genome~variable, value.var='value')
row.names(heat.wide) <- heat.wide$genome
heat.wide <- heat.wide[,-1]
heat.wide <- apply(heat.wide, 2, as.numeric)
heat.wide[is.na(heat.wide)] <- 0
col.order <- hclust(dist(heat.wide))$order
row.order <- hclust(dist(t(heat.wide)))$order
heat.wide <- heat.wide[row.order, col.order]

heat.long$genome <- factor(heat.long$genome, levels = colnames(heat.wide))
heat.long$variable <- factor(heat.long$variable, levels = colnames(heat.wide))

ggplot(heat.long, aes(x=genome, y=variable, fill=value)) +
  geom_tile() +
  #labs(x='',y='') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title=element_blank()) +
  scale_fill_manual(values=c('blue','red','purple'), na.value="white")
  
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
tmp <- as.matrix(get_upper_tri(heat.wide))
tmp$genome <- colnames(tmp)
heat.long <- melt(tmp, id.vars = c('genome'))
heat.long <- heat.long[!heat.long$genome == heat.long$variable,]
heat.long$value <- as.factor(heat.long$value)
heat.long$genome <- factor(heat.long$genome, levels = colnames(heat.wide))
heat.long$variable <- factor(heat.long$variable, levels = colnames(heat.wide))



###prep seqs for CRISPR search
crisprIndex <- read_xlsx('raw_data/blastdb/HoneyBee-Virome-2020/pnas.2000228117.sd04.xlsx', sheet = 1)
crisprSeqs <- read_xlsx('raw_data/blastdb/HoneyBee-Virome-2020/pnas.2000228117.sd04.xlsx', sheet = 2)
for(i in 1:nrow(crisprSeqs)){
  #write.fasta(crisprSeqs[i,]$sequence, names = crisprSeqs[i,]$sequenceID, file.out = 'raw_data/blastdb/HoneyBee-Virome-2020/crisprseqs.fasta', open = 'a')
}
#blastn -db crisprseqs.fasta -query ../../../analysis/contigs_plasmids_phage.fa -ungapped -dust no -soft_masking false -outmft 6

hits <- read.table('raw_data/blastdb/HoneyBee-Virome-2020/contigs_plasmids_phage_vs_crispr.blastn')
hits$host <- 'na'
for(i in 1:nrow(hits)){
  find <- hits[i,]$V2
  fetchRefContig <- crisprSeqs[crisprSeqs$sequenceID == find,]$genomeID
  fetchRefName <- unique(crisprIndex[which(str_detect(crisprIndex$Contig, fetchRefContig)),]$Strain)
  if(length(fetchRefName) > 1){
    print('found more than one bacterial genome for crisprseq')
  }else{
    hits[i,]$host <- fetchRefName
  }
}
write.csv(hits, 'analysis/crisprBlast.csv', quote = F, row.names = F)



###################compare to DASH crispr database
#blastn -db raw_data/crispr/SpacersDB.fasta -query analysis/contigs_plasmids_phage.fa -ungapped -dust no -soft_masking false -perc_identity 100 -outfmt 6 -num_threads 4 > analysis/crisprDash.blastn

##################search for 16S using phyloflash
#When running phyloFlash, please provide the location of
#the databases with the following option:
#-dbhome /mnt/ceph/jvanleuven/138.1

#################use sortmeRNA to find 16s
#sortmerna will not run with reads >30000
#seqkit seq -M 300000 -g pacbio.fasta > pacbio30000.fasta
#sortmerna --ref ~/138.1/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.fasta --reads ~/bees/pacbio30000.fasta

#################used metaQUAST to compare assemblies, ended up being somewhat useful for identifiying assembly of bacteria
#open up the results from pacbio-noWGA and pull out bacterial contigs
bac_contigs <- read.table("raw_data/metaQUAST/alignments_flye.tsv", fill = T)
#drop ecoli contamination
bac_contigs <- bac_contigs[!str_detect(bac_contigs$V1, 'Escherichia_coli_'),]
bac_contigs$V1 <- str_split(bac_contigs$V1, "_",simplify = T)[,1]
bac_contigs_t <- data.frame(t(bac_contigs[]))
names(bac_contigs_t) <- bac_contigs_t[1,]
bac_contigs_t <- bac_contigs_t[-1,]
bac_contigs_t$dummy <- 1:nrow(bac_contigs_t)
dat <- data.frame(matrix(nrow=1, ncol=2))
for(i in 1:ncol(bac_contigs_t)){
  contigs <- na.omit(bac_contigs_t[,i])
  dat <- rbind(dat, cbind(rep(names(bac_contigs_t)[i], length(bac_contigs_t[,i])), bac_contigs_t[,i]))
}

