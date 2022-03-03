library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(ggsci)
library(dplyr)
library(UpSetR)
library(ggforce)
library(dplyr)
display.brewer.all(colorblindFriendly = TRUE)


comb_tab <- read.csv("analysis/phaster_summary2.csv", sep = "\t")
plot_comb <- comb_tab
plot_comb$type <- str_replace(plot_comb$type,",\\w+","")
plot_comb$cat <- paste(plot_comb$prediction, plot_comb$type, sep=",")
types <- table(plot_comb$cat)
plot_comb_short <- plot_comb[plot_comb$cat %in% names(types[types > 10]),]

ggplot(plot_comb_short, aes(x=length, y=cov, colour=match, shape=circ)) +
  geom_point(alpha=0.75, size=2) +
  facet_wrap(~ cat) +
  theme_bw() +
  scale_color_nejm()


ggplot(comb_tab, aes(x=length, y=cov, fill=prediction, colour=prediction)) +
  geom_point()

tmp <- comb_tab
tmp$prediction <- str_replace(tmp$prediction, 'organism,organism', 'organism')
tmp$prediction <- str_replace(tmp$prediction, 'organism,organism,organism', 'organism')
tmp$prediction <- str_replace(tmp$prediction, 'organism,organism,organism,organism', 'organism')
tmp$prediction <- str_replace(tmp$prediction, 'organism,organism,organism,plasmid,organism,organism', 'organism')
tmp$prediction <- str_replace(tmp$prediction, 'organism,plasmid', 'organism')
tmp[is.na(tmp$prediction),]$prediction <- 'unknown'
tmp[!tmp$prediction %in% c('unknown', 'virus', 'organism', 'plasmid'),] <- 'unknown'
tmp2 <- tmp %>%
  count(prediction)
ggplot(tmp2,aes(x=prediction, y=n))+
  geom_bar(stat='identity')



###simple lengths distribution plots 
comb_tab <- read.table("analysis/genome_combine_summary.tab", sep="\t", header=T)             
comb_tab.phage <- comb_tab[comb_tab$IsPhage,]
comb_tab.phage$Vibrant.type <- str_replace(comb_tab.phage$Vibrant.type, 'lytic,lytic', 'lytic')
comb_tab.phage$Flye.length <- comb_tab.phage$Flye.length/1000
comb_tab.phage[which(is.na(comb_tab.phage$Vibrant.type)),]$Vibrant.type <- 'unidentified'
comb_tab.phage$Vibrant.type <- factor(comb_tab.phage$Vibrant.type, , levels=c('lytic','lysogenic','unidentified'))
ggplot(comb_tab.phage, aes(y=Flye.length, x=Flye.cov, color=Vibrant.type)) +
  geom_point(size=2, alpha=0.5, shape=16) +
  theme_classic() +
  facet_wrap(~Vibrant.type, labeller=labeller(Vibrant.type=c("lytic" = "virulent", "lysogenic"="temperate", "unidentified"="unidentified"))) +
  labs(x='Read coverage', y='Contig length (kb)') +
  theme(legend.position = 'none', strip.background.x=element_rect(color = NA), strip.text = element_text(size=12)) +
  scale_color_brewer(palette='Dark2')
#ggsave(filename = "plots/cov_length.pdf", width = 6.5, height= 3, units = "in")
#ggsave(filename = "plots/cov_length.png", width = 6.5, height= 3, units = "in")


##compare completeness scores
drop <- c('\\d+','\\(','\\)')
comb_tab.phage$Phaster.COMPLETENESS.score <- str_remove_all(comb_tab.phage$Phaster.COMPLETENESS.score, paste(drop, collapse = '|'))
table(comb_tab.phage$Phaster.COMPLETENESS.score)
comb_tab.phage$Phaster.COMPLETENESS.score <- str_replace(comb_tab.phage$Phaster.COMPLETENESS.score, ',incomplete', '')
table(comb_tab.phage$Phaster.COMPLETENESS.score)
comb_tab.phage$Vibrant.Quality <- str_remove(comb_tab.phage$Vibrant.Quality, ',complete circular')

join <- paste(comb_tab.phage$Phaster.COMPLETENESS.score, comb_tab.phage$Vibrant.Quality, sep=',')
complete.tab <- data.frame(phaster=NA, vibrant=NA, table(join))
complete.tab$phaster <- str_split(complete.tab$join, ',', simplify = T)[,1]
complete.tab$vibrant <- str_split(complete.tab$join, ',', simplify = T)[,2]
complete.tab$Freq <- as.numeric(complete.tab$Freq) 
complete.tab$phaster <- factor(complete.tab$phaster, 
                               levels = c('NA',  'questionable','incomplete', 'intact'))
complete.tab$vibrant <- factor(complete.tab$vibrant, 
                               levels = c('NA',  'low quality draft','medium quality draft', 'high quality draft'))
ggplot(complete.tab, aes(x=vibrant, y=phaster, fill=Freq)) +
  geom_tile() +
  theme_cowplot() +
  geom_text(aes(label = round(Freq, 1)), size=5) +
  scale_fill_gradient2(low='#dadaeb', high='#54278f') +
  labs(x='Vibrant', y='Phaster') +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.1), legend.position = 'none', plot.margin = unit(c(0.1,0.35,0.1,0.1), 'in'))
#ggsave(filename = "plots/quality.pdf", width = 5, height= 5, units = "in")
#ggsave(filename = "plots/quality.png", width = 5, height= 5, units = "in")


##build new upsetR plot with phaster results
deepvir <- read.csv("raw_data/identified_contigs_by_tools/deepvirfinder.txt", header = F)
metaphinder <- read.csv("raw_data/identified_contigs_by_tools/metaphinder.txt", header = F)
PPRmeta <- read.csv("raw_data/identified_contigs_by_tools/PPRmeta.txt", header = F)
vibrant <- read.csv("raw_data/identified_contigs_by_tools/vibrant.txt", header = F)
virfinder <- read.csv("raw_data/identified_contigs_by_tools/virfinder.txt", header = F)
virsorter_vir <- read.csv("raw_data/identified_contigs_by_tools/virsorter-virome.txt", header = F)
virsort <- read.csv("raw_data/identified_contigs_by_tools/virsorter.txt", header = F)
vibrant_vir <- read.csv("raw_data/identified_contigs_by_tools/vibrant-virome.txt", header = F)
comb_tab <- read.csv("analysis/phaster_summary2.csv", sep="\t")
phaster <- comb_tab[!is.na(comb_tab$COMPLETENESS.score.),]$scaffold

find_phage <- list(deepvir=deepvir$V1, 
                   PPRmeta=PPRmeta$V1, 
                   vibrant=vibrant$V1, 
                   virfinder=virfinder$V1,
                   virsorter_virome=virsorter_vir$V1,
                   virsorter=virsort$V1,
                   vibrant_virome=vibrant_vir$V1,
                   phaster=phaster)

upset(fromList(find_phage), order.by="freq", nsets = 9)

devtools::unload('UpSetR')
library(ComplexUpset)
upset_tab <- comb_tab
upset_tab$deepvir <- upset_tab$scaffold %in% deepvir$V1
upset_tab$metaphinder <- upset_tab$scaffold %in% metaphinder$V1
upset_tab$PPRmeta <- upset_tab$scaffold %in% PPRmeta$V1
upset_tab$vibrant <- upset_tab$scaffold %in% vibrant$V1
upset_tab$virfinder <- upset_tab$scaffold %in% virfinder$V1
upset_tab$virsorter <- upset_tab$scaffold %in% virsort$V1
upset_tab$virsorter_virome <- upset_tab$scaffold %in% virsorter_vir$V1
upset_tab$vibrant_virome <- upset_tab$scaffold %in% vibrant_vir$V1
upset_tab$phaster <- upset_tab$scaffold %in% phaster
upset_tab[is.na(upset_tab$Quality),]$Quality <- "NA"

methods <- names(find_phage)
get <- which(apply(upset_tab[,methods], 1, any)) 
upset_tab <- upset_tab[get,]
upset_tab$Quality <- str_replace(upset_tab$Quality, "low quality draft,low quality draft", "low quality draft")
upset_tab$Quality <- str_replace(upset_tab$Quality, ",complete circular", "")
upset_tab$Quality <- factor(upset_tab$Quality, levels = c("high quality draft", "medium quality draft", "low quality draft", "NA"))

upset(
  data=upset_tab,
  intersect=methods,
  min_size=5,
  keep_empty_groups=FALSE,
)

upset(
  upset_tab,
  methods,
  min_size=5,
  dot_size = 1.5,
  keep_empty_groups=FALSE,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=F
    )
  ),
  annotations = list(
    'Genome Quality'=list(
      aes=aes(x=intersection, fill=Quality),
      geom=list(
        geom_bar(stat='count', position='fill'),
        scale_y_continuous(labels=scales::percent_format()),
        scale_fill_manual(values=c("grey10", "grey40", "grey70", "grey100"))
      )
    ),
    'Circular'=list(
      aes=aes(x=intersection, fill=circ),
      geom=list(
        geom_bar(stat='count', position='fill'),
        scale_y_continuous(labels=scales::percent_format()),
        scale_fill_manual(values=c("grey80","grey20"))
      )
    )
  ),
  width_ratio=0.1
)
#ggsave(filename = "plots/upsetR.pdf", width = 8.5, height= 11, units = "in")
#write.table(upset_tab, "analysis/upset_table.tab", quote = F, row.names = F, sep="\t")



###make pie charts showing overlapping clusters
library(packcircles)
library(scatterpie)
library(ggthemes)

genome_comp <- read.csv('analysis/genome_combine_summary.tab', sep='\t')
table(genome_comp$VC.Status)
cluster_venn <- genome_comp %>%
  filter(VC.Status=='Clustered') %>%
  filter(IsPhage) %>%
  select(VC, VC.members, VC.Size) %>%
  distinct()

cluster_venn$belgCnt <- str_count(cluster_venn$VC.members, 'MN85')
cluster_venn$swissCnt <- str_count(cluster_venn$VC.members, 'JAAOBB')
cluster_venn$texCnt <- str_count(cluster_venn$VC.members, 'contig_|scaffold_|circular')
cluster_venn$refCnt <- cluster_venn$VC.Size - rowSums(cluster_venn[,c('swissCnt', 'belgCnt', 'texCnt')])
cluster_venn$group <- factor(1:nrow(cluster_venn))
  
###use to calculate location of pies
packing <- circleProgressiveLayout(cluster_venn$VC.Size, sizetype='area')
cluster_venn <- cbind(cluster_venn, packing)
dat.gg2 <- circleLayoutVertices(packing, npoints = 100)
ggplot() + 
  geom_polygon(data = dat.gg2, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  geom_text(data = cluster_venn, aes(x, y, size=VC.Size, label = VC)) +
  scale_size_continuous(range = c(1,4)) +
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()


ggplot() +
  geom_scatterpie(data=cluster_venn, 
                  mapping=aes(x=x, y=y, group=group, r=radius),
                  cols=c('refCnt', 'swissCnt', 'belgCnt', 'texCnt'),
                  colour='white') + 
    coord_equal() +
    theme_void() +
    #scale_fill_brewer(palette = 'Dark2') +
    scale_fill_colorblind() +
    geom_scatterpie_legend(cluster_venn$radius, x=-22, y=-10)
#ggsave(filename = "plots/vc_pies.pdf", width = 8.5, height= 8.5, units = "in")
#ggsave(filename = "plots/vc_pies.png", width = 8.5, height= 8.5, units = "in")


######################################look closer at 13 clusters that have all all samples in them (conserved phages)
overlap <- cluster_venn %>%
  filter(swissCnt > 0 & texCnt > 0 & belgCnt > 0)
#swiss_seqs <- read.fasta("raw_data/blastdb/JAAOBB01.1.engel.fasta")
#belg_seqs <- read.fasta("raw_data/blastdb/matthijnssen.fasta")
#tex_seqs <- read.fasta("raw_data/assembly_flye_combined.fasta")
phage_seqs <- read.fasta("analysis/all_phage.fasta")
outdir <- "analysis/clusters"
if (!dir.exists(outdir)) {dir.create(outdir)}
for(i in overlap$VC){
  outdir <- paste("analysis/clusters/", i, sep='')
  if (!dir.exists(outdir)) {dir.create(outdir)}
  tigs <- str_split(cluster_venn[cluster_venn$VC == i,]$VC.members, ',', simplify = T)
  get <- which(names(phage_seqs) %in% tigs)
  for(j in 1:length(get)){
    write.fasta(phage_seqs[get[j]], 
                names = names(phage_seqs)[get[j]], 
                file.out = paste(outdir, '/', names(phage_seqs)[get[j]], '.fasta', sep=''))
  }
  write.fasta(phage_seqs[get], 
              names = names(phage_seqs)[get], 
              file.out = paste('analysis/clusters/', i, '.fasta', sep=''))
}

#######okay. fasta files written. tested Prokka and Clinker and works. Loop through clusters and make clinker visualizations for everything
###could not run in R because of conda. see shell scripts in sh_scripts

##okay, so here I'm seeing 11 shared between the three bee samples and 2 with refseq. 
##figure out why simple venn is showing 9 and 3

##simple venn
###combined viral contigs identified by >3 programs (includes plasmids) with european phage contigs.
##ran all though vContact2 using kbase
#import output file and compare
#seem to have run 2 different vConact2 analyses. Maybe didn't use the same one for both plots. 

#eurocomp <- read.csv("analysis/vContact2_kbase_euro.csv")
eurocomp <- read.csv("raw_data/vContact2/genome_by_genome_overview.csv")
eurocomp <- eurocomp[,-1]

eurocomp$source <- NA
eurocomp[str_detect(eurocomp$Genome, "JAAOBB"),]$source <- "Swiss"
eurocomp[str_detect(eurocomp$Genome, "MN85"),]$source <- "Belgian"
eurocomp[str_detect(eurocomp$Genome, "contig_|circular|scaffold_"),]$source <- "Texas"
eurocomp[is.na(eurocomp$source),]$source <- "RefSeq"
print(nrow(eurocomp))  ##total of 3507 contigs from all datasets 
table(eurocomp$source)

counts <- eurocomp %>%
    group_by(source) %>%
    summarise(contigs=n(), 
              contigs_clustered=sum(!VC == "NaN"),
              clusters=n_distinct(VC)) %>%
    arrange(-contigs)
counts.freq <- counts
counts.freq[,2:4] <- apply(counts.freq[,2:4], 2, function(x) x/sum(x))
counts.f.l <- melt(counts.freq)
counts.l <- melt(counts)

ggplot(counts.f.l, aes(y=variable, x=value, fill=source)) +
  geom_bar(position='stack', stat='identity')

ggplot(counts.l, aes(y=variable, x=value, fill=source, group=value)) +
  geom_bar(position='stack', stat='identity') +
  theme_minimal() +
  scale_y_discrete(limits=c('clusters', 'contigs_clustered', 'contigs')) +
  scale_fill_manual(values = c( '#56b4e9', 'grey','#e69f00','#009e73')) +
  labs(x="", y="")
#ggsave(filename = "plots/cluter_comp.pdf", width = 6.5, height= 1.85, units = "in")
#ggsave(filename = "plots/cluster_comp.png", width = 6.5, height= 1.85, units = "in")


eurocomp_bee <- eurocomp %>%
  filter(source %in% c('Swiss', 'Belgian', 'Texas'))
print(nrow(eurocomp_bee))  ##total of 1203 contigs from bee datasets 

venn_list <- list(Swiss = eurocomp[eurocomp$source == "Swiss",]$VC,
                  Belgian = eurocomp[eurocomp$source == "Belgian",]$VC,
                  Texan = eurocomp[eurocomp$source == "Texas",]$VC,
                  RefSeq = eurocomp[eurocomp$source == "RefSeq",]$VC)

##remove contigs that didn't cluster
table(venn_list$Swiss=="NaN")
table(venn_list$Texan=="NaN")
table(venn_list$Belgian=="NaN")
table(venn_list$RefSeq=="NaN")

venn_list$Swiss <- venn_list$Swiss[!venn_list$Swiss=="NaN"]
venn_list$Belgian <- venn_list$Belgian[!venn_list$Belgian=="NaN"]
venn_list$Texan <- venn_list$Texan[!venn_list$Texan=="NaN"]
venn_list$RefSeq <- venn_list$RefSeq[!venn_list$RefSeq=="NaN"]

venn_list$Swiss <- venn_list$Swiss[!venn_list$Swiss==""]
venn_list$Belgian <- venn_list$Belgian[!venn_list$Belgian==""]
venn_list$Texan <- venn_list$Texan[!venn_list$Texan==""]
venn_list$RefSeq <- venn_list$RefSeq[!venn_list$RefSeq==""]

length(unique(venn_list$Swiss))
length(unique(venn_list$Texan))
length(unique(venn_list$Belgian))
length(unique(venn_list$RefSeq))

##find shared
Reduce(intersect, venn_list)
find <- Reduce(intersect, venn_list)
#View(eurocomp %>%
  filter(VC %in% find))

library(VennDiagram)
#swiss #e69f00
#belg #56b4e9
#tex #009e73

venn.diagram(venn_list,
             fill= c('#e69f00', '#56b4e9', '#009e73', 'grey'),
             lwd = 2,
             col = 'white',
             cex = 2,
             cat.cex=2,
             cat.fontface='bold',
             filename = "plots/venn.tif")

venn <- venn.diagram(venn_list,
             fill= c('#e69f00', '#56b4e9', '#009e73', 'grey'),
             lwd = 2,
             col = 'white',
             cex = 2,
             cat.cex=2,
             cat.fontface='bold',
             filename = NULL)
library(grDevices)
pdf(file='plots/venn.pdf')
  grid.draw(venn)
dev.off()

########taxonomy contig count plots
{
sumTab <- readxl::read_excel('analysis/genome_combine_summary_edit.xlsx')
sumTab$virusTax <- str_replace(sumTab$virusTax, 'Gokushovirinae', 'Microviridae')
sumTab$virusTax <- str_replace(sumTab$virusTax, 'Lactobacillus Firm-5', 'Lactobacillus nr. melliventris')
sumTab$virusTax <- str_replace(sumTab$virusTax, 'Lactobacillus Firm-4', 'Bombilactobacillus spp.')
  
hostTab <- sumTab %>%
  filter(IsPhage, !jtHost=='NA') %>%
  group_by(jtHost) %>%
  summarise(n = n(), freq = sum(Flye.cov)) %>%
  arrange(n)

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nrow(hostTab))
hostTab$jtHost <- factor(hostTab$jtHost, levels=hostTab$jtHost)

hostTax.p <- ggplot(hostTab, aes(y=jtHost, x=n, fill=jtHost)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=n), vjust=0.5, hjust=-0.5, color="grey30", size=3.5) +
  #geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25)
  theme_classic() +
  xlim(0, max(hostTab$n)+7) +
  scale_fill_manual(values=mycolors) +
  theme(legend.position = 'null')+
  labs(x="", y="")
hostTax.p
#ggsave(filename = "plots/host_count.pdf", width = 4.5, height= 3.75, units = "in")
#ggsave(filename = "plots/host_count.png", width = 4.5, height= 3.75, units = "in")

vir.count <- sumTab %>%
  filter(IsPhage, !virusTax=='NA') %>%
  group_by(virusTax) %>%
  summarise(n = n(), freq = sum(Flye.cov)) %>%
  arrange(n)

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nrow(vir.count))
vir.count$virusTax <- factor(vir.count$virusTax, levels=vir.count$virusTax)

virusTax.p <- ggplot(vir.count, aes(y=virusTax, x=n, fill=virusTax)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=n), vjust=0.5, hjust=-0.5, color="grey30", size=3.5) +
  #geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25)
  theme_classic() +
  xlim(0, max(vir.count$n)+2) +
  scale_fill_manual(values=mycolors) +
  theme(legend.position = 'null') +
  labs(x="", y="")
virusTax.p 
#ggsave(filename = "plots/virustax_count.pdf", width = 4.5, height= 3.75, units = "in")
#ggsave(filename = "plots/virustax_count.png", width = 4.5, height= 3.75, units = "in")
  
  
p.tax <- plot_grid(hostTax.p, virusTax.p + labs(x="Count"), nrow=2, rel_heights = c(1,0.55), align = 'v', labels=c("A", "D"))  
p.tax
#ggsave(filename = "plots/taxCombo_count.pdf", width = 3.5, height= 6, units = "in")
#ggsave(filename = "plots/taxCombo_count.png", width = 3.5, height= 6, units = "in")


######tax coverage plots
host.cov <- sumTab %>%
  filter(IsPhage, !jtHost=='NA')
host.cov$jtHost <- factor(host.cov$jtHost, levels=hostTab$jtHost)

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nrow(hostTab))
hostTaxCov.p <- ggplot(host.cov, aes(y=jtHost, x=Flye.cov, fill=jtHost, color=jtHost)) +
#hostTaxCov.p <- ggplot(host.cov, aes(y=jtHost, x=Flye.cov)) +
  geom_point(position = position_jitter(height = 0.2)) +
  theme_classic() +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values=mycolors) +
  theme(legend.position = 'null', 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(x="", y="")
hostTaxCov.p


vir.cov <- sumTab %>%
  filter(IsPhage, !virusTax=='NA')
vir.cov$virusTax <- factor(vir.cov$virusTax, levels=vir.count$virusTax)

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nrow(vir.count))
virusTaxCov.p <- ggplot(vir.cov, aes(y=virusTax, x=Flye.cov, fill=virusTax, color=virusTax)) +
  geom_point(position = position_jitter(height = 0.2)) +
  theme_classic() +
  theme(legend.position = 'null', 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(x="", y="") +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values=mycolors)
virusTaxCov.p 

p.covs <- plot_grid(hostTaxCov.p, virusTaxCov.p+labs(x='Read coverage'), nrow=2, rel_heights = c(1,0.55), align = 'v', labels=c("B", "E"))  
p.covs
#ggsave(filename = "plots/taxCombo_cov.pdf", width = 3.5, height= 6, units = "in")
#ggsave(filename = "plots/taxCombo_cov.png", width = 3.5, height= 6, units = "in")


plot_grid(p.tax, p.covs, nrow=1, align = 'v', rel_widths = c(1,0.7))
#ggsave(filename = "plots/taxCombo_both.pdf", width = 6, height= 6, units = "in")
#ggsave(filename = "plots/taxCombo_both.png", width = 6, height= 6, units = "in")

##stacked bar charts
hostTab$freq <- hostTab$freq/sum(hostTab$freq)
hostTab$jtHost <- factor(hostTab$jtHost, levels=hostTab$jtHost)
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nrow(hostTab))
hostTaxPerc.p <- ggplot(hostTab, aes(y=freq, x=1, fill=jtHost, group=-freq)) +
  geom_bar(stat='identity', position="stack") +
  #geom_text(aes(label=n), vjust=0.5, hjust=-0.5, color="grey30", size=3.5) +
  #geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25)
  theme_classic() +
  scale_fill_manual(values=mycolors) +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(t = 5, r = 0, b = 0, l = 5, unit = "pt")) +
  labs(x="", y="")
hostTaxPerc.p


vir.count$freq <- vir.count$freq/sum(vir.count$freq)
vir.count$virusTax <- factor(vir.count$virusTax, levels=vir.count$virusTax)
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nrow(vir.count))
virTaxPerc.p <- ggplot(vir.count, aes(y=freq, x=1, fill=virusTax, group=-freq)) +
  geom_bar(stat='identity', position="stack") +
  #geom_text(aes(label=n), vjust=0.5, hjust=-0.5, color="grey30", size=3.5) +
  #geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25)
  theme_classic() +
  scale_fill_manual(values=mycolors) +
  theme(legend.position = 'none',
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(t = 0, r = 0, b = 15, l = 5, unit = "pt")) +
  labs(x="Proportion", y="")
virTaxPerc.p



p.freqs <- plot_grid(hostTaxPerc.p, virTaxPerc.p, nrow=2, rel_heights = c(1,1), align = 'hv', labels=c("C", "F"))  
p.freqs
plot_grid(p.tax, p.covs,  p.freqs, nrow=1, axis='bt', align = 'v', rel_widths = c(1,0.7, 0.4))
#ggsave(filename = "plots/taxCombo_all.pdf", width = 8, height= 6, units = "in")
#ggsave(filename = "plots/taxCombo_all.png", width =8, height= 6, units = "in")

}

###plot silva taxonomy from sortmeRNA 
##used pacbio reads >30kb in length (>90% of the reads)
##used SILVAnr96 database
##not sure if this is the best way. 
blast <- read.csv('raw_data/sortmerna/pacbio_SILVAnr96_sortmerna.blast', sep='\t', header = F)
names(blast) <- c("seqid", 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
taxtab <- read.csv('raw_data/sortmerna/SILVA_SSU.noLSU.masked.trimmed.NR96.fixed.headers', header = F)
taxtab <- data.frame(str_split_fixed(taxtab$V1, " ", n = 2))
taxtab$X1 <- str_remove(taxtab$X1, ">")
names(taxtab) <- c('sseqid', 'tax')
#taxtab <- cbind(taxtab, data.frame(str_split_fixed(taxtab$X2, ";", n=9)))

ggplot(blast, aes(x=length, y=pident) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept=1400, linetype="dashed", color = "red") +
  geom_vline(xintercept=1600, linetype="dashed", color = "red") +
  theme_bw()

#ggsave(filename = "plots/silva_matches_pacbio.pdf", width = 8, height= 6, units = "in")
#ggsave(filename = "plots/silva_matches_pacbio.png", width =8, height= 6, units = "in")


###selected 13 vContact2 clusters that contained a representitve phage from every location and looked at them
###more closely. Wrote a short shell script to run prokka on all the contigs and clinker to create visualizations
###clinker outputs pairwise blast comparisons between all genes in all the contigs.
###want to look at these a little more closely here.



cnttab <- blast %>%
  filter(length>1400) %>%
  filter(length<1700) %>%
  group_by(sseqid) %>%
  summarise(nreads=n(),
            minpident=min(pident),
            maxpident=max(pident),
            avgpident=mean(pident)) %>%
  arrange(-nreads)
cnttab <- merge(cnttab, taxtab, by='sseqid')
cnttab <- cnttab %>%
  arrange(-nreads)
write.csv(cnttab, file='analysis/sortmerna_counts.csv', quote = F, row.names = F)
cnttab.b <- cnttab[str_detect(cnttab$tax, 'Bacteria'), ]
tmp <- data.frame(str_split(cnttab.b$tax, ';', simplify = T))
cnttab.b <- cbind(cnttab.b, tmp)
cnttab.b <- cnttab.b[c('sseqid', 'nreads', 'minpident', 'maxpident', 'avgpident', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7')]
names(cnttab.b) <- c('sseqid', 'nreads', 'minpident', 'maxpident', 'avgpident', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

sum(cnttab.b$nreads)
hist(cnttab.b$nreads, breaks=50)

genus.tab <- cnttab.b %>%
  group_by(Genus) %>%
  summarise(nreads=sum(nreads)) %>%
  arrange(-nreads)

genus.tab$freq <- genus.tab$nreads/sum(genus.tab$nreads)

##trying to sort out the SILVIA hits
#first column are the SILVIA identifiers
#CP031513.1051213.1052778 = Bombilactobacillus bombi strain BI-2.5 from Glycophosate paper. This is Firm-4 (CP031513), but don't understand why it is from bumblebees
#AEKZ01001951.393.1940 = not sure. probably Lactobacillus Firm-5
#MH842189.1.1534 = Firm-5, CP009531 = strain wkB8 or possibly something slightly diverged
#MH842189.1.1534 = Firm-5, CP009531 = strain wkB8
#AYWN01000015.724.2340 = 	Mesorhizobium
#CP007799.4461895.4463619 = E. coli = contamination from other pooled library
clusts <- dir('analysis/clusters/')
clusts <- clusts[!str_detect(clusts,'fasta')]
clusts <- clusts[!str_detect(clusts,'redo')]
clusts <- clusts[str_detect(clusts,'_')]

for(i in clusts){
  wd <- paste('analysis/clusters/', i, '/', sep="")
  bfile <- list.files(wd)
  bfile <- bfile[str_detect(bfile, 'clinkerout')]
  blast <- read.table(paste(wd,bfile,sep=''), header = F, fill = T)
  blast$ref <- NA
  blast$query <- NA
  for(j in 1:nrow(blast)){
    tmp <- blast[j,1]
    if(str_detect(tmp, '-------')){
      ref <- blast[j-1,1]
      query <- blast[j-1,3]
    }
    blast[j,]$ref <- ref
    blast[j,]$query <- query
  }
  names(blast)[1:4] <- blast[3,1:4] 
  blast <- blast[str_detect(blast$Target,'_'),]
  blast$comp <- paste(blast$ref, blast$query, sep="-")
  blast$Identity <- as.numeric(blast$Identity)
  blast$Similarity <- as.numeric(blast$Similarity)
  order <- unique(c(blast$ref, blast$query))
  
  p <- ggplot(blast, aes(Similarity)) +
    geom_histogram(stat = "bin", bins = 10) +
    theme_bw() +
    facet_grid(vars(ref), vars(query), drop = T)
  w <- length(order) * 0.85
  #ggsave(p, filename = paste("plots/clusters/", i, '.pdf', sep=''), width = w, height= w, units = "in")
  #ggsave(p, filename = paste("plots/clusters/", i, '.png', sep=''), width =w, height= w, units = "in")
  
  sum <- blast %>% 
    group_by(comp) %>%
    summarise(mean_ident=mean(Identity), mean_sim=mean(Similarity))
  sum$cluster <- i
  write.csv(sum, paste("plots/clusters/", i, '.csv', sep=''), quote = F, row.names = F)
}
