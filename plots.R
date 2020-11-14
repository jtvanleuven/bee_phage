library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(ape)
library(ggsci)
library(UpSetR)


#functions
{
  set_size = function(w, h, factor=1.5) {
    s = 1 * factor
    options(
      repr.plot.width=w * s,
      repr.plot.height=h * s,
      repr.plot.res=100 / factor,
      jupyter.plot_mimetypes='image/png',
      jupyter.plot_scale=1
    )
  }
}

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
                   metaphinder=metaphinder$V1, 
                   PPRmeta=PPRmeta$V1, 
                   vibrant=vibrant$V1, 
                   virfinder=virfinder$V1,
                   virsorter_virome=virsorter_vir$V1,
                   virsorter=virsort$V1,
                   vibrant_virome=vibrant_vir$V1,
                   phaster=phaster)

upset(fromList(find_phage), order.by="freq", nsets = 9)

library(ComplexUpset)
upset_tab <- comb_tab
upset_tab$deepvir <- upset_tab$scaffold %in% deepvir$V1
upset_tab$metaphinder <- upset_tab$scaffold %in% metaphinder$V1
upset_tab$PPRmeta <- upset_tab$scaffold %in% metaphinder$V1
upset_tab$vibrant <- upset_tab$scaffold %in% vibrant$V1
upset_tab$virfinder <- upset_tab$scaffold %in% virfinder$V1
upset_tab$virsorter <- upset_tab$scaffold %in% virsorter$V1
upset_tab$virsorter_virome <- upset_tab$scaffold %in% virsorter_vir$V1
upset_tab$vibrant_virome <- upset_tab$scaffold %in% vibrant_vir$V1
upset_tab$phaster <- upset_tab$scaffold %in% phaster
upset_tab[is.na(upset_tab$Quality),]$Quality <- "NA"

methods <- names(find_phage)

upset_tab <- tibble(upset_tab)

set_size()
upset(
  upset_tab,
  methods,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=F,
      aes=aes(fill=Quality)
    )
  ),
  width_ratio=0.1
)

