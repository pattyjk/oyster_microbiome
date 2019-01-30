## R analysis of oyster 16S rRNA data
```
setwd('oysters/')

##read map and OTU table
oyster_table<-read.delim('oyster_dna_hpcc1000_rare.txt', row.names=1)
oyster_map<-read.delim('Oyster_Map.txt')
oyster_tax<-oyster_table$taxonomy
oyster_table<-oyster_table[,-190]
oyster_table<-as.data.frame(t(oyster_table))
oyster_table$SampleID<-row.names(oyster_table)

##merge meta data to OTU table
oyster_table<-merge(oyster_table, oyster_map[,c(1,4)], by='SampleID', all=T)
#write.csv(oyster_table, 'oyster_table.csv')
#write.table(t(oyster_table), 'oyster_table2.txt')

##subset table to include only oyster tissue
oyster_tissue<-subset(oyster_table, oyster_table$Sample_Type == 'Stomach' | oyster_table$Sample_Type == 'Gills')

oyster_spl
#write.csv(oyster_tissue, 'oyster_tissue.csv')

##remove extra stuff
oyster_tissue<-oyster_tissue[,-1]
oyster_tissue<-as.data.frame(oyster_tissue[,-grep('Sample_Type', colnames(oyster_tissue))])

##take log10 column means
oyster_otu_abund<-colMeans(oyster_tissue)

###########Plot data
per_ocu<-read.delim('data/per_ocu.txt', header=T)
library(ggplot2)
ggplot(per_ocu, aes(per_ocu$Log_av_abun, per_ocu$Perc_occu, colour=as.factor(per_ocu$No_Sites)))+
  geom_point(aes(size=1.2))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  scale_colour_manual(values=c('blue', 'grey', 'black', 'orange', 'red'))+
  ylab("Percent Occupany")+
  xlab("Log10 Mean OTU Abundance")+
  ggtitle('By Habitat')

ggplot(per_ocu, aes(per_ocu$Log_av_abun, per_ocu$Perc_occu, colour=as.factor(per_ocu$No_Sites), shape=per_ocu$Site))+
  geom_point(aes(size=1.2))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  scale_colour_manual(values=c('grey', 'black', 'orange', 'green'))+
  ylab("Percent Occupany")+
  xlab("Log10 Mean OTU Abundance")+
  ggtitle('By site')+
  scale_shape_manual(values=c(19, 21,22,23,1,2,3))

library(gridExtra)
grid.arrange(ocu_hab, ocu_site, ncol=2)


ggplot(per_ocu, aes(per_ocu$Log_av_abun, per_ocu$Perc_occu, colour=Of_interest))+
  geom_point(aes(size=1.2))+
  coord_cartesian(ylim=c(0,1))+
    ylab("Percent Occupany")+
  xlab("Log10 Mean OTU Abundance")+
  theme_bw()
  

########Filter core OTUs in QIIME
filter_otus_from_otu_table.py -i oyster_dna_hpcc1000_rare.biom -o core_otus.biom -e core_oyster_otus.txt --negate_ids_to_exclude
##########
core_taxa<-read.delim('core_taxa_tot_abund.txt', header=T)
library(plyr)
core_sum<-ddply(core_taxa, c('Site', 'Habitat'), summarise, mean=mean(Core_Counts), n=length(Core_Counts), sd=sd(Core_Counts), se=sd/n)

library(ggplot2)
ggplot(core_sum, aes(Site, Habitat))+
  geom_tile(aes(fill=mean))+
  theme_bw()+
  scale_fill_gradient(high = "red", low = "steelblue")+
  ggtitle('Mean total core reads per 1000 reads')

#Core OTU abudance
core_otu_table<-read.delim('core_otu_table.txt', header=T)
library(reshape2)
core_otu_table_m<-melt(core_otu_table)
library(plyr)
core_tab_sum<-ddply(core_otu_table_m, c('Habitat', 'variable', 'Site'), summarise, mean=mean(value), sd=sd(value), n=length(value), se=sd/n, sum=sum(value), rel=mean/sum)


ggplot(core_tab_sum, aes(variable, rel, fill=Habitat))+
  geom_bar(stat='identity')+
   #geom_errorbar(aes(ymax=mean+se, ymin=mean, width=0.2, colour=Habitat), stat="identity")+
  theme_bw()+
  ylab("Mean Relative Abundance")+
  xlab('')+
  facet_wrap(~Site, scales='free')+
  coord_flip()
  


ggplot(core_tab_sum, aes(variable, Habitat))+
  geom_tile(aes(fill=mean))+
  theme_bw()+
    scale_fill_gradient(high = "red", low = "steelblue")+
  facet_wrap(~Site, ncol=1)+
  xlab('')+
  ylab('')+
  ggtitle('Mean reads per 1000 reads')+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))

#add taxonomy
core_otu_tax<-read.delim('core_otu_tax.txt', header=T)

core_otus<-merge(core_tab_sum, core_otu_tax, by='variable')

#split taxonomy into groups
library(tidyr)
core_otus<-separate(core_otus, col=taxonomy, sep=":", into=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'))

pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
ggplot(core_otus, aes(variable, mean, fill=Class))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values = pal)+
  ylab('Mean reads per 1000 reads')+
  ggtitle('Mean reads per 1000 reads')+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))


#################################################################
###############Core OTUs, proportion of reads####################
#################################################################
oyster_rel<-read.delim('core_rel.csv', header=T)
library(reshape2)
library(ggplot2)
rel_m<-melt(oyster_rel)

pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
ggplot(rel_m, aes(variable, value, fill=Habitat))+
  geom_bar(stat='identity')+
  facet_wrap(~Site)+
  theme_bw()+
  coord_flip()+
  xlab('')+
  scale_fill_manual(values = pal)+
  ylab("Proportion of reads")+
  theme(text = element_text(size=14))

