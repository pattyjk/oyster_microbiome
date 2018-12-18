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

length(which())

###########Plot data
per_ocu<-read.delim('per_ocu.txt', header=T)
library(ggplot2)
ocu_hab<-ggplot(per_ocu, aes(per_ocu$Log_av_abun, per_ocu$Perc_occu, colour=as.factor(per_ocu$No_Sites)))+
  geom_point(aes(size=1.2))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  scale_colour_manual(values=c('blue', 'grey', 'black', 'orange', 'red'))+
  ylab("Percent Occupany")+
  xlab("Log10 Mean OTU Abundance")+
  ggtitle('By Habitat')

ocu_site<-ggplot(per_ocu, aes(per_ocu$Log_av_abun, per_ocu$Perc_occu, colour=as.factor(per_ocu$No_Sites), shape=per_ocu$Site))+
  geom_point(aes(size=1.2))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  scale_colour_manual(values=c('grey', 'black', 'orange'))+
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
###############QIIME Code
beta_diversity.py -m bray_curtis -i duxbury1000_table.biom -o dux_beta
beta_diversity.py -m bray_curtis -i WF1000_table.biom -o WF_beta
beta_diversity.py -m bray_curtis -i RI1000_table.biom -o ri_beta

make_distance_boxplots.py -m Oyster_Map.txt -o WF_beta/wf_dis -d WF_beta/bray_curtis_WF1000_table.txt --save_raw_data -f 'Sample_Type'
make_distance_boxplots.py -m Oyster_Map.txt -o dux_beta/dux_dis -d dux_beta/bray_curtis_duxbury1000_table.txt --save_raw_data -f 'Sample_Type'
make_distance_boxplots.py -m Oyster_Map.txt -o ri_beta/ri_dis -d ri_beta/bray_curtis_RI1000_table.txt --save_raw_data -f 'Sample_Type'


###########################################################################
##################Habitat WU boxplots######################################
###########################################################################
oyst_dis_bray<-read.delim('oyster_dist_bray2.txt', header=T)
library(ggplot2)
#library(reshape2)
#oyst_dis_m<-melt(oyst_dis_bray)
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
ggplot(oyst_dis_bray, aes(Habitats, Sim, fill=Habitats))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Site)+
  xlab('')+
  ylab('Weighted UniFrac Similarity')+
   theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = pal)+
  coord_cartesian(ylim=c(0,1))

dist_aov<-aov(Sim ~Habitats + Site, data=oyst_dis_bray)

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

####################################################################
####################Random Forest###################################
####################################################################
setwd("/home/pattyjk/Dropbox/R/oysters/")
library("randomForest")
library("plyr")
library("rfUtilities")
library("caret")

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#read in OTU table
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
dim(s16)
#236396 OTUs by 202 samples

#remove taxonomy
s16<-s16[,-202]
dim(s16)
#236396 OTUs by 201 samples

#how many nonzero counts?
otu_nonzero_counts<-apply(s16, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

#remove rare taxa
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=s16, cutoff_pro=0.1)
otu_table_rare_removeda<-otu_table_rare_removed
dim(otu_table_rare_removed)
#6932 OTUs by 201 samples

#scale data
otu_table_scaled <- scale(otu_table_rare_removed, center = TRUE, scale = TRUE)  

#add meta
otu_table_scaled_treatment <- data.frame(t(otu_table_scaled))
otu_table_scaled_treatment$SampleID<-row.names(otu_table_scaled_treatment)
meta_sub<-as.data.frame(meta[,c(1,4)])
names(meta_sub)<-c("SampleID", "Sample_Type")
otu_table_scaled_treatment<-merge(otu_table_scaled_treatment, meta_sub, by=c("SampleID"))
dim(otu_table_scaled_treatment)
#199 by 6934
otu_table_scaled_treatment<-otu_table_scaled_treatment[,-1]
dim(otu_table_scaled_treatment)
##199 by 6933

#run RF model
set.seed(501)
RF_treatment_classify<-randomForest(x=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], y=otu_table_scaled_treatment[ , ncol(otu_table_scaled_treatment)] , ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
#RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=1000 , ntree=501)

#identifying important features
RF_state_classify_imp <- as.data.frame(RF_treatment_classify$importance)
RF_state_classify_imp$features <- rownames( RF_state_classify_imp )
RF_state_classify_imp_sorted <- arrange( RF_state_classify_imp  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#top 50 features
barplot(RF_state_classify_imp_sorted[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF") 

top50_feat<-as.data.frame(RF_state_classify_imp_sorted$features[1:50])
names(top50_feat)<-c("OTU")
str(top50_feat)

#read in OTU table, extract taxonomy
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
tax<-as.data.frame(s16$taxonomy)
names(tax)<-c("taxonomy")
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)
str(s16_tax)
s16<-s16[,-202]

#change to relative abundance
library(vegan)
s16<-decostand(s16, method = 'total')
#check to make sure rel abund
#rowSums(s16[,-52])

#extract top 50 from OTU table
s16.top50<-s16[rownames(s16) %in% top50_feat$OTU,]
dim(s16.top50)  
s16.top50$OTU<-row.names(s16.top50)
names(s16.top50)
rownames(s16.top50)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m<-melt(s16.top50)
names(top50_m)<-c("OTU", "SampleID", 'Rel_abund')
library(plyr)
top50_m<-merge(top50_m, meta, by='SampleID')
top50_sum<-ddply(top50_m, c('OTU', "Sample_Type"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum<-merge(top50_sum, s16_tax, by='OTU')

#split taxonomy into groups
library(tidyr)
library(stringi)
top50_sum<-separate(top50_sum, taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#plot the data
colour_tax<-rainbow(41, s=1, v=1)[sample(1:41,41)]
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

library(ggplot2)
ggplot(top50_sum, aes(OTU, mean, fill=Class))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  facet_wrap(~Sample_Type, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 OTUs in distinguishing categories- across all sites")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")
  ```
