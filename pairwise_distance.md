## Examining community similairty with pairwise distances
### Calculate similarity comparisons in QIIME
```
#QIIME 1.9.1
beta_diversity.py -m bray_curtis -i duxbury1000_table.biom -o dux_beta
beta_diversity.py -m bray_curtis -i WF1000_table.biom -o WF_beta
beta_diversity.py -m bray_curtis -i RI1000_table.biom -o ri_beta

make_distance_boxplots.py -m Oyster_Map.txt -o WF_beta/wf_dis -d WF_beta/bray_curtis_WF1000_table.txt --save_raw_data -f 'Sample_Type'
make_distance_boxplots.py -m Oyster_Map.txt -o dux_beta/dux_dis -d dux_beta/bray_curtis_duxbury1000_table.txt --save_raw_data -f 'Sample_Type'
make_distance_boxplots.py -m Oyster_Map.txt -o ri_beta/ri_dis -d ri_beta/bray_curtis_RI1000_table.txt --save_raw_data -f 'Sample_Type'
```

### Plot similarities in R
```
#R . 3.4.1
library(ggplot2)
library(reshape2)
library(plyr)

#read in data
oyst_dis_bray<-read.delim('oyster_dist_bray2.txt', header=T)

#calculate summary stats of WU similarities
oyst_dis_m<-melt(oyst_dis_bray)
oyst_dis_sum<-ddply(oyst_dis_m, c('Habitats', 'Site'), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#write summary to file
write.table(oyst_dis_sum, "oyst_dis_sum.txt", quote=F, row.names=F, sep="\t")

#plot with boxplots
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

#analyze sig. with ANOVA/TukeyHSD
dist_aov<-aov(Sim ~Habitats + Site, data=oyst_dis_bray)
tukeyHSD(dist_aov)
```
