## Source Tracker Analysis- full dataset
```
setwd("oyster_microbiome")
# load sample metadata
metadata <- read.table('data/map_source2.txt', sep='\t',h=T, row.names = 1)

# load OTU table
otus <- read.table('data/oyster_dna_hpcc1000_rare.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')


otus <- read.table('data/source_table2.txt',sep='\t', header=T, row.names=1)
#remove taxonomy
otus<-otus[,-190]

#transpose OTU table
otus <- as.data.frame(t(otus))

# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}
#its good

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='Source')
test.ix <- which(metadata$SourceSink=='Sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
#source('/Users/patty/OneDrive/Desktop/sourcetracker-1.0.1/sourcetracker-1.0.1/src/SourceTracker.r')
source('/home/pattyjk/Dropbox/R/sourcetracker-1.0.1/src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

#data was written to data files for subsequent use
saveRDS(results, file = "s_tracked_results.rds")
saveRDS(results.train, file = "s_tracker_train.rds",)
```

## Plot data- full dataset
```
results<-readRDS(file = "s_tracker_results.rds")
results.train<-readRDS(file='s_tracker_train.rds')

#extract proportions
proportions<-as.data.frame(results$proportions)
proportions$SampleID<-row.names(proportions)

library(reshape2)
library(ggplot2)
library(plyr)

#reshape data
proportions_m<-melt(proportions)

#tack on metadata
metadata <- read.table('data/map_source2.txt', sep='\t', header=T)
proportions_m<-merge(proportions_m, metadata, by='SampleID')
proportions_m$variable<-gsub("Sed", "stuff", proportions_m$variable)
proportions_m$variable<-gsub("WC", "Sed", proportions_m$variable)
proportions_m$variable<-gsub("stuff", "WC", proportions_m$variable)
proportions_m$Sample_Type<-gsub("Sed", "Stomach", proportions_m$Sample_Type)
proportions_m$Sample_Type<-gsub("WC", "Gills", proportions_m$Sample_Type)

#summarize data
proportion_sum<-ddply(proportions_m, c("Location", "variable", "Sample_Type"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#plot data as barplot
ggplot(proportion_sum, aes(Location, mean, fill=variable))+
geom_bar(stat='identity')+
theme_bw()+
ylab("Proportion of Community")+
xlab("")+
facet_wrap(~Sample_Type)

#plot data as dotplot
ggplot(proportion_sum, aes(Location, mean, colour=variable))+
geom_point(position=position_dodge(0.5), aes(size=2))+
theme_bw()+
theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))+
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.5), width=.2)+
ylab("Mean Proportion of the community")+
xlab("")+
coord_flip()+
facet_wrap(~Sample_Type)
```

## Source tracker- core microbiome
```
#calculate source for core microbes

setwd("oyster_microbiome")
# load sample metadata
metadata <- read.table('data/map_source2.txt', sep='\t',h=T, row.names = 1)

# load core microbiome OTU table
otus <- read.table('data/core_otu_table.txt',sep='\t', header=T,row.names=1,check=F,comment='')

#remove metadata
otus<-otus[,c(-1,-2)]

otus <- read.table('data/source_table2.txt',sep='\t', header=T, row.names=1)
#remove taxonomy
otus<-otus[,-190]

#transpose OTU table
otus <- as.data.frame(t(otus))

# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}
#its good

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='Source')
test.ix <- which(metadata$SourceSink=='Sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('/home/pattyjk/Dropbox/GitHub/oyster_microbiome/sourcetracker-1.0.1/src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st2 <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions
results2 <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train2 <- predict(st2, alpha1=alpha1, alpha2=alpha2)

#data was written to data files for subsequent use
saveRDS(results2, file = "s_tracked_results_core.rds")
saveRDS(results.train2, file = "s_tracker_train_core.rds")

##plot data
library(reshape2)
library(ggplot2)
library(plyr)

#extract proportions
proportions2<-as.data.frame(results2$proportions)
proportions2$SampleID<-row.names(proportions2)

#reshape data
proportions_m2<-melt(proportions2)

#tack on metadata
metadata <- read.table('data/map_source2.txt', sep='\t', header=T)
proportions_m2<-merge(proportions_m2, metadata, by='SampleID')
proportions_m2$variable<-gsub("WC", "stuff", proportions_m2$variable)
proportions_m2$variable<-gsub("Unknown", "WC", proportions_m2$variable)
proportions_m2$variable<-gsub("stuff", "Unknown", proportions_m2$variable)


#summarize data
proportion_sum2<-ddply(proportions_m2, c("Location", "variable"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#plot data as barplot
ggplot(proportion_sum2, aes(Location, mean, fill=variable))+
geom_bar(stat='identity')+
theme_bw()+
ylab("Proportion of Community")+
xlab("")


#plot data as dotplot
ggplot(proportion_sum2, aes(Location, mean, colour=variable))+
geom_point(position=position_dodge(0.5), aes(size=2))+
theme_bw()+
theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14))+
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.5), width=.2)+
ylab("Mean Proportion of the community")+
xlab("")+
coord_flip()
```
