## Source Tracker Analysis
```
################oyster source tracker

# load sample metadata

#windows
#metadata <- read.table('/Users/patty/Dropbox/R/oysters/map_source.txt', sep='\t',h=T,row.names=1,check=F,comment='')
#linux
metadata <- read.table('/home/pattyjk/Dropbox/R/oysters/map_source2.txt', sep='\t',h=T, row.names = 1)

# load OTU table
otus <- read.table('/Users/patty/Dropbox/R/oysters/oyster_dna_hpcc1000_rare.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')

#linux
otus <- read.table('/home/pattyjk/Dropbox/R/oysters/source_table2.txt',sep='\t', header=T, row.names=1)
#remove taxonomy
otus<-otus[,-190]
otus <- as.data.frame(t(otus))

# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}
#is clean son

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

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
plot(results, labels[test.ix], type='bar')
plot(results, labels[test.ix], type='dist')
plot(results.train, labels[train.ix], type='pie')
plot(results.train, labels[train.ix], type='bar')
plot(results.train, labels[train.ix], type='dist')



# plot results with legend
# plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))


#####################Ternary plots
setwd("./oysters/")
library(ggtern)
oysters_global_tern<-read.delim("oyster_global_tern.txt", header=T)
oysters_global_tern$mean<-rowMeans(oysters_global_tern[,-1])
  
ggtern(oysters_global_tern, aes(x=WC, y=Sed, z=Stomach, colour=Class, size=mean))+
         geom_point()+
  theme_bw()+
  guides(colour=F)+
  ggtitle("All Sites")

ggtern(oysters_global_tern, aes(x=WC, y=Sed, z=Gills, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=F)+
  ggtitle("All Sites")

##DU1
DU1<-read.delim("DU1_tern.txt", header=T)
DU1$mean<-rowMeans(DU1[,-1])
ggtern(DU1, aes(x=WC, y=Sed, z=Stomach, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=F)+
  ggtitle("DU1")

ggtern(DU1, aes(x=WC, y=Sed, z=Gills, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=F)+
  ggtitle("DU1")

###DU2
DU2<-read.delim("DU2_tern.txt", header=T)
DU2$mean<-rowMeans(DU2[,-1])
ggtern(DU2, aes(x=WC, y=Sed, z=Stomach, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  ggtitle("DU2")+
  guides(colour=F)

ggtern(DU2, aes(x=WC, y=Sed, z=Gills, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  ggtitle("DU2")+
  guides(colour=F)

##WF
WF<-read.delim("WF_tern.txt", header=T)
WF$mean<-rowMeans(WF[,-1])
ggtern(WF, aes(x=WC, y=Sed, z=Stomach, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=FALSE)+
  ggtitle("WF")

ggtern(WF, aes(x=WC, y=Sed, z=Gills, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=FALSE)+
  ggtitle("WF")

###RI
RI<-read.delim("RI_tern.txt", header=T)
RI$mean<-rowMeans(RI[,-1])
ggtern(RI, aes(x=WC, y=Sed, z=Stomach, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=FALSE)+
  ggtitle("RI")

ggtern(RI, aes(x=WC, y=Sed, z=Gills, colour=Class, size=mean))+
  geom_point()+
  theme_bw()+
  guides(colour=FALSE)+
  ggtitle("RI")
  ```
