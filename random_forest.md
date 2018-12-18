## Random forest analysis
```
#######################Random Forest Per Site##########################

####################################################################
####################DU1#############################################
####################################################################
setwd("/home/pattyjk/Dropbox/R/oysters/")
library("randomForest")
library("plyr")
library("rfUtilities")
library("caret")

#read in OTU table
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
dim(s16)
#236396 OTUs by 202 samples

#remove taxonomy
s16<-s16[,-202]
dim(s16)
#236396 OTUs by 201 samples

#how many nonzero counts?
otu_nonzero_counts_DU1<-apply(s16_DU1, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts_DU1, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

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

otu_table_rare_removed_DU1 <- remove_rare(table=s16_DU1, cutoff_pro=0.1)
dim(otu_table_rare_removed_DU1)
#3765 OTUs by 62 samples

#scale data
otu_table_scaled_DU1 <- scale(otu_table_rare_removed_DU1, center = TRUE, scale = TRUE)  

#add meta
otu_table_scaled_treatment_DU1 <- data.frame(t(otu_table_scaled_DU1))
otu_table_scaled_treatment_DU1$SampleID<-row.names(otu_table_scaled_treatment_DU1)
meta_sub<-as.data.frame(meta[,c(1,4)])
names(meta_sub)<-c("SampleID", "Sample_Type")
otu_table_scaled_treatment_DU1<-merge(otu_table_scaled_treatment_DU1, meta_sub, by=c("SampleID"))
dim(otu_table_scaled_treatment_DU1)
#62 by 3767
otu_table_scaled_treatment_DU1<-otu_table_scaled_treatment_DU1[,-1]
dim(otu_table_scaled_treatment_DU1)
#62 by 3766

#run RF model
set.seed(501)
RF_treatment_classify_DU1<-randomForest(x=otu_table_scaled_treatment_DU1[,1:(ncol(otu_table_scaled_treatment_DU1)-1)], y=otu_table_scaled_treatment_DU1[ , ncol(otu_table_scaled_treatment_DU1)] , ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=1000 , ntree=501)

#identifying important features
RF_state_classify_imp_DU1 <- as.data.frame(RF_treatment_classify_DU1$importance)
RF_state_classify_imp_DU1$features <- rownames( RF_state_classify_imp_DU1 )
RF_state_classify_imp_sorted_DU1 <- arrange( RF_state_classify_imp_DU1  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted_DU1$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#top 50 features
barplot(RF_state_classify_imp_sorted_DU1[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted_DU1[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF") 

top50_feat_DU1<-as.data.frame(RF_state_classify_imp_sorted_DU1$features[1:50])
names(top50_feat_DU1)<-c("OTU")
str(top50_feat_DU1)

#read in OTU table, extract taxonomy
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
tax<-as.data.frame(s16$taxonomy)
names(tax)<-c("taxonomy")
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)
str(s16_tax)
s16<-s16[,-202]

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_DU1<-s16[,colnames(s16) %in% meta_split$DU1$SampleID]
dim(s16_DU1)
#236396 by 62

s16_DU2<-s16[,colnames(s16) %in% meta_split$DU2$SampleID]
dim(s16_DU2)
#236396 by 43

s16_RI<-s16[,colnames(s16) %in% meta_split$RI$SampleID]
dim(s16_RI)
#236396 by 58

s16_WF<-s16[,colnames(s16) %in% meta_split$WF$SampleID]
dim(s16_WF)
#236396 by 36

#change to relative abundance
library(vegan)
s16_DU1<-decostand(s16_DU1, method = 'total')

#extract top 50 from OTU table
s16.top50_DU1<-s16_DU1[rownames(s16_DU1) %in% top50_feat_DU1$OTU,]
dim(s16.top50_DU1)  
s16.top50_DU1$OTU<-row.names(s16.top50_DU1)
names(s16.top50_DU1)
rownames(s16.top50_DU1)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m_DU1<-melt(s16.top50_DU1)
names(top50_m_DU1)<-c("OTU", "SampleID", 'Rel_abund')
library(plyr)
top50_m_DU1<-merge(top50_m_DU1, meta, by='SampleID')
top50_sum_DU1<-ddply(top50_m_DU1, c('OTU', "Sample_Type"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum_DU1<-merge(top50_sum_DU1, s16_tax, by='OTU')

#split taxonomy into groups
library(tidyr)
library(stringi)
top50_sum_DU1<-separate(top50_sum_DU1, taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#plot the data
colour_tax<-rainbow(41, s=1, v=1)[sample(1:41,41)]
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

library(ggplot2)
ggplot(top50_sum_DU1, aes(OTU, mean, fill=Class))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  facet_wrap(~Sample_Type, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 OTUs in distinguishing categories- across all sites- DU1")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")

####################################################################
####################DU2#############################################
####################################################################
#read in OTU table
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
dim(s16)
#236396 OTUs by 202 samples

#remove taxonomy
s16<-s16[,-202]
dim(s16)
#236396 OTUs by 201 samples

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_DU2<-s16[,colnames(s16) %in% meta_split$DU2$SampleID]
dim(s16_DU2)
#236396 by 43

#nonzero counts
otu_nonzero_counts_DU2<-apply(s16_DU2, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts_DU2, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

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

otu_table_rare_removed_DU2 <- remove_rare(table=s16_DU2, cutoff_pro=0.1)
dim(otu_table_rare_removed_DU2)
#10537 OTUs by 43 samples

#scale data
otu_table_scaled_DU2 <- scale(otu_table_rare_removed_DU2, center = TRUE, scale = TRUE)  

#add meta
otu_table_scaled_treatment_DU2 <- data.frame(t(otu_table_scaled_DU2))
otu_table_scaled_treatment_DU2$SampleID<-row.names(otu_table_scaled_treatment_DU2)
meta_sub<-as.data.frame(meta[,c(1,4)])
names(meta_sub)<-c("SampleID", "Sample_Type")
otu_table_scaled_treatment_DU2<-merge(otu_table_scaled_treatment_DU2, meta_sub, by=c("SampleID"))
dim(otu_table_scaled_treatment_DU2)
#43 by 10539
otu_table_scaled_treatment_DU2<-otu_table_scaled_treatment_DU2[,-1]
dim(otu_table_scaled_treatment_DU2)
#43 by 10538

#run RF model
set.seed(501)
RF_treatment_classify_DU2<-randomForest(x=otu_table_scaled_treatment_DU2[,1:(ncol(otu_table_scaled_treatment_DU2)-1)], y=otu_table_scaled_treatment_DU2[ , ncol(otu_table_scaled_treatment_DU2)] , ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
#RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=1000 , ntree=501)

#identifying important features
RF_state_classify_imp_DU2 <- as.data.frame(RF_treatment_classify_DU2$importance)
RF_state_classify_imp_DU2$features <- rownames( RF_state_classify_imp_DU2 )
RF_state_classify_imp_sorted_DU2 <- arrange( RF_state_classify_imp_DU2  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted_DU2$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#top 50 features
barplot(RF_state_classify_imp_sorted_DU2[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted_DU2[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF") 
top50_feat_DU2<-as.data.frame(RF_state_classify_imp_sorted_DU2$features[1:50])
names(top50_feat_DU2)<-c("OTU")
str(top50_feat_DU2)

#read in OTU table, extract taxonomy
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
tax<-as.data.frame(s16$taxonomy)
names(tax)<-c("taxonomy")
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)
str(s16_tax)
s16<-s16[,-202]

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_DU2<-s16[,colnames(s16) %in% meta_split$DU2$SampleID]
dim(s16_DU2)
#236396 by 43

#change to relative abundance
library(vegan)
s16_DU2<-decostand(s16_DU2, method = 'total')

#extract top 50 from OTU table
s16.top50_DU2<-s16_DU2[rownames(s16_DU2) %in% top50_feat_DU2$OTU,]
dim(s16.top50_DU2)  
#50 by 43
s16.top50_DU2$OTU<-row.names(s16.top50_DU2)
names(s16.top50_DU2)
rownames(s16.top50_DU2)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m_DU2<-melt(s16.top50_DU2)
names(top50_m_DU2)<-c("OTU", "SampleID", 'Rel_abund')
library(plyr)
top50_m_DU2<-merge(top50_m_DU2, meta, by='SampleID')
top50_sum_DU2<-ddply(top50_m_DU2, c('OTU', "Sample_Type"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum_DU2<-merge(top50_sum_DU2, s16_tax, by='OTU')

#split taxonomy into groups
library(tidyr)
library(stringi)
top50_sum_DU2<-separate(top50_sum_DU2, taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#plot the data
colour_tax<-rainbow(41, s=1, v=1)[sample(1:41,41)]
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

library(ggplot2)
ggplot(top50_sum_DU2, aes(OTU, mean, fill=Class))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  facet_wrap(~Sample_Type, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 OTUs in distinguishing categories- across all sites- DU2")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+sd, ymin=mean, width=0.2, colour=Class), stat="identity")


####################################################################
####################RI##############################################
####################################################################
#read in OTU table
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
dim(s16)
#236396 OTUs by 202 samples

#remove taxonomy
s16<-s16[,-202]
dim(s16)
#236396 OTUs by 201 samples

otu_nonzero_counts_RI<-apply(s16_RI, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts_RI, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

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

otu_table_rare_removed_RI <- remove_rare(table=s16_RI, cutoff_pro=0.1)
dim(otu_table_rare_removed_RI)
#5611 OTUs by 58 samples

#scale data
otu_table_scaled_RI <- scale(otu_table_rare_removed_RI, center = TRUE, scale = TRUE) 

#add meta
otu_table_scaled_treatment_RI <- data.frame(t(otu_table_scaled_RI))
otu_table_scaled_treatment_RI$SampleID<-row.names(otu_table_scaled_treatment_RI)
meta_sub<-as.data.frame(meta[,c(1,4)])
names(meta_sub)<-c("SampleID", "Sample_Type")
otu_table_scaled_treatment_RI<-merge(otu_table_scaled_treatment_RI, meta_sub, by=c("SampleID"))
dim(otu_table_scaled_treatment_RI)
#58 by 5613
otu_table_scaled_treatment_RI<-otu_table_scaled_treatment_RI[,-1]
dim(otu_table_scaled_treatment_RI)
#58 by 5612

#run RF model
set.seed(501)
RF_treatment_classify_RI<-randomForest(x=otu_table_scaled_treatment_RI[,1:(ncol(otu_table_scaled_treatment_RI)-1)], y=otu_table_scaled_treatment_RI[ , ncol(otu_table_scaled_treatment_RI)] , ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
#RF_treatment_classify_sig_RI<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=1000 , ntree=501)

#identifying important features
RF_state_classify_imp_RI <- as.data.frame(RF_treatment_classify_RI$importance)
RF_state_classify_imp_RI$features <- rownames( RF_state_classify_imp_RI )
RF_state_classify_imp_sorted_RI <- arrange( RF_state_classify_imp_RI  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted_RI$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#top 50 features
barplot(RF_state_classify_imp_sorted_RI[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted_RI[1:10,"features"] , ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF") 

top50_feat_RI<-as.data.frame(RF_state_classify_imp_sorted_RI$features[1:50])
names(top50_feat_RI)<-c("OTU")
str(top50_feat_RI)

#read in OTU table, extract taxonomy
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
tax<-as.data.frame(s16$taxonomy)
names(tax)<-c("taxonomy")
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)
str(s16_tax)
s16<-s16[,-202]

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_RI<-s16[,colnames(s16) %in% meta_split$RI$SampleID]
dim(s16_RI)
#236396 by 58

#change to relative abundance
library(vegan)
s16_RI<-decostand(s16_RI, method = 'total')

#extract top 50 from OTU table
s16.top50_RI<-s16_RI[rownames(s16_RI) %in% top50_feat_RI$OTU,]
dim(s16.top50_RI)  
#50 by 58
s16.top50_RI$OTU<-row.names(s16.top50_RI)
names(s16.top50_RI)
rownames(s16.top50_RI)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m_RI<-melt(s16.top50_RI)
names(top50_m_RI)<-c("OTU", "SampleID", 'Rel_abund')
library(plyr)
top50_m_RI<-merge(top50_m_RI, meta, by='SampleID')
top50_sum_RI<-ddply(top50_m_RI, c('OTU', "Sample_Type"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum_RI<-merge(top50_sum_RI, s16_tax, by='OTU')

#split taxonomy into groups
library(tidyr)
library(stringi)
top50_sum_RI<-separate(top50_sum_RI, taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#plot the data
colour_tax<-rainbow(41, s=1, v=1)[sample(1:41,41)]
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

library(ggplot2)
ggplot(top50_sum_RI, aes(OTU, mean, fill=Class))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  facet_wrap(~Sample_Type, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 OTUs in distinguishing categories- across all sites- RI")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")


####################################################################
####################WF##############################################
####################################################################
#read in OTU table
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
dim(s16)
#236396 OTUs by 202 samples

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_WF<-s16[,colnames(s16) %in% meta_split$WF$SampleID]
dim(s16_WF)
#236396 by 36

#non zero counts
otu_nonzero_counts_WF<-apply(s16_WF, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts_WF, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

#remove taxonomy
s16<-s16[,-202]
dim(s16)
#236396 OTUs by 201 samples

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

otu_table_rare_removed_WF <- remove_rare(table=s16_WF, cutoff_pro=0.2)
dim(otu_table_rare_removed_WF)
#4506 OTUs by 36 samples

#scale data
otu_table_scaled_WF <- scale(otu_table_rare_removed_WF, center = TRUE, scale = TRUE)  

#add meta
otu_table_scaled_treatment_WF <- data.frame(t(otu_table_scaled_WF))
otu_table_scaled_treatment_WF$SampleID<-row.names(otu_table_scaled_treatment_WF)
meta_sub<-as.data.frame(meta[,c(1,4)])
names(meta_sub)<-c("SampleID", "Sample_Type")
otu_table_scaled_treatment_WF<-merge(otu_table_scaled_treatment_WF, meta_sub, by=c("SampleID"))
dim(otu_table_scaled_treatment_WF)
#36 by 9173
otu_table_scaled_treatment_WF<-otu_table_scaled_treatment_WF[,-1]
dim(otu_table_scaled_treatment_WF)
#36 by 9172

#run RF model
set.seed(501)
RF_treatment_classify_WF<-randomForest(x=otu_table_scaled_treatment_WF[,1:(ncol(otu_table_scaled_treatment_WF)-1)], y=otu_table_scaled_treatment_WF[ , ncol(otu_table_scaled_treatment_WF)] , ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
#RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify, xdata=otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)], nperm=1000 , ntree=501)

#identifying important features
RF_state_classify_imp_WF <- as.data.frame(RF_treatment_classify_WF$importance)
RF_state_classify_imp_WF$features <- rownames( RF_state_classify_imp_WF )
RF_state_classify_imp_sorted_WF <- arrange( RF_state_classify_imp_WF  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted_WF$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (VaWFable Importance)", main="RF Classification VaWFable Importance DistWFbution")

#top 50 features
barplot(RF_state_classify_imp_sorted_WF[1:50,"MeanDecreaseAccuracy"], las=2, names.arg=RF_state_classify_imp_sorted_WF[1:10,"features"] , ylab="Mean Decrease in Accuracy (VaWFable Importance)", main="Classification RF") 
top50_feat_WF<-as.data.frame(RF_state_classify_imp_sorted_WF$features[1:50])
names(top50_feat_WF)<-c("OTU")
str(top50_feat_WF)

#read in OTU table, extract taxonomy
s16<-read.delim("oyster_dna.txt", header=T, row.names=1)
tax<-as.data.frame(s16$taxonomy)
names(tax)<-c("taxonomy")
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)
str(s16_tax)
s16<-s16[,-202]

#read in meta data
meta<-read.delim("Oyster_Map-n.txt", header=T)

#split meta by site
meta_split<-split(meta, meta$Location)

#split OTU table by site
s16_WF<-s16[,colnames(s16) %in% meta_split$WF$SampleID]
dim(s16_WF)
#236396 by 62

#change to relative abundance
library(vegan)
s16_WF<-decostand(s16_WF, method = 'total')

#extract top 50 from OTU table
s16.top50_WF<-s16_WF[rownames(s16_WF) %in% top50_feat_WF$OTU,]
dim(s16.top50_WF)  
#50 by 37
s16.top50_WF$OTU<-row.names(s16.top50_WF)
names(s16.top50_WF)
rownames(s16.top50_WF)

#get mean relative abundance of each OTU in the 4 cats
library(reshape)
top50_m_WF<-melt(s16.top50_WF)
names(top50_m_WF)<-c("OTU", "SampleID", 'Rel_abund')
library(plyr)
top50_m_WF<-merge(top50_m_WF, meta, by='SampleID')
top50_sum_WF<-ddply(top50_m_WF, c('OTU', "Sample_Type"), summarize, mean=mean(Rel_abund), sd=sd(Rel_abund), n=length(Rel_abund), se=sd/n)
top50_sum_WF<-merge(top50_sum_WF, s16_tax, by='OTU')

#split taxonomy into groups
library(tidyr)
library(stringi)
top50_sum_WF<-separate(top50_sum_WF, taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#plot the data
colour_tax<-rainbow(41, s=1, v=1)[sample(1:41,41)]
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

library(ggplot2)
ggplot(top50_sum_WF, aes(OTU, mean, fill=Class))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  facet_wrap(~Sample_Type, ncol = 4)+
  theme_bw()+
  coord_flip()+
  ggtitle("Top 50 OTUs in distinguishing categories- across all sites- WF")+
  guides(fill=guide_legend(ncol=1))+
  ylab("Mean Relative Abundance")+
  xlab("")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")
```
