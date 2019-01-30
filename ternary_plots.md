#####################Ternary plots
```
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
