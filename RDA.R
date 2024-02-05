##############packages################
library(vegan)
library(ggplot2)
library(ggrepel)

##############RDA################
fc=read.csv("factor.csv",header = T,row.names = 1)
sp=read.csv("sp.csv",header = T,row.names = 1)
spp=decostand(sp,method = "hellinger")
fcc=log10(fc)
uu=rda(spp~.,fcc)
ii=summary(uu,scaling =2) 
sp=as.data.frame(ii$species[,1:2])/1.5
st=as.data.frame(ii$sites[,1:2])/1.5
yz=as.data.frame(ii$biplot[,1:2])/4
grp=as.data.frame(c(rep("a",4),rep("b",4),rep("c",4),rep("d",4)))
colnames(grp)="group"

ggplot() +
  geom_point(data = st,aes(RDA1,RDA2,shape=grp$group,fill=grp$group),size=4)+
  scale_shape_manual(values = c(21:25))+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "red")+
  geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)))+
  labs(x=paste("RDA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())

anova.cca(uu)
envfit(uu,fcc)