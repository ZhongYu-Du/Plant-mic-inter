##############packages################
library(vegan)
library(ggplot2)
library(ggrepel)

##############CCA################
fc=read.csv("factor.csv",header = T,row.names = 1)
sp=read.csv("sp.csv",header = T,row.names = 1)
spp=decostand(sp,method = "hellinger")
fcc=log10(fc)
uu=cca(spp~.,fcc)
ii=summary(uu)
sp=as.data.frame(ii$species[,1:2])*20
st=as.data.frame(ii$constraints[,1:2])
yz=as.data.frame(ii$biplot[,1:2])*5
fenzu=as.data.frame(c(rep("a",4),rep("b",4),rep("c",4),rep("d",4)))
colnames(fenzu)="group"

ggplot() +
  geom_point(data = st,aes(CCA1,CCA2,shape=fenzu$group,fill=fenzu$group),size=4)+
  scale_shape_manual(values = c(21:25))+
  geom_text_repel(data = sp,aes(CCA1,CCA2,label=row.names(sp)))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "black")+
  geom_text_repel(data = yz,aes(CCA1,CCA2,label=row.names(yz)))+
  labs(x=paste("CCA 1 (", format(100 *ii$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("CCA 2 (", format(100 *ii$cont[[1]][2,2], digits=4), "%)", sep=""))+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  stat_ellipse(geom = "polygon",data = st,aes(CCA1,CCA2,fill=fenzu$group),
               alpha=0.5)+theme_bw()+theme(panel.grid=element_blank())