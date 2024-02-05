library(ggplot2)
library(tidyverse)
library(agricolae)
library(car)
library(reshape2)

##Multiple comparisons
df <- read.table("data.txt",header = T, row.names = 1)

shapiro.test(df$values)#normal distribution
bartlett.test(df$values ~ df$group, data=df)#Homogeneity Of Variance

variance<-aov(values ~ group, data=df)
MC <- LSD.test(variance,"group", p.adj="none")#结果显示：标记字母法out$group
mark <- data.frame(MC$groups)
mark$group = row.names(mark)

ggplot(df,aes(group, value, color=group))+
  geom_boxplot()+
  geom_jitter()+
  scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781","#355783"))+
  geom_text(data=mark,
            aes(x=group,y=yield+2,label=groups),
            color="black",
            size = 5,
            fontface="bold")+#添加字母标记
  xlab("")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14,face = 'bold'))

##Nonparametric test
df<-read.csv("data2.csv",header = T, row.names = 1)

shapiro.test(df$value)#Non-normal distribution
leveneTest(df$value ~ df$condition, data=df)#Non-homogeneity of variance

le.dun <- dunnTest(value ~ condition,method="bonferroni",df)
le.le <- cldList(P.adj~Comparison,
                 data = le.dun$res[order(le.dun$res$Z,decreasing = T),],
                 threshold = 0.05 )

LE.sd<-aggregate(df$value,by=list(df$condition),FUN=sd)
LE.mn<-aggregate(df$value,by=list(df$condition),FUN=mean)
LE<-merge(LE.mn,LE.sd)
LE<-merge(LE,le.le)
colnames(LE)<-c('condition','mean','sd','letter')

ggplot(data=LE,aes(condition, mean))+
  labs(y="Leaf length")+
  geom_bar(stat="identity",position="dodge",aes(fill=condition),width =0.8)+
  scale_fill_brewer(palette="Paired")+
  geom_errorbar(aes(ymax=mean+sd,ymin=mean-sd),position=position_dodge(0.9),width=0.15)+
  geom_text(data=LE,
            aes(x=condition,y=mean+sd+0.3,label=letter))+
  theme_classic()+
  theme(panel.background = element_rect(color='black'))
