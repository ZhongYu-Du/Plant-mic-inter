############Packages############
library(NST)
library(ggpubr)
library(ggpubr)
library(ggsignif)
library(ggplot2)

############Bgroup############
comm<-read.csv("MBBOTU.csv", header = T,row.names = 1)
group<-read.csv("Bgroup.csv", header = T,row.names = 1)
comm<-t(comm)

set.seed(123)
tnst <- tNST(comm = comm, group = group, dist.method = 'jaccard', null.model = 'PF', 
             rand = 1000, nworker = 4)
nst_group <- tnst$index.pair.grp
nst_group$group<-factor(nst_group$group,levels=c("CA","MA","CB","MB"))
options(scipen=200)
my_comparisons <- list(c("CA", "MA"), c("CB", "MB"),c("CA", "CB"), c("MA", "MB"))

p_NSTB<-ggboxplot(data = nst_group, x = 'group', y = 'NST.ij.ruzicka', color = 'group') +
  stat_compare_means(method = 'anova',label.y = 1.5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  geom_jitter(aes(color=group),position = position_jitter(0.2),
              satckdir="center", dotsize=0.9,stroke = 1,size=1,alpha=1)+
  coord_cartesian(ylim = c(0,1.5))+
  theme_bw()+theme(panel.grid=element_blank())

############Fgroup############
comm<-read.csv("FOTU.csv", header = T,row.names = 1)
group<-read.csv("Fgroup.csv", header = T,row.names = 1)
comm<-t(comm)

set.seed(123)
tnst <- tNST(comm = comm, group = group, dist.method = 'jaccard', null.model = 'PF', 
             rand = 1000, nworker = 4)
nst_group <- tnst$index.pair.grp
nst_group$group<-factor(nst_group$group,levels=c("CA","MA","CB","MB"))
options(scipen=200)
my_comparisons <- list(c("CA", "MA"), c("CB", "MB"),c("CA", "CB"), c("MA", "MB"))

p_NSTF<-ggboxplot(data = nst_group, x = 'group', y = 'NST.ij.ruzicka', color = 'group') +
  stat_compare_means(method = 'anova',label.y = 1.5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  geom_jitter(aes(color=group),position = position_jitter(0.2),
              satckdir="center", dotsize=0.9,stroke = 1,size=1,alpha=1)+
  coord_cartesian(ylim = c(0,1.5))+
  theme_bw()+theme(panel.grid=element_blank())
