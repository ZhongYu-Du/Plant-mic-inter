############Packages############
library(ggplot2)
library(corrplot)
library(Hmisc)
library(dplyr)
library(igraph)
library(linkET)
############Correlation############
##Probiotic correlation
b<-read.csv("Probiotic correlation.csv", row.names = 1)

b1<-cor(b)
b2<-b1[1:6,7:13]
b3<-t(b2)

colb3 <- colorRampPalette(c("#8A2BE2","#EE6A50")) 
resb1 <- cor.mtest(b1, conf.level = .95)
resb2<-resb1$p[1:6,7:13]

p_b<-corrplot(b2,is.corr = FALSE,method = "color",
            p.mat = resb2,insig = "label_sig", sig.level = c(.001, .05), pch.cex = 1.8, 
            pch.col = "black", col = colb3(100),
            addgrid.col="black",tl.srt =90,tl.col="black",cl.cex=0.8,tl.cex=0.8)

############Mantel analysis############
############Bacteria############
soi<-read.csv("Soil.csv", row.names = 1)
cor<-rcorr(as.matrix(soi), type = "pearson")
data1<-read.csv("Bment.csv", row.names = 1)

mantel <- mantel_test(data1, soi,
                      spec_select = list(Keystone = 1:12,
                                         Biomarker = 13:33)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
         dd = cut(r, breaks = c(-Inf,0,Inf),
                  labels = c("r < Negative", "r > Positive")))

P_bment<-qcorrplot(cor, type = "upper", show.diag =FALSE,
          corr.test = T, cluster.type = "all") +
  geom_square() +
  scale_fill_gradient2(high = '#CC0000', low = '#0088A8')+
  geom_mark(r = NA,
            only_mark = T,
            color = "white", 
            size=4)+ 
  geom_couple(aes(colour = pd, 
                  size = rd, 
                  linetype = dd),
              data = mantel,
              curvature = nice_curvature()) + 
  scale_size_manual(values = c(1, 1.5, 2.5))+
  scale_colour_manual(values = c("#D2691E","#006400","grey")) +
  scale_linetype_manual(values = c("dotted", "solid"))+
  guides(fill = guide_colourbar(title = "corr", order = 1),
         colour = guide_legend(title = "Mantel's p", order = 2),
         size = guide_legend(title = "Mantel's r", order = 3),
         linetype = guide_legend("", order = 4))


############Fungi############
data2<-read.csv("Fment.csv", row.names = 1)

mantel <- mantel_test(data2, soi,
                      spec_select = list(Keystone = 1:7,
                                         Biomarker = 8:22)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
         dd = cut(r, breaks = c(-Inf,0,Inf),
                  labels = c("r < Negative", "r > Positive")))


P_fment<-qcorrplot(cor, type = "upper", show.diag =FALSE,
                   corr.test = T, cluster.type = "all") +
  geom_square() +
  scale_fill_gradient2(high = '#CC0000', low = '#0088A8')+
  geom_mark(r = NA,
            only_mark = T,
            color = "white", 
            size=4)+ 
  geom_couple(aes(colour = pd, 
                  size = rd, 
                  linetype = dd),
              data = mantel,
              curvature = nice_curvature()) + 
  scale_size_manual(values = c(1, 1.5, 2.5))+
  scale_colour_manual(values = c("#D2691E","#006400","grey")) +
  scale_linetype_manual(values = c("dotted", "solid"))+
  guides(fill = guide_colourbar(title = "corr", order = 1),
         colour = guide_legend(title = "Mantel's p", order = 2),
         size = guide_legend(title = "Mantel's r", order = 3),
         linetype = guide_legend("", order = 4))

