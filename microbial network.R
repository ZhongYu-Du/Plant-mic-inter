setwd()
##############packages############
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(Biostrings)

##############Bacteria################
##############B_M################
taxonomy = read.csv("B_tax.csv", row.names=1,header=T)
otu_BM = read.csv("B_M.csv", row.names=1,header=T)
ps_BM = phyloseq(otu_table(as.matrix(otu_BM), taxa_are_rows=TRUE),tax_table(as.matrix(taxonomy)))

###############network by phylum###############
result_netClu_BM = corMicro (ps = ps_BM,
                               N = 300,
                               method.scale = "TMM",
                               r.threshold=0.60,
                               p.threshold=0.01,
                               method = "pearson")

ps_net_BM = result_netClu_BM[[3]]
asv_table_BM = as.data.frame(t(vegan_otu(ps_net_BM)))
tax_BM = as.data.frame(vegan_tax(ps_net_BM))
cor_netClu_BM = result_netClu_BM[1][[1]]
cor_netClu_BM[is.na(cor_netClu_BM)]=0
result2_netClu_BM <- model_Gephi.2(cor = cor_netClu_BM,
                                     method = "cluster_fast_greedy",
                                     seed = 12)
node_netClu_BM = result2_netClu_BM[[1]]

#---node annotation#-----------
nodes_netClu_BM = nodeadd(plotcord =node_netClu_BM,otu_table = asv_table_BM,tax_table = tax_BM)
netClu_BM<-result2_netClu_BM[[3]]
row.names(netClu_BM) = netClu_BM$ID
nodeG_BM = merge(nodes_netClu_BM,netClu_BM,by = "row.names",all  =FALSE)

#-----edge#--------
edge_netClu_BM = edgeBuild(cor = cor_netClu_BM,node = node_netClu_BM)
colnames(edge_netClu_BM)[8] = "cor"
edge_netClu_BM$weight<-abs(edge_netClu_BM$weight)
edge_netClu_BM<-edge_netClu_BM[,-1:-2]
names(edge_netClu_BM)[1:2]<-c('Source','Target')

#-----remove isolated points#--------
s_t_BM<-as.data.frame(c(edge_netClu_BM$Source,edge_netClu_BM$Target))
s_t_BM<-distinct(s_t_BM)
names(s_t_BM)<-'ID'
nd_BM<-merge(s_t_BM, nodeG_BM,  by= 'ID' , all.x = TRUE )

write.csv(edge_netClu_BM, 'edge_BM.csv', row.names = F)
write.csv(nd_BM, 'node_BM.csv', row.names = F)

##############zipi analysis################
result_zipi_BM = corMicro (ps = ps_BM,
                             N = 300,
                             r.threshold=0.60,
                             p.threshold=0.01,
                             method = "pearson")

cor_zipi_BM = result_zipi_BM[[1]]
result2_zipi_BM = nodeEdge(cor = cor_zipi_BM)
edge_zipi_BM = result2_zipi_BM[[1]]
node_zipi_BM = result2_zipi_BM[[2]]
igraph_zipi_BM  = igraph::graph_from_data_frame(edge_zipi_BM, directed = FALSE, vertices = node_zipi_BM)
res_zipi_BM <- zipi_text(igraph = igraph_zipi_BM,method = "cluster_fast_greedy")
res_zipi_data_BM<-merge(res_zipi_BM,taxonomy,by='row.names')
row.names(res_zipi_data_BM)<-res_zipi_data_BM[,1]
res_zipi_data_BM<-res_zipi_data_BM[,-1]
zipi_BM_p <- zipi_plot(res_zipi_data_BM) + ggrepel::geom_text_repel(data = res_zipi_data_BM, 
                                                                                       aes(x = p, y = z, color = Phylum, label = label), 
                                                                                       size = 4)

##############network topology################
dat_np_BM = net_properties.4(igraph_zipi_BM,n.hub = T)
result_rdn_BM = random_Net_compate(igraph = igraph_zipi_BM, type = "gnm", step = 10000, netName = layout)
p_rdn_BM<-result_rdn_BM[[1]]
rdn_BM = result_rdn_BM[[4]]
rdn_BM<-as.data.frame(rdn_BM)
rdn_BM<-rownames_to_column(rdn_BM)
colnames(rdn_BM)<-c('properties','BM_net','z_mean_BM','z_sd_BM')


##############B_C################
otu_BC = read.csv("B_C.csv", row.names=1,header=T)
ps_BC = phyloseq(otu_table(as.matrix(otu_BC), taxa_are_rows=TRUE),tax_table(as.matrix(taxonomy)))

###############network by phylum###############
result_netClu_BC = corMicro (ps = ps_BC,
                             N = 300,
                             method.scale = "TMM",
                             r.threshold=0.60,
                             p.threshold=0.01,
                             method = "pearson")

ps_net_BC = result_netClu_BC[[3]]
asv_table_BC = as.data.frame(t(vegan_otu(ps_net_BC)))
tax_BC = as.data.frame(vegan_tax(ps_net_BC))
cor_netClu_BC = result_netClu_BC[1][[1]]
cor_netClu_BC[is.na(cor_netClu_BC)]=0
result2_netClu_BC <- model_Gephi.2(cor = cor_netClu_BC,
                                   method = "cluster_fast_greedy",
                                   seed = 12)
node_netClu_BC = result2_netClu_BC[[1]]

#---node annotation#-----------
nodes_netClu_BC = nodeadd(plotcord =node_netClu_BC,otu_table = asv_table_BC,tax_table = tax_BC)
netClu_BC<-result2_netClu_BC[[3]]
row.names(netClu_BC) = netClu_BC$ID
nodeG_BC = merge(nodes_netClu_BC,netClu_BC,by = "row.names",all  =FALSE)

#-----edge#--------
edge_netClu_BC = edgeBuild(cor = cor_netClu_BC,node = node_netClu_BC)
colnames(edge_netClu_BC)[8] = "cor"
edge_netClu_BC$weight<-abs(edge_netClu_BC$weight)
edge_netClu_BC<-edge_netClu_BC[,-1:-2]
names(edge_netClu_BC)[1:2]<-c('Source','Target')

#-----Remove isolated points#--------
s_t_BC<-as.data.frame(c(edge_netClu_BC$Source,edge_netClu_BC$Target))
s_t_BC<-distinct(s_t_BC)
names(s_t_BC)<-'ID'
nd_BC<-merge(s_t_BC, nodeG_BC,  by= 'ID' , all.x = TRUE )

write.csv(edge_netClu_BC, 'edge_BC.csv', row.names = F)
write.csv(nd_BC, 'node_BC.csv', row.names = F)

##############zipi analysis################
result_zipi_BC = corMicro (ps = ps_BC,
                           N = 300,
                           r.threshold=0.60,
                           p.threshold=0.01,
                           method = "pearson")

cor_zipi_BC = result_zipi_BC[[1]]
result2_zipi_BC = nodeEdge(cor = cor_zipi_BC)
edge_zipi_BC = result2_zipi_BC[[1]]
node_zipi_BC = result2_zipi_BC[[2]]
igraph_zipi_BC  = igraph::graph_from_data_frame(edge_zipi_BC, directed = FALSE, vertices = node_zipi_BC)
res_zipi_BC <- zipi_text(igraph = igraph_zipi_BC,method = "cluster_fast_greedy")
res_zipi_data_BC<-merge(res_zipi_BC,taxonomy,by='row.names')
row.names(res_zipi_data_BC)<-res_zipi_data_BC[,1]
res_zipi_data_BC<-res_zipi_data_BC[,-1]
zipi_BC_p <- zipi_plot(res_zipi_data_BC) + ggrepel::geom_text_repel(data = res_zipi_data_BC, 
                                                                    aes(x = p, y = z, color = Phylum, label = label), 
                                                                    size = 4)

##############network topology################
dat_np_BC = net_properties.4(igraph_zipi_BC,n.hub = T)
result_rdn_BC = random_Net_compate(igraph = igraph_zipi_BC, type = "gnm", step = 10000, netName = layout)
p_rdn_BC<-result_rdn_BC[[1]]
rdn_BC = result_rdn_BC[[4]]
rdn_BC<-as.data.frame(rdn_BC)
rdn_BC<-rownames_to_column(rdn_BC)
colnames(rdn_BC)<-c('properties','BC_net','z_mean_BC','z_sd_BC')

##############Fungi################
##############F_M################
taxonomy = read.csv("F_tax.csv", row.names=1,header=T)
otu_FM = read.csv("F_M.csv", row.names=1,header=T)
ps_FM = phyloseq(otu_table(as.matrix(otu_FM), taxa_are_rows=TRUE),tax_table(as.matrix(taxonomy)))

###############network by phylum###############
result_netClu_FM = corMicro (ps = ps_FM,
                             N = 300,
                             method.scale = "TMM",
                             r.threshold=0.60,
                             p.threshold=0.01,
                             method = "pearson")

ps_net_FM = result_netClu_FM[[3]]
asv_table_FM = as.data.frame(t(vegan_otu(ps_net_FM)))
tax_FM = as.data.frame(vegan_tax(ps_net_FM))
cor_netClu_FM = result_netClu_FM[1][[1]]
cor_netClu_FM[is.na(cor_netClu_FM)]=0
result2_netClu_FM <- model_Gephi.2(cor = cor_netClu_FM,
                                   method = "cluster_fast_greedy",
                                   seed = 12)
node_netClu_FM = result2_netClu_FM[[1]]

#---node annotation#-----------
nodes_netClu_FM = nodeadd(plotcord =node_netClu_FM,otu_table = asv_table_FM,tax_table = tax_FM)
netClu_FM<-result2_netClu_FM[[3]]
row.names(netClu_FM) = netClu_FM$ID
nodeG_FM = merge(nodes_netClu_FM,netClu_FM,by = "row.names",all  =FALSE)

#-----edge#--------
edge_netClu_FM = edgeBuild(cor = cor_netClu_FM,node = node_netClu_FM)
colnames(edge_netClu_FM)[8] = "cor"
edge_netClu_FM$weight<-abs(edge_netClu_FM$weight)
edge_netClu_FM<-edge_netClu_FM[,-1:-2]
names(edge_netClu_FM)[1:2]<-c('Source','Target')

#-----remove isolated points#--------
s_t_FM<-as.data.frame(c(edge_netClu_FM$Source,edge_netClu_FM$Target))
s_t_FM<-distinct(s_t_FM)
names(s_t_FM)<-'ID'
nd_FM<-merge(s_t_FM, nodeG_FM,  by= 'ID' , all.x = TRUE )

write.csv(edge_netClu_FM, 'edge_FM.csv', row.names = F)
write.csv(nd_FM, 'node_FM.csv', row.names = F)

##############zipi anlysis################
result_zipi_FM = corMicro (ps = ps_FM,
                           N = 300,
                           r.threshold=0.60,
                           p.threshold=0.01,
                           method = "pearson")

cor_zipi_FM = result_zipi_FM[[1]]
result2_zipi_FM = nodeEdge(cor = cor_zipi_FM)
edge_zipi_FM = result2_zipi_FM[[1]]
node_zipi_FM = result2_zipi_FM[[2]]
igraph_zipi_FM  = igraph::graph_from_data_frame(edge_zipi_FM, directed = FALSE, vertices = node_zipi_FM)
res_zipi_FM <- zipi_text(igraph = igraph_zipi_FM,method = "cluster_fast_greedy")
res_zipi_data_FM<-merge(res_zipi_FM,taxonomy,by='row.names')
row.names(res_zipi_data_FM)<-res_zipi_data_FM[,1]
res_zipi_data_FM<-res_zipi_data_FM[,-1]
zipi_FM_p <- zipi_plot(res_zipi_data_FM) + ggrepel::geom_text_repel(data = res_zipi_data_FM, 
                                                                    aes(x = p, y = z, color = Phylum, label = label), 
                                                                    size = 4)

##############network topology################
dat_np_FM = net_properties.4(igraph_zipi_FM,n.hub = T)
result_rdn_FM = random_Net_compate(igraph = igraph_zipi_FM, type = "gnm", step = 10000, netName = layout)
p_rdn_FM<-result_rdn_FM[[1]]
rdn_FM = result_rdn_FM[[4]]
rdn_FM<-as.data.frame(rdn_FM)
rdn_FM<-rownames_to_column(rdn_FM)
colnames(rdn_FM)<-c('properties','FM_net','z_mean_FM','z_sd_FM')

##############F_C################
otu_FC = read.csv("B_C.csv", row.names=1,header=T)
ps_FC = phyloseq(otu_table(as.matrix(otu_FC), taxa_are_rows=TRUE),tax_table(as.matrix(taxonomy)))

###############network by phylum###############
result_netClu_FC = corMicro (ps = ps_FC,
                             N = 300,
                             method.scale = "TMM",
                             r.threshold=0.60,
                             p.threshold=0.01,
                             method = "pearson")

ps_net_FC = result_netClu_FC[[3]]
asv_table_FC = as.data.frame(t(vegan_otu(ps_net_FC)))
tax_FC = as.data.frame(vegan_tax(ps_net_FC))
cor_netClu_FC = result_netClu_FC[1][[1]]
cor_netClu_FC[is.na(cor_netClu_FC)]=0
result2_netClu_FC <- model_Gephi.2(cor = cor_netClu_FC,
                                   method = "cluster_fast_greedy",
                                   seed = 12)
node_netClu_FC = result2_netClu_FC[[1]]

#---node annotation#-----------
nodes_netClu_FC = nodeadd(plotcord =node_netClu_FC,otu_table = asv_table_FC,tax_table = tax_FC)
netClu_FC<-result2_netClu_FC[[3]]
row.names(netClu_FC) = netClu_FC$ID
nodeG_FC = merge(nodes_netClu_FC,netClu_FC,by = "row.names",all  =FALSE)

#-----edge#--------
edge_netClu_FC = edgeBuild(cor = cor_netClu_FC,node = node_netClu_FC)
colnames(edge_netClu_FC)[8] = "cor"
edge_netClu_FC$weight<-abs(edge_netClu_FC$weight)
edge_netClu_FC<-edge_netClu_FC[,-1:-2]
names(edge_netClu_FC)[1:2]<-c('Source','Target')

#-----Remove isolated points#--------
s_t_FC<-as.data.frame(c(edge_netClu_FC$Source,edge_netClu_FC$Target))
s_t_FC<-distinct(s_t_FC)
names(s_t_FC)<-'ID'
nd_FC<-merge(s_t_FC, nodeG_FC,  by= 'ID' , all.x = TRUE )

write.csv(edge_netClu_FC, 'edge_FC.csv', row.names = F)
write.csv(nd_FC, 'node_FC.csv', row.names = F)

##############zipi anlysis################
result_zipi_FC = corMicro (ps = ps_FC,
                           N = 300,
                           r.threshold=0.60,
                           p.threshold=0.01,
                           method = "pearson")

cor_zipi_FC = result_zipi_FC[[1]]
result2_zipi_FC = nodeEdge(cor = cor_zipi_FC)
edge_zipi_FC = result2_zipi_FC[[1]]
node_zipi_FC = result2_zipi_FC[[2]]
igraph_zipi_FC  = igraph::graph_from_data_frame(edge_zipi_FC, directed = FALSE, vertices = node_zipi_FC)
res_zipi_FC <- zipi_text(igraph = igraph_zipi_FC,method = "cluster_fast_greedy")
res_zipi_data_FC<-merge(res_zipi_FC,taxonomy,by='row.names')
row.names(res_zipi_data_FC)<-res_zipi_data_FC[,1]
res_zipi_data_FC<-res_zipi_data_FC[,-1]
zipi_FC_p <- zipi_plot(res_zipi_data_FC) + ggrepel::geom_text_repel(data = res_zipi_data_FC, 
                                                                    aes(x = p, y = z, color = Phylum, label = label), 
                                                                    size = 4)

##############network topology################
dat_np_FC = net_properties.4(igraph_zipi_FC,n.hub = T)
result_rdn_FC = random_Net_compate(igraph = igraph_zipi_FC, type = "gnm", step = 10000, netName = layout)
p_rdn_FC<-result_rdn_FC[[1]]
rdn_FC = result_rdn_FC[[4]]
rdn_FC<-as.data.frame(rdn_FC)
rdn_FC<-rownames_to_column(rdn_FC)
colnames(rdn_FC)<-c('properties','FC_net','z_mean_FC','z_sd_FC')

############complexity############
multimerge<-function(dat=list(),...){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  dat[[1]]<-NULL
  for(i in dat){
    mergedat<-merge(mergedat,i,...)
  }
  return(mergedat)
}
mul_mg<-multimerge(list(rdn_BM,rdn_BC,rdn_FM,rdn_FC))
write.csv(mul_mg,"microbial network\\mul_mg.csv")

############stability############
library(patchwork)
library(tidyverse)
library(pulsar)
##############Bacteria################
medata = read.csv("group.csv",header = T, row.names=1)
otutab = read.csv("Botu.csv",header = T, row.names=1)
taxonomy = read.csv("B_tax.csv",header = T, row.names=1)

ps = phyloseq(sample_data(medata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))

tab.r = network.pip(
  ps = ps,
  N = 300,
  layout_net = "model_maptree2",
  r.threshold = 0.60,
  p.threshold = 0.01,
  method = "spearman",
  clu_method = "cluster_fast_greedy")
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab 

##############robustness################
##robustness1
res_rn <- Robustness.Targeted.removal(ps = ps,
                                      corg = cortab,
                                      degree = TRUE,
                                      zipi = FALSE)
res_rn_data <- res_rn[[2]]

p_rn = ggplot(res_rn_data[res_rn_data$weighted == "unweighted", ], aes(x = Number.hub.removed, 
                                                                       y = remain.mean, color = Group)) + 
  geom_line()+ 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 0.2) + 
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+theme_classic()+
  theme(panel.background = element_rect(color='black'))

##robustness2
res_rn2 = Robustness.Random.removal(ps = ps,
                                    corg = cortab,
                                    Top = 0)
res_rn_data2 = res_rn2[[2]]

p_rn2 = ggplot(res_rn_data2[res_rn_data2$weighted == "unweighted", ], aes(x = Proportion.removed, 
                                                                          y = remain.mean, color = Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 0.2) + 
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+theme_classic()+
  theme(panel.background = element_rect(color='black'))

##############survivability################
res_kh = natural.con.microp (
  ps = ps,
  corg = cortab,
  norm = TRUE,
  end = 200,
  start = 0)

res_kh_data = res_kh[[2]]

p_kh = ggplot(res_kh_data, aes(Num.of.remove.nodes, Natural.connectivity, 
                               color = Group)) + geom_point(alpha = 0.1) + 
  geom_smooth(level = F) + 
  theme_set(theme_bw(base_size = 25)) + theme(panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) + 
  theme_classic()+
  theme(panel.background = element_rect(color='black'))


##############Fungi################
medata = read.csv("group.csv",header = T, row.names=1)
otutab = read.csv("Fotu.csv",header = T, row.names=1)
taxonomy = read.csv("F_tax.csv",header = T, row.names=1)

ps = phyloseq(sample_data(medata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))

tab.r = network.pip(
  ps = ps,
  N = 300,
  layout_net = "model_maptree2",
  r.threshold = 0.60,
  p.threshold = 0.01,
  method = "spearman",
  clu_method = "cluster_fast_greedy")
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab 

##############robustness################
##robustness1
res_rn <- Robustness.Targeted.removal(ps = ps,
                                      corg = cortab,
                                      degree = TRUE,
                                      zipi = FALSE)
res_rn_data <- res_rn[[2]]

p_rn = ggplot(res_rn_data[res_rn_data$weighted == "unweighted", ], aes(x = Number.hub.removed, 
                                                                       y = remain.mean, color = Group)) + 
  geom_line()+ 
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 0.2) + 
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+theme_classic()+
  theme(panel.background = element_rect(color='black'))

##robustness2
res_rn2 = Robustness.Random.removal(ps = ps,
                                    corg = cortab,
                                    Top = 0)
res_rn_data2 = res_rn2[[2]]

p_rn2 = ggplot(res_rn_data2[res_rn_data2$weighted == "unweighted", ], aes(x = Proportion.removed, 
                                                                          y = remain.mean, color = Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin = remain.mean - remain.sd, 
                      ymax = remain.mean + remain.sd), size = 0.2) + 
  xlab("Proportion of species removed") + 
  ylab("Proportion of species remained")+theme_classic()+
  theme(panel.background = element_rect(color='black'))

##############survivability################
res_kh = natural.con.microp (
  ps = ps,
  corg = cortab,
  norm = TRUE,
  end = 200,
  start = 0)

res_kh_data = res_kh[[2]]

p_kh = ggplot(res_kh_data, aes(Num.of.remove.nodes, Natural.connectivity, 
                               color = Group)) + geom_point(alpha = 0.1) + 
  geom_smooth(level = F) + 
  theme_set(theme_bw(base_size = 25)) + theme(panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) + 
  theme_classic()+
  theme(panel.background = element_rect(color='black'))
