###################Packages#################

library(psych)
library(reshape2)
library(ggplot2) 
library(randomForest)
library(patchwork)
library(vegan)

###################CSF#################
myro <- as.data.frame(read.csv("random forest\\CSF.csv", header=TRUE,row.names = 1))
spearman <- corr.test(myro[,1:8], myro[,9:13], method = 'spearman', adjust = 'none')

r <- data.frame(spearman$r)  
r$myro <- rownames(r)
r <- melt(r, id = 'myro')
spearman <- cbind(r)
spearman$myro<- factor(spearman$myro, c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP"))

p <- ggplot() +
  theme_bw() +
  geom_tile(data = spearman, aes(x = variable, y = myro, fill = value)) +
  scale_fill_gradientn(colors = c('#2D6DB1', 'white', '#DC1623'), limit = c(-1, 1)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black'), legend.key = element_blank(), 
        axis.text.x = element_text(color = 'black', angle =45, hjust = 1, vjust = 1), axis.text.y = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', x = '', fill = 'Correlation')

set.seed(123)
Ascomycota_forest <- randomForest(Ascomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Basidiomycota_forest <- randomForest(Basidiomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Mortierellomycota_forest <- randomForest(Mortierellomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
unclassified_forest <- randomForest(unclassified~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Rozellomycota_forest <- randomForest(Rozellomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)

Ascomycota<- data.frame(importance(Ascomycota_forest, scale = TRUE), check.names = FALSE)
Basidiomycota <- data.frame(importance(Basidiomycota_forest, scale = TRUE), check.names = FALSE)
Mortierellomycota <- data.frame(importance(Mortierellomycota_forest, scale = TRUE), check.names = FALSE)
unclassified <- data.frame(importance(unclassified_forest, scale = TRUE), check.names = FALSE)
Rozellomycota <- data.frame(importance(Rozellomycota_forest, scale = TRUE), check.names = FALSE)

IMportance_t <- data.frame(cbind(Ascomycota$`%IncMSE`,Basidiomycota$`%IncMSE`,Mortierellomycota$`%IncMSE`,
                                 unclassified$`%IncMSE`, Rozellomycota$`%IncMSE`))
colnames(IMportance_t) <- c("Ascomycota", "Basidiomycota", "Mortierellomycota",
                            "unclassified","Rozellomycota")
rownames(IMportance_t) <- c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP")
IMportance_t[IMportance_t<0] <- 0
IMportance_t[is.na(IMportance_t)] <- 0
write.csv(IMportance_t,file = "random forest\\IMportance_t1.CSV")

impdata <- read.csv("random forest\\IMportance_t1.CSV",header = TRUE) 
measure_name=setdiff(colnames(impdata),c('Items'))

data1=melt(impdata,           
           id.vars='Items',
           measure.vars=measure_name,
           variable.name = "sample",
           value.name = "expr")

p1 <- p+geom_point(data = data1, aes(x = sample, y = Items, size = expr*10), shape = 1) +
  scale_size_continuous(range = c(0,10)) + 
  labs(size = 'Importance (%)')

exp <- read.csv("random forest\\exp.CSV",header = TRUE)

exp$OTUs <- factor(exp$OTUs, levels = c("Ascomycota", "Basidiomycota", "Mortierellomycota",
                                                  "unclassified","Rozellomycota"))

p2 <- ggplot(exp, aes(OTUs, Values)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+coord_cartesian(ylim = c(0, 60))+
  ylab("Explained variation(%)")

##P3
aa <- randomForest(shannon ~ DOC+AK+AP+TAs+TSb+ASb+pH+SMC, data=myro, importance=TRUE,ntree=500)

oo=importance(aa,ntree=500)
barplot(oo)
write.csv(oo,file = "RF.csv")

RF <- read.csv("RF.CSV",header = TRUE)
randomForest(formula = shannon ~ ., data =myro, importance = TRUE,ntree = 500) 


RF$X <- factor(RF$X, c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP"))
p3 <- ggplot(RF, aes(X, X.IncMSE)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+coord_cartesian(xlim = c(0, 10))+
  ylab("Explained variation(%)")+
  coord_flip()

p2/p1|p2/p3


###################CSB#################
myro <- as.data.frame(read.csv("random forest\\CSB.csv", header=TRUE,row.names = 1))
spearman <- corr.test(myro[,1:8], myro[,9:13], method = 'spearman', adjust = 'none')

r <- data.frame(spearman$r)  
r$myro <- rownames(r)
r <- melt(r, id = 'myro')
spearman <- cbind(r)
spearman$myro<- factor(spearman$myro, c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP"))

p <- ggplot() +
  theme_bw() +
  geom_tile(data = spearman, aes(x = variable, y = myro, fill = value)) +
  scale_fill_gradientn(colors = c('#2D6DB1', 'white', '#DC1623'), limit = c(-1, 1)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black'), legend.key = element_blank(), 
        axis.text.x = element_text(color = 'black', angle =45, hjust = 1, vjust = 1), axis.text.y = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', x = '', fill = 'Correlation')

set.seed(123)
Ascomycota_forest <- randomForest(Ascomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Basidiomycota_forest <- randomForest(Basidiomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Mortierellomycota_forest <- randomForest(Mortierellomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
unclassified_forest <- randomForest(unclassified~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)
Rozellomycota_forest <- randomForest(Rozellomycota~DOC+AP+AK+TAs+TSb+ASb+pH+SMC, data =myro, importance = TRUE, ntree = 500)

Ascomycota<- data.frame(importance(Ascomycota_forest, scale = TRUE), check.names = FALSE)
Basidiomycota <- data.frame(importance(Basidiomycota_forest, scale = TRUE), check.names = FALSE)
Mortierellomycota <- data.frame(importance(Mortierellomycota_forest, scale = TRUE), check.names = FALSE)
unclassified <- data.frame(importance(unclassified_forest, scale = TRUE), check.names = FALSE)
Rozellomycota <- data.frame(importance(Rozellomycota_forest, scale = TRUE), check.names = FALSE)

IMportance_t <- data.frame(cbind(Ascomycota$`%IncMSE`,Basidiomycota$`%IncMSE`,Mortierellomycota$`%IncMSE`,
                                 unclassified$`%IncMSE`, Rozellomycota$`%IncMSE`))
colnames(IMportance_t) <- c("Ascomycota", "Basidiomycota", "Mortierellomycota",
                            "unclassified","Rozellomycota")
rownames(IMportance_t) <- c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP")
IMportance_t[IMportance_t<0] <- 0
IMportance_t[is.na(IMportance_t)] <- 0
write.csv(IMportance_t,file = "random forest\\IMportance_t1.CSV")

impdata <- read.csv("random forest\\IMportance_t1.CSV",header = TRUE) 
measure_name=setdiff(colnames(impdata),c('Items'))

data1=melt(impdata,           
           id.vars='Items',
           measure.vars=measure_name,
           variable.name = "sample",
           value.name = "expr")

p1 <- p+geom_point(data = data1, aes(x = sample, y = Items, size = expr*10), shape = 1) +
  scale_size_continuous(range = c(0,10)) + 
  labs(size = 'Importance (%)')

exp <- read.csv("random forest\\exp.CSV",header = TRUE)

exp$OTUs <- factor(exp$OTUs, levels = c("Ascomycota", "Basidiomycota", "Mortierellomycota",
                                        "unclassified","Rozellomycota"))

p2 <- ggplot(exp, aes(OTUs, Values)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+coord_cartesian(ylim = c(0, 60))+
  ylab("Explained variation(%)")

##P3
aa <- randomForest(shannon ~ DOC+AK+AP+TAs+TSb+ASb+pH+SMC, data=myro, importance=TRUE,ntree=500)

oo=importance(aa,ntree=500)
barplot(oo)
write.csv(oo,file = "RF.csv")

RF <- read.csv("RF.CSV",header = TRUE)
randomForest(formula = shannon ~ ., data =myro, importance = TRUE,ntree = 500) 


RF$X <- factor(RF$X, c("pH", "AK", "TAs","SMC","ASb","TSb","DOC","AP"))
p3 <- ggplot(RF, aes(X, X.IncMSE)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+coord_cartesian(xlim = c(0, 10))+
  ylab("Explained variation(%)")+
  coord_flip()

p2/p1|p2/p3
