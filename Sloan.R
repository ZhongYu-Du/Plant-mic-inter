############packages############
library(Hmisc)
library(minpack.lm)
library(stats4)
library(dplyr)

############tools############
sloan<-function(df){
  spp<-df
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  out_sloan<-data.frame(p,freq,freq.pred,pred.ci[,2:3])
  label<-paste("R^2 =",round(Rsqr,4),"\n","Nm =", round(coef(m.fit)*N,4) ,"\n","m =", round(coef(m.fit),4))
  sloan_label<-list(out_sloan,label)
  return(sloan_label)}

############SCB############
data_SCB<- read.csv('SCB.csv', row.names = 1)
data_SCB<-t(data_SCB)

##picture
slo_SCB<-sloan(data_SCB)[[1]]
slo_SCB = mutate(slo_SCB)
head(slo_SCB)
slo_SCB$group[slo_SCB[,2]<slo_SCB[,4]]="Low"
slo_SCB$group[slo_SCB[,2]>slo_SCB[,4]]="High"
slo_SCB$group[(slo_SCB[,2]>=slo_SCB[,4])&(slo_SCB[,2]<=slo_SCB[,5])]="Med"
slo_SCB$group <- factor(slo_SCB$group)

p_SCB<-ggplot() +
  geom_line(data = slo_SCB,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SCB,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SCB,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency") +
  geom_point(data = slo_SCB,aes(x=log(p),y=freq, color = group, fill = group),
             size = 2)+
  scale_color_manual(values = c("chocolate","darkorange4","blueviolet"))+
  annotate("text", x=-12,y=1, label=sloan(data_SCB)[[2]])+
  theme_bw()+
  geom_line(data = slo_SCB,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SCB,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SCB,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")

############SMB############
data_SMB<- read.csv('SMB.csv', row.names = 1)
data_SMB<-t(data_SMB)

##picture
slo_SMB<-sloan(data_SMB)[[1]]
slo_SMB = mutate(slo_SMB)
head(slo_SMB)
slo_SMB$group[slo_SMB[,2]<slo_SMB[,4]]="Low"
slo_SMB$group[slo_SMB[,2]>slo_SMB[,4]]="High"
slo_SMB$group[(slo_SMB[,2]>=slo_SMB[,4])&(slo_SMB[,2]<=slo_SMB[,5])]="Med"
slo_SMB$group <- factor(slo_SMB$group)

p_SMB<-ggplot() +
  geom_line(data = slo_SMB,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SMB,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SMB,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency") +
  geom_point(data = slo_SMB,aes(x=log(p),y=freq, color = group, fill = group),
             size = 2)+
  scale_color_manual(values = c("chocolate","darkorange4","blueviolet"))+
  annotate("text", x=-12,y=1, label=sloan(data_SMB)[[2]])+
  theme_bw()+
  geom_line(data = slo_SMB,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SMB,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SMB,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")

############SCF############
data_SCF<- read.csv('SCF.csv', row.names = 1)
data_SCF<-t(data_SCF)

##picture
slo_SCF<-sloan(data_SCF)[[1]]
slo_SCF = mutate(slo_SCF)
head(slo_SCF)
slo_SCF$group[slo_SCF[,2]<slo_SCF[,4]]="Low"
slo_SCF$group[slo_SCF[,2]>slo_SCF[,4]]="High"
slo_SCF$group[(slo_SCF[,2]>=slo_SCF[,4])&(slo_SCF[,2]<=slo_SCF[,5])]="Med"
slo_SCF$group <- factor(slo_SCF$group)

p_SCF<-ggplot() +
  geom_line(data = slo_SCF,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SCF,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SCF,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency") +
  geom_point(data = slo_SCF,aes(x=log(p),y=freq, color = group, fill = group),
             size = 2)+
  scale_color_manual(values = c("chocolate","darkorange4","blueviolet"))+
  annotate("text", x=-12,y=1, label=sloan(data_SCF)[[2]])+
  theme_bw()+
  geom_line(data = slo_SCF,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SCF,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SCF,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")

############SMF############
data_SMF<- read.csv('SMF.csv', row.names = 1)
data_SMF<-t(data_SMF)

##picture
slo_SMF<-sloan(data_SMF)[[1]]
slo_SMF = mutate(slo_SMF)
head(slo_SMF)
slo_SMF$group[slo_SMF[,2]<slo_SMF[,4]]="Low"
slo_SMF$group[slo_SMF[,2]>slo_SMF[,4]]="High"
slo_SMF$group[(slo_SMF[,2]>=slo_SMF[,4])&(slo_SMF[,2]<=slo_SMF[,5])]="Med"
slo_SMF$group <- factor(slo_SMF$group)

p_SMF<-ggplot() +
  geom_line(data = slo_SMF,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SMF,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SMF,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency") +
  geom_point(data = slo_SMF,aes(x=log(p),y=freq, color = group, fill = group),
             size = 2)+
  scale_color_manual(values = c("chocolate","darkorange4","blueviolet"))+
  annotate("text", x=-12,y=1, label=sloan(data_SMF)[[2]])+
  theme_bw()+
  geom_line(data = slo_SMF,aes(x=log(p),y=freq.pred),size = 1.5,linetype = 1)+
  geom_line(data = slo_SMF,aes(x=log(p),y=Lower),size = 1.5,linetype = 2, color="darkgreen")+
  geom_line(data = slo_SMF,aes(x=log(p),y=Upper),size = 1.5,linetype = 2, color="darkgreen")+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")




