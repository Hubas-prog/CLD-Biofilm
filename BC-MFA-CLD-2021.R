########################### SCRIPT ################################
# Chlordecone-contaminated epipelic biofilms show exopolymeric secretions with increased adsorption capacities in Guadeloupe tropical rivers
# by : Hubas C* et al. 
# Submitted to STOTEN
# *Corresponding author : cedric.hubas@mnhn.fr
# Credits : Script = C.H., Data = C.H. and all other co-authors
########################### SCRIPT ################################

# In house R functions
# source the following function : https://github.com/Hubas-prog/BC-MFA/blob/main/BC-MFA_custom_function.R

source("https://github.com/Hubas-prog/BC-MFA/blob/main/BC-MFA_custom_function.R")
#source("/Users/cedric.hubas/BC-MFA/BC-MFA_custom_function.R")

# Packages
library(ggplot2)
library(scales)
library(corrplot)

# aesthetics
my.palette <- colorRampPalette(c("red3","orange","yellow","green3","royalblue"))
my.theme<-theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                axis.text=element_text(size=14),
                axis.title=element_text(size=14))
row.pal<-colorRampPalette(c("green3","orange","red3"))
col.pal<-colorRampPalette(c("#FF0000","#310000","#FF9300"))

# biochemical data upload (aggregated)
mean.PERCENT.PIG<-read.csv("mean.aggregate.pigment.csv")
mean.PERCENT.MONO<-read.csv("mean.aggregate.monosaccharides.csv")
mean.PERCENT.FATTY<-read.csv("mean.aggregate.fatty.csv",row.names=1)

# Checking paired biochemical data
dim(mean.PERCENT.PIG)
names(mean.PERCENT.PIG)
dim(mean.PERCENT.MONO)
names(mean.PERCENT.MONO)
dim(mean.PERCENT.FATTY)
names(mean.PERCENT.FATTY)
sitecheck<-data.frame(mean.PERCENT.PIG[,1:3],mean.PERCENT.MONO[,1:3],mean.PERCENT.FATTY[-36,1:3])
sitecheck

# merge biochemical datasets
names(mean.PERCENT.PIG)
names(mean.PERCENT.MONO)
names(mean.PERCENT.FATTY)
colnames(mean.PERCENT.FATTY)<-gsub(".n","-n",colnames(mean.PERCENT.FATTY),fix=T)
colnames(mean.PERCENT.FATTY)<-gsub("X","",colnames(mean.PERCENT.FATTY),fix=T)

mean.PERCENT.PIG$group<-paste(mean.PERCENT.PIG$Group.1,
                              mean.PERCENT.PIG$Group.2,
                              mean.PERCENT.PIG$Group.3)
mean.PERCENT.MONO$group<-paste(mean.PERCENT.MONO$Group.1,
                               mean.PERCENT.MONO$Group.2,
                               mean.PERCENT.MONO$Group.3)
mean.PERCENT.FATTY$group<-paste(mean.PERCENT.FATTY$Group.1,
                                mean.PERCENT.FATTY$Group.2,
                                mean.PERCENT.FATTY$Group.3)

data.pig.mono<-merge(mean.PERCENT.PIG[,4:27],mean.PERCENT.MONO[,4:13],by="group")
all.data<-merge(data.pig.mono,mean.PERCENT.FATTY[,4:44],by="group")
names(all.data)
all.data$river<-substr(all.data$group,1,3)
all.data$site<-substr(all.data$group,5,6)
all.data$time<-substr(all.data$group,8,10)

# Data upload CLD biochemical
pol.data.all<-read.csv("cld.data.biochemical.csv",sep=",")
data.frame(all.data$group,
           pol.data.all$group)

# biochemical data upload (list) 
newdata.pig.list<-read.csv("data.pig.list.csv")
newdata.mono.list<-read.csv("data.mono.list.csv")
newdata.fatty.list<-read.csv("data.fatty.list.csv",row.names = 1)

# T-RFLP data upload
archeaL<-read.csv("Arch.csv",h=T,sep=";")
bacteriaL<-read.csv("Bact.csv",h=T,sep=";") 
eukaryaL<-read.csv("Euk.csv",h=T,sep=";")

dim(archeaL) ; dim(bacteriaL) ; dim(eukaryaL)
data.frame(archeaL[,1:3],bacteriaL[,1:3],eukaryaL[,1:3])

ARC.L<-archeaL[,6:173]/(apply(archeaL[,6:173],1,sum))*100
BAC.L<-bacteriaL[,6:184]/(apply(bacteriaL[,6:184],1,sum))*100
EUK.L<-eukaryaL[,6:95]

colnames(ARC.L)<-gsub("X","A",colnames(ARC.L))
colnames(BAC.L)<-gsub("X","B",colnames(BAC.L))
colnames(EUK.L)<-gsub("X","E",colnames(EUK.L))

vector.arc<-apply(ARC.L,2,sum)!=0
vector.bac<-apply(BAC.L,2,sum)!=0
vector.euk<-apply(EUK.L,2,sum)!=0

newdata.trflp<-data.frame(bacteriaL[,2:4],
                          ARC.L[,vector.arc],
                          BAC.L[,vector.bac],
                          EUK.L[,vector.euk])
newdata.trflp$River<-gsub("RCA","GRC",newdata.trflp$River,fix=T)
newdata.trflp$River<-gsub("GRCR","RCA",newdata.trflp$River,fix=T)

mean.PERCENT.TRFLP<-aggregate(newdata.trflp[,4:386],
                              by=list(newdata.trflp$River,
                                      newdata.trflp$Position,
                                      newdata.trflp$Time),
                              mean)

# Data upload CLD molecular
pol.data.molec<-read.csv("cld.data.molecular.csv",sep=",")
data.frame(mean.PERCENT.TRFLP[,1:3],
           pol.data.molec[,1:3])

# Correlations
give.cor<-function(cor.data,alpha.error) {
  testRes<-cor.mtest(cor.data, conf.level = 0.95)
  data.frame(variables=rownames(cor(cor.data)),
             pearson=data.frame(cor(cor.data))[,dim(cor.data)[2]],
             pvalue=data.frame(testRes)[,dim(cor.data)[2]])
}

cor.data.pig<-data.frame(all.data[all.data$site=="AV",c(2:24)],
           CLD=pol.data.all$CLD[all.data$site=="AV"])

cor.data.mono<-data.frame(all.data[all.data$site=="AV",c(25:33)],
                         CLD=pol.data.all$CLD[all.data$site=="AV"])

cor.data.fatty<-data.frame(all.data[all.data$site=="AV",c(34:73)],
                         CLD=pol.data.all$CLD[all.data$site=="AV"])

give.cor(cor.data.pig,0.95)
give.cor(cor.data.mono,0.95)
give.cor(cor.data.fatty,0.95)

# Plots carbohydrates
plot1<-ggplot(newdata.mono.list,
              aes(y=percent,x=carbohydrate))+
  geom_boxplot()+
  ylab("Relative concentration (%)")+
  xlab("")+
  theme_bw()+
  my.theme


M<-ggplot(all.data[all.data$site=="AV",],
          aes(x=log(Myo.inositol),
              y=pol.data.all$CLD[pol.data.all$site=="AV"],
              label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("log(Myo-inositol relative concentration (%))")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

all.data$ratio1<-all.data$Myo.inositol/all.data$Scyllo.inositol

MS<-ggplot(all.data[all.data$site=="AV",],
           aes(x=log(ratio1),
               y=pol.data.all$CLD[pol.data.all$site=="AV"],
               label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("log(Myo-/Scyllo-inositol ratio)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

all.data$ratio2<-all.data$Myo.inositol/all.data$Glucuronic.acid

MG<-ggplot(all.data[all.data$site=="AV",],
           aes(x=log(ratio2),
               y=pol.data.all$CLD[pol.data.all$site=="AV"],
               label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("log(Myo-inositol/Glucuronic acid ratio)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

plot_grid(plot1,M,MS,MG,ncol=2,align = "v", axis="lr",labels = c("a","b","c","d"))

# Plots pigments
plot2<-ggplot(newdata.pig.list,
              aes(y=percent,x=pigments))+
  geom_boxplot()+
  ylab("Relative concentration (%)")+
  xlab("")+
  theme_bw()+
  my.theme

CHlA<-ggplot(all.data[all.data$site=="AV",],
             aes(x=Chl.a,
                 y=pol.data.all$CLD[pol.data.all$site=="AV"],
                 label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("Chlorophyll a relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme


FUCO<-ggplot(all.data[all.data$site=="AV",],
             aes(x=Fuco,
                 y=pol.data.all$CLD[pol.data.all$site=="AV"],
                 label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("Fucoxanthin relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

C3<-ggplot(all.data[all.data$site=="AV",],
           aes(x=Chl.c3,
               y=pol.data.all$CLD[pol.data.all$site=="AV"],
               label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("Chlorophyll c3 relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

plot_grid(plot2,CHlA,FUCO,C3,ncol=2,align = "v", axis="lr",labels = c("a","b","c","d"))

# Plots fatty acids
plot3<-ggplot(newdata.fatty.list,
              aes(y=percent,x=fatty))+
  geom_boxplot()+
  ylab("Relative concentration (%)")+
  xlab("")+
  theme_bw()+
  my.theme

EPA<-ggplot(all.data[all.data$site=="AV",],
            aes(x=all.data[all.data$site=="AV",]$"20.5-n3",
                y=pol.data.all$CLD[pol.data.all$site=="AV"],
                label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("EPA relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

ARA<-ggplot(all.data[all.data$site=="AV",],
            aes(x=all.data[all.data$site=="AV",]$"20.4-n6",
                y=pol.data.all$CLD[pol.data.all$site=="AV"],
                label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("ARA relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  my.theme

SAT<-ggplot(all.data[all.data$site=="AV",],
            aes(x=all.data[all.data$site=="AV",]$"12.0"+
                  all.data[all.data$site=="AV",]$"13.0"+
                  all.data[all.data$site=="AV",]$"14.0"+
                  all.data[all.data$site=="AV",]$"15.0"+
                  all.data[all.data$site=="AV",]$"16.0"+
                  all.data[all.data$site=="AV",]$"17.0"+
                  all.data[all.data$site=="AV",]$"18.0"+
                  all.data[all.data$site=="AV",]$"19.0"+
                  all.data[all.data$site=="AV",]$"20.0"+
                  all.data[all.data$site=="AV",]$"22.0"+
                  all.data[all.data$site=="AV",]$"24.0",
                y=pol.data.all$CLD[pol.data.all$site=="AV"],
                label=time))+
  geom_text(aes(color=river))+
  scale_color_manual(values=my.palette(6))+
  geom_smooth(method="lm")+
  xlab("SFA relative concentration (%)")+
  ylab(expression(paste("CLD (µg.",L^-1,")")))+
  ylim(c(-0.5,2))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14))

plot_grid(plot3,EPA,ARA,SAT,ncol=2,align = "v", axis="lr",labels = c("a","b","c","d"))

# BC-MFA biochemical
blo.biochem<-c(dim(mean.PERCENT.PIG[,4:26])[2],
               dim(mean.PERCENT.MONO[,4:12])[2],
               dim(mean.PERCENT.FATTY[,4:43])[2])

bc.mfa(all.data[all.data$site=="AV",2:73],
       bloc=blo.biochem,
       fac=factor(pol.data.all$Context[all.data$site=="AV"],
                  level=c("UnP","Med","Ext")),
       CLD.concentration=sqrt(pol.data.all$CLD[all.data$site=="AV"]/pi),
       spcos=0.25)

# BC-MFA molecular

DATA.bcmfa<-mean.PERCENT.TRFLP[,apply(mean.PERCENT.TRFLP[,4:386],2,sum)!=0]
names(DATA.bcmfa)

bc.mfa(DATA.bcmfa[DATA.bcmfa$"Group.2"=="av",4:386],
             bloc=c(length(4:155),length(156:334),length(335:386)),
             fac=factor(pol.data.molec$Context[pol.data.molec$"Group.2"=="av"],
                        level=c("UnP","Med","Ext")),
             CLD.concentration=sqrt(pol.data.molec$CLD[pol.data.molec$"Group.2"=="av"]/pi),
             spcos=0)
