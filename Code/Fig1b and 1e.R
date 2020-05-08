library(ggplot2)
library(cowplot)
dat<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/mO2dat_for_model_fitting.csv',header=T)
mod<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/O2_model_preds_and_CI.csv',sep=";")

## Fig 1B
ggplot()+geom_jitter(data=dat,aes(x= DO_conc,y= mo2,fill=factor(devStage)),width=.2,alpha=.6,shape=21,size=2)+facet_wrap(~temp,ncol=4)+
	geom_line(data=mod,aes(x=DO,y=pred,group=factor(devStage)),size=1,alpha=.8,color='black')+
	geom_ribbon(data=mod,aes(x=DO,ymax=upper,ymin=lower,fill=factor(devStage)),alpha=.5)+
	#	geom_point(data=mod,aes(x=DO,y= meanObs,color=factor(devStage)),size=3)+theme_cowplot()
   scale_color_viridis(option='viridis',discrete=TRUE,end = .9,direction = 1)+
     scale_fill_viridis(option='viridis',discrete=TRUE,end = .9,direction = 1)+
     theme_cowplot()+
     xlab("Dissolved oxygen concentration μg/ml")+
     ylab("Metabolic rate (μgO2/hr/egg)")+theme(
  strip.background = element_blank(),
  strip.text = element_blank())

sumdat<-aggregate(dat$mo2, by = list(temp=dat$temp,devDay=dat$devDay, DOsat =dat$DOsat), FUN = "mean")
names(sumdat)[4]<-"B"
mass = exp(-12.630 + 2.692 *log(sumdat$devDay)+ 2.838* log(12))
sumdat$D  =  0.2283*mass^1*exp(0.0792* sumdat$temp)
sumdat$DO_satconc<-9.7
sumdat$DO_satconc[sumdat$temp==14.5]<-10.2
sumdat$DO_satconc[sumdat$temp==12]<-10.8
sumdat$DO_satconc[sumdat$temp==8]<-11.8
sumdat$DO_conc<- sumdat$DOsat*sumdat$DO_satconc/100
D=sumdat$D  
k=1.47
G=1.76
Co=sumdat$DO_conc
sumdat$Ci=(-k*G - D + G*Co+sqrt(k^2*G^2+(D-G*Co)^2+2*k*G*(D+G*Co)))/(2*G)
K=1.47 #estimated K
sumdat$pred<-sumdat$Ci/(K+ sumdat$Ci)
sumdat$minpred<-sumdat$Ci/(1.11+ sumdat$Ci) #lower CI of K
sumdat$maxpred<-sumdat$Ci/(1.95+ sumdat$Ci) #upper CI of K

### Fig 1E
ggplot(sumdat,aes(x=Ci,y=pred))+geom_line(aes(x=Ci,y=pred),size=1,alpha=0)+ geom_ribbon(aes(x=Ci,ymin=minpred,ymax=maxpred),alpha=.3)+ geom_point(aes(x= Ci,y= B/D, group=interaction(temp, devDay),fill=D),alpha=.6,shape=21,size=2)+
theme_cowplot()+ylim(0,1)+xlim(0,10)+ scale_fill_viridis(option='inferno',begin = .2,end=.9,direction = 1)+theme(legend.position=c(.7,.3)) +ylab("")+theme(plot.margin=margin(.6,.8,.6,.2, "cm"))+geom_line(aes(x=Ci,y=pred),size=1,alpha=.6)+
annotate("segment", x = 1.47, xend = 1.47, y = 0, yend = .5,colour = "black",linetype="dashed")+
annotate("segment", x = -0, xend = 1.47, y = .5, yend = .5,colour = "black",linetype="dashed")+theme(legend.title=element_blank())

