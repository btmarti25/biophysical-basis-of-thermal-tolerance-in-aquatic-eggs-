library(reshape2)
library(ggplot2)
library(lme4)
library(gridExtra)
library(MASS)
library(cowplot)
library(dplyr)
library(Hmisc)
library(viridis)


#write.csv(sumdat,'sumdat.csv')
sumdat<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/Final/Data/SurvivalToHatch.csv',header=T)
TreatmentID <-sumdat$TreatmentID
splitdat<-split(sumdat, TreatmentID)
ID<-0
Treat_mean<-0
Treat_se<-0
for (i in 1:length(splitdat)){

model<-glmer(cbind(Survived, morts)~ +(1| RepNumber),family = binomial(link=logit),control = glmerControl(optimizer="Nelder_Mead",optCtrl = list(maxfun = 5000000)),data= splitdat[[i]]) 

ID[i]<-	i
splitdat[[i]]$Treat_mean<-as.numeric(fixef(model))
mm = unique(model.matrix(model))
Tpred<-mm%*%fixef(model)
mm*fixef(model)*t(mm)
splitdat[[i]]$Treat_se <- sqrt(diag(mm %*% tcrossprod(vcov(model),mm)))
}
sumdat<-unsplit(splitdat, TreatmentID)

sumdat$upper<-sumdat$Treat_mean+ sumdat$Treat_se * 1.96
sumdat$lower<-sumdat$Treat_mean- sumdat$Treat_se * 1.96

sumdat$Treat_mean <-exp(sumdat$Treat_mean)/(exp(sumdat$Treat_mean)+1) 
sumdat$lower<-exp(sumdat$lower)/(exp(sumdat$lower)+1)
sumdat$upper <-exp(sumdat$upper)/(exp(sumdat$upper)+1)

sumdat$PeriodNames<-"Early"
sumdat$PeriodNames[sumdat$DevPeriod==2]<-"Mid"
sumdat$PeriodNames[sumdat$DevPeriod==3]<-"Late"

sumdat$PeriodNames <- factor(sumdat$PeriodNames, levels = c("Early","Mid", "Late"))

sumdat <-sumdat[order(sumdat $PeriodNames), ]

cols<-c("#FE9261","#4AB78E","#5D1774")

	p2<-ggplot(sumdat,aes(x=factor(PeriodNames),y=Survived/25))+geom_jitter(aes(fill=factor(O2)),height=.0,width=.1,alpha=.7,size=1.5,shape=21,color='black')+facet_wrap(~Temp,ncol=1,scales='free')+
	geom_errorbar(aes(x=factor(PeriodNames),ymax= upper,ymin=lower,color=factor(O2)),alpha=.15,width=.15,size=.6,color='black')+
	geom_point(aes(x=factor(PeriodNames),y= Treat_mean,fill=factor(O2)),size=3,alpha=1,shape=21,color='black')+	
	theme_cowplot()+
	xlab("Exposure stage")+ylab("Proportion surviving to hatch")+ ylim(0,1)+
	theme(
	  strip.background = element_blank(),
	  strip.text.x = element_blank()
	)+ theme(legend.position = "none")+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)

	

dat<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/Final/Data/survivaldata_QAQCed.csv')
mdat<-subset(dat,Category=="Mort")
TreatmentN=mdat$RepNumber
splitDat<-split(mdat,TreatmentN)

length(splitDat)
for (i in 1:length(splitDat)){
n = 25
N=25
	for (t in 1:nrow(splitDat[[i]])){
		n= n - splitDat[[i]]$value[t]
		N[t]=n		
	}
	splitDat[[i]]$S=N
}
mdat<-unsplit(splitDat,TreatmentN)
library(grid)

p1<-ggplot(mdat,aes(x= devDay,y=S/25,group= RepNumber))+geom_line(aes(color=factor(O2)),alpha=.8,size=.8)+facet_wrap(~ Temp + DevPeriod, ncol = 3, scales = "free_x")+
 theme_cowplot()+
theme(strip.background = element_blank(),strip.text = element_blank())+
  ylab("Proportion surviving")+xlab("Development day")+ theme(legend.position = c(0.15, 0.45))+xlim(13,45)+
    scale_color_manual(values= cols)+ theme(legend.title = element_blank())


grid.arrange(p1, p2,ncol=2,widths=c(2,.75))


###
parms<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/Final/Data/BootStrappedParms.csv',header=F)

k=mean(parms[,1])
G=mean(parms[,2])
b0=mean(parms[,3])
b1=mean(parms[,4])
mass = exp(-12.630 + 2.692 *log(sumdat$devDay)+ 2.838* log(12))
D  =  0.2283*mass^1*exp(0.0792* sumdat$Temp)
Co= sumdat$O2
Ci=(-k*G - D + G*Co+sqrt(k^2*G^2+(D-G*Co)^2+2*k*G*(D+G*Co)))/(2*G)
f=Ci/(k+Ci)


sumdat$O2lim=f

#fit survival model as function of O2lim
calc_likelihood <- function(x) {
s=x[1]
b=x[2]
b2=x[3]
#pred=  s* exp(-b* Odif  b2 * Odif^2)
pred=s*exp(-exp(b+ sumdat$O2lim*b2))
lprob=dbinom(sumdat$Survived,size=25,prob=pred,log=T)
error= -sum(lprob)
}

#initial parameters
x0<-c(.9,4,-37)
mod<-optim(x0, calc_likelihood)

## use fitted parameters to make prediction
s=mod$par[1]
b=mod$par[2]
b2=mod$par[3]
O2lim<-seq(0,1,.001)
pred=  s*exp(-exp(b+O2lim*b2))

treatdat<-sumdat[!duplicated(sumdat$TreatmentID),]
treatdat$treat_pred=s*exp(-exp(b+ treatdat $O2lim*b2))
ggplot(treatdat,aes(x= treat_pred,y= Treat_mean))+geom_point(shape=21,size=2,fill='coral2',alpha=.7)+theme_cowplot()+geom_line(aes(x=Treat_mean,y= Treat_mean))+ scale_x_continuous(breaks=c(0, 1))+ scale_y_continuous(breaks=c(0, 1))+xlab("")+ylab("")

preddf<-data.frame(O2lim , pred)


#plot 1c
ggplot()+geom_line(data=preddf,aes(x=1-O2lim,y= pred),size=2,color='#FE9261')+
	geom_jitter(data=sumdat,aes(x=1-O2lim,y= Survived/25),width=.004,alpha=.5)+
	theme_cowplot()+
	ylab("Proportion surviving")+
	xlab("Oxygen limitation 1-B/D")

