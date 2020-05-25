hatchDat =read.csv(file="~/Google Drive/Projects/ChinookEggExperiment/Final/Data/HatchDayData.csv",sep=";")
meanhatchDat <-aggregate(hatchDat, by=list(hatchDat$RepNumber),FUN=mean, na.rm=TRUE)
meanhatchDat$PeriodNames<-"Early"
meanhatchDat$PeriodNames[meanhatchDat $DevPeriod==2]<-"Mid"
meanhatchDat$PeriodNames[meanhatchDat $DevPeriod==3]<-"Late"

meanhatchDat $PeriodNames <- factor(meanhatchDat $PeriodNames, levels = c("Early","Mid", "Late"))
meanhatchDat <-meanhatchDat[order(meanhatchDat $PeriodNames), ]


treatMeanhatchDat <-aggregate(meanhatchDat$HatchDay, by=list(meanhatchDat$Temp,meanhatchDat$O2,meanhatchDat$PeriodNames),
 FUN=function(x) c(mean=mean(x), se=std.error(x)))
treatMeanhatchDat <- do.call(data.frame, treatMeanhatchDat)
names(treatMeanhatchDat)<-c("Temp","O2","PeriodNames","HatchMean","HatchSE")

treatMeanhatchDat$upper= treatMeanhatchDat$HatchMean + 1.96* treatMeanhatchDat$HatchSE
treatMeanhatchDat$lower= treatMeanhatchDat$HatchMean - 1.96* treatMeanhatchDat$HatchSE


cols<-c("#FE9261","#4AB78E","#5D1774")
ggplot()+geom_jitter(data= meanhatchDat ,aes(x=Temp,y= HatchDay,fill=factor(O2)),height=.0,width=.1,alpha=.4,size=1.5,shape=21,color='black')+facet_wrap(~ PeriodNames,ncol=3)+
	#geom_errorbar(data= treatMeanhatchDat ,aes(x=factor(Temp),ymax= upper,ymin=lower,color=factor(O2)),width=.15,size=.6,color='black')+
	geom_line(data= treatMeanhatchDat,aes(x=Temp,y= HatchMean,color=factor(O2)),size=1)+
		geom_ribbon(data= treatMeanhatchDat ,aes(x=Temp,ymax= upper,ymin=lower,fill=factor(O2)),alpha=.4)+
	geom_point(data= treatMeanhatchDat,aes(x=Temp,y= HatchMean,fill=factor(O2)),size=3,shape=21,color='black')+
    theme_cowplot()+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    xlab("Temperature °C")+
    ylab("Days to hatch ")+ 
    theme(legend.title = element_blank())+
    scale_x_continuous( breaks =c(12,14.5,17),limits=c(11.5,17.5))+ 
    theme(legend.position = c(0.9, 0.72))
