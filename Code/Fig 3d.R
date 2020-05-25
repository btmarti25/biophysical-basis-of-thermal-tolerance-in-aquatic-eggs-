library(ggplot2)
library(viridis)
Olim<-read.csv(file="~/Google Drive/Projects/ChinookEggExperiment/Final/Data/CFD_Olim.csv")

aggOlim <- aggregate(Olim, by = list(Olim $temp,Olim$velocity), FUN = mean, na.rm = TRUE)


ggplot()+
	geom_line(data=aggOlim,aes(x= velocity,y= Olim,group=factor(temp),color=factor(temp)),size=2,alpha=.5)+
	geom_jitter(data=Olim,aes(x= velocity,y= Olim,group=factor(temp),fill=factor(temp)), height=.00,width=.0015,alpha=.15, shape = 21, size = 2)+
	scale_color_viridis(option = "plasma", discrete = TRUE, end = 0.7, direction = 1)+
	scale_fill_viridis(option = "plasma", discrete = TRUE, end = 0.7, direction = 1) +theme_cowplot()+xlim(0,.1)+
	geom_hline(yintercept=.80,size=1.5,linetype='solid')+
	ylab("Oxygen limitation 1-B/D")+xlab("Flow velocity (cm/s)")+
	theme(legend.position = c(0.7, 0.9)) + 
	theme(legend.title = element_blank()) 

