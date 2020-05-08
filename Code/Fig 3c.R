
dat<-read.csv(file="~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/DO_distribution_data.csv")

# calc mean O2 measurment for each sensory across time points
aggdat<-aggregate(dat, by = list(Level=dat$Level, Location=dat $Location,in_out=dat$in_out,Tube=dat$Tube,temp=dat$temp), FUN = mean, na.rm = TRUE)

# dummy variable for pairwise analysis so that sensors in different tubes have large distances, and we only evaluate pairwise differences between sensors in the same tube
aggdat$zDummy<- aggdat$temp * 1000 + aggdat$Tube*100 + aggdat$z

aggdat<-data.frame(Level=aggdat$Level, Location =aggdat$Location, in_out = aggdat$in_out, Tube =aggdat$Tube, temp =aggdat$temp, z=aggdat$z, x=aggdat$x, y=aggdat$y,zDummy=aggdat$zDummy, oxygen.mg.L..ppm.=aggdat$oxygen.mg.L..ppm.)

# calc pairwise differences in O2 meansured between sensors
O2dif<-as.matrix(dist(aggdat$oxygen.mg.L..ppm.)) 
dists <- as.matrix(dist(aggdat[,7:9])) 
O2dif <- as.vector(O2dif)
dists <- as.vector(dists)

comb<-data.frame(O2dif, dists)
comb <-subset(comb,dists<10 & dists>0)
comb$Color <-"Color"
ggplot(comb,aes(x=dists,y=O2dif))+geom_jitter(width=.03,alpha=.1)+theme_cowplot()+stat_smooth(aes(color= Color))+scale_color_viridis(option = "plasma", discrete = TRUE, begin = 0.5, direction = 1)+theme(legend.position = "none") +xlab("Distance (cm)")+ylab(expression(paste("Pairwise absolute DO \n    difference (ug/ml)")))+  theme(plot.margin = margin(1, 1, 1, 1, "cm"))


ggplot(aggdat,aes(x=oxygen.mg.L..ppm.,group=factor(temp)))+geom_density(aes(fill=factor(temp)))+	scale_fill_viridis(option = "plasma", discrete = TRUE,begin=.3, end = 0.8, direction = 1,alpha=.5)+
	theme(legend.position = "none") +theme_cowplot()+	theme(legend.position = c(0.7, 0.8)) +theme(legend.title = element_blank()) +
	xlab("Dissolved oxygen ugO2/l")


