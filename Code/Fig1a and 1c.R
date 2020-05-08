library(reshape2)
library(ggplot2)
library(chron)
library(gridExtra)
library(plotrix)
#setwd('/Users/benjamin')
#setwd("/Volumes/LaCie/")
#load in csv with trial info
df <- read.csv(file = "~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/RespTrialData.csv")
names(df)[1]<-"FileName"

alldat = 0

#loop through the 15 trials to process respirometry data and compile into one dataframe
for (j in 1:15) {
	
	#load in an individual respirometry trial data file
	trial_N <- j
	filename = toString(df$FileName[trial_N])
	x <- toString(df$Keeper.eggs[trial_N])
	eggs <- as.list(strsplit(x, ",")[[1]])
	x <- toString(df$blanks[trial_N])
	blanks <- as.list(strsplit(x, ",")[[1]])
	#setwd("/Users/benjamin/Google Drive/Projects/ChinookEggExperiment/Respirometry/MicroResp Data For analysis/")
	setwd("~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/MicroResp Data For analysis/")
	dat <- read.csv(file = filename, header = T)
	dat <- dat[, c(1:2, 4, 53:ncol(dat))]
	names(dat)[1:ncol(dat)] <- c("date", "time", "relTime", "A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2", "C3", "C4", "C5", "C6", "D1", 
		"D2", "D3", "D4", "D5", "D6")

	mdat <- melt(data = dat, id.vars = c("date", "time", "relTime"))
	
	# insert trial info (temp, time, egg/blank)
	mdat$temp <- df$Temp[j]
	mdat$devStage <- df$DevStage[j]
	mdat$devDay <- df$DevDay[j]
	mdat$relTime <- chron(times = mdat$relTime)
	mdat$hours <- as.numeric(mdat$relTime * 24) #24 hours per day
	eggdat <- subset(mdat, variable %in% eggs)
	eggdat$type <- "egg"
	blankdat <- subset(mdat, variable %in% blanks)
	blankdat$type <- "blank"
	cdat <- rbind(eggdat, blankdat)	
	cdat <- droplevels(cdat)
	levs = levels(cdat$variable)
	cdat <- subset(cdat, hours < 10)
	cdat <- cdat[order(cdat$variable), ]
	
	# run fit a spline to each DO trace from a trial
	dataOut = 0
	for (i in 1:length(levs)) {
		sdat = subset(cdat, variable == levs[i])
		splineX <- smooth.spline(sdat$hours, sdat$value, nknots = 5)
		predXd <- predict(splineX, deriv = 1)
		predX <- predict(splineX)
		s_O2 <- as.vector(predX$y)
		s_dO2 <- as.vector(predXd$y)
		time = sdat$hours
		#plot(time ,-s_dO2)
		tempout = data.frame(s_O2, s_dO2, time)
		dataOut = rbind(dataOut, tempout)
	}
	dataOut <- dataOut[2:nrow(dataOut), 1:2]
	combdat = cbind(cdat, dataOut)
	
	# calculate the mean DO depletion in the blanks to use for correcting
	aggdat <- aggregate(list(mean_dO2 = combdat$s_dO2, mean_O2 = combdat$s_O2), by = list(combdat$devDay, combdat$type, combdat$hours), FUN = mean, na.rm = TRUE)
	names(aggdat)[1:3] <- c("devDay", "type", "hours")
	b_agg <- subset(aggdat, type == "blank")
	b_agg <- b_agg[, c(1, 3:5)]
	combdat = merge(combdat, b_agg, by = c("devDay", "hours"))
	alldat = rbind(alldat, combdat)
}

alldat <- alldat[-c(1), ]
## 14.5 trial 1 correction
#we used factory calibration setting that missestimated O2 saturation. This code corrects for this by rescaling O2 level and mO2 rate by the ratio of the blanks between 14.1 devPeriod 1 and the remaining 14.5 trials
cordat1 <- subset(alldat, temp == 14.5 & devStage == 1 & type == "blank")
cordat2 <- subset(alldat, temp == 14.5 & devStage != 1 & type == "blank")
corFact_14.5 = mean(cordat2$s_O2)/mean(cordat1$s_O2)
alldat$s_O2[alldat$temp == 14.5 & alldat$devStage == 1] <- alldat$s_O2[alldat$temp == 14.5 & alldat$devStage == 1] * corFact_14.5
alldat$s_dO2[alldat$temp == 14.5 & alldat$devStage == 1] <- alldat$s_dO2[alldat$temp == 14.5 & alldat$devStage == 1] * corFact_14.5
## 

# correct mO2 (in % O2 sat/min) for blank depletion  where correction factor is O2_focal_chamber * mean(mO2_blanks)/mean(O2_blanks)
blank_depetion_correction <- alldat$s_O2  * alldat$mean_dO2 / alldat$mean_O2
alldat$corr_dO2 <- -1 * (alldat$s_dO2 - alldat$s_O2  * alldat$mean_dO2 / alldat$mean_O2)

alldat$DO_sat <- 11.8
alldat$DO_sat[alldat$temp == 12] <- 10.8
alldat$DO_sat[alldat$temp == 14.5] <- 10.2
alldat$DO_sat[alldat$temp == 17] <- 9.7
alldat$DO_conc<-alldat$s_O2 * alldat$DO_sat/100


EggVol = 4/3 * pi * 0.37^3
Vol = 1.7 - EggVol
Vol_L = Vol
alldat$corrFac = (alldat$DO_sat * Vol_L)/100 # ugO2/hr
alldat$corr_dO2 = alldat$corr_dO2 * alldat$corrFac


## calc MO2 in normoxia
highdat <- subset(alldat, DO_conc > 7.5 & DO_conc < 10 & type == "egg")

#replicate means
high_aggdat <- aggregate(list(mean_corr_dO2 = highdat$corr_dO2), by = list(devStage = highdat$devStage, devDay = highdat$devDay, temp = highdat$temp, variable = highdat$variable), 
	FUN = mean, na.rm = TRUE)
	
#treatment means
high_groupdat <- aggregate(list(mean_corr_dO2 = high_aggdat$mean_corr_dO2), by = list(devStage = high_aggdat$devStage, devDay = high_aggdat$devDay, temp = high_aggdat$temp), FUN = function(x) c(mean = mean(x), 
	se = std.error(x)))
high_groupdat <- do.call(data.frame, high_groupdat)
names(high_groupdat)[4] <- "meanMo2"
names(high_groupdat)[5] <- "seMo2"
high_groupdat $upper <- high_groupdat $meanMo2 + high_groupdat $seMo2 * 1.96
high_groupdat $lower <- high_groupdat $meanMo2 - high_groupdat $seMo2 * 1.96

## plot fig1a
ggplot() + geom_jitter(data = high_aggdat, aes(x = devDay, y = mean_corr_dO2, fill = factor(temp)), height=0,width = 0.25, alpha = 0.6,size=1.8,shape=21) + 
	geom_point(data = high_groupdat, aes(x = devDay, y = meanMo2, fill = factor(temp)), size = 3,alpha=0.6,shape=21) + 
	geom_line(data = high_groupdat, aes(x = devDay, y = meanMo2, color = factor(temp)), size = 1)+ylab("MO2 ug/egg/hr")+ylim(0,18)+
		geom_ribbon(data = high_groupdat, aes(x = devDay, ymin = lower, ymax = upper, fill = factor(temp)), alpha = 0.3) + 
		theme_cowplot()+
     scale_color_viridis(option='inferno',discrete=TRUE,end = .8,direction = 1)+
     scale_fill_viridis(option='inferno',discrete=TRUE,end = .8,direction = 1)+
     xlab("Days post fertilization")+
     ylab("Metabolic rate (ugO2/hr/egg)")+    theme(legend.position = c(0.05, 0.85))+ 
    theme(legend.title = element_blank())


### Figure 1C
parms<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/MatlabModel/BootStrappedParms.csv',header=F)
datout=0
devDayvec <-15:37
tempvec<-c(8,12,14.5,17)
mean_mo2<-0
upper_mo2 <-0
lower_mo2 <-0
temp_mo2<-0
devDay_mo2<-0
for (j in 1:length(tempvec)){
temp= tempvec[j]	
for (i in 1:length(devDayvec)){
devDay= devDayvec[i]
mass = exp(-12.630 + 2.692 *log(devDay)+ 2.838* log(12));
mo2  =  parms[,3]*mass^1*exp(parms[,4]* temp);
mean_mo2[i] =  as.numeric(mean(mo2))
upper_mo2[i] = as.numeric(quantile(mo2,.975))
lower_mo2[i] = as.numeric(quantile(mo2,.025))
temp_mo2[i]<-temp
devDay_mo2[i]<-devDay
tempdat<-data.frame(temp=temp_mo2, devDay=devDay_mo2, mean_mo2, upper_mo2, lower_mo2)
}
datout =rbind(datout,tempdat)
}

datout <- datout[-c(1), ]
plotx<-ggplot()+geom_line(data=datout,aes(x= devDay,y= mean_mo2,group= temp,color= factor(temp)),size=1.5)+geom_ribbon(data=datout,aes(x= devDay,ymax= upper_mo2,ymin=lower_mo2,fill= factor(temp)),alpha=.2)+xlab("Day post fetilization")+ylab("MO2 (ug/hr/egg)")+ labs( color = "Temperature",fill="Temperature")+ 
geom_jitter(data = high_aggdat, aes(x = devDay, y = mean_corr_dO2, fill = factor(temp)), height=0,width = 0.25, alpha = 0.6,size=1.8,shape=21)+
theme_cowplot()+
     scale_color_viridis(option='inferno',discrete=TRUE,end = .8,direction = 1)+
     scale_fill_viridis(option='inferno',discrete=TRUE,end = .8,direction = 1)+
     xlab("Days post fertilization")+
     ylab("Metabolic rate (ugO2/hr/egg)")+theme(legend.position = "none") #theme(legend.position = c(0.05, 0.85))
   # theme(legend.title = element_blank())+ geom_jitter(data = high_aggdat, aes(x = devDay, y = mean_corr_dO2, fill = factor(temp)), height=0,width = 0.25, alpha = 0.6,size=1.8,shape=21) +l

comb<-merge(high_groupdat ,datout,by=c("temp","devDay"))
comb$dummy<-"B"

	inset<-ggplot(comb,aes(x= mean_mo2, y= meanMo2))+geom_point(shape=21,fill='black',alpha=.0)+ annotate("segment", x = 0, xend = 15, y = 0, yend = 15,colour = "black")+
stat_smooth(aes(col=dummy))+theme_cowplot()+
ylab("B")+xlab("D")+xlim(0,50)+ylim(0,15)+theme(legend.position="none")+
geom_point(aes(x= mean_mo2, y= meanMo2),shape=21,fill='black',alpha=.5)	+
 annotate("text", x = 7, y = 14.5, label = "1:1")+
	theme(text = element_text(size=10)) + theme(axis.text = element_text(size = 10))


	
plotx + 
  annotation_custom(
    ggplotGrob(inset), 
    xmin = 14, xmax = 27, ymin = 30, ymax = 64
  )

	