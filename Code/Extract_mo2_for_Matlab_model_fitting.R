library(reshape2)
library(ggplot2)
library(chron)
library(signal)
library(gridExtra)
library(plotrix)
#setwd('/Users/benjamin')
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
	cdat <- subset(cdat, hours < 12)
	cdat <- cdat[order(cdat$variable), ]
	# temp control not established for ~30 mins during trial 1
	if (cdat$devDay[1] == 15){
	cdat <- subset(cdat,hours > .5 )
	}
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

alldat<-subset(alldat, s_O2>5)
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

alldat$corrFac = (alldat$DO_sat * Vol_L)/100 
alldat$corr_dO2 = alldat$corr_dO2 * alldat$corrFac # ugO2/hr


eggdat<- subset(alldat,type=='egg')
eggdat <- eggdat[order(eggdat$variable),]
eggdat$ID <- cumsum(!duplicated(eggdat[c(1,6)]))

DO_points = c(10,20,30,40,50,60,70,80,90,100)

splitdat<-split(eggdat, eggdat$ID)

datout<-data.frame()
for (i in 1:length(splitdat)){
	
	for (j in 1:length(DO_points)){
	DOlev = DO_points[j]
	if (min(abs(DOlev-splitdat[[i]]$s_O2)) < 2){ 
	rowid<-which(abs(DOlev-splitdat[[i]]$s_O2)== min(abs(DOlev-splitdat[[i]]$s_O2)))
	temp<-splitdat[[i]][rowid,]
	temp$s_O2 <-DOlev
	datout <-rbind(datout,temp)
	}
	}
}

datout$treatID <- cumsum(!duplicated(datout[c(1)]))

sdatout<-data.frame(devDay=datout$devDay,temp=datout$temp,devStage=datout$devStage, wellID = datout$variable, treatID = datout$treatID, DOsat = datout$s_O2,DO_conc=datout$DO_conc,mo2=datout$corr_dO2,DO_satconc = datout$DO_sat)

#write.csv(sdatout,'mO2dat_for_model_fitting.csv')