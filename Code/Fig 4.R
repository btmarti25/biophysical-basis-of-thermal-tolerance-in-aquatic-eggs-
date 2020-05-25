library(reshape2)
library(ggplot2)
library(lme4)
library(gridExtra)
library(MASS)
library(cowplot)
library(dplyr)
library(tidyverse)
library(cowplot)
GravelSurvivalData <- read.csv(file = "~/Google Drive/Projects/ChinookEggExperiment/Final/Data/GravelSurvivalData.csv", header = T)

## fit mixed model for survial as function of flow and temperature
model <- glmer(cbind(survivors, morts) ~ factor(temp) * factor(flow) + (1 | Tube), family = binomial(link = logit), 
	control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1e+06)), data = GravelSurvivalData)
summary(model)

## get Predictions and CI from mixed model
predictionData <- expand.grid(temp = c(11, 15), flow = c(0.01, 0.04, 0.12), distance = 0)
predictionData$pred <- predict(model, predictionData, re.form = NA)
mm = unique(model.matrix(model))
Tpred <- mm %*% fixef(model)
mm * fixef(model) * t(mm)
CI <- sqrt(diag(mm %*% tcrossprod(vcov(model), mm))) * 1.96
predictionData$CI <- CI
predictionData$sPred <- exp(predictionData$pred)/(1 + exp(predictionData$pred))
predictionData$upr <- exp(predictionData$pred + predictionData$CI)/(1 + exp(predictionData$pred + predictionData$CI))
predictionData$lwr <- exp(predictionData$pred - predictionData$CI)/(1 + exp(predictionData$pred - predictionData$CI))

## Get mean of observed survival of top and bottom egg pockets
aggregateSurvival <- aggregate(GravelSurvivalData, by = list(GravelSurvivalData$Tube), FUN = mean, na.rm = TRUE)

##plot observed and predicted survival
pSurvival <- ggplot() + 
	geom_line(data = predictionData, aes(x = flow, y = sPred, group = factor(temp), 
	color = factor(temp)), size = 1) + 
	geom_ribbon(data = predictionData, aes(x = flow, ymin = lwr, 
	ymax = upr, fill = factor(temp)), alpha = 0.3, width = 0.004, size = 1) + 
	geom_jitter(data = aggregateSurvival, aes(x = flow, y = propS, fill = factor(temp)), width = 0.003, alpha = 	0.6, shape = 21, size = 2) +
	theme_cowplot()+
	theme(legend.position = c(0.7, 0.2)) + 
	scale_color_viridis(option = "plasma", discrete = TRUE, end = 0.65, direction = 1)+
	scale_fill_viridis(option = "plasma", discrete = TRUE, end = 0.65, direction = 1) + 
	xlab("Flow (cm/s)") + ylab("Proportion surviving") + ylim(0, 1) + xlim(0, 0.126) +
	theme(legend.title = element_blank())+ theme(legend.position = c(0.7, 0.2))  
	


### sublethal plot
dat <- read.csv(file = "~/Google Drive/Projects/ChinookEggExperiment/Final/Data/Alevin_lengths.csv")
aggLengthData <- aggregate(dat, by = list(dat$Tube), FUN = mean, na.rm = TRUE)
tube <- read.csv(file = "~/Google Drive/Projects/ChinookEggExperiment/Final/Data/gravel_tudeIds.csv")
comb <- merge(dat, tube, by = "Tube")
combAgg <- merge(aggLengthData, tube, by = "Tube")


meanMod <- lm(length ~ flow * temp, data = combAgg)
summary(meanMod)
DatComb <- cbind(combAgg, predict(meanMod, combAgg, interval = "predict"))

plength <- ggplot(DatComb, aes(x = flow, y = fit, group = temp)) +
	geom_line(aes(x = flow, y = fit,	color = factor(temp)), size = 1) + 
	geom_jitter(aes(x = flow, y = length, fill = factor(temp)), width = 0.003, alpha = 0.6, shape = 21, 
		size = 2)+
	stat_smooth(aes(x = flow, y = length, color = factor(temp), fill = factor(temp)), method = "lm") +
	theme_cowplot()+
	scale_color_viridis(option = "plasma", discrete = TRUE, end = 0.65, direction = 1) + 
	scale_fill_viridis(option = "plasma", discrete = TRUE, end = 0.65, direction = 1) + 
	ylim(18, 23) + xlab("Flow (cm/s)") + ylab("Length (mm)") + 
	theme(legend.position = c(0.7, 0.2)) + 
	theme(legend.title = element_blank()) + xlim(0, 0.126)



#####
GravelSurvivalData$logitS=log(GravelSurvivalData$propS/(1-GravelSurvivalData$propS))
backPocket<-subset(GravelSurvivalData,Location=="Top")
names(backPocket)[10]<-"logitSback"
frontPocket<-subset(GravelSurvivalData,Location=="Bottom")
names(frontPocket)[10]<-"logitSfront"

comb<-merge(frontPocket ,backPocket,by=c("temp","flow","rep"))

t.test(comb$logitSfront, comb$logitSback, paired = TRUE, alternative = "two.sided")

pairedDat<-data.frame(Front=comb$logitSfront, Back=comb$logitSback)
pairedDat$rep<-1:29
pairedDat =melt(pairedDat,id="rep")

pd<-pairedDat %>%
    tidyr::spread(variable, value) %>%
    dplyr::mutate(is_increasing = Front < Back) %>%
    tidyr::gather("variable", "value", 2:3)
pd <- pd %>% mutate(variable = as.factor(variable))
pd$variable <- factor(pd$variable,levels(pd$variable)[c(2,1)])

pFrontBack<-    ggplot(pd,aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = variable), alpha = 0.3, col = "grey") +
    geom_jitter(aes(fill=variable),height=0,width=.03,alpha = 0.6, shape = 21, size = 2) +   
    geom_line(aes(group = rep, col = is_increasing),alpha=.5) +
    theme_cowplot()+
    scale_color_viridis(option = "plasma", discrete = TRUE, begin = 0.35,end=.85,direction = 1) + 
	scale_fill_viridis(option = "plasma", discrete = TRUE, begin = 0.35, end=.85,direction = 1) +
		theme(legend.position = "none") + xlab("Egg cluster postion")+ylab("Logit survival")


grid.arrange(pSurvival, pFrontBack, plength, ncol = 3)
