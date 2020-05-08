library(ggplot2)
library(cowplot)
library(gridExtra)
parms<-read.csv(file='~/Google Drive/Projects/ChinookEggExperiment/PNAS/Data/BootStrappedParms.csv',header=F)

Co <-11
k<-mean(parms[,1])
G<-mean(parms[,2])
D<-seq(1,50,1)

dBdD=(sqrt(2 *D *G *(k - Co) + G^2 *(Co + k)^2 + D^2) + Co *G - D - G *k)/(2 *sqrt(2* D *G *(k - Co) + G^2 *(Co + k)^2 + D^2))
dBdCo=(2 * D *G^2 *k)/(sqrt(2 *D *G *(k - Co) + G^2 * (Co + k)^2 + D^2)* (sqrt(2 * D * G *(k - Co) + G^2 *(Co + k)^2 + D^2) + Co * G - D + G * k))
data<-data.frame(D, dBdD, dBdCo)


ggplot(data,aes(x=D,y= dBdCo/G))+geom_line(size=2,col="#EDA13D",alpha=.85)+ theme_cowplot()+xlab("Metabolic demand")+ylab("Supply/Demand sensitivity")+
geom_line(aes(x=D,y=dBdD),col='#681E64',size=2,alpha=.85)+ylim(0,1)        
        
