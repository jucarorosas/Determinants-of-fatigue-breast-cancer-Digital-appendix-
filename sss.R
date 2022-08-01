################# Prediction from linear latent class model
library(lcmm)

data(data_hlme)
## fitted model
m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
        subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
                                             -14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))
## newdata for predictions plot
newdata<-data.frame(Time=seq(0,5,length=100),
                    X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
plot(predictL(m,newdata,var.time = "Time"),legend.loc="right",bty="l")
## data from the first subject for predictions plo
firstdata<-data_hlme[1:3,]
plot.predict(m,firstdata,"Time","right",bty="l")

plot.