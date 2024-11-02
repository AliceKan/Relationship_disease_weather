library(skimr)
library(foreign)
library(ggplot2)
library(gnm)
library(splines)
library(tsModel)
library(Epi)
library(dlnm)
library(dplyr)
library(corrplot)

rm(list=ls())
load(file="NMMAPS3cities.RData")
ls()
skim(chic)
class(chic$date)
chic$date <- as.Date(chic$date, format=("%Y %m %d"))

#time series plots of resp death and two pollutants
layout(1:3)
plot(chic$date, chic$resp, ylab="Resp Deaths", xlab="Date", cex=0.6, 
     main="Time series of daily respiratory deaths")
plot(chic$date, chic$o3, ylab="Ozone levels", xlab="Date", cex=0.6, col=2,
     main="Time series of daily ozone level")
plot(chic$date, chic$pm10, ylab="pm10 levels", xlab="Date", cex=0.2, col=2,
     main="Time series of daily pm10 level")
layout(1)

# for pm10, ggplot to clearly see the wave
ggplot(chic, aes(x=date, y=pm10)) + 
  geom_point() +
  geom_smooth(method="loess", span=.1, col="red") +
  labs(title="Time series of daily deaths with a loess line")

# temperature
layout(1:4)
plot(chic$date, chic$tmpd, ylab="tmpd", xlab="Date", cex=0.2, 
     main="Time series of daily mean temperature")
plot(chic$date, chic$tmin, ylab="tmin", xlab="Date", cex=0.2,
     main="Time series of daily min temperature")
plot(chic$date, chic$tmax, ylab="tmax", xlab="Date", cex=0.2,
     main="Time series of daily max temperature")
plot(chic$date, chic$dptp, ylab="dptp", xlab="Date", cex=0.2,
     main="Time series of daily dew point temperature")
layout(1)

#correlation
sub <- select(chic, tmpd:resp)
corr <- cor( sub[complete.cases(sub), ] )
par(mfrow=c(1,2))
corrplot(corr, method="number", tl.col="black", tl.srt=45)
corrplot(corr, type="upper", tl.col="black", tl.srt=45)

options(na.action="na.exclude")
# seasonal and long-term trend
#spline
spltime <- bs(chic$date, df=7*14)
mod.a <- glm(chic$resp ~ spltime, data=chic, family=quasipoisson)
pred.a <- predict(mod.a, type="response")
summary(mod.a)
plot(chic$date, chic$resp, ylab="Resp Deaths", xlab="date", cex=0.6,col=grey(0.6),
     main="Modelling seasonality and long-term trend (df=7*14)")
lines(chic$date, pred.a, col="red")


# other confounders------------------------------------------------------------
#scale
chic$ozone10 <- chic$o3/10
chic$pm10ten <- chic$pm10/10
#naive
mod.b <- glm(chic$resp ~ ozone10, chic, family=quasipoisson)
summary(mod.b)
ci.exp(mod.b, subset="ozone10")

#adjust for season and long term trend, and dow
mod.c <- update(mod.b, .~. + spltime +dow)
summary(mod.c)
ci.exp(mod.c, subset="ozone10")

#adjust for temperature 
#mean
spltmpd <- ns(chic$tmpd, df=3) # df=3, qaic=26105, <df=4(26107) or df=2(26110) 
mod.d1 <- update(mod.c, .~. + spltmpd)
summary(mod.d1)
ci.exp(mod.d1, subset="ozone10")
#min
spltmin <- ns(chic$tmin, df=3) 
mod.d2 <- update(mod.c, .~. + spltmin)
summary(mod.d2)
ci.exp(mod.d2, subset="ozone10")
#max
spltmax <- ns(chic$tmax, df=3)  
mod.d3 <- update(mod.c, .~. + spltmax)
summary(mod.d3)
ci.exp(mod.d3, subset="ozone10")
#dew point
spldptp <- ns(chic$dptp, df=3)   
mod.d4 <- update(mod.c, .~. + spldptp)
summary(mod.d4)
ci.exp(mod.d4, subset="ozone10")
#all
mod.d5 <- update(mod.c, .~. + spldptp+ spltmin+ spltmax+ spltmpd)
summary(mod.d5)
ci.exp(mod.d5, subset="ozone10")

QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois( model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}
#lmodels <- list(mod.d1)
lmodels <- list(mod.d1,mod.d2,mod.d3,mod.d4,mod.d5)
sapply(lmodels, QAIC)#d5 has least QAIC

#plus pm10
mod.e <- update(mod.d5, .~. +pm10ten)
summary(mod.e)
ci.exp(mod.e,subset="ozone10")

lmodels <- list(mod.e,mod.d5)
sapply(lmodels, QAIC)# e got lower QAIC
#anova(mod.d5, mod.e)

#AIC(mod.a,mod.b,mod.c)
#anova(mod.a, mod.b, test="Chisq")
#summary(mod.a)$dispersion


# delay (lags)
# single lag------------------------------------------------------------
o3lag <- Lag(chic$ozone10, 0:7)
colnames(o3lag) <- c(sprintf("o3lag%d", 0:7))
chic <- cbind(chic, o3lag)
pm10lag <- Lag(chic$pm10ten, 0:7)
colnames(pm10lag) <- c(sprintf("pm10lag%d", 0:7))
chic <- cbind(chic, pm10lag)

N <- ncol(o3lag)
M <- matrix(numeric(N*3), ncol=3, byrow=TRUE, 
            dimnames=list(colnames(o3lag), c("rr", "cil", "ciu")))

for(i in 1:N){
  cat(i," ")
  
  slag <- chic[, colnames(o3lag)[i]]
  pmlag <- chic[, colnames(pm10lag)[i]]
  slmod.loop <- glm(resp ~ slag +spldptp+spltmin+spltmax+spltmpd +spltime +dow+pmlag, 
                    chic, family=quasipoisson)
  
  print(QAIC(slmod.loop))
  # save estimates
  M[i,] <- ci.exp(slmod.loop, subset="slag")
}
M <- data.frame(M)

# Plot
# windows(10,6)
par(mfrow=c(1,2))
rng <- c(0.960,1.050)
lagday <- 0:7
plot(M$rr ~ lagday, ylim=rng, type="n", 
     xlab="Lag day", ylab="RR for Ozone", main="Single lag model")
abline(h=1, lty=2, col="grey")
segments(lagday, M$cil, lagday, M$ciu)
points(lagday, M$rr, pch=19, col=2)

# unconstrained DLM-----------------------------------------------------------

mod.g1 <- glm(chic$resp ~ Lag(chic$ozone10, 0:7)+Lag(chic$pm10ten, 0:7)
             +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)
summary(mod.g1)
ci.exp(mod.g1, subset="ozone10")
mod.g2 <- glm(chic$resp ~ Lag(chic$ozone10, 0:7)+bs(chic$pm10ten)
               +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)
summary(mod.g2)
ci.exp(mod.g2, subset="ozone10")
QAIC(mod.g1)
QAIC(mod.g2)
M2.dlm <- data.frame( ci.exp(mod.g1, subset="ozone10"))
colnames(M2.dlm) <- c("rr", "cil", "ciu")
M2.dlm

summary(mod.g1)$coef[2:9,]
exp(sum(summary(mod.g1)$coef[2:9,1]))

summary(mod.g2)$coef[2:9,]
exp(sum(summary(mod.g2)$coef[2:9,1]))

# Plot
rng <- c(0.960,1.050)
lagday <- 0:7
plot(M2.dlm$rr ~ lagday, ylim=rng, type="n", 
     xlab="Lag day", ylab="RR for Ozone", main="Distributed lag model")
abline(h=1, lty=2, col="grey")
segments(lagday, M2.dlm$cil, lagday, M2.dlm$ciu)
points(lagday, M2.dlm$rr, pch=19, col=2)


# constrained DLM-------------------------------------------------------------
cbo3spl <- crossbasis(chic$o3, lag=7, argvar=list(fun="lin"), # constrained DLM
                      arglag=list(fun="ns",knots=c(1,3)))
cbo3int <- crossbasis(chic$o3, lag=7, argvar=list(fun="lin"), # unconstrained DLM
                      arglag=list(fun="integer"))
cbpm10spl <- crossbasis(chic$pm10, lag=7, argvar=list(fun="lin"), # constrained DLM
                      arglag=list(fun="ns",knots=c(1,3)))
cbpm10int <- crossbasis(chic$pm10, lag=7, argvar=list(fun="lin"), # unconstrained DLM
                      arglag=list(fun="integer"))

mod.spl1 <- glm(chic$resp ~ cbo3spl +cbpm10spl +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)
mod.int1 <- glm(chic$resp ~ cbo3int +cbpm10int +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)
# ci.exp(mod.spl1, subset="cbo3")
# ci.exp(mod.int1, subset="cbo3")
QAIC(mod.spl1)
QAIC(mod.int1)
#mod.spl2 <- glm(chic$resp ~ cbo3spl +ns(chic$pm10) +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)
mod.int2 <- glm(chic$resp ~ cbo3int +ns(chic$pm10) +spldptp+spltmin+spltmax+spltmpd +spltime +dow, chic, family=quasipoisson)

# Prediction
pred.o3spl <- crosspred(cbo3spl, mod.spl1,cen=0, by=1, bylag=0.2)
pred.o3int <- crosspred(cbo3int, mod.int1, cen=0, by=1)
#pred.o3int <- crosspred(cbo3int, mod.int2, cen=0, by=1)

summary(pred.o3spl)
summary(pred.o3int)

plot(pred.o3spl, "slices", var=10, type="l", col=4, lwd=2, 
     ylab="RR for 10 unit increase", ylim=c(0.950,1.070), main="Distributed Lag Models")
points(pred.o3int, "slices", var=10, ci="bars", pch=19, col=2)     
legend("topright", c("Unconstrained DLM","Spline DLM"),col=c(2,4),lwd=2,ncol=1,cex=0.8)

# Plot of overall cumulative association
plot(pred.o3spl, "overall", col=4, ci="area", lwd=2,
     ylab="RR (95% CI)", xlab="Ozone", main="Overall cumulative association")
lines(pred.o3int, "overall", ci="lines", col=2, lty=2)
legend("topleft", c("Unconstrained DLM","Spline DLM"),col=c(2,4),lwd=2,ncol=1,inset=0.04,cex=0.8)

with(pred.o3spl, cbind(allRRfit, allRRlow, allRRhigh))["10",]
with(pred.o3int, cbind(allRRfit, allRRlow, allRRhigh))["10",] 

#TSCCO---------------------------
for (num in 0:7)
{
  tmpds <- paste("tmpd_l", num, sep="")       
  assign(tmpds, dplyr::lag(chic$tmpd, num)) # assign the lagged values in turn
  tmins <- paste("tmin_l", num, sep="")
  assign(tmins, dplyr::lag(chic$tmin, num))
  tmaxs <- paste("tmax_l", num, sep="")
  assign(tmaxs, dplyr::lag(chic$tmax, num))
  dptps <- paste("dptp_l", num, sep="")
  assign(dptps, dplyr::lag(chic$dptp, num))
}

chic$l07tmpd <- (tmpd_l0+tmpd_l1+tmpd_l2+tmpd_l3+tmpd_l4+tmpd_l5+tmpd_l6+tmpd_l7)/8  
chic$l07tmin <- (tmin_l0+tmin_l1+tmin_l2+tmin_l3+tmin_l4+tmin_l5+tmin_l6+tmin_l7)/8 
chic$l07tmax <- (tmax_l0+tmax_l1+tmax_l2+tmax_l3+tmax_l4+tmax_l5+tmax_l6+tmpd_l7)/8
chic$l07dptp <- (dptp_l0+dptp_l1+dptp_l2+dptp_l3+dptp_l4+dptp_l5+dptp_l6+dptp_l7)/8

chic$year  = as.numeric(format(chic$date,"%Y"))
chic$month = as.numeric(format(chic$date,"%m"))
chic$month <- as.factor(chic$month)
chic$year  <- as.factor(chic$year)
chic$dow   <- as.factor(chic$dow)
chic$stratum <- with(chic, as.factor(year:month:dow))

m1 <- cbind(chic$pm10lag0,chic$pm10lag1,chic$pm10lag2,chic$pm10lag3,chic$pm10lag4,chic$pm10lag5,chic$pm10lag6)
chic$lagpm10 <- apply(m1,1,mean,na.rm=TRUE)
m2 <- cbind(chic$o3lag0,chic$o3lag1,chic$o3lag2,chic$o3lag3,chic$o3lag4,chic$o3lag5,chic$o3lag6)
chic$lagozone <- apply(m2,1,mean,na.rm=TRUE)

mod.h <- gnm(resp ~ lagozone+lagpm10+ns(chic$l07dptp, df=3)+ns(chic$l07tmin, df=3)+ns(chic$l07tmax, df=3)+ns(chic$l07tmpd, df=3), 
                data=chic, family=quasipoisson, eliminate=stratum)
ci.exp(mod.h,subset="lagozone")
summary(mod.h)
QAIC(mod.h)

# res plot----------------------------------------------------------------
layout(1)
res <- residuals(mod.g1, type="response")
plot(chic$date, res ,ylim=c(-20,60),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals for distributed lag model controlling for PM10 by splines",
     ylab="Residuals",xlab="Date")
abline(h=0,lty=2,lwd=1)
# res for no adj
mod.noadj <- glm(resp ~ 1, data=chic, family=quasipoisson)
res.noadj <- residuals(mod.noadj, type="response") # residuals from the null model
plot(chic$date, res.noadj ,ylim=c(-20,60),pch=19,cex=0.4,col=grey(0.6),
     main="Residuals over time (w/o adjustment)",
     ylab="Residuals (observed-fitted)",xlab="Date")
abline(h=0,lty=2,lwd=1)

# interaction model
chic$season = ifelse(chic$month %in% 4:9, 2, 1)
chic$season = factor(chic$season, levels=c(1,2), labels=c("cold","warm"))
mod.interact <- glm(chic$resp ~ ozone10+ozone10:season+ spltime +dow+ spldptp+ spltmin+ spltmax+ spltmpd+pm10, chic, family=quasipoisson)
ci.exp(mod.interact,subset="ozone10")