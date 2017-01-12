require(survival)
require(MortalitySmooth)
require(parfm)
require(rms)
require(msm)
require(DEoptim)
require(numDeriv)
require(eha)
require(car)
require(fitdistrplus)
require(logspline)
require(EWGoF)

source('LIB_LT.R')
source('my.parfm.r')
source('KS_WEIBULL_BOOT.R')
source('my.plot.cox.zph.r')

#Reading and preparing data
dscdane1=try(read.csv('dane1.csv'),silent=T)
if(class(dscdane1)=='try-error') dscdane1=read.csv('C:/Users/Maciek/Documents/R-PRJ/DariuszMalek/dane1.csv')

names(dscdane1)[1]='mother'
dscdane1=dscdane1[,-c(6,7,10,11)] # remove not used columns
dscdane1$sex=as.factor(c('Males','Females')[dscdane1$sex+1]) #originaly Males=0, Females=1
dscdane1$tr=as.factor(c('Virgin','Reproducing')[dscdane1$tr+1]) #originaly Virgin=0, Reproducing=1
dscdane1$trsex=paste(dscdane1$tr,tolower(dscdane1$sex),sep=' ')
dscdane1$bean1=as.numeric(dscdane1$bean1)
dscdane1$mass1=as.numeric(dscdane1$mass1)
U=unique(dscdane1$trsex)

#resolution of figures
res=300 #dpi

#Data description
head(dscdane1)
#mother: id of the mother (random effect/grouping variable/shared frailty)
#bean1: mass of the seed of host plant
#mass1: adult body mass
#sex: gender males and females
#tr: reproduction treatment, virgin and reproducing
#timesurvived: age at death
#status: 1-exact obs. 0-random censoring

options(contrasts = c("contr.treatment","contr.poly"))
#please change to: options(contrasts = c("contr.sum","contr.poly"))
#if you want to operate on orthogonal contrasts

########################################
#Testing proportional hazard assumption
########################################

###########################
#1) log-log plots for categorical variables
###########################

csf=npsurv(Surv(timesurvived,status) ~ sex + tr, data = dscdane1)
plot(csf,col=1:4)
legend('bottomleft',legend=U,col=1:4,lty=1,lwd=2,bty='n',y.intersp = 1.1)


tiff(filename='./results/log-log_plot.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(mar=c(4,4,1,1))
survplot(csf,col=1:4,loglog=T,logt=T,xlim=c(1.5,3.5),
         xlab='log Survival time',ylab='log(-log Survivorship)',conf='bands',
         label.curves=F,lty=1,lwd=2)
legend('bottomleft',legend=U,col=1:4,lty=1,lwd=2,bty='n',y.intersp = 1.1)
box()
dev.off()

#no clear evidence for non-proportional hazard

###########################
#2) Hazard fits and plots for categorical variables
###########################

LT=list()
for (u in U){
  ind=dscdane1$trsex==u
  LT=c(LT,list(BuildLifeTable(dscdane1$timesurvived[ind],dscdane1$status[ind],time.units=2)))
}
names(LT)=U
length(LT)

ind1=cumsum(LT[[1]]$Deaths)>0
ind2=cumsum(LT[[2]]$Deaths)>0
ind3=cumsum(LT[[3]]$Deaths)>0
ind4=cumsum(LT[[4]]$Deaths)>0

M1=Mort1Dsmooth(x=LT[[1]]$x2,y=LT[[1]]$Deaths,offset=log(LT[[1]]$Exposures),overdispersion=T,ndx=5,pord=3)
M2=Mort1Dsmooth(x=LT[[2]]$x2,y=LT[[2]]$Deaths,offset=log(LT[[2]]$Exposures),overdispersion=T,ndx=5,pord=3)
M3=Mort1Dsmooth(x=LT[[3]]$x2,y=LT[[3]]$Deaths,offset=log(LT[[3]]$Exposures),overdispersion=T,ndx=5,pord=3)
M4=Mort1Dsmooth(x=LT[[4]]$x2,y=LT[[4]]$Deaths,offset=log(LT[[4]]$Exposures),overdispersion=T,ndx=3,pord=3)

tiff(filename='./results/mortality_plot.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(mar=c(4,4,1,1))
plot(LT[[1]]$x2,log(LT[[1]]$Hazar),xlim=c(0,32),ylim=c(-5,1),ylab='log mortality',xlab='Age')
lines(LT[[2]]$x2,log(LT[[2]]$Hazar),type='p',col=2)
lines(LT[[3]]$x2,log(LT[[3]]$Hazar),type='p',col=3)
lines(LT[[4]]$x2,log(LT[[4]]$Hazar),type='p',col=4)

lines(M1$x[ind1],M1$logmortality[ind1],type='l',col=1,lwd=2)
lines(M2$x[ind2],M2$logmortality[ind2],type='l',col=2,lwd=2)
lines(M3$x[ind3],M3$logmortality[ind3],type='l',col=3,lwd=2)
lines(M4$x[ind4],M4$logmortality[ind4],type='l',col=4,lwd=2)
box()
legend('bottomright',U,col=1:4,lwd=2,lty=1,bty='n')
dev.off()


###########################
#Accessing proportionality test via cox regression
###########################

#Hazard plots suggest interaction between sex:tr
#we include this interaction into cox model
ModCox_<-coxph(Surv(timesurvived,status) ~ sex + 
                 tr + bean1 + mass1 + sex:tr+
                 frailty(mother,frailty='gamma',eps=1e-11),
               outer.max=50,
               data = dscdane1)

G=cox.zph(ModCox_)
G

G.=as.matrix(round(G$table,4))
G.=cbind(c('Sex (Males)','Treatment (Virgin)','Bean size','Adult body mass','Interaction (sex:treatment)','GLOBAL'),G)

write.csv(G.,'./results/cox.zph.test.csv',row.names=F)

#proportional assumption hold for all terms and globally

####################################
#Schoenfeld residuals
graphics.off()
par(mfrow=c(3,2))
my.plot.cox.zph(G)
par(mfrow=c(1,1))

#ploting exemplary Schoenfeld residuals for continous variables and interaction
tiff(filename='./results/Examplary_Schoenfeld.tiff',width=res*4.5,height=res*6,compression ='lzw',res=res,units='px')
par(mar=c(1,4,1,1))
par(oma=c(3,0,1,1))
par(mfrow=c(3,1))
my.plot.cox.zph(G,var=3,ylab='Beta(t) for bean size',xlab=''); #omit first two
my.plot.cox.zph(G,var=4,ylab='Beta(t) for adult body mass',xlab=''); #omit first two
my.plot.cox.zph(G,var=5,ylab='Beta(t) for interaction',xlab=''); #omit first two
mtext('Age',1,line=2.5,cex=0.7)
dev.off()

#ploting exemplary Schoenfeld residuals for continous variables and interaction
tiff(filename='./results/Whole_Schoenfeld.tiff',width=res*6,height=res*8,compression ='lzw',res=res,units='px')
par(mar=c(1,4,1,1))
par(oma=c(3,0,1,1))
par(mfrow=c(3,2))
my.plot.cox.zph(G,var=1,ylab='Beta(t) for sex (Males)',xlab=''); #omit first two
my.plot.cox.zph(G,var=2,ylab='Beta(t) for treatment (Virgin)',xlab=''); #omit first two
my.plot.cox.zph(G,var=3,ylab='Beta(t) for bean size',xlab=''); #omit first two
my.plot.cox.zph(G,var=4,ylab='Beta(t) for adult body mass',xlab=''); #omit first two
my.plot.cox.zph(G,var=5,ylab='Beta(t) for interaction',xlab=''); #omit first two
mtext('Age',1,line=2.5,cex=0.7)
dev.off()
#No clear departures from linearity, ph assumption holds

###########################
#CONCLUSINS
###########################

#Proportionality of hazards is not heavily violated
#It seems that there is interaction between treatment and sex. 
#This interaction will be tested and then included in all models


#########################################################
# Testing colinearity of parameters by mean of cox model
#########################################################
#Computes variance inflation factors from the covariance matrix of parameter estimates, 
#using the method of Davis et al. (1986), which is based on the correlation matrix 
#from the information matrix

rms:::vif(ModCox_)
rms:::vif(ModCox_)<5

#No strong colinearity

#############################################################################################
#############################################################################################
#############################################################################################

#Applying proportional hazard model
#Initial model: all main effects, one interaction tr:sex

##################################
#initial model selection gompertz vs. weibull
##################################

p1<-select.parfm(Surv(timesurvived,status) ~ sex + 
                   tr + bean1 + mass1 + tr: sex,
                 cluster='mother',  
                 frailty=c("none",'gamma'),
                 dist=c("gompertz","weibull",'exponential'),
                 data = dscdane1)
p1
p1.AIC=cbind(c('Gompertz','Weibull','Exonential'),round(as.matrix(p1$AIC),2))
write.csv(p1.AIC,'./results/baseline.gamma.aic.csv',row.names=F)

parfm:::plot.select.parfm(p1)

#1)  frailty model  improve the fit according to AIC. 
#AIC select models with gamma frailty, BIC without
#We decided to include gamma frailty even if its effect may be tiny.
#2) Weibull fits clearly better, however Gompertz "is more biologically relevant":
#more people use it for biological data,
#Notice: Weibull assumes that mortality starts from zero, which maybe not very realistic for biological data

#Decission: use Weibull baseline 
distr='weibull'

#First, test if data follows Weibull distribution
#Because we have multiple predictors the test is not decisive, but only informative:
#The mixture of different Weibull distributions is not necessarily a Weibull distribution

#We don't have censored observations so it is easier to test the distribution
sum(dscdane1$status)-length(dscdane1$status)

sub.data=subset(dscdane1$timesurvived,(dscdane1$sex=='Males')&(dscdane1$tr=='Virgin'))
FIT=fitdist(sub.data, "weibull")
tiff(filename='./results/WeibullTest_VirginMales.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(oma=c(0,0,0,0))
par(mar=c(2.5,2,1.5,0.5))
par(cex=0.9)
fitdistrplus:::plot.fitdist(FIT)
dev.off()
VM_1=ks.weibull.test(sub.data) #method 1
VM_2=WLK.test(sub.data,nsim=1e4) #method 2
fitdistrplus:::plot.fitdist(FIT)

sub.data=subset(dscdane1$timesurvived,(dscdane1$sex=='Females')&(dscdane1$tr=='Virgin'))
FIT=fitdist(sub.data, "weibull")
tiff(filename='./results/WeibullTest_VirginFemales.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(oma=c(0,0,0,0))
par(mar=c(2.5,2,1.5,0.5))
par(cex=0.9)
fitdistrplus:::plot.fitdist(FIT)
dev.off()
VF_1=ks.weibull.test(sub.data) #method 1
VF_2=WLK.test(sub.data,nsim=1e4) #method 2
fitdistrplus:::plot.fitdist(FIT)

sub.data=subset(dscdane1$timesurvived,(dscdane1$sex=='Males')&(dscdane1$tr=='Reproducing'))
FIT=fitdist(sub.data, "weibull")
tiff(filename='./results/WeibullTest_ReproducingMales.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(oma=c(0,0,0,0))
par(mar=c(2.5,2,1.5,0.5))
par(cex=0.9)
fitdistrplus:::plot.fitdist(FIT)
dev.off()
RM_1=ks.weibull.test(sub.data) #method 1
RM_2=WLK.test(sub.data,nsim=1e4) #method 2   -> significant deviation from Weibull
fitdistrplus:::plot.fitdist(FIT)

sub.data=subset(dscdane1$timesurvived,(dscdane1$sex=='Females')&(dscdane1$tr=='Reproducing'))
FIT=fitdist(sub.data, "weibull")
tiff(filename='./results/WeibullTest_ReproducingFemales.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(oma=c(0,0,0,0))
par(mar=c(2.5,2,1.5,0.5))
par(cex=0.9)
fitdistrplus:::plot.fitdist(FIT)
dev.off()
RF_1=ks.weibull.test(sub.data) #method 1
RF_2=WLK.test(sub.data,nsim=1e4) #method 2
fitdistrplus:::plot.fitdist(FIT)

WeiTe=cbind("KS Bootstrap*"=c(RF_1,VF_1,RM_1,VM_1),
"Generalized Gamma**"=c(RF_2$p.value,VF_2$p.value,RM_2$p.value,VM_2$p.value))
rownames(WeiTe)=U
WeiTe
write.csv(cbind(U,round(WeiTe,4)),'./results/Weibull_test.csv',row.names=F)

#There are some deviations, but using Weibull is not unreasonable


#################################################################################
#################################################################################
#FITTING WEIBULL PROPORTIONAL HAZARD MODEL ######################################
#################################################################################
#################################################################################


#I strongly recommend to use my modification to parfm routines
source('my.parfm.r')
#Using them:
#1:you can avoid common problems of lack of convergence
#2:Hessian calculation is much more precise (negative variance problem is less frequent)
#3:Optim calculations are more precise -> better estimates
#However still far from perfect for very complicated models, 
#for them selecting initial parameters is crucial

################
#The simplest possible model including all main factors and the suggested interaction
#The original parfm function is used
OrgMod_<-parfm(Surv(timesurvived,status) ~ sex + 
                 tr + bean1 + mass1 + sex:tr,
               cluster='mother',  
               frailty='gamma',
               dist=distr,
               data = dscdane1)

#The same result are obtained for when using my function. 
#But my.parfm  is more trustable for more complicated models, when estimation problems can occur
Mod_<-my.parfm(Surv(timesurvived,status) ~ sex + 
                 tr + bean1 + mass1 + sex:tr,
               cluster='mother',
               frailty='gamma',
               dist=distr,
               data = dscdane1)

#yes, the same results, even slightly more precise
attributes(Mod_)$loglik-attributes(OrgMod_)$loglik

#First Let's test if interaction sex:tr is significant via LRT
#the most powerful test is almost always based on LRT (only bootstrap can be better in some cases)

NoInter<-my.parfm(Surv(timesurvived,status) ~ sex + 
                    tr + bean1 + mass1 ,
                  cluster='mother',  
                  frailty='gamma',
                  dist=distr,
                  data = dscdane1)

my.LRT(obj1=Mod_,obj2=NoInter)
#truly significant!
#Mod_ is now our basic model

rownames(Mod_)
Z=c('Theta (frailty par.)','Rho (Weibull shape par.)','Lambda (Weibull scale par.)','Sex (Males)',
    'Treatment (Virgin)','Bean size','Addult body mass','Interaction (sex:treatment)')
write.csv(cbind(Z,round(Mod_,4)),'./results/BasicModel_i_test.csv',row.names=F)


#parameters of the model will be used as starting values for a "one-step" more complex models
iniP=c(Mod_[,1][-1],0)
iniFP=Mod_[1,1]


##################################
# Here we test all additional interactions
# w add an single interaction to the basic model Mod_ and check its significance via LRT
# The procedure is equivalent to add1() function from stats package
# We cannot do exhaustive search for a most parsimonious model, because computations would last too long. But see Cox model at the end of this file
##################################

##################################
#1) testing tr:bean1 interaction
##################################

Mod_tr_bean1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                         tr + bean1 + mass1 + sex:tr +tr:bean1,
                       inip=iniP,
                       iniFpar=iniFP,
                       cluster='mother',
                       frailty='gamma',
                       dist=distr,
                       data = dscdane1)

LRT_tr_bean=my.LRT(obj1=Mod_tr_bean1,obj2=Mod_)
LRT_tr_bean
#tr:bean1 is not significant, we don’t have to include this interaction in the model

##################################
#2) testing sex:bean1 interaction
##################################

Mod_sex_bean1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                          tr + bean1 + mass1 + sex:tr +sex:bean1,
                        inip=iniP,
                        iniFpar=iniFP,
                        cluster='mother',
                        frailty='gamma',
                        dist=distr,
                        data = dscdane1)
LRT_sex_bean=my.LRT(obj1=Mod_sex_bean1,obj2=Mod_)
LRT_sex_bean
#sex:bean1 is not significant, we don't have to include this interaction in the model


##################################
#3) testing sex:mass1 interaction
##################################

Mod_sex_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                          tr + bean1 + mass1 + sex:tr +sex:mass1,
                        inip=iniP,
                        iniFpar=iniFP,
                        cluster='mother',
                        frailty='gamma',
                        dist=distr,
                        data = dscdane1)
LRT_sex_mass=my.LRT(obj1=Mod_sex_mass1,obj2=Mod_)
LRT_sex_mass
#sex:mass1 is not significant, we don't have to include this interaction in the model


##################################
#4) testing tr:mass1 interaction
##################################

Mod_tr_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                         tr + bean1 + mass1 + sex:tr +tr:mass1,
                       inip=iniP,
                       iniFpar=iniFP,
                       cluster='mother',
                       frailty='gamma',
                       dist=distr,
                       data = dscdane1)
LRT_tr_mass=my.LRT(obj1=Mod_tr_mass1,obj2=Mod_)
LRT_tr_mass
#tr:mass1 is not significant, we don't have to include this interaction in the model


##################################
#5) testing bean1:mass1 interaction
##################################

Mod_bean1_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                            tr + bean1 + mass1 + sex:tr +bean1:mass1,
                          inip=iniP,
                          iniFpar=iniFP,
                          cluster='mother',
                          frailty='gamma',
                          dist=distr,
                          data = dscdane1)
LRT_bean_mass=my.LRT(obj1=Mod_bean1_mass1,obj2=Mod_)
LRT_bean_mass
#bean1:mass1 is not significant, we don't have to include this interaction in the model
#This model would not converge if standard parfm were used.

##################################
#6) testing sex:bean1:mass1 interaction
##################################

Mod_sex_bean1_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                                tr + bean1 + mass1 + sex:tr +sex:bean1:mass1,
                              inip=c(iniP,0),
                              iniFpar=iniFP,
                              cluster='mother',
                              frailty='gamma',
                              dist=distr,
                              data = dscdane1)
LRT_sex_bean_mass=my.LRT(obj1=Mod_sex_bean1_mass1,obj2=Mod_)
LRT_sex_bean_mass

#sex:bean1:mass1 is not significant, we don't have to include this interaction in the model

##################################
#7) testing tr:bean1:mass1 interaction
##################################

Mod_tr_bean1_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                                tr + bean1 + mass1 + sex:tr +tr:bean1:mass1,
                              inip=c(iniP,0),
                              iniFpar=iniFP,
                              cluster='mother',
                              frailty='gamma',
                              dist=distr,
                              data = dscdane1)
LRT_tr_bean_mass=my.LRT(obj1=Mod_tr_bean1_mass1,obj2=Mod_)
LRT_tr_bean_mass
#tr:bean1:mass1 is not significant, we don't have to include this interaction in the model


##################################
#8) testing tr:sex:mass1 interaction
##################################

Formula=Surv(timesurvived,status) ~ sex + tr + bean1 + mass1 + sex:tr +sex:mass1:tr
Mod_sex_tr_mass1<-my.parfm(Formula,
                           inip=c(iniP,0,0,0),
                           iniFpar=iniFP,
                           cluster='mother',
                           frailty='gamma',
                           dist=distr,
                           data = dscdane1)
LRT_sex_tr_mass=my.LRT(obj1=Mod_sex_tr_mass1,obj2=Mod_)
LRT_sex_tr_mass

#sex:tr:mass1 is not significant
#however there is a problem with calculation of covariance matrix, which is indeed more generic
#(and independent of frailty)
testmod <- try(eha:::coxreg(Formula, data = dscdane1,center = FALSE),silent=F)
testmod
testmod <- try(eha:::phreg(formula = Formula, data = dscdane1, dist = distr, shape = 0, center = FALSE),silent=F)
testmod
#sex:tr:mass1 will be not included into the main model

##################################
#9) testing tr:sex:bean1 interaction
##################################

Formula=Surv(timesurvived,status) ~ sex + tr + bean1 + mass1 + sex:tr +sex:bean1:tr
Mod_sex_tr_bean1<-my.parfm(Formula,
                           inip=c(iniP,0,0,0),
                           iniFpar=iniFP,
                           cluster='mother',
                           frailty='gamma',
                           dist=distr,
                           data = dscdane1)
LRT_sex_tr_bean=my.LRT(obj1=Mod_sex_tr_bean1,obj2=Mod_)
LRT_sex_tr_bean
#sex:bean1:tr is not significant and doesn't have to be included into the model

LRT_TAB=rbind(LRT_bean_mass,
LRT_sex_mass,LRT_sex_bean,
LRT_tr_mass,LRT_tr_bean,
LRT_sex_tr_bean,LRT_sex_tr_mass,
LRT_sex_bean_mass,LRT_tr_bean_mass)
rnam=c('Bean size : adult body mass',
       'Sex : adult body mass',
       'Sex : bean size',
       'Treatment : adult body mass',
       'Treatment : bean size',
       'Sex : treatment : bean size',
       'Sex : treatment : adult body mass',
       'Sex : bean size : adult body mass',
       'Treatment : bean size : adult body mass')
LRT_TAB
n=c('Added interaction',colnames(LRT_TAB))
LRT_TAB=cbind(rnam,round(LRT_TAB[,1],2),LRT_TAB[,2],round(LRT_TAB[,3],4))
colnames(LRT_TAB)=n

write.csv(LRT_TAB,'./results/Interactions_LRT.csv',row.names=F)

#####################################################################
#All likelihood ratio tests show that we don't have to extend further the model 
#by including chosen interactions.

#The basic summary of the model 
Mod_

#theta is frailty parameter
Mod_[1,1]
#rho and lambda are Weibull baseline parameters
Mod_[2,1]
Mod_[3,1]

#the results can be summarized in the plot
plot.parfm(Mod_)  #orthogonal coding can improve the plot
#but better visible as text
ci.parfm(Mod_,level=0.05) 

#predicted frailty values
u<-predict.parfm(Mod_)
D=as.matrix(sort(u))
colnames(D)='Frailty'
D


#Addiding frailty to data

##################################################################################################
##################################################################################################

###############################################
#Analysis of Deviance for parfm model
###############################################

#The p-values returned in parfm/my.parfm output are calculated from MLE of parameters and 
#their assyptotic standar errors ("logL hessian method"). 
#They are based on Wald t (or t ratio statistic)

#More "reliable" p-values are based on LRT
#We can easily manually perform Analysis of Deviance with marginality assumption

Mod_excl_sex<-my.parfm(Surv(timesurvived,status) ~  
                         tr + bean1 + mass1 , #both sex and interaction is out
                       cluster='mother',  
                       frailty='gamma',
                       dist=distr,
                       data = dscdane1)

Mod_excl_tr<-my.parfm(Surv(timesurvived,status) ~ sex + 
                        bean1 + mass1, #both tr and interaction is out
                      cluster='mother',  
                      frailty='gamma',
                      dist=distr,
                      data = dscdane1)

Mod_excl_bean1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                           tr +  mass1 + sex:tr, #bean1 is out
                         cluster='mother',  
                         frailty='gamma',
                         dist=distr,
                         data = dscdane1)

Mod_excl_mass1<-my.parfm(Surv(timesurvived,status) ~ sex + 
                           tr + bean1  + sex:tr, #mass1 is out
                         cluster='mother',  
                         frailty='gamma',
                         dist=distr,
                         data = dscdane1)

Mod_excl_sex_tr<-my.parfm(Surv(timesurvived,status) ~ sex + 
                            tr + bean1  + mass1, #interaction is out
                          cluster='mother',  
                          frailty='gamma',
                          dist=distr,
                          data = dscdane1)

cat('Analysis of deviance type II\n')
Tab=rbind(my.LRT(obj2=Mod_excl_sex,obj1=Mod_),
          my.LRT(obj2=Mod_excl_tr,obj1=Mod_),
          my.LRT(obj2=Mod_excl_bean1,obj1=Mod_),
          my.LRT(obj2=Mod_excl_mass1,obj1=Mod_),
          my.LRT(obj2=Mod_excl_sex_tr,obj1=Mod_))
stati=rbind(matrix(NA,3,2),Tab[,-2])
Tab=round(cbind(estimate=Mod_[,1],stati),6)
Tab

#However interaction is significant and type III is typically recommended
#read: http://goanna.cs.rmit.edu.au/~fscholer/anova.php
#read: https://www.mail-archive.com/r-help@stat.math.ethz.ch/msg69781.html
#However,
#The Type II analysis (Yates's method of fitting constants) is usually not
#preferred because of the underlying assumption of no interactions. This argument
#is, however, also founded on unrealistic models. Furthermore, by considering the
#power of the two methods, it is clear that Type II is preferable.
#read: Oyvind Langsrud. "ANOVA for unbalanced data: Use Type II instead of Type III sums of squares", Statistics and Computing, Volume 13, Number 2, pp. 163-167, 2003. 
#      https://www.researchgate.net/profile/Oyvind_Langsrud/publication/220286726_ANOVA_for_unbalanced_data_Use_type_II_instead_of_type_III_sums_of_squares/links/556f3ffd08aefcb861dd622f.pdf
#Coding type III is more complicated than coding type II, and I will need more time for that
#>>>Under construction

#####################################################3

#comparison of parameter estimates for models with frailty and without frailty
Formula=Surv(timesurvived,status) ~ sex + tr + bean1 + mass1 + sex:tr 
#phreg of eha package
noFr=eha:::phreg(formula = Formula, data = dscdane1, dist = distr, shape = 0, center = FALSE)
#parfm with no frailty of parfm package
Mod_noFr<-my.parfm(Formula,frailty='none', dist=distr,     data = dscdane1)

#eha package has different coding of Weibull parameters, let's convert them
L=length(noFr$coef)
logshape <- as.numeric(noFr$coef[substr(names(noFr$coef), 5, 9) == "shape"])
logscale <- as.numeric(noFr$coef[substr(names(noFr$coef), 5, 9) == "scale"])
we=exp(c(logshape, -exp(logshape) * logscale))

#tale with comparison
cbind("Gamma frailty (parfm)"=Mod_[,1][-1],"No frailty (parfm)"=Mod_noFr[,1],"No frailty (phreg)"=c(we,coef(noFr)[c(-L,-L+1)]))

#********
#the predict() function of parfm package returns only estimated frailties for different subjects
#we will use these frailties for further analysis

Pr=predict.parfm(Mod_)
FrailtyVec=sapply(dscdane1$mother, function(k) Pr[names(Pr)==k])

####################################################################
####################################################################
#Analysis of the interaction, frailty parameter excluded or included
####################################################################
####################################################################

graphics.off()

include.frailty=T
if (include.frailty) Model=Mod_ else Model=Mod_noFr

#function to calculate predicted hazards and marginal hazards for Mod_noFr
my.predict<-function(Data,Model,sex=NULL,treatment=NULL,bean=mean(Data$bean1),mass=mean(Data$mass1),
                     x=0:max(Data$timesurvived)){
  
  include.frailty=class(try(Model["theta",1],silent=T))!="try-error"
  if(include.frailty) cat('Frailty included\n')
  
  rho=Model["rho",1]
  lambda=Model["lambda",1]
  baseline.haz=lambda*rho*x^(rho-1)
  mm=model.matrix(Formula,data=Data)
  coef.val=as.list(Model[,1])
  
  if ((length(sex)>0)&&(length(treatment)>0)&&(length(bean)==1)&&(length(mass)==1)) {
    #Single estimates, 
    id=(Data$tr==treatment)&(Data$sex==sex)    
    if (bean=='mean') bean=mean(Data$bean1[id])
    if (mass=='mean') mass=mean(Data$mass1[id])
    sex.c=as.data.frame(mm)$sexMales[id][1]
    tr.c=as.data.frame(mm)$trVirgin[id][1]
    int.c=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[id][1]
    y=baseline.haz*exp(coef.val$mass1*mass+
                         coef.val$sexMales*sex.c+
                         coef.val$bean1*bean+
                         coef.val$trVirgin*tr.c+
                         unname(unlist(coef.val["sexMales:trVirgin"]))*int.c)
    if (include.frailty==T) y=y*mean(FrailtyVec)
  } else if ((length(sex)>0)&&(length(treatment)>0)&&(length(bean)==0)&&(length(mass)==1)) {
    
    id=(Data$tr==treatment)&(Data$sex==sex)    
    if (mass=='mean') mass=mean(Data$mass1[id])
    
    haz=list()
    for (beani in  seq_along(Data$bean1)){
      beanv=Data$bean1[beani]
      id=(Data$tr==treatment)&(Data$sex==sex)&(Data$bean1==beanv)    
      if (sum(id)>0){
        sex.c=as.data.frame(mm)$sexMales[id][1]
        tr.c=as.data.frame(mm)$trVirgin[id][1]
        int.c=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[id][1]
        if (include.frailty==T) ui=FrailtyVec[id] else ui=1
        haz=c(haz,list(ui*baseline.haz*exp(coef.val$mass1*mass+
                                          coef.val$sexMales*sex.c+
                                          coef.val$bean1*beanv+
                                          coef.val$trVirgin*tr.c+
                                          unname(unlist(coef.val["sexMales:trVirgin"]))*int.c)))
      }
    }
    hazMat=sapply(haz,c)
    survMat=exp(-sapply(haz,cumsum))
    y=rowSums(hazMat*survMat)/rowSums(survMat)
    
  } else if ((length(sex)>0)&&(length(treatment)>0)&&(length(bean)==1)&&(length(mass)==0)) {
    
    id=(Data$tr==treatment)&(Data$sex==sex)    
    if (bean=='mean') bean=mean(Data$bean1[id])
    
    haz=list()
    for (massi in  seq_along(Data$mass1)){
      massv=Data$mass1[massi]
      id=(Data$tr==treatment)&(Data$sex==sex)&(Data$mass1==massv)    
      if (sum(id)>0){
        sex.c=as.data.frame(mm)$sexMales[id][1]
        tr.c=as.data.frame(mm)$trVirgin[id][1]
        int.c=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[id][1]
        if (include.frailty==T) ui=FrailtyVec[id] else ui=1
        haz=c(haz,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                          coef.val$sexMales*sex.c+
                                          coef.val$bean1*bean+
                                          coef.val$trVirgin*tr.c+
                                          unname(unlist(coef.val["sexMales:trVirgin"]))*int.c)))
      }
    }
    hazMat=sapply(haz,c)
    survMat=exp(-sapply(haz,cumsum))
    y=rowSums(hazMat*survMat)/rowSums(survMat)
  } else if ((length(sex)==0)&&(length(treatment)==0)&&(length(bean)==0)&&(length(mass)==0)) {
    #marginal model for whole population
    
    haz=list()
    for (massi in  seq_along(Data$mass1))
      for (beani in  seq_along(Data$bean1)){
        massv=Data$mass1[massi]
        beanv=Data$bean1[beani]
        id=(Data$mass1==massv)&(Data$bean1==beanv)    
        if (sum(id)>0){
          sex.c=as.data.frame(mm)$sexMales[id][1]
          tr.c=as.data.frame(mm)$trVirgin[id][1]
          int.c=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[id][1]
          if (include.frailty==T) ui=FrailtyVec[id] else ui=1
          haz=c(haz,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                            coef.val$sexMales*sex.c+
                                            coef.val$bean1*beanv+
                                            coef.val$trVirgin*tr.c+
                                            unname(unlist(coef.val["sexMales:trVirgin"]))*int.c)))
        }
      }
    hazMat=sapply(haz,c)
    survMat=exp(-sapply(haz,cumsum))
    y=rowSums(hazMat*survMat)/rowSums(survMat)
  } else if ((length(sex)>0)&&(length(treatment)>0)&&(length(bean)==0)&&(length(mass)==0)) {
    
    haz=list()
    for (massi in  seq_along(Data$mass1))
      for (beani in  seq_along(Data$bean1)){
        massv=Data$mass1[massi]
        beanv=Data$bean1[beani]
        id=(Data$tr==treatment)&(Data$sex==sex)&(Data$mass1==massv)&(Data$bean1==beanv)    
        if (sum(id)>0){
          sex.c=as.data.frame(mm)$sexMales[id][1]
          tr.c=as.data.frame(mm)$trVirgin[id][1]
          int.c=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[id][1]
          if (include.frailty==T) ui=FrailtyVec[id] else ui=1
          haz=c(haz,list(baseline.haz*ui*exp(coef.val$mass1*massv+
                                            coef.val$sexMales*sex.c+
                                            coef.val$bean1*beanv+
                                            coef.val$trVirgin*tr.c+
                                            unname(unlist(coef.val["sexMales:trVirgin"]))*int.c)))
        }
      }
    hazMat=sapply(haz,c)
    survMat=exp(-sapply(haz,cumsum))
    y=rowSums(hazMat*survMat)/rowSums(survMat)
    

  } else if ((length(sex)==0)&&(length(treatment)==0)&&(length(bean)==1)&&(length(mass)==1)) {
    #Calculate mixture of the hazards
    
    idVF=(Data$sex=='Females')&(Data$tr=='Virgin')
    idRF=(Data$sex=='Females')&(Data$tr=='Reproducing')
    idVM=(Data$sex=='Males')&(Data$tr=='Virgin')
    idRM=(Data$sex=='Males')&(Data$tr=='Reproducing')
    if (bean=='mean') bean.VF=mean(Data$bean1[idVF]) else bean.VF=bean #sorry for copy-paste ;)
    if (bean=='mean') bean.VM=mean(Data$bean1[idVM]) else bean.VM=bean #I will improve this later
    if (mass=='mean') mass.VF=mean(Data$mass1[idVF]) else mass.VF=mass
    if (mass=='mean') mass.VM=mean(Data$mass1[idVM]) else mass.VM=mass
    if (bean=='mean') bean.RF=mean(Data$bean1[idRF]) else bean.RF=bean
    if (bean=='mean') bean.RM=mean(Data$bean1[idRM]) else bean.RM=bean
    if (mass=='mean') mass.RF=mean(Data$mass1[idRF]) else mass.RF=mass
    if (mass=='mean') mass.RM=mean(Data$mass1[idRM]) else mass.RM=mass
    
    sex.c.VF=as.data.frame(mm)$sexMales[idVF][1]
    tr.c.VF=as.data.frame(mm)$trVirgin[idVF][1]
    int.c.VF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVF][1]
    sex.c.RF=as.data.frame(mm)$sexMales[idRF][1]
    tr.c.RF=as.data.frame(mm)$trVirgin[idRF][1]
    int.c.RF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRF][1]
    
    sex.c.VM=as.data.frame(mm)$sexMales[idVM][1]
    tr.c.VM=as.data.frame(mm)$trVirgin[idVM][1]
    int.c.VM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVM][1]
    sex.c.RM=as.data.frame(mm)$sexMales[idRM][1]
    tr.c.RM=as.data.frame(mm)$trVirgin[idRM][1]
    int.c.RM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRM][1]
    
    
    haz.VF=baseline.haz*exp(coef.val$mass1*mass.VF+
                              coef.val$sexMales*sex.c.VF+
                              coef.val$bean1*bean.VF+
                              coef.val$trVirgin*tr.c.VF+
                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VF)
    haz.RF=baseline.haz*exp(coef.val$mass1*mass.RF+
                              coef.val$sexMales*sex.c.RF+
                              coef.val$bean1*bean.RF+
                              coef.val$trVirgin*tr.c.RF+
                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RF)
    haz.VM=baseline.haz*exp(coef.val$mass1*mass.VM+
                              coef.val$sexMales*sex.c.VM+
                              coef.val$bean1*bean.VM+
                              coef.val$trVirgin*tr.c.VM+
                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VM)
    haz.RM=baseline.haz*exp(coef.val$mass1*mass.RM+
                              coef.val$sexMales*sex.c.RM+
                              coef.val$bean1*bean.RM+
                              coef.val$trVirgin*tr.c.RM+
                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RM)
    if (include.frailty==T) {
      ui=mean(FrailtyVec[id]) 
      warning('Frailty not implemented')
    } else ui=1
    LVF=sum(idVF); LRF=sum(idRF); LVM=sum(idVM); LRM=sum(idRM)
    S=LVF+LVM+LRF+LRM
    LVF=LVF/S; LRF=LRF/S; LVM=LVM/S; LRM=LRM/S
    
    #marginal hazard
    y=ui*(haz.VF*LVF*exp(-cumsum(haz.VF))+haz.RF*LRF*exp(-cumsum(haz.RF))+
         haz.VM*LVM*exp(-cumsum(haz.VM))+haz.RM*LRM*exp(-cumsum(haz.RM)))/
      (LVM*exp(-cumsum(haz.VM))+LRM*exp(-cumsum(haz.RM))+
         LVF*exp(-cumsum(haz.VF))+LRF*exp(-cumsum(haz.RF)))
  } else if ((length(sex)==0)&&(length(treatment)==0)&&(length(bean)==1)&&(length(mass)==0)) {
    #Calculate mixture of the hazards
    
    idVF=(Data$sex=='Females')&(Data$tr=='Virgin')
    idRF=(Data$sex=='Females')&(Data$tr=='Reproducing')
    idVM=(Data$sex=='Males')&(Data$tr=='Virgin')
    idRM=(Data$sex=='Males')&(Data$tr=='Reproducing')
    if (bean=='mean') bean.VF=mean(Data$bean1[idVF]) else bean.VF=bean
    if (bean=='mean') bean.VM=mean(Data$bean1[idVM]) else bean.VM=bean
    if (bean=='mean') bean.RF=mean(Data$bean1[idRF]) else bean.RF=bean
    if (bean=='mean') bean.RM=mean(Data$bean1[idRM]) else bean.RM=bean
    
    sex.c.VF=as.data.frame(mm)$sexMales[idVF][1]
    tr.c.VF=as.data.frame(mm)$trVirgin[idVF][1]
    int.c.VF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVF][1]
    sex.c.RF=as.data.frame(mm)$sexMales[idRF][1]
    tr.c.RF=as.data.frame(mm)$trVirgin[idRF][1]
    int.c.RF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRF][1]
    
    sex.c.VM=as.data.frame(mm)$sexMales[idVM][1]
    tr.c.VM=as.data.frame(mm)$trVirgin[idVM][1]
    int.c.VM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVM][1]
    sex.c.RM=as.data.frame(mm)$sexMales[idRM][1]
    tr.c.RM=as.data.frame(mm)$trVirgin[idRM][1]
    int.c.RM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRM][1]
    
    haz.VF=list()
    for (massi in  seq_along(Data$mass1[idVF])){
      massv=Data$mass1[idVF][massi]
      if (include.frailty==T) ui=FrailtyVec[idVF][massi] else ui=1
      haz.VF=c(haz.VF,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                              coef.val$sexMales*sex.c.VF+
                                              coef.val$bean1*bean.VF+
                                              coef.val$trVirgin*tr.c.VF+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VF)))
    }
    hazMat.VF=sapply(haz.VF,c)
    survMat.VF=exp(-sapply(haz.VF,cumsum))
    
    haz.VM=list()
    for (massi in  seq_along(Data$mass1[idVM])){
      massv=Data$mass1[idVM][massi]
      if (include.frailty==T) ui=FrailtyVec[idVM][massi] else ui=1
      haz.VM=c(haz.VM,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                              coef.val$sexMales*sex.c.VM+
                                              coef.val$bean1*bean.VM+
                                              coef.val$trVirgin*tr.c.VM+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VM)))
    }
    hazMat.VM=sapply(haz.VM,c)
    survMat.VM=exp(-sapply(haz.VM,cumsum))
    
    haz.RF=list()
    for (massi in  seq_along(Data$mass1[idRF])){
      massv=Data$mass1[idRF][massi]
      if (include.frailty==T) ui=FrailtyVec[idRF][massi] else ui=1
      haz.RF=c(haz.RF,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                              coef.val$sexMales*sex.c.RF+
                                              coef.val$bean1*bean.RF+
                                              coef.val$trVirgin*tr.c.RF+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RF)))
    }
    hazMat.RF=sapply(haz.RF,c)
    survMat.RF=exp(-sapply(haz.RF,cumsum))
    
    haz.RM=list()
    for (massi in  seq_along(Data$mass1[idRM])){
      massv=Data$mass1[idRM][massi]
      if (include.frailty==T) ui=FrailtyVec[idRM][massi] else ui=1
      haz.RM=c(haz.RM,list(ui*baseline.haz*exp(coef.val$mass1*massv+
                                              coef.val$sexMales*sex.c.RM+
                                              coef.val$bean1*bean.RM+
                                              coef.val$trVirgin*tr.c.RM+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RM)))
    }
    hazMat.RM=sapply(haz.RM,c)
    survMat.RM=exp(-sapply(haz.RM,cumsum))
    
    y=rowSums(cbind(hazMat.RM*survMat.RM,hazMat.VM*survMat.VM,hazMat.RF*survMat.RF,hazMat.VF*survMat.VF))/
      rowSums(cbind(survMat.RM,survMat.VM,survMat.RF,survMat.VF))
    
  } else if ((length(sex)==0)&&(length(treatment)==0)&&(length(bean)==0)&&(length(mass)==1)) {
    #Calculate mixture of the hazards
    
    idVF=(Data$sex=='Females')&(Data$tr=='Virgin')
    idRF=(Data$sex=='Females')&(Data$tr=='Reproducing')
    idVM=(Data$sex=='Males')&(Data$tr=='Virgin')
    idRM=(Data$sex=='Males')&(Data$tr=='Reproducing')
    if (mass=='mean') mass.VF=mean(Data$mass1[idVF]) else mass.VF=mass
    if (mass=='mean') mass.VM=mean(Data$mass1[idVM]) else mass.VM=mass
    if (mass=='mean') mass.RF=mean(Data$mass1[idRF]) else mass.RF=mass
    if (mass=='mean') mass.RM=mean(Data$mass1[idRM]) else mass.RM=mass
    
    sex.c.VF=as.data.frame(mm)$sexMales[idVF][1]
    tr.c.VF=as.data.frame(mm)$trVirgin[idVF][1]
    int.c.VF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVF][1]
    
    sex.c.RF=as.data.frame(mm)$sexMales[idRF][1]
    tr.c.RF=as.data.frame(mm)$trVirgin[idRF][1]
    int.c.RF=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRF][1]
    
    sex.c.VM=as.data.frame(mm)$sexMales[idVM][1]
    tr.c.VM=as.data.frame(mm)$trVirgin[idVM][1]
    int.c.VM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idVM][1]
    
    sex.c.RM=as.data.frame(mm)$sexMales[idRM][1]
    tr.c.RM=as.data.frame(mm)$trVirgin[idRM][1]
    int.c.RM=unname(unlist(as.data.frame(mm)["sexMales:trVirgin"]))[idRM][1]
    
    haz.VF=list()
    for (beani in  seq_along(Data$bean1[idVF])){
      beanv=Data$bean1[idVF][beani]
      if (include.frailty==T) ui=FrailtyVec[idVF][beani] else ui=1
      haz.VF=c(haz.VF,list(ui*baseline.haz*exp(coef.val$mass1*mass.VF+
                                              coef.val$sexMales*sex.c.VF+
                                              coef.val$bean1*beanv+
                                              coef.val$trVirgin*tr.c.VF+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VF)))
    }
    hazMat.VF=sapply(haz.VF,c)
    survMat.VF=exp(-sapply(haz.VF,cumsum))
    
    haz.VM=list()
    for (beani in  seq_along(Data$bean1[idVM])){
      beanv=Data$bean1[idVM][beani]
      if (include.frailty==T) ui=FrailtyVec[idVM][beani] else ui=1
      haz.VM=c(haz.VM,list(ui*baseline.haz*exp(coef.val$mass1*mass.VM+
                                              coef.val$sexMales*sex.c.VM+
                                              coef.val$bean1*beanv+
                                              coef.val$trVirgin*tr.c.VM+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.VM)))
    }
    hazMat.VM=sapply(haz.VM,c)
    survMat.VM=exp(-sapply(haz.VM,cumsum))
    
    haz.RF=list()
    for (beani in  seq_along(Data$bean1[idRF])){
      beanv=Data$bean1[idRF][beani]
      if (include.frailty==T) ui=FrailtyVec[idRF][beani] else ui=1
      haz.RF=c(haz.RF,list(ui*baseline.haz*exp(coef.val$mass1*mass.RF+
                                              coef.val$sexMales*sex.c.RF+
                                              coef.val$bean1*beanv+
                                              coef.val$trVirgin*tr.c.RF+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RF)))
    }
    hazMat.RF=sapply(haz.RF,c)
    survMat.RF=exp(-sapply(haz.RF,cumsum))
    
    haz.RM=list()
    for (beani in  seq_along(Data$bean1[idRM])){
      beanv=Data$bean1[idRM][beani]
      if (include.frailty==T) ui=FrailtyVec[idRM][beani] else ui=1
      haz.RM=c(haz.RM,list(ui*baseline.haz*exp(coef.val$mass1*mass.RM+
                                              coef.val$sexMales*sex.c.RM+
                                              coef.val$bean1*beanv+
                                              coef.val$trVirgin*tr.c.RM+
                                              unname(unlist(coef.val["sexMales:trVirgin"]))*int.c.RM)))
    }
    hazMat.RM=sapply(haz.RM,c)
    survMat.RM=exp(-sapply(haz.RM,cumsum))
    
    y=rowSums(cbind(hazMat.RM*survMat.RM,hazMat.VM*survMat.VM,hazMat.RF*survMat.RF,hazMat.VF*survMat.VF))/
      rowSums(cbind(survMat.RM,survMat.VM,survMat.RF,survMat.VF))
    
  } else stop('Combination of argumnts is not implemented')
  y
}

#############################################################################################
#predicted hazards for different sex and treatments at global mean values for mass1 and bean1
#This prediction is clearly missleading. Similarly misleading is using global mean values for plotting
#survivorships as Darek did with survreg (?)
#############################################################################################
#NOT FOR PUBLICATION
L2=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Virgin')
L4=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Virgin')
L1=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Reproducing')
L3=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Reproducing')
ylim=range(c(L1[-1],L2[-1],L3[-1],L4[-1]))
plot(x,log(L1),type='l',ylim=log(ylim),ylab='log fitted hazard',xlab='Age',main='Hazard at global mean values of bean and mass ')
lines(x,log(L2),col=2)
lines(x,log(L3),col=3)
lines(x,log(L4),col=4)
legend('bottomright',legend=U,col=1:4,lty=1,bty='n')

#############################################################################################
#predicted hazards for different sex and treatments at local (sex/tr - specific) 
#mean values for mass1 and bean1
#############################################################################################
#Better than previous, but these are just curves for hypothetical single individuals, we need curve for sub-populations
L2=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Virgin',mass='mean',bean='mean')
L4=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Virgin',mass='mean',bean='mean')
L1=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Reproducing',mass='mean',bean='mean')
L3=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Reproducing',mass='mean',bean='mean')
ylim=range(c(L1[-1],L2[-1],L3[-1],L4[-1]))
plot(x,log(L1),type='l',ylim=log(ylim),ylab='log fitted hazard',xlab='Age',main='Hazard at sex nad treatment specific mean values of bean and mass ')
lines(x,log(L2),col=2)
lines(x,log(L3),col=3)
lines(x,log(L4),col=4)
legend('bottomright',legend=U,col=1:4,lty=1,bty='n')

#############################################################################################
#predicted marginal hazards for different sex and treatments 
#calculated as marginal hazards for mass1 and bean1
#############################################################################################
#Marginal hazards for bean and mass as a function of sex and treatment, suggested plot for publication
L0=my.predict(Data=dscdane1,Model=Model,sex=NULL,treatment=NULL,mass=NULL,bean=NULL)
L2=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Virgin',mass=NULL,bean=NULL)
L4=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Virgin',mass=NULL,bean=NULL)
L1=my.predict(Data=dscdane1,Model=Model,sex='Females',treatment='Reproducing',mass=NULL,bean=NULL)
L3=my.predict(Data=dscdane1,Model=Model,sex='Males',treatment='Reproducing',mass=NULL,bean=NULL)
ylim=range(c(L1[-1],L2[-1],L3[-1],L4[-1]))

tiff(filename='./results/marginal_hazard_1_plot.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(mar=c(4,4,1,1))
plot(x,log(L1),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard for bean, mass, and frailty')
lines(x,log(L2),col=2);lines(x,log(L3),col=3);lines(x,log(L4),col=4);lines(x,log(L0),col='gold',lwd=2,lty=2)
legend('bottomright',legend=c(U,'Whole population'),col=c(1:4,'gold'),lty=c(1,1,1,1,2),lwd=c(1,1,1,1,2),bty='n',cex=0.8)
dev.off()

plot(x,log(L1),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard for bean, mass and frailty')
lines(x,log(L2),col=2);lines(x,log(L3),col=3);lines(x,log(L4),col=4);lines(x,log(L0),col='gold',lwd=2,lty=2)
legend('bottomright',legend=c(U,'Whole population'),col=c(1:4,'gold'),lty=c(1,1,1,1,2),lwd=c(1,1,1,1,2),bty='n')

#############################################################################################
#predicted marginal hazard for different values of bean size (within range of the data)
#calculated as marginal hazards for mass1 and sex and treatment
#############################################################################################
#Marginal hazards for sex, treatment and mass as function of bean size, suggested plot for publication
N=10
colpal=rev(sapply(seq(1,0.2,len=N),function(k) adjustcolor(col='white',green.f=0,alpha.f=k,red.f=k,blue.f=1-k^2)))
V=(seq(min(dscdane1$bean1),max(dscdane1$bean1), len=N))
L=sapply(V,function(my.bean) my.predict(Data=dscdane1,Model=Model,sex=NULL,treatment=NULL,mass=NULL,bean=my.bean))
ylim=range(L[-1,])

tiff(filename='./results/marginal_hazard_2_plot.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(mar=c(4,4,1,1))
plot(x,log(L[,1]),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard as function of bean size',col='white')
for (j in 2:N) lines(x,log(L[,j]),col=colpal[j])
lines(x,log(L0),col='gold',lwd=2,lty=2)
legend('bottomright',legend=c('Bean size:',round(V,1),'Whole population'),col=c(NA,colpal,'gold'),lty=c(NA,rep(1,N),2),lwd=c(NA,rep(1,N),2),bty='n',cex=0.8)
dev.off()

plot(x,log(L[,1]),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard as function of bean size',col='white')
for (j in 2:N) lines(x,log(L[,j]),col=colpal[j])
lines(x,log(L0),col='gold',lwd=2,lty=2)
legend('bottomright',legend=c('Bean size:',round(V,1),'Whole population'),col=c(NA,colpal,'gold'),lty=c(NA,rep(1,N),2),lwd=c(NA,rep(1,N),2),bty='n')

#############################################################################################
#predicted marginal hazard for different values of adult body size (within range of the data)
#calcualted as marginal hazards for mass1 and bean1
#############################################################################################
#Marginal hazards for bean size, sex and treatment as a function of adult body size, suggested plot for publication
N=10
V=seq(min(dscdane1$mass1),max(dscdane1$mass1), len=N)
colpal=rev(sapply(seq(1,0.2,len=N),function(k) adjustcolor(col='white',green.f=0,alpha.f=k,red.f=k,blue.f=1-k^2)))
L=sapply(V,function(my.mass) my.predict(Data=dscdane1,Model=Model,sex=NULL,treatment=NULL,bean=NULL,mass=my.mass))
ylim=range(L[-1,])

tiff(filename='./results/marginal_hazard_3_plot.tiff',width=res*6,height=res*4,compression ='lzw',res=res,units='px')
par(mar=c(4,4,1,1))
plot(x,log(L[,1]),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard as function of adult body mass',col='white')
for (j in 1:N) lines(x,log(L[,j]),col=colpal[j])
lines(x,log(L0),col='green3',lwd=2,lty=1)
legend('bottomright',legend=c('Adult body mass:',round(V,1),'Whole population'),col=c(NA,colpal,'green3'),lty=c(NA,rep(1,N),1),lwd=c(NA,rep(1,N),2),bty='n',cex=0.8)
dev.off()

plot(x,log(L[,1]),type='l',ylim=log(ylim),ylab='log marginal hazard of fitted model',xlab='Age',main='Marginal hazard as function of adult body mass',col='white')
for (j in 1:N) lines(x,log(L[,j]),col=colpal[j])
lines(x,log(L0),col='green3',lwd=2,lty=1)
legend('bottomright',legend=c('Adult body mass:',round(V,1),'Whole population'),col=c(NA,colpal,'green3'),lty=c(NA,rep(1,N),1),lwd=c(NA,rep(1,N),2),bty='n')


#____________________________________________________________________________________________________
#_____________________________________BONUS__________________________________________________________
#____________________________________________________________________________________________________

#################################################################################
#################################################################################
#FITTING COX PROPORTIONAL HAZARD MODEL ##########################################
#################################################################################
#################################################################################

#The alternative approach is semi-parametric cox ph model with gamma frailty
#read intro to cox: http://socserv.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf

#Parametric vs. semiparametric approach

#read: https://krex.k-state.edu/dspace/bitstream/handle/2097/8787/angelacrumer2011.pdf?sequence=3
#"When the shape parameter is known, the Weibull model far out performs the 
#Cox proportional hazards model, but when the shape parameter is unknown, the 
#Cox proportional hazards model and the Weibull model give comparable results."

#sometimes parametric is better
#read: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4545306/
#read: https://www.ncbi.nlm.nih.gov/pubmed/18159979
#read: http://cemsiis.meduniwien.ac.at/fileadmin/msi_akim/CeMSIIS/KB/volltexte/Nardi_Schemper_2003_Statistics_in_Medicine.pdf

#sometimes cox is better
#rad: https://stats.stackexchange.com/questions/64739/in-survival-analysis-why-do-we-use-semi-parametric-models-cox-proportional-haz


ModCox_<-coxph(Surv(timesurvived,status) ~ sex + 
                 tr + bean1 + mass1 + sex:tr+
                 frailty(mother,frailty='gamma',eps=1e-11),
               outer.max=50,
               data = dscdane1)

summary(ModCox_)

#more clear to see
as.matrix(round(summary(ModCox_)$coefficients,6))

#qualitatively the same results as for parfm Mod_
#as only bean1 is not significant

#################################################
#Exhaustive AIC model selection for coxph model
#Main terms always present
#################################################

dscdane1<-within(dscdane1, SurvVec<-Surv(timesurvived,status))
dependent='SurvVec'
basicterms=paste('sex','tr','bean1','mass1',sep=' + ')

SelectCoxModel<-function(terms){
  modelslist<-lapply(seq_along(terms),function(k) 
    paste(dependent,'~',basicterms,'+',apply(X=combn(terms,k),MARGIN=2,paste, collapse='+'),"+frailty(mother,frailty='gamma',eps=1e-11)",sep=''))
  modelsvec=unlist(modelslist)
  fitslist<-lapply(modelsvec,function(k){
    cat('.')
    fitAIC=extractAIC(coxph(as.formula(k),data=dscdane1,outer.max=50))
    data.frame(n.terms=fitAIC[1],AIC=fitAIC[2],formula=k)
  })
  
  TAB=(do.call(rbind,fitslist))
  o=order(TAB$AIC)
  TAB=TAB[o,]
  TAB
}

terms=c('sex:tr','mass1:tr','bean1:tr','sex:mass1','sex:bean1','mass1:bean1','sex:tr:mass1','sex:tr:bean1',
        'tr:mass1:bean1','sex:mass1:bean1')
TAB=SelectCoxModel(terms)
head(TAB)

bestmodel=coxph(as.formula(as.character(TAB$formula[1])),data=dscdane1,outer.max=50)
W=summary(bestmodel)
D=as.matrix(round(W$coefficients[,c(1,4:6)],6))
colnames(D)=c('coef','Chisq','df','p-val')
D

#It is commonly known that AIC selects too complicated models
#Most effects of significant interactions are tiny
#The interpretation of such a model is extremely difficult
#The model has clearly too many parameters, we can try to reduce it by eliminating 3-way interactions
#3-way interactions would be anyway hard to explain

terms=c('sex:tr','mass1:tr','bean1:tr','sex:mass1','sex:bean1','mass1:bean1')
TAB2=SelectCoxModel(terms)
head(TAB2)

bestmodel=coxph(as.formula(as.character(TAB2$formula[1])),data=dscdane1,outer.max=50)
W=summary(bestmodel)
D=as.matrix(round(W$coefficients[,c(1,4:6)],6))
colnames(D)=c('coef','Chisq','df','p-val')
D

#Much better than above and comparable to the results obtained by parfm
#let's do LRT selection for the last interaction tr:bean as it seems to be non-significant

Z=drop1(bestmodel,test='Chisq')
as.matrix(Z[5,])
#it is clear that the last interaction tr:bean can be dropped from the model
#as it is not significant

bestmodel = ModCox_ #not used later, just informative

###################################################################
#Checking model coxph for proportionality assumption
###################################################################
#***this was already tested at the beginning, just to remind
G=cox.zph(ModCox_)
G
#proportional assumption hold for all terms and globally
########################
#Schoenfeld residuals
par(mfrow=c(3,2));plot(G);par(mfrow=c(1,1))
#No clear departures from linearity

###################################################################
#Checking model coxph for influential observations
###################################################################
dfbeta <- residuals(ModCox_, type="dfbeta")
par(mfrow=c(3, 2))
for (j in 1:5) {
  plot(dfbeta[, j], ylab=names(coef(ModCox_))[j])
  abline(h=0, lty=2)
}
par(mfrow=c(1, 1))
#No striking influential observations

###################################################################
#Checking model coxph for linearity via martingale residuals
#only relevant for continuous variables: mass1 and bean1
###################################################################

res <- residuals(ModCox_, type="martingale")
X <- as.matrix(dscdane1[, c("mass1", "bean1")]) # matrix of covariates
par(mfrow=c(2, 2))
for (j in 1:2) {  # residual plots
  plot(X[, j], res, xlab=c("mass1", "bean1")[j], ylab="residuals")
  abline(h=0, lty=2)
  lines(lowess(X[, j], res, iter=0))
}
b <- coef(ModCox_)[c(2,3)]  # regression coefficients
for (j in 1:2) {  # component-plus-residual plots
  plot(X[, j], b[j]*X[, j] + res, xlab=c("mass1", "bean1")[j],
       ylab="component+residual")
  abline(lm(b[j]*X[, j] + res ~ X[, j]), lty=2)
  lines(lowess(X[, j], b[j]*X[, j] + res, iter=0))
}
par(mfrow=c(1, 1))

#linearity hold almost perfectly for both predictors

###################################################################
# Analysis of deviance of the cox ph model
###################################################################

# This is more powerful test so it can detect even tiny effect of 
#bean1 and frailty not detected in coxph summary
options(contrasts = c("contr.treatment","contr.poly")) #default options
Anova(ModCox_,type=2)

#even more clear if more appropriate type III is used
options(contrasts = c("contr.sum","contr.poly")) #set orthogonal contrasts
Anova(ModCox_,type=3)
options(contrasts = c("contr.treatment","contr.poly")) #back to default options

#compared to wald test of coxph()
as.matrix(round(summary(ModCox_)$coefficients,6))

#Models that contain a frailty term are a special case: 
#due to the technical difficulty, when there is a newdata argument (in predict.coxph)
#the predictions will always be for a random effect of zero. 
#prediction of cox model are under construction now
