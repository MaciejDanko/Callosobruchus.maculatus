#Constructing lifetable for individual data with censored observations
#EventTime -> Event time
#Censoring -> 1-exact observation, 0-censoring
#time.units -> output age interval
BuildLifeTable<-function (EventTime,Censoring,time.start=0,time.end=max(EventTime),time.units=7,rm.zero.ex=T){
  Censoring=as.numeric(Censoring)
  if (any(Censoring>1)|any(Censoring<0)){
    print("Error")
    stop
    break
  }
  N.time=(floor((time.end-time.start)/time.units)+2)*time.units+time.start
  X=seq(from=time.start,by=time.units,to=N.time)
  Haz=rep(NA,length(X)-1)
  Exposures=rep(NA,length(X)-1)
  Exposures2=rep(NA,length(X)-1)
  Deaths=rep(NA,length(X)-1)
  NRisk=rep(NA,length(X)-1)
  NXg=length(Censoring)
  NX=NXg
  for (xi in 1:length(X)-1){
    index=(EventTime>X[xi])&(EventTime<=X[xi+1])
    index2=(EventTime>X[xi])
    DX=sum(Censoring[index])
    NX=length(Censoring[index2])
    NXm=(NX-length(Censoring[index]))
    NRisk[xi]=NX
    Haz[xi]=-log(1-DX/NX)
    Exposures[xi]=(NXm+sum((EventTime[index]-X[xi])))
    Deaths[xi]=DX
    
  }
  Mat=NULL
  X=X[1:length(X)-1]
  if (rm.zero.ex) ind=Exposures>0 else ind=1:length(Exposures)
  Mat$x=X[ind]
  Mat$x2=(X+time.units/2)[ind]
  Mat$Hazard=Haz[ind] 
  Mat$Deaths=Deaths[ind]
  Mat$Exposures=Exposures[ind]
  Mat$TimeUnit=time.units
  Mat$NRisk=NRisk[ind]
  return(Mat)
}

#Wrapper of previous function
#D = Event time
#C= Censoring -> 1-exact observation, 0-censoring
#time.units -> output age interval
Get_lx<-function(D,C,time.units=0.25){
  LT=BuildLifeTable(D,C,time.units=time.units)
  lx=exp(-cumsum(c(0,LT$Hazard[-length(LT$Hazard)])))*length(D)
  lx
}

#some basic "cohort" operations when there is no censoring
dx_2_lx<-function(dx) sum(dx)-c(0,cumsum(dx))[1:length(dx)]
lx_2_dx<-function(lx) c(-diff(lx),lx[length(lx)])
dx_2_Times_No_C<-function(dx,x) unlist(sapply(1:(length(x)),function (kk) rep(my.step(x)[kk],dx[kk])))
Times_2_dx_No_C<-function(times,x) (hist(x=times,plot=F,right=F,breaks=c(x,x[length(x)]+1000))$counts)
