#Kolmogorov-Smirnov test simulation
#read: https://stats.stackexchange.com/questions/126539/testing-whether-data-follows-t-distribution/126552#126552
#read: https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best 
#read: https://stats.stackexchange.com/questions/60511/weibull-distribution-parameters-k-and-c-for-wind-speed-data/60530#60530
#read: Jogesh Babu, G., and C. R. Rao. "Goodness-of-fit tests when parameters are estimated." Sankhya: The Indian Journal of Statistics 66 (2004): 63-74.

require(fitdistrplus)
require(logspline)

ks.weibull.test<-function(data, boot.null=1e4){
  
  #fit weibull distribution to the data
  Fit=fitdist(data, "weibull")
  
  #simulation of null distribution of the KS statistics
  KS.null.simulation = replicate(boot.null, {      
    #draw random sample form weibull distribution
    r = rweibull(n = length(data), shape= Fit$estimate["shape"], scale = Fit$estimate["scale"])
    #calculate statisics
    z = ks.test(r, "pweibull", shape= Fit$estimate["shape"], scale = Fit$estimate["scale"])$statistic     
  })
  
  #Fits a logspline density using splines to approximate the log-density
  lspline.fit = logspline(KS.null.simulation)
  
  #Calculate KS statistic for the data vs. Weibull cdf
  KS.stat=ks.test(data, "pweibull", shape= Fit$estimate["shape"], 
                  scale = Fit$estimate["scale"])$statistic
  
  #calculate p-value using the simulated null distribution of the KS-statistics
  pval = 1 - plogspline(KS.stat, lspline.fit) #1-cdf
    
  c(pval=pval)
}

#UNDER CONSTRUCTION
# #Use the non-parametric bootstrap to construct pointwise confidence intervals around the PDF and CDF of the estimated Weibull distribution
# ks.weibull.boot<-function(data, data.range=seq(min(data),max(data), len=1000), n.boot=1e3, what=c('pdf','cdf')){
#   
#   if (!what %in% c('pdf','cdf')) stop('Wrong what argumnet, use pdf or cdf')
#   what=what[1]
#   if (what==pdf) bfun=dweibull else bfun=pweibull
#   
#   #fit weibull distribution to the data
#   Fit=fitdist(data, "weibull")
#   
#   #Theor <- rweibull(1e6, shape= Fit$estimate["shape"], scale = Fit$estimate["scale"])
#   
#   wboot <- sapply(1:n.boot, function(i) {
#     #random sample 
#     xi <- sample(data, size=length(data), replace=TRUE)
#     #fit weibool with Maximum likelihood method
#     MLE.est <- suppressWarnings(fitdist(xi, distr="weibull"))  
#     #use estimate weibull MLE 
#     bfun(data.range, shape=MLE.est$estimate["shape"],  scale = MLE.est$estimate["scale"])
#   })
#   
#   par(bg="white", las=1, cex=1.2)
#   gcol=adjustcolor(1,alpha.f=0.5)
#   plot(data.range, wboot[, 1], type="l", col=gcol, ylim=range(wboot),
#        xlab="data", ylab=what)
#   for(i in 2:ncol(wboot)) lines(data.range, wboot[, i], col=gcol)
#   
#   # Add pointwise confidence bands
#   
#   quants <- apply(wboot, 1, quantile, c(0.025, 0.5, 0.975))
#   lines(data.range, quants[1, ], col=4, lwd=1.5, lty=2)
#   lines(data.range, quants[3, ], col=4, lwd=1.5, lty=2)
#   lines(data.range, quants[2, ], col=4, lwd=2)   
#   
#   
# }