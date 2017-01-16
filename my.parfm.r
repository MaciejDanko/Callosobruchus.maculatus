require(DEoptim)
require(numDeriv)
require(eha)
#substitution for optim function

#Function performing Likelihood ratio test for parfm objects
my.LRT<-function(obj1,obj2){
  #obj1 is alternative model
  #obj2 is nested model
  l1=attributes(obj1)$loglik  
  l2=attributes(obj2)$loglik  
  df1=dim(obj1)[1]
  df2=dim(obj2)[1]
  Z=c(chisq=2*(l1-l2),df=df1-df2,pval=1 - pchisq(2*(l1-l2), df1-df2))
  Z
}

#extracting variance covariance matrix, it works only with my.parfm() it does not work with parfm()
vcov.parfm<-function(fit) attributes(fit)$varcov

#Variance Inflation Factor routines based on rms:::vif() function
#remove should contain name of equivalent of the intercept
#for Weibull it is "lmbda" parameter
vif.parfm<-function (fit,remove=c('lambda')) {
  #read: Davis CE, Hyde JE, Bangdiwala SI, Nelson JJ: 
  #An example of dependencies among variables in a conditional 
  #logistic regression. In Modern Statistical Methods in Chronic Disease Epidemiology, 
  #Eds SH Moolgavkar and RL Prentice, pp. 140–147. New York: Wiley; 1986.
  #read: http://www.how2stats.net/2011/09/variance-inflation-factor-vif.html
  require(Hmisc)
  v <- vcov.parfm(fit)
  nam <- dimnames(fit)[[1]]
  dimnames(v)=list(nam,nam)
  
  #dropping equivalent of intercept given in remove par
  pos=sapply(remove,function(k) which(nam==k))
  if (length(pos)>0) {
    v=v[-pos,-pos]
    nam=nam[-pos]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

#function to split string (txt) accoring to given character (char)
splitstr<-function(char,txt){
  (sapply(seq_along(txt),function(j){
    char=substr(char,1,1)
    z=gregexpr(char,txt[j])
    if(z<0){
      txt[j]
    } else {
      st=c(1,z[[1]]+1)
      en=c(z[[1]],nchar(txt[j])+1)-1
      gsub(' ','',sapply(seq_along(en),function(k) substr(txt[j],st[k],en[k])),fixed=T)
    }
  }))  
}

rowProds<-function(M) apply(M,1,prod)
colProds<-function(M) apply(M,2,prod)

#function predicting marginal hazard or linear predictor
#Var.Name = variable name to be substituted
#Value = new value of this variable
#Subset = Subset of Data
#x = an age/time vector
#Model = model ftted by my.parfm

#rm(list=c('f1','ZZ','mm','j','rho','lambda','Terms','Terms2','T1','T2','ind','k','coefi','new_mm','y'))
#l.p - linear predictors not tested yet, 
my.predict.fit.parfm<-function(Model,Data,max.x=NULL,Var.Name=NULL,Value=NULL,Type=c('marginal','l.p'),Subset=NULL){
  Type=Type[1]
  if (length(Subset)>0) Data=subset(Data,Subset)
  if (missing(Data)) stop('Data must be given')
  if(!grepl('Surv',attr(Model,'formula'))) {
    warning('formula should be specified explicitely, trying recovering it, but errors can occur.')
    f1=eval(parse(text=attr(Model,'formula')))
    formula=paste(f1[2],f1[3],sep='~') 
  } else formula=attr(Model,'formula')
  if (length(max.x)==0) {
    ZZ=substr(formula,1,regexpr('~',formula,fixed=T)-1)
    S=eval(parse(text=paste('with(data = Data,',ZZ,')',sep='')))
    x=sort(unique(as.numeric(gsub(' ','',gsub('+','',S,fixed=T)))))
    x=seq(0,max(x),1)
  } else x=seq(0,max.x,1)
  rho=Model["rho",1]
  lambda=Model["lambda",1]
  baseline.haz=lambda*rho*x^(rho-1)
  mm=model.matrix.default(as.formula(attributes(Model)$formula),data=Data)[,-1]
  Terms=colnames(mm)
  Terms2=attr(Model,'terms')
  TermsX=Terms2[!grepl(':',Terms2,fixed=T)]
  if (length(Terms)!=length(Terms2)) stop('Something wrong')
  T2=splitstr(':',Terms2); T1=splitstr(':',Terms)
  if (!all(sapply(seq_along(unlist(T1)),function(k) grepl(unlist(T2)[k],unlist(T1)[k],fixed=T)))) stop('Something wrong')
  if (length(Var.Name)>0) {
    if (length(Value)>1) {
      warning('Only first element of Value was considered')
    }
    if (length(Value)==0) stop(paste('Please give Value for',Var.Name))
    ind.nam=(Terms2==Var.Name)
    if (sum(ind.nam)==0) stop (paste('Unknown Var.Name, use:  ',paste(TermsX,collapse=', ')))
    mm[,ind.nam]=Value
  }
  #ind.m=sapply(T1,length)==1
  new_mm=mm*NA
  for(j in seq_along(T1))
    if (length(T1[[j]])==1) new_mm[,j]=mm[,j] else new_mm[,j]=rowProds(sapply(seq_along(T1[[j]]),function(k) mm[,Terms%in%T1[[j]][k]]))
  
  coefi=Model[,1]
  ind=names(coefi)%in%colnames(new_mm)
  coefi=coefi[ind]
  M=new_mm*(t(matrix(coefi,length(coefi),dim(new_mm)[1])))
  
  if (Type=='marginal') {  
    Pr=predict.parfm(Model)
    Fr=eval(parse(text=paste('Data','$',attributes(Pr)$clustname,sep='')))
    ui=sapply(Fr, function(k) Pr[names(Pr)==k])
    hazMat=(ui*exp(rowSums(M)))%*%t(baseline.haz)
    cumhazMat=t(apply(hazMat,1,cumsum))
    survMat=exp(-cumhazMat)
    y=colSums(hazMat*survMat)/colSums(survMat)
  } else if (Type=='l.p'){
    y=sum(M)
  } else stop('Unknown type')
  y
}

#Calculate Martingale residuals for the fit obtained from my.parfm()
#Function under construction
# - Calcualtes only martingale now
# - Only Weibull model! 
# - "data" must be included in function call
#Fit = model ftted by my.parfm
my.residuals.parfm<-function(fit, type='martingale'){
  if (attributes(fit)$dist!='weibull') stop('Only weibull implemented so far')
  
  if(!grepl('Surv',attr(fit,'formula'))) {
    warning('formula should be specified explicitely, trying recovering it, but errors can occur.')
    f1=eval(parse(text=attr(fit,'formula')))
    formula=paste(f1[2],f1[3],sep='~') 
  } else formula=attr(fit,'formula')
  
  datan=(attributes(fit)$call$data)
  Z=substr(formula,1,regexpr('~',formula,fixed=T)-1)
  if ((length(datan)!=0)&&(nchar(datan)>0)) {
    S=eval(parse(text=paste('with(data = ',datan,',',Z,')',sep='')))
  } else{
    S=eval(parse(text=paste(Z)))
    stop('Data must be given in function call')
  }  
  Time=as.numeric(gsub(' ','',gsub('+','',S,fixed=T)))
  Status=(!grepl('+',S,fixed=T))*1
  include.frailty=attributes(fit)$frailty!='none'
  if (include.frailty){
    Pr=predict.parfm(fit)
    Fr=eval(parse(text=paste(datan,'$',attributes(Pr)$clustname,sep='')))
    FrailtyVec=sapply(Fr, function(k) Pr[names(Pr)==k])
  } else FrailtyVec=1
  rho=fit["rho",1]
  lambda=fit["lambda",1]
  cum.baseline.haz=lambda*Time^rho
  data=eval(datan)
  mm=model.matrix(as.formula(formula),data=data)[,-1]# remove intercept
  coefi=fit[,1]
  ind=names(coefi)%in%colnames(mm)
  coefi=coefi[ind]
  M=exp(rowSums(mm*(t(matrix(coefi,length(coefi),dim(mm)[1])))))
  resid=Status-cum.baseline.haz*M*FrailtyVec
  resid
}

#Improved optim function that can solve some convergence problems
#Under construction. In next version an option for controling negative vaiance problem will be added
toptim<-function(par,control,method,...,max.times=10,min.times=3) {
  cat(method,': ')
  orgmethod=method
  control$reltol=(.Machine$double.eps)^0.575
  cat(1,'(',control$maxit,') ')
  R=stats:::optim(par=par,control=control,method=method,...)
  for (j in 1:(max.times-1)) {
    control$maxit=control$maxit*2
    cat(j+1,'(',control$maxit,') ')
    if (R$convergence==10) {
      cat("\nForced switch of the optim method to BFGS\n")
      method="BFGS" 
    } else method=orgmethod
    R=stats:::optim(par=R$par,control=control,method=method,...)
    if ((R$convergence==0) && (j>=(min.times-1))) break
  }
  cat('Anticipated convergence = ',R$convergence==0,'\n')
  R
}

#My optimization function for my.parfm() (uses toptim())
my.optim<-function(par,fn,obs,dist,frailty,correct,lower,upper,do.DE.optim=FALSE,method='BFGS',hessian=TRUE,control=list()){
  #run several methods of optimization
  #hessian is set to true to predict some errors (will be removed in next version)
  Res1=suppressWarnings(try(toptim(par=par,fn=parfm:::Mloglikelihood,gr=NULL,method="Nelder-Mead",hessian=TRUE,control=control,lower=-Inf,upper=Inf,
                    obs=obs,dist=dist,frailty=frailty,correct=correct),silent=T))
  tmpRes1=Res1
  if (class(Res1)=="try-error") Res1=NULL 
  Res2=suppressWarnings(try(toptim(par=par,fn=parfm:::Mloglikelihood,gr=NULL,method="BFGS",hessian=TRUE,control=control,lower=-Inf,upper=Inf,
                     obs=obs,dist=dist,frailty=frailty,correct=correct),silent=T))
  tmpRes2=Res2
  if (class(Res2)=="try-error") Res2=NULL 
  if (do.DE.optim){ #currently not used Differential Evolution Algorithm
    ctrl=DEoptim.control(itermax = max(50,control$maxit), trace = FALSE)
    deo=DEoptim(fn=fn, lower=lower,upper=upper,control=ctrl,
                obs=obs,dist=dist,frailty=frailty,correct=correct)
    Res3=suppressWarnings(toptim(par=deo$optim$bestmem,fn=parfm:::Mloglikelihood,gr=NULL,method=method,hessian=TRUE,control=control,lower=-Inf,upper=Inf,
             obs=obs,dist=dist,frailty=frailty,correct=correct))
    if (class(Res3)=="try-error") Res3=NULL 
  } else Res3=NULL
  #selecting best method
  V=list(Res1,Res2,Res3)
  V=V[!sapply(V,is.null)]
  if (length(V)==0) stop('Error: no convergence')
  ind=which.min(sapply(1:length(V),function(k) V[[k]]$value))[1]
  cat('  Best method: ',c('Nelder-Mead','BFGS','DE+opt')[ind])
  Res=V[[ind]]
  #much more accurate calcualtion of hessian than used in optim()
  Res$hessian=try(numDeriv:::hessian(fn,Res$par,obs=obs,dist=dist,frailty=frailty,correct=correct),silent=T)
  Res
}

#rewritten from parfm package
weibull<-function (pars, t, what) {
  if (what == "H") 
    return(pars[2] * t^(pars[1]))
  else if (what == "lh") 
    return(log(pars[1]) + log(pars[2]) + ((pars[1] - 1) * 
                                            log(t)))
}

#rewritten from parfm package
gompertz<-function (pars, t, what) {
  if (what == "H") 
    return(pars[2]/pars[1] * (exp(pars[1] * t) - 1))
  else if (what == "lh") 
    return(log(pars[2]) + pars[1] * t)
}

#Maximum likelihood function rewrited from parfm package
#Modifications: 
#1) exchange NaNs from the output into 1e100 values: necessary for DE algorithm
my.Mloglikelihood<-function (p, obs, dist, frailty, correct) {
  #print(head(p))
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(p[1])
  }
  else if (frailty == "lognormal") {
    sigma <- exp(p[1])
  }
  else if (frailty == "possta") {
    nu <- exp(-exp(p[1]))
    D <- max(obs$dqi)
    Omega <- Omega(D, correct = correct, nu = nu)
  }
  if (dist %in% c("weibull", "inweibull")) {
    pars <- cbind(rho = exp(p[obs$nFpar + 1:obs$nstr]), lambda = exp(p[obs$nFpar + 
                                                                         obs$nstr + 1:obs$nstr]))
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  }
  else if (dist == "exponential") {
    pars <- cbind(lambda = exp(p[obs$nFpar + 1:obs$nstr]))
    beta <- p[-(1:(obs$nFpar + obs$nstr))]
  }
  else if (dist == "gompertz") {
    pars <- cbind(gamma = exp(p[obs$nFpar + 1:obs$nstr]), 
                  lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  }
  else if (dist == "lognormal") {
    pars <- cbind(mu = p[obs$nFpar + 1:obs$nstr], sigma = exp(p[obs$nFpar + 
                                                                  obs$nstr + 1:obs$nstr]))
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  }
  else if (dist == "loglogistic") {
    pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr], kappa = exp(p[obs$nFpar + 
                                                                     obs$nstr + 1:obs$nstr]))
    beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
  }
  rownames(pars) <- levels(as.factor(obs$strata))
  dist <- eval(parse(text = dist))
  cumhaz <- NULL
  if (frailty != "none") {
    cumhaz <- matrix(unlist(sapply(levels(as.factor(obs$strata)), 
                                   function(x) {
                                     t(cbind(dist(pars[x, ], obs$time[obs$strata == x], what = "H") * 
                                               exp(as.matrix(obs$x)[obs$strata == x, -1, drop = FALSE] %*% as.matrix(beta)), 
                                             obs$cluster[obs$strata == x]))
                                   })), ncol = 2, byrow = TRUE)
    cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                        FUN = sum)[, 2, drop = FALSE]
    if (!is.null(obs$trunc)) {
      cumhazT <- matrix(unlist(sapply(levels(as.factor(obs$strata)), 
                                      function(x) {
                                        t(cbind(dist(pars[x, ], obs$trunc[obs$strata == 
                                                                            x], what = "H") * exp(as.matrix(obs$x)[obs$strata == 
                                                                                                                     x, -1, drop = FALSE] %*% as.matrix(beta)), 
                                                obs$cluster[obs$strata == x]))
                                      })), ncol = 2, byrow = TRUE)
      cumhazT <- aggregate(cumhazT[, 1], by = list(cumhazT[, 
                                                           2]), FUN = sum)[, 2, drop = FALSE]
    }
  }  else {
    cumhaz <- sum(apply(cbind(rownames(pars), pars), 1, function(x) {
      sum(dist(as.numeric(x[-1]), obs$time[obs$strata == 
                                             x[1]], what = "H") * exp(as.matrix(obs$x[obs$strata == 
                                                                                        x[1], -1, drop = FALSE]) %*% as.matrix(beta)))
    }))
    if (!is.null(obs$trunc)) {
      cumhazT <- sum(apply(cbind(rownames(pars), pars), 
                           1, function(x) {
                             sum(dist(as.numeric(x[-1]), obs$trunc[obs$strata == 
                                                                     x[1]], what = "H") * exp(as.matrix(obs$x[obs$strata == 
                                                                                                                x[1], -1, drop = FALSE]) %*% as.matrix(beta)))
                           }))
    }
  }
  loghaz <- NULL
  if (frailty != "none") {
    loghaz <- matrix(unlist(sapply(levels(as.factor(obs$strata)), 
                                   function(x) {
                                     t(cbind(obs$event[obs$strata == x] * (dist(pars[x, 
                                                                                     ], obs$time[obs$strata == x], what = "lh") + 
                                                                             as.matrix(obs$x)[obs$strata == x, -1, drop = FALSE] %*% 
                                                                             as.matrix(beta)), obs$cluster[obs$strata == 
                                                                                                             x]))
                                   })), ncol = 2, byrow = TRUE)
    loghaz <- aggregate(loghaz[, 1], by = list(loghaz[, 2]), 
                        FUN = sum)[, 2, drop = FALSE]
  } else {
    loghaz <- sum(apply(cbind(rownames(pars), pars), 1, function(x) {
      sum(obs$event[obs$strata == x[1]] * (dist(as.numeric(x[-1]), 
                                                obs$time[obs$strata == x[1]], what = "lh") + 
                                             as.matrix(obs$x[obs$strata == x[1], -1, drop = FALSE]) %*% 
                                             as.matrix(beta)))
    }))
  }
  logSurv <- NULL
  if (frailty == "gamma") {
    logSurv <- mapply(parfm:::fr.gamma, k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), what = "logLT")
  }
  else if (frailty == "ingau") {
    logSurv <- mapply(parfm:::fr.ingau, k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      theta = rep(theta, obs$ncl), what = "logLT")
  }
  else if (frailty == "possta") {
    logSurv <- sapply(1:obs$ncl, function(x) parfm:::fr.possta(k = obs$di[x], 
                                                       s = as.numeric(cumhaz[[1]])[x], nu = nu, Omega = Omega, 
                                                       what = "logLT", correct = correct))
  }
  else if (frailty == "lognormal") {
    logSurv <- mapply(parfm:::fr.lognormal, k = obs$di, s = as.numeric(cumhaz[[1]]), 
                      sigma = rep(sigma, obs$ncl), what = "logLT")
  }
  else if (frailty == "none") {
    logSurv <- mapply(parfm:::fr.none, s = cumhaz, what = "logLT")
  }
  if (!is.null(obs$trunc)) {
    logSurvT <- NULL
    if (frailty == "gamma") {
      logSurvT <- mapply(parfm:::fr.gamma, k = 0, s = as.numeric(cumhazT[[1]]), 
                         theta = rep(theta, obs$ncl), what = "logLT")
    }
    else if (frailty == "ingau") {
      logSurvT <- mapply(parfm:::fr.ingau, k = 0, s = as.numeric(cumhazT[[1]]), 
                         theta = rep(theta, obs$ncl), what = "logLT")
    }
    else if (frailty == "possta") {
      logSurvT <- sapply(1:obs$ncl, function(x) parfm:::fr.possta(k = 0, 
                                                          s = as.numeric(cumhazT[[1]])[x], nu = nu, Omega = Omega, 
                                                          what = "logLT", correct = correct))
    }
    else if (frailty == "lognormal") {
      logSurvT <- mapply(parfm:::fr.lognormal, k = 0, s = as.numeric(cumhazT[[1]]), 
                         sigma = rep(sigma, obs$ncl), what = "logLT")
    }
    else if (frailty == "none") {
      logSurvT <- mapply(parfm:::fr.none, s = cumhazT, what = "logLT")
    }
  }
  Mloglik <- -sum(as.numeric(loghaz[[1]]) + logSurv)
  if (!is.null(obs$trunc)) {
    Mloglik <- Mloglik + sum(logSurvT)
  }
  attr(Mloglik, "cumhaz") <- as.numeric(cumhaz[[1]])
  if (!is.null(obs$trunc)) {
    attr(Mloglik, "cumhazT") <- as.numeric(cumhazT[[1]])
  }
  else {
    attr(Mloglik, "cumhazT") <- NULL
  }
  attr(Mloglik, "loghaz") <- as.numeric(loghaz[[1]])
  attr(Mloglik, "logSurv") <- (logSurv)
  if (!is.null(obs$trunc)) {
    attr(Mloglik, "logSurvT") <- (logSurvT)
  }
  Mloglik[is.na(Mloglik)]=1e100
  return(Mloglik)
}

#parfm function copied and modified from parfm package
#Modifications: 
#0) use my.optim as default
#1) adding posiblity to change optimization procedures, init: optim=my.optim 
#  all necessary parameters icncluded in optim() call
#2) automatic calcualtion of lower and upper bounds for parameters, needed by DE optimization
#   currently not used, init: do.DE.optim=FALSE
my.parfm<-function (formula, cluster = NULL, strata = NULL, data, inip = NULL, lower=NULL,upper=NULL,
          iniFpar = NULL, dist = "weibull", frailty = "none", method = "BFGS", do.DE.optim=FALSE,
          maxit = 2500, Fparscale = 1, showtime = TRUE, correct = 0, optim=my.optim)   {
  varcov='Not calculated'
  if (missing(data)) {
    data <- eval(parse(text = paste("data.frame(", paste(all.vars(formula), 
                                                         collapse = ", "), ")")))
  }
  if (!(dist %in% c("exponential", "weibull", "gompertz", "loglogistic", 
                    "lognormal"))) {
    stop("invalid baseline hazard")
  }
  if (!(frailty %in% c("none", "gamma", "ingau", "possta", 
                       "lognormal"))) {
    stop("invalid frailty distribution")
  }
  if (frailty == "none" && !is.null(cluster)) {
    warning(paste("With frailty='none' the cluster variable '", 
                  cluster, "' is not used!", sep = ""))
  }
  if (frailty == "none" && !is.null(iniFpar)) {
    warning("With frailty='none' the argument 'iniFpar' is not used!")
  }
  if (frailty == "possta") {
    if (10^correct == Inf || 10^-correct == 0) {
      stop("'correct' is too large!")
    }
    if (10^correct == 0 || 10^-correct == Inf) {
      stop("'correct' is too small!")
    }
  }
  else if (correct != 0) {
    warning(paste("'correct' has no effect when 'frailty = ", 
                  frailty, "'", sep = ""))
  }
  obsdata <- NULL
  if (length(formula[[2]]) == 3) {
    obsdata$time <- eval(formula[[2]][[2]], envir = data)
    obsdata$event <- eval(formula[[2]][[3]], envir = data)
  }
  else if (length(formula[[2]]) == 4) {
    obsdata$trunc <- eval(formula[[2]][[2]], envir = data)
    obsdata$time <- eval(formula[[2]][[3]], envir = data)
    obsdata$event <- eval(formula[[2]][[4]], envir = data)
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object", 
               "in the left-hand side of the formula object", "must be either 0 (no event) or 1 (event)."))
  }
  obsdata$x <- as.data.frame(model.matrix(formula, data = data))
  if (is.null(cluster)) {
    if (frailty != "none") {
      stop(paste("if you specify a frailty distribution,\n", 
                 "then you have to specify the cluster variable as well"))
    }
    else {
      obsdata$cluster <- rep(1, nrow(data))
    }
    obsdata$ncl <- 1
    obsdata$di <- sum(obsdata$event)
  }  else {
    if (!cluster %in% names(data)) {
      stop(paste("object '", cluster, "' not found", sep = ""))
    }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster), 
                            FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
  }
  if (is.null(strata)) {
    obsdata$strata <- rep(1, length(obsdata$time))
    obsdata$nstr <- 1
    obsdata$dq <- sum(obsdata$event)
  }  else {
    if (!strata %in% names(data)) {
      stop(paste("object '", strata, "' not found", sep = ""))
    }
    obsdata$strata <- eval(as.name(strata), envir = data)
    obsdata$nstr <- length(levels(as.factor(obsdata$strata)))
    obsdata$dq <- aggregate(obsdata$event, by = list(obsdata$strata), 
                            FUN = sum)[, , drop = FALSE]
    snames <- obsdata$dq[, 1]
    obsdata$dq <- as.vector(obsdata$dq[, 2])
    names(obsdata$dq) <- snames
  }
  
  if (!is.null(cluster) && !is.null(strata)) {
    obsdata$dqi <- xtabs(x ~ Group.1 + Group.2, data = aggregate(obsdata$event, 
                                                                 by = list(obsdata$cluster, obsdata$strata), FUN = sum))
    dimnames(obsdata$dqi) <- list(cluster = dimnames(obsdata$dqi)[[1]], 
                                  strata = dimnames(obsdata$dqi)[[2]])
  } else if (!is.null(cluster)) {
    obsdata$dqi <- obsdata$di
  } else if (!is.null(strata)) {
    obsdata$dqi <- obsdata$dq
  } else {
    obsdata$dqi <- sum(obsdata$event)
  }
  if (frailty == "none") {
    nFpar <- 0
  } else if (frailty %in% c("gamma", "ingau", "possta", "lognormal")) {
    nFpar <- 1
  }
  obsdata$nFpar <- nFpar
  if (dist == "exponential") {
    nBpar <- 1
  } else if (dist %in% c("weibull", "gompertz", "lognormal", 
                       "loglogistic")) {
    nBpar <- 2
  }
  obsdata$nBpar <- nBpar
  nRpar <- ncol(obsdata$x) - 1
  obsdata$nRpar <- nRpar
  if (!is.null(inip)) {
    if (length(inip) != nBpar * obsdata$nstr + nRpar) {
      stop(paste("number of initial parameters 'inip' must be", 
                 nBpar * obsdata$nstr + nRpar))
    }
    p.init <- inip
    if (dist %in% c("exponential", "weibull", "gompertz")) {
      if (any(p.init[1:obsdata$nstr] <= 0)) {
        stop(paste("with that baseline, the 1st parameter has to be > 0"))
      }
      p.init[1:obsdata$nstr] <- log(p.init[1:obsdata$nstr])
    }
    if (dist %in% c("weibull", "gompertz", "lognormal", "loglogistic")) {
      if (any(p.init[obsdata$nstr + 1:obsdata$nstr] <= 
                0)) {
        stop(paste("with that baseline, the 2nd parameter has to be > 0"))
      }
      p.init[obsdata$nstr + 1:obsdata$nstr] <- log(p.init[obsdata$nstr + 
                                                            1:obsdata$nstr])
    }
  } else {
    coxformula <- formula
    if (!is.null(strata)) {
      if (dist != "exponential") {
        coxformula <- eval(parse(text = paste("update(formula, .~.+strata(", 
                                              strata, "))", sep = "")))
      } else {
        coxformula <- eval(parse(text = paste("update(formula, .~.+", 
                                              strata, ")", sep = "")))
      }
    }
    shape <- 0
    d <- dist
    if (d == "exponential") {
      d <- "weibull"
      shape <- 1
    }
    coxMod <- eha:::phreg(formula = coxformula, data = data, dist = d, 
                    shape = shape, center = FALSE, control = list(maxiter = maxit))
    logshape <- as.numeric(coxMod$coef[substr(names(coxMod$coef), 5, 9) == "shape"])
    logscale <- as.numeric(coxMod$coef[substr(names(coxMod$coef), 5, 9) == "scale"])
    
    if (!is.null(strata) && dist == "exponential") {
      logscale <- logscale - c(0, coxMod$coef[substr(names(coxMod$coef), 
                                                     1, nchar(strata)) == strata])
    }
    if (dist == "exponential") {
      p.init <- -logscale
    }
    else if (dist == "weibull") {
      p.init <- c(logshape, -exp(logshape) * logscale)
    }
    else if (dist == "gompertz") {
      p.init <- c(-logscale, logshape)
    }
    else if (dist == "lognormal") {
      p.init <- c(logscale, -logshape)
    }
    else if (dist == "loglogistic") {
      p.init <- c(-exp(logshape) * logscale, logshape)
    }
    if (nRpar > 0) {
      p.init <- c(p.init, as.numeric(coxMod$coef[setdiff(names(coxMod$coef), 
                                                         c("(Intercept)", "log(scale)", "log(shape)"))]))
    }
  }
  if (frailty == "none") {
    pars <- NULL
  } else if (frailty %in% c("gamma", "ingau", "lognormal")) {
    if (is.null(iniFpar)) {
      iniFpar <- 1
    }
    else if (iniFpar <= 0) {
      stop("initial heterogeneity parameter (theta) has to be > 0")
    }
    pars <- log(iniFpar)
  }  else if (frailty == "possta") {
    if (is.null(iniFpar)) {
      iniFpar <- 0.5
    } else if (iniFpar <= 0 || iniFpar >= 1) {
      stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
    }
    pars <- log(-log(iniFpar))
  }
  pars <- c(pars, p.init)
  res <- NULL
  if ((frailty == "none") && is.null(inip)) {
    todo <- expression({
      res <- list(par = pars)
    })
    if (showtime) {
      extime <- system.time(eval(todo))[1]
    }
    else {
      eval(todo)
      extime <- NULL
    }
    it <- NULL
    lL <- coxMod$loglik[2]
  } else {
    if ((length(upper)==0) || (length(lower)==0)){
      lower=-25*abs(pars)-10
      upper=25*abs(pars)+10
    }
    
    todo <- expression({
      res <- optim(par = pars, fn = my.Mloglikelihood, method = method, upper=upper,lower=lower,
                   obs = obsdata, dist = dist, frailty = frailty, do.DE.optim=do.DE.optim,
                   correct = correct, hessian = TRUE, 
                   control = list(maxit = maxit, parscale = c(rep(Fparscale, nFpar), 
                                                              rep(1,nBpar * obsdata$nstr + nRpar))))
    })
    if (showtime) {
      extime <- system.time(eval(todo))[1]
    }
    else {
      eval(todo)
      extime <- NULL
    }
    if (res$convergence > 0) {
      warning("optimisation procedure did not converge,\n              conv = ", 
              bquote(.(res$convergence)), ": see optim() for details")
    }
    it <- res$counts[1]
    lL <- as.double(-res$value)
    if (frailty == "possta") {
      lL <- lL + correct * log(10) * obsdata$ncl
    }
  }
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(res$par[1:nFpar])
    sigma2 <- NULL
    nu <- NULL
  }
  else if (frailty == "lognormal") {
    theta <- NULL
    sigma2 <- exp(res$par[1:nFpar])
    nu <- NULL
  }
  else if (frailty == "possta") {
    theta <- NULL
    sigma2 <- NULL
    nu <- exp(-exp(res$par[1:nFpar]))
  }
  else if (frailty == "none") {
    theta <- NULL
    sigma2 <- NULL
    nu <- NULL
  }
  if (dist == "exponential") {
    lambda <- exp(res$par[nFpar + 1:obsdata$nstr])
    ESTIMATE <- c(lambda = lambda)
  }
  else if (dist %in% c("weibull")) {
    rho <- exp(res$par[nFpar + 1:obsdata$nstr])
    lambda <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(rho = rho, lambda = lambda)
  }
  else if (dist == "gompertz") {
    gamma <- exp(res$par[nFpar + 1:obsdata$nstr])
    lambda <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(gamma = gamma, lambda = lambda)
  }
  else if (dist == "lognormal") {
    mu <- res$par[nFpar + 1:obsdata$nstr]
    sigma <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(mu = mu, sigma = sigma)
  }
  else if (dist == "loglogistic") {
    alpha <- res$par[nFpar + 1:obsdata$nstr]
    kappa <- exp(res$par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(alpha = alpha, kappa = kappa)
  }
  if (nRpar == 0) {
    beta <- NULL
  }
  else {
    beta <- res$par[-(1:(nFpar + nBpar * obsdata$nstr))]
    names(beta) <- paste("beta", names(obsdata$x), sep = ".")[-1]
  }
  ESTIMATE <- c(theta = theta, sigma2 = sigma2, nu = nu, ESTIMATE, 
                beta = beta)
  if ((frailty == "none") && is.null(inip)) {
    var <- coxMod$var
    if (nRpar == 0) {
      seBeta <- NULL
    } else {
      seBeta <- sqrt(diag(var)[setdiff(names(coxMod$coef), 
                                       c("(Intercept)", "log(scale)", "log(shape)"))])
      PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                2 * pnorm(q = -abs(beta/seBeta)))
    }
    if (obsdata$nstr == 1) {
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(-x1), mean = logscale[x], 
                      cov = var["log(scale)", "log(scale)"], ses = TRUE)
        })
        STDERR <- c(seLambda = seLambda)
      }
      else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(x1), mean = logshape[x], 
                      cov = var["log(shape)", "log(shape)"], ses = TRUE)
        })
        seLambda <- deltamethod(g = ~exp(-exp(x2) * x1), 
                                mean = c(logscale, logshape), cov = var[c("log(scale)", 
                                                                          "log(shape)"), c("log(scale)", "log(shape)")], 
                                ses = TRUE)
        STDERR <- c(seRho = seRho, seLambda = seLambda)
      }
      else if (dist == "gompertz") {
        seGamma <- deltamethod(g = ~exp(-x1), mean = logscale, 
                               cov = var["log(scale)", "log(scale)"], ses = TRUE)
        seLambda <- deltamethod(g = ~exp(x1), mean = logshape, 
                                cov = var["log(shape)", "log(shape)"], ses = TRUE)
        STDERR <- c(seGamma = seGamma, seLambda = seLambda)
      }
      else if (dist == "lognormal") {
        seMu <- sqrt(var["log(scale)", "log(scale)"])
        seSigma <- deltamethod(g = ~exp(-x1), mean = logshape, 
                               cov = var["log(shape)", "log(shape)"], ses = TRUE)
        STDERR <- c(seMu = seMu, seSigma = seSigma)
      }
      else if (dist == "loglogistic") {
        seAlpha <- deltamethod(g = ~-exp(x2) * x1, mean = c(logscale, 
                                                            logshape), cov = var[c("log(scale)", "log(shape)"), 
                                                                                 c("log(scale)", "log(shape)")], ses = TRUE)
        seKappa <- deltamethod(g = ~exp(x1), mean = logshape, 
                               cov = var["log(shape)", "log(shape)"], ses = TRUE)
        STDERR <- c(seAlpha = seAlpha, seKappa = seKappa)
      }
    } else {
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(-x1), mean = logscale[x], 
                      cov = var[paste("log(scale)", x, sep = ":"), 
                                paste("log(scale)", x, sep = ":")], ses = TRUE)
        })
        STDERR <- c(seLambda = seLambda)
      }
      else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(x1), mean = logshape[x], 
                      cov = var[paste("log(shape)", x, sep = ":"), 
                                paste("log(shape)", x, sep = ":")], ses = TRUE)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(-exp(x2) * x1), mean = c(logscale[x], 
                                                        logshape[x]), cov = var[paste(c("log(scale)", 
                                                                                        "log(shape)"), x, sep = ":"), paste(c("log(scale)", 
                                                                                                                              "log(shape)"), x, sep = ":")], ses = TRUE)
        })
        STDERR <- c(seRho = seRho, seLambda = seLambda)
      }
      else if (dist == "gompertz") {
        seGamma <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(-x1), mean = logscale[x], 
                      cov = var[paste("log(scale)", x, sep = ":"), 
                                paste("log(scale)", x, sep = ":")], ses = TRUE)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(x1), mean = logshape, 
                      cov = var[paste("log(shape)", x, sep = ":"), 
                                paste("log(shape)", x, sep = ":")], ses = TRUE)
        })
        STDERR <- c(seGamma = seGamma, seLambda = seLambda)
      }
      else if (dist == "lognormal") {
        seMu <- sqrt(diag(var[substr(rownames(var), 5, 
                                     9) == "scale", substr(rownames(var), 5, 9) == 
                                "scale"]))
        seSigma <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(-x1), mean = logshape, 
                      cov = var[paste("log(shape)", x, sep = ":"), 
                                paste("log(shape)", x, sep = ":")], ses = TRUE)
        })
        STDERR <- c(seMu = seMu, seSigma = seSigma)
      }
      else if (dist == "loglogistic") {
        seAlpha <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~-exp(x2) * x1, mean = c(logscale, 
                                                   logshape), cov = var[paste(c("log(scale)", 
                                                                                "log(shape)"), x, sep = ":"), paste(c("log(scale)", 
                                                                                                                      "log(shape)"), x, sep = ":")], ses = TRUE)
        })
        seKappa <- sapply(1:obsdata$nstr, function(x) {
          deltamethod(g = ~exp(x1), mean = logshape, 
                      cov = var[paste("log(shape)", x, sep = ":"), 
                                paste("log(shape)", x, sep = ":")], ses = TRUE)
        })
        STDERR <- c(seAlpha = seAlpha, seKappa = seKappa)
      }
    }
    STDERR <- c(STDERR, se.beta = seBeta)
  }  else {
    varcov = solve(res$hessian)
    var <- try(diag(varcov), silent = TRUE)
    if (class(var) == "try-error") {
      warning(var[1])
      STDERR <- rep(NA, nFpar + nBpar * obsdata$nstr + 
                      nRpar)
      PVAL <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
    } else {
      if (any(var <= 0)) {
        warning(paste("negative variances have been replaced by NAs\n", 
                      "Please, try other initial values", "or another optimisation method"))
      }
      if (frailty %in% c("gamma", "ingau")) {
        seTheta <- sapply(1:nFpar, function(x) {
          ifelse(var[x] > 0, sqrt(var[x] * theta[x]^2), 
                 NA)
        })
        seSigma2 <- seNu <- NULL
      }
      else if (frailty == "lognormal") {
        seSigma2 <- sapply(1:nFpar, function(x) {
          ifelse(var[x] > 0, sqrt(var[x] * sigma2[x]^2), 
                 NA)
        })
        seTheta <- seNu <- NULL
      }
      else if (frailty == "possta") {
        seNu <- sapply(1:nFpar, function(x) {
          ifelse(var[x] > 0, sqrt(var[x] * (nu * log(nu))^2), 
                 NA)
        })
        seTheta <- seSigma2 <- NULL
      }
      if (dist == "exponential") {
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + 
                                                x] * lambda[x]^2), NA)
        })
        STDERR <- c(seLambda = seLambda)
      }
      else if (dist %in% c("weibull")) {
        seRho <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + 
                                                x] * rho[x]^2), NA)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                               obsdata$nstr + x] * lambda[x]^2), NA)
        })
        STDERR <- c(seRho = seRho, seLambda = seLambda)
      }
      else if (dist == "gompertz") {
        seGamma <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + 
                                                x] * gamma[x]^2), NA)
        })
        seLambda <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                               obsdata$nstr + x] * lambda[x]^2), NA)
        })
        STDERR <- c(seGamma = seGamma, seLambda = seLambda)
      }
      else if (dist == "lognormal") {
        seMu <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + 
                                                x]), NA)
        })
        seSigma <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                               obsdata$nstr + x] * sigma[x]^2), NA)
        })
        STDERR <- c(seMu = seMu, seSigma = seSigma)
      }
      else if (dist == "loglogistic") {
        seAlpha <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + 
                                                x]), NA)
        })
        seKappa <- sapply(1:obsdata$nstr, function(x) {
          ifelse(var[nFpar + obsdata$nstr + x] > 0, sqrt(var[nFpar + 
                                                               obsdata$nstr + x] * kappa[x]^2), NA)
        })
        STDERR <- c(seAlpha = seAlpha, seKappa = seKappa)
      }
      if (nRpar == 0) {
        seBeta <- NULL
      }
      else {
        seBeta <- numeric(nRpar)
        varBeta <- var[-(1:(nFpar + nBpar * obsdata$nstr))]
        for (i in 1:nRpar) {
          seBeta[i] <- ifelse(varBeta[i] > 0, sqrt(varBeta[i]), 
                              NA)
        }
        PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                  2 * pnorm(q = -abs(beta/seBeta)))
      }
      STDERR <- c(STDERR, se.beta = seBeta)
      if (frailty != "none") {
        STDERR <- c(se.theta = seTheta, se.sigma2 = seSigma2, 
                    se.nu = seNu, STDERR)
      }
    }
  }
  resmodel <- cbind(ESTIMATE = ESTIMATE, SE = STDERR)
  rownames(resmodel) <- gsub("beta.", "", rownames(resmodel))
  if (nRpar > 0) {
    resmodel <- cbind(resmodel, `p-val` = PVAL)
  }
  class(resmodel) <- c("parfm", class(resmodel))
  Call <- match.call()
  if (!match("formula", names(Call), nomatch = 0)) 
    stop("A formula argument is required")
  Terms <- terms(formula, data = data)
  attributes(resmodel) <- c(attributes(resmodel), list(call = Call, 
                                                       convergence = res$convergence, 
                                                       it = it, 
                                                       extime = extime, 
                                                       nobs = nrow(data), 
                                                       shared = (nrow(data) > obsdata$ncl), 
                                                       loglik = lL, 
                                                       dist = dist, 
                                                       varcov = varcov, #varince covaariance matrix added
                                                       cumhaz = attributes(my.Mloglikelihood(p = res$par,
                                                                           obs = obsdata, dist = dist, frailty = frailty, correct = correct))$cumhaz, 
                                                       cumhazT = attributes(my.Mloglikelihood(p = res$par, obs = obsdata, 
                                                                                           dist = dist, frailty = frailty, correct = correct))$cumhazT, 
                                                       di = obsdata$di, dq = obsdata$dq, dqi = obsdata$dqi, 
                                                       frailty = frailty, clustname = cluster, stratname = strata, 
                                                       correct = correct, formula = as.character(Call[match("formula", 
                                                                                                            names(Call), nomatch = 0)]), terms = attr(Terms, 
                                                                                                                                                      "term.labels")))
  if (frailty != "none") {
    names(attr(resmodel, "cumhaz")) <- names(attr(resmodel, 
                                                  "di")) <- unique(obsdata$cluster)
  }
  if (showtime) {
    cat("\nExecution time:", extime, "second(s) \n")
  }
  return(resmodel)
}


#Improved results display
my.anova.parfm<-function (object, ...) {
  if (length(list(object, ...)) > 1) 
    return(anova.parfmlist(list(object, ...)))
  else {
    if (!inherits(object, "parfm")) 
      stop(paste("The argument must be a parametric frailty model,", 
                 "object of class 'parfm'"))
    termlist <- attr(object, "terms")
    nmodels <- length(termlist)
    df <- integer(nmodels + 1)
    loglik <- double(nmodels + 1)
    df[1] <- 0
    mycall <- attributes(object)$call
    mycall[["showtime"]] <- FALSE
    cat('Intercept model: ')
    mycall[[2]] <- update.formula(mycall[[2]], . ~ 1)
    loglik[1] <- attr(eval(mycall), "loglik")
    cat('\n')
    for (i in 1:nmodels) {
      df[1 + i] <- df[i] + 1
      cat(paste(paste(as.character(mycall[[2]])[c(2,1,3)])),'+',termlist[i])
      Forma=parse(text = paste("update.formula(mycall[[2]], . ~ . +",termlist[i], ")"))
      mycall[[2]] <- eval(Forma)
      loglik[1 + i] <- attr(eval(mycall), "loglik")
      cat('\n')
    }
    #diff means sequential
    table <- data.frame(loglik = loglik, Chisq = c(NA, 2 * diff(loglik)), Df = c(NA, diff(df)))
    table[["Pr(>|Chi|)"]] <- 1 - pchisq(table$Chisq, table$Df)
    row.names(table) <- c("NULL", termlist)
    title <- paste("Sequential Analysis of Deviance Table\n",
                   '\nRemember that in R the order in which you enter the variables in the formula affects the anova results. By default, anova fits the variables sequentially ("type I sum of squares"), while at the same time respecting hierarchy. Some other statistics programs such as SAS and JMP use marginal fitting of terms instead ("type III sum of squares") instead. ',
                   'for this kind of models drop1 is more preferable\n',
                   'see: https://stats.stackexchange.com/questions/37805/why-do-anova-and-drop1-provided-different-answers-for-glmms\n',
                   'you can use my.drop1.parfm() function\n',
                   "\nSpecifying a single object gives a sequential analysis of deviance table for that fit. That is, the reductions in the model log-likelihood as each term of the formula is added in turn are given in as the rows of a table, plus the log-likelihoods themselves.\n",
                   "\nIf more than one object is specified, the table has a row for the degrees of freedom and loglikelihood for each model. For all but the first model, the change in degrees of freedom and loglik is also given. (This only make statistical sense if the models are nested.) It is conventional to list the models from smallest to largest, but this is up to the user.\n",
                   "Parametric frailty model: response is ", 
                   deparse(attr(object, "call")[match("formula", names(attr(object, 
                                                                            "call")))][[1]][[2]]), "\nTerms added sequentially (first to last)\n", 
                   sep = "")
    structure(table, heading = title, class = c("anova", 
                                                "data.frame"))
  }
}

#UNDER CONSTRUCTION
# my.drop1.parfm<-function (object) {
#     if (!inherits(object, "parfm")) 
#       stop(paste("The argument must be a parametric frailty model,", 
#                  "object of class 'parfm'"))
#     termlist <- attr(object, "terms")
#     nmodels <- length(termlist)
#     df <- integer(nmodels)
#     loglik <- double(nmodels)
#     
#     fullmycall <- attributes(object)$call
#     fullmycall[["showtime"]] <- FALSE
#     mycall=fullmycall 
#     fullloglik <- attr(object, "loglik")
#     fulldf<- length(termlist)
#     cat('\n')
#     for (i in  1:nmodels) {
#       fullmyform=update.formula(fullmycall, ~.)
#       Forma1=update.formula(paste(paste(paste(as.character(fullmyform)[c(2,1,3)]),collapse=''),'-',termlist[i]),~.)
#       cat(as.character(Forma1[c(2,1,3)]))
#       Forma=parse(text = paste("update.formula(fullmycall[[2]], . ~ . -",termlist[i], ")"))
#       mycall[[2]] <- eval(Forma)
#       obj2=eval(mycall)
#       cat('\nPerformed function call: '); print(attributes(obj2)$call)
#       loglik[i] <- attr(obj2, "loglik")
#       cat('\n')
#     }
#     table <- data.frame(loglik = loglik, ratio=(fullloglik-loglik),Chisq = c( 2 * (fullloglik-loglik)), Df = 1)
#     table[["Pr(>|Chi|)"]] <- round(1 - pchisq(table$Chisq, table$Df),6)
#     row.names(table) <- c(termlist)
#     title <- paste("Analysis of Deviance Table by dropping terms\n",
#                    "Marginality assumption is not fullfilled (just drop of single terms)\n",
#                    "Parametric frailty model: response is ", 
#                    deparse(attr(object, "call")
#                            [match("formula",names(attr(object,"call")))][[1]][[2]]), 
#                    "\nTerms added sequentially (first to last)\n", 
#                    sep = "")
#     structure(table, heading = title, class = c("anova", 
#                                                 "data.frame"))
#   
# }
