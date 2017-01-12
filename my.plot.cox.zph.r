#labels correction for survival:::plot.cox.zph
my.plot.cox.zph<-function (x, resid = TRUE, se = TRUE, df = 4, nsmo = 40, var, xlab = "Age",ylab=NULL,...) {
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (qmat$rank < df) 
    stop("Spline fit is singular, try a smaller degrees of freedom")
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
  }
  if (length(ylab)==0) ylab <- paste("Beta(t) for", dimnames(yy)[[2]]) else ylab <- rep(ylab,nvar)
  if (missing(var)) 
    var <- 1:nvar
  else {
    if (is.character(var)) 
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) < 
          1) 
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx), 
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  for (i in var) {
    y <- yy[, i]
    yhat <- pmat %*% qr.coef(qmat, y)
    if (resid) 
      yr <- range(yhat, y)
    else yr <- range(yhat)
    if (se) {
      temp <- 2 * sqrt(x$var[i, i] * seval)
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") 
      plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i], 
           ...)
    else if (x$transform == "log") 
      plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i], 
           log = "x", ...)
    else {
      plot(range(xx), yr, type = "n", xlab = xlab, ylab = ylab[i], 
           axes = FALSE, ...)
      axis(1, xaxisval, xaxislab)
      axis(2)
      box()
    }
    if (resid) 
      points(xx, y)
    lines(pred.x, yhat)
    if (se) {
      lines(pred.x, yup, lty = 2)
      lines(pred.x, ylow, lty = 2)
    }
  }
}