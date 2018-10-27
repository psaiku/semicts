library(truncnorm)
library(scatterplot3d)
library(truncdist)
# Semi-continuous
#
# This package contains functions useful to work with semi-continuous data
# with point mass at 0 and continuous on the positive real line.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Generates a random sample from a semi-continuous distribution. Currently,
#' truncated normal (truncnorm), log-normal (lnorm), and gamma distributions are
#' supported for the continuous part of the distribution.
#'
#' @param n Number of random variables to generate
#' @param pzero Point mass at 0
#' @param cts.density Name of a continuous density with support on the positive real line. Supported values: truncnorm (default), lnorm, and gamma
#' @param  cts.params An array containing the parameters for cts.density (default: c(1,1) for mean, standard deviation of truncated normal). For
#' log-normal, it should be an array containing meanlog, and sdlog of the distribution. For gamma,
#' an array of shape, and rate values must be supplied.
#' @return An array of semi-continuous random variables.
#' @examples
#' rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' @import truncnorm
rsemicts.old <- function(n, pzero = 0.5, cts.density = "truncnorm", cts.params = c(1, 1)) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("truncnorm package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(pzero==0) {
    stop("Pr(Y=0) cannot be 0!")
  }
  u <- runif(n)
  x <- u
  x <- (u>pzero)*u

  idx <- which(x>0)

  if(length(idx)>0) {
    switch(cts.density,
           truncnorm={
             x[idx] <- qtruncnorm((u[idx]-pzero)/(1-pzero), a=0, b=Inf, mean=cts.params[1], sd=cts.params[2])
            },
           lnorm={
             x[idx] <- qlnorm((u[idx]-pzero)/(1-pzero), meanlog=cts.params[1], sdlog=cts.params[2])
           },
           gamma={
             x[idx] <- qgamma((u[idx]-pzero)/(1-pzero), shape=cts.params[1], rate=cts.params[2])
           })
  }
  class(x) <- "semicts"
  return(x)
}

#' Generates a random sample from a semi-continuous distribution. Currently,
#' truncated normal (truncnorm), log-normal (lnorm), and gamma distributions are
#' supported for the continuous part of the distribution.
#'
#' @param n Number of random variables to generate
#' @param pzero Point mass at 0
#' @param cts.density Name of a continuous density with support on the positive real line. Supported values: truncnorm (default), lnorm, and gamma
#' @param  cts.param An array containing the parameters for cts.density (default: c(1,1) for mean, standard deviation of truncated normal). For
#' log-normal, it should be an array containing meanlog, and sdlog of the distribution. For gamma,
#' an array of shape, and rate values must be supplied.
#' @return An array of semi-continuous random variables.
#' @examples
#' rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' @export
#' @import truncnorm
rsemicts <- function(n, pzero = 0.5, r.func = NA, cts.density = "truncnorm",
                      cts.param=list(), left.args=list(), right.args=list()) {
  if(pzero==0 & is.na(r.func)) {
    stop("Pr(Y=0) cannot be 0!")
  }
  #if(!missing(FUN)) {
  if(!is.na(r.func)) {
    #pzero = FUN(cts.param)
    r.func <- match.fun(r.func)
    #pzero = do.call(r.func, as.list(cts.param))
    pzero = r.func(cts.param)
  }
  u <- runif(n)
  x <- u
  x <- (u>pzero)*u

  idx <- which(x>0)

  fapp <- function(t, f, argss) {
    return(do.call(f, append((u[t]-pzero)/(1-pzero), argss)))
  }

  if(length(idx)>0) {
    f <- match.fun(paste0("q", cts.density))
    x[idx] <- apply(as.matrix(idx), 1, fapp, f, append(append(left.args, cts.param), right.args))
  }
  class(x) <- "semicts"
  return(x)
}

#' @export
dsemicts <- function(y, pzero = 0, r.func=NA, cts.density="truncnorm", cts.param,
                    left.args=c(), right.args=c()) {
  cts.param = as.numeric(cts.param)
  left.args = as.numeric(left.args)
  right.args = as.numeric(right.args)
  f <- match.fun(paste0("d", cts.density))
  if(pzero==0 & is.na(r.func)) {
    stop("Pr(Y=0) cannot be 0!")
  }

  #if(!missing(FUN)) {
  if(!is.na(r.func)) {
    #pzero = FUN(cts.param)
    pzero = do.call(r.func, as.list(cts.param))
  }

  fapp <- function(y1, pzero, f, argss) {
    dlt <- (y1==0)
    return(ifelse(dlt, log(pzero), log((1-pzero)*do.call(f, as.list(c(y1*(1-dlt), argss))))))
  }

  ll <- sum(apply(as.matrix(y), 1, fapp, pzero, f, c(left.args, cts.param, right.args)))

  return(ll)
}

#' @export
mle.semicts <- function(y, pzero = 0, r.func=NA, cts.density="truncnorm", cts.param,
                        left.args=c(), right.args=c()) {
  ll <- function(param) {
    return(-1*dsemicts(y, pzero, r.func, cts.density, param, left.args, right.args))
  }

  gamma.hat <- 0.0
  cts.param.hat <- rep(0, length(cts.param))

  # Two-part estimation
  if(is.na(r.func) & pzero > 0) {
    gamma.hat <- sum(y==0)/length(y)
    y <- y[which(y>0)]
    fit0 <- optim(cts.param, ll)
    cts.param.hat <- fit0$par
  } else {
    fit0 <- optim(cts.param, ll)
    cts.param.hat <- fit0$par
    gamma.hat <- FUN(cts.param.hat)
  }

  return(list(gamma.hat=gamma.hat, cts.param.hat = cts.param.hat))
}

#' @export
plotll <- function(y, pzero = 0, r.func, cts.density="truncnorm", cts.param,
                    left.args=c(), right.args=c(), par1, par2, xlab, ylab,
                    zlab, main) {
  n.grid <- length(par1)
  ff <- function(par) {
    return(dsemicts(y, pzero, r.func, cts.density, par, left.args, right.args))
  }
  if(is.na(par2)) {
    ll <- apply(as.matrix(par1), 1, ff)
    plot(pars, ll, type="l", xlab=xlab, ylab=ylab, main=main)
  } else {
    d <- expand.grid(par1, par2)
    ll <- apply(d, 1, ff)
    #persp(x = as.numeric(d[,1]), y = as.numeric(d[,2]), z = ll, main=main, xlab=xlab, ylab=ylab,
     #      zlab=zlab)
    plot_ly(x = as.numeric(d[,1]), y = as.numeric(d[,2]), z = ll, main=main, xlab=xlab,
            ylab=ylab, zlab=zlab)
  }
}

#' Generate a random sample from a latent normal distribution (tobit)
#'
#' @param n Number of random variables to generate
#' @param mean mean of the latent normal distribution
#' @param  sd standard deviation of the latent normal distribution
#' @return An array of semi-continuous random variables.
#' @examples
#' rtobit(100)
#' @export
rtobit <- function(n, mean=1, sd=1) {
  pzero <- 1-pnorm(mean/sd)
  return(rsemicts(n, pzero=pzero, cts.density="truncnorm", cts.param=c(mean,sd)))
}

#' Generates a histogram for the semi-continuous data provided
#'
#' @param obj an array of semi-continuous data
#' @param xlab label for x-axis (has a default)
#' @param ylab label for y-axis (has a default)
#' @param main title (has a default)
#' @param cols two colors as an array: one for the proportion of zeroes, and the other for the continuous data (has a default)
#' @param legends names for the legend as an array of two elements (has a default)
#' @return histogram object
#' @examples
#' x <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' hist(x)
#' @export
hist1.semicts <- function(obj, xlab="x", ylab="Density", main="",
                        cols=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)),
                        legends=c("Intensity", "Proportion of Zeros")) {
  par(mar=c(5, 4, 4, 6) + 0.1)
  obs.zero.prop <- length(obj[obj == 0]) / length(obj)
  obs.hist <- hist.default(obj[which(obj>0)], breaks = seq(0, max(obj)+1, by = 0.5), plot = FALSE)
  obs.hist$counts <- obs.hist$counts / sum(obs.hist$counts)
  plot(obs.hist, col = cols[1], xlab="", ylab="", main = main, axes=FALSE)
  mtext(ylab, side=2, line=2.5)
  axis(2, ylim=c(0,1), col="black", las=1)
  par(new=TRUE)
  plot(c(max(obj)+1), c(obs.zero.prop), xlab="", ylab="", ylim=c(0,1), axes=FALSE)
  mtext("Point mass", side=4, col="black", line=2.5)
  axis(4, ylim=c(0,1), col="black", col.axis="black", las=1)
  points(x = c(max(obj)+1), y = c(obs.zero.prop), col = cols[2], pch = 19, cex = 1.5)
  axis(1, pretty(range(obj),10))
  mtext("x", side=1, col="black", line=2.5)
  if(!is.null(legends)) {
    legend("topright", legends, cex = 1, lwd = c(10, NA), col = cols, pch = c(NA, 19))
  }
}

#' Generates a histogram for the semi-continuous data provided
#'
#' @param obj an array of semi-continuous data
#' @param xlab label for x-axis (has a default)
#' @param ylab label for y-axis (has a default)
#' @param main title (has a default)
#' @param cols two colors as an array: one for the proportion of zeroes, and the other for the continuous data (has a default)
#' @param legends names for the legend as an array of two elements (has a default)
#' @return histogram object
#' @examples
#' x <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.params=c(1,1))
#' hist(x)
hist.semicts <- function(obj, xlab="X", ylab="Density", main="",
                         cols=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), legends=c("Intensity", "Proportion of Zeros")) {
  obs.zero.prop <- length(obj[obj == 0]) / length(obj)
  obs.hist <- hist.default(obj, breaks = seq(0, max(obj)+1, by = 0.5), plot = FALSE)
  obs.hist$counts <- obs.hist$counts / sum(obs.hist$counts)
  #plot(obs.hist, col = rgb(0, 0, 1, 0.5), xlab = "Precipitation (mm/day)", ylab = "Density", main = "Observed daily precipitation - 1949-2000")
  plot(obs.hist, col = cols[1], xlab = xlab, ylab = ylab, main = main)
  points(x = c(0), y = c(obs.zero.prop), col = cols[2], pch = 19, cex = 1.5)
  legend("topright", legends, cex = 1.5, lwd = c(10, NA), col = cols, pch = c(NA, 19))
}


#' Prints the proportion of zeroes, and summary of the positive data in the
#' semicts object supplied.
#'
#' @param obj A semicts object (for ex. returned from the rsemicts function)
#' @examples
#' x <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' summary(x)
#' @export
summary.semicts <- function(obj,...) {
  propz <- length(which(obj==0))/length(obj)
  writeLines(sprintf("Proportion of zeroes: %s", paste(round(100*propz, 2), "%", sep="")))
  summary.default(obj[which(obj>0)])
}

#' Returns MPSE, MAD, average of positive predictions when the true value is zero
#' (and vice-versa), proportion of matched zeroes, and proportion of matched positives.
#'
#' @param y_pred A semicts object (for ex. returned from the rsemicts function)
#' @param y A semicts object
#' @examples
#' x <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' y <- rsemicts(100, pzero=0.4, cts.density="lnorm", cts.param=c(1,1))
#' pred.perf(x, y)
#' @export
pred.perf <- function(y_pred, y) {
  n_11 <- ifelse(y>0 & y_pred >0, 1, 0)
  n_01 <- ifelse(y==0 & y_pred>0, 1, 0)
  n_10 <- ifelse(y>0 & y_pred==0, 1, 0)
  n_00 <- ifelse(y==0 & y_pred==0, 1, 0)
  n_z <- sum(y==0)
  n <- length(y)
  mse11 <- mean((y-y_pred)^2*n_11)
  lad11 <- mean(abs(y-y_pred)*n_11)
  m01 <- mean(y_pred*n_01)
  m10 <- mean(y*n_10)
  p.dd <- ifelse(n_z==0, 1, sum(n_00)/n_z)
  p.rr <- sum(n_11)/(n-n_z)
  return(list("mse.rr"=mse11, "lad.rr"=lad11, "mse.dr"=m01, "mse.rd"=m10, "p.dd"=p.dd, "p.rr"=p.rr))
}
