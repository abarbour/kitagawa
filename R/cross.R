#' Calculate the cross-spectrum of two timeseries
#' @param x numeric; timeseries
#' @param y numeric; timeseries. if missing, assumed to be column no. 2 in \code{x}
#' @param k integer; the number of sine tapers, unless this is \code{NULL}; in the latter case
#' a Welch-based spectrum is calculated rather than a multitaper spectrum. There are distinct
#' advantages and disadvantages to either of these.
#' @param samp numeric; the sampling rate (e.g., \code{\link{deltat}}) of the data; must be the same for \code{x} and \code{y}
#' @param q numeric; the probability quantile [0,1] to calculate coherence significance levels; if missing, a 
#' pre-specified sequence is included. This is will be ignored for Welch-based spectra (see \code{k}).
#' @param verbose logical; should messages be printed?
#' @param ... additional arguments to \code{\link[sapa]{SDF}}
#' 
#' @export
#' 
#' @examples
#' require(stats)
#' 
#' n <- 1000
#' ramp <- seq_len(n)
#' parab <- ramp^2
#' 
#' set.seed(1255)
#' X <- ts(rnorm(n) + ramp/2)
#' Y <- ts(rnorm(n) + ramp/10 + parab/100)
#' 
#' # Calculate the cross spectrum, multitaper if K is not NULL
#' csd <- cross_spectrum(X, Y, k=20)
#' 
#' with(csd,{
#'   x <- Frequency
#'   px <- c(x, rev(x))
#'   
#'   y1 <- y1s <- Coherence
#'   cpr <- attr(csd, 'mtcsd')[['Coh.probs']]
#'   sig99 <- max(cpr$coh.sig) # 99% confidence
#'   insig <- y1 < sig99
#'   y1s[insig] <- NA
#'   plot(x, y1, col=NA, main='Coherence')
#'   lines(x, y1, col='grey')
#'   lines(x, y1s)
#'   
#'   y2 <- y2s <- Admittance
#'   y2e <- 2*Admittance.stderr
#'   y2s[insig] <- NA
#'   py <- c(y2 + y2e, rev(y2 - y2e))
#'   plot(x, y2, col=NA, ylim=range(py), main='Admittance')
#'   polygon(px, py, col='lightgrey', border=NA)
#'   lines(x, y2, lty=3)
#'   lines(x, y2s)
#'   
#'   y3 <- y3s <- Phase * 180/pi
#'   y3s[insig] <- NA
#'   plot(x, y3, col=NA, ylim=90*c(-1,1), main='Phase')
#'   lines(x, y3, col='grey')
#'   lines(x, y3s)
#' })
#' 
#' # turn off k to get Welch's Overlapping
#' wcsd <- cross_spectrum(X, Y, k=NULL)
#' 
cross_spectrum <- function(x, ...) UseMethod('cross_spectrum')

#' @rdname cross_spectrum
#' @export
cross_spectrum.mts <- function(x, ...){
    xsamp <- stats::deltat(x)
    cross_spectrum.default(x=unclass(x), samp=xsamp, ...)
}

#' @rdname cross_spectrum
#' @export
cross_spectrum.default <- function(x, y, k=10, samp=1, q, verbose=TRUE, ...){
  
  XY <- if (missing(y)){
    as.matrix(x)
  } else {
    cbind(as.double(x), as.double(y))
  }
  if (any(is.na(XY))) stop("Neither series can contain NA values")
  
  # Calculate sine-mt CS
  do.mt <- !is.null(k)
  cs <- if (do.mt){
    if (verbose) message("calculating sine-multitaper spectra...")
    sapa::SDF(XY, method='multitaper', n.taper=k, sampling.interval=samp, ...)
  } else {
    if (verbose) message("calculating Welch spectra...")
    sapa::SDF(XY, method='wosa', sampling.interval=samp, ...)
  }
  csa <- attributes(cs)
  if (do.mt){
    if (!all.equal(k, csa[['n.taper']])) warning('discrepancy in no. of tapers')
  }
  
  nfrq <- which(names(csa) == 'frequency')
  mtcsa <- csa[-nfrq]
  freq <- as.double(csa[[nfrq]])
  period <- 1/freq
  
  # Spectral content (complex)
  S <- as.matrix(cs)
  colnames(S) <- attr(cs, 'labels')
  S11 <- S[,'S11']
  S12 <- S[,'S12']
  S22 <- S[,'S22']
  
  # Coherence
  Coh <- abs(Mod(S12)^2 / (S11 * S22))
  
  # Admittance (gain)
  G <- abs(sqrt(Coh * S22 / S11))
  G.err <- sqrt((1 - Coh) / ifelse(do.mt,k,1))
  
  # Phase
  Phi <- atan2(x = Re(S12), y = Im(S12))
  Phi2 <- Arg(S12)
  stopifnot(all.equal(Phi, Phi2))
  
  # Assemble results
  csd. <- dplyr::data_frame(Frequency=freq, Period=1/freq, 
                           Coherence=Coh, 
                           Admittance=G, Admittance.stderr=G.err, 
                           Phase=Phi,
                           S11, S12, S22)

  # Calculate minimum coherence significance levels
  mtcsa[['Coh.probs']] <- if (do.mt){
    if (missing(q)) q <- c(0.01, seq(0.05, 0.95, by=0.05), 0.99)
    data.frame(q, coh.sig=stats::pf(q, 2, 4*k))
  } else {
    NA
  }
  attr(csd., 'mtcsd') <- mtcsa
  
  class(csd.) <- c(ifelse(do.mt,'mtcsd','wosacsd'), class(csd.))
  
  return(csd.)
}

# simple threshold unwrapper to remove lower wrapping in phase
unwrap.phase.lower <- function(ang, thresh=0) {
  while (any(ang <= thresh)){
    inds <- ang <= thresh
    ang[inds] <- ang[inds] + 360
  }
  return(ang)
}

#' Logarithmic smoothing with loess
#'
#' @param x numeric; the index series (cannot contain \code{NA})
#' @param y numeric; the series of values associated with x
#' @param x.is.log logical; determines whether the series in \code{x} has
#' been log-transformed already. If \code{FALSE} then \code{log10} is used.
#' @param ... additional parameters (e.g., \code{span}) passed to \code{\link[stats]{loess.smooth}}
#'
#' @seealso \code{\link[stats]{loess.smooth}} and \code{\link[stats]{approxfun}}
#' 
#' @references Barbour, A. J., and D. C. Agnew (2011), Noise levels on 
#' Plate Boundary Observatory borehole strainmeters in southern California,
#' Bulletin of the Seismological Society of America, 101(5), 2453-2466,
#' doi: 10.1785/0120110062
#' 
#' @return The result of \code{\link[stats]{loess.smooth}}
#' @export
#'
#' @examples
#' set.seed(11133)
#' n <- 101
#' lx <- seq(-1,1,length.out=n)
#' y <- rnorm(n) + cumsum(rnorm(n))
#' plot(lx, y, col='grey')
#' lines(logsmoo(lx, y, x.is.log=TRUE))
#' 
logsmoo <- function(x, y, x.is.log=FALSE, ...){
  if (any(is.na(x))) stop('x vector cannot contain NA values')
  fx <- if (x.is.log){
    x
  } else {
    log10(x)
  }
  fy <- y
  nx <- length(fx)
  COHFUN <- stats::approxfun(fx, fy, rule=2)
  fxnew <- seq(min(fx), max(fx), length.out=nx*2)
  fynew <- COHFUN(fxnew)
  lsmoo <- stats::loess.smooth(fxnew, fynew, family='gaussian', evaluation=nx, ...)
  if (!x.is.log) lsmoo[['x']] <- 10 ** lsmoo[['x']]
  return(lsmoo)
}