#' @title Generic methods for objects with class \code{'wrsp'}.
#' 
#' @description
#' An object with class 'wrsp' is a list containing the
#' response information, and the
#' mechanical, hydraulic, and material properties used to
#' generate the response.
#' 
#' @details
#' The response information is a
#' matrix with three columns: frequency, amplitude, and phase 
#' [\eqn{\omega}, \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be as they were input,
#' \eqn{A_\alpha (\omega)} is in meters per strain, 
#' and \eqn{\Phi_\alpha (\omega)} is in radians.
#' 
#' @name wrsp-methods
#' @aliases wrsp
#' @rdname wrsp-methods
#' @docType methods
#'
#' @seealso 
#' \code{\link{well_response}}
#' 
#' \code{\link{kitagawa-package}}
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' 
#' @param object wrsp object
#' @param series character; the series to plot (amplitude or phase)
#' @param pch point character, as in \code{\link{par}}
#' @param xlims limits for x-axis (applies to both amp and phs frames)
#' @param ylims optional list of limits for y-axis (i.e., \code{list(amp=c(..),phs=c(...))})
#' @param ... optional arguments
#' 
#' @examples
#' W <- well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)
#' str(W)
#' print(W)
#' print(summary(W))
#' #
#' # Plot the response
#' plot(rnorm(10), xlim=c(-1,11), ylim=c(-2,2))
#' lines(W)
#' lines(W, "phs", col="red")
#' points(W)
#' points(W, "phs")
#' #
#' Wdf <- as.data.frame(W)
#' plot(Amp. ~ omega, Wdf)
#' plot(Phs. ~ omega, Wdf)
#' #
#' # or use the builtin method
#' plot(W)
NULL

#' @rdname wrsp-methods
#' @alias as.data.frame.wrsp
#' @method as.data.frame wrsp
#' @S3method as.data.frame wrsp
as.data.frame.wrsp <- function(object, ...){
  WR <- object[["Response"]]
  df <- as.data.frame(WR)
  #names(df) <- "n.wrsp"
  return(df)
}
#' @rdname wrsp-methods
#' @alias data.frame.wrsp
#' @method data.frame wrsp
#' @S3method data.frame wrsp
data.frame.wrsp <- as.data.frame.wrsp

#' @rdname wrsp-methods
#' @aliases print.wrsp
#' @method print wrsp
#' @S3method print wrsp
print.wrsp <- function(object, ...){
  stopifnot(is.wrsp(object))
  message("Sealed well-response:")
  WR <- object[["Response"]]
  WR <- as.data.frame(WR)
  print(head(WR,3))
  message("\t...")
  print(tail(WR,3))
}

#' @rdname wrsp-methods
#' @aliases summary.wrsp
#' @method summary wrsp
#' @S3method summary wrsp
summary.wrsp <- function(object, ...){
  stopifnot(is.wrsp(object))
  WR <- object[["Response"]]
  toret <- summary.default(WR)
  class(toret) <- "summary.wrsp"
  return(toret)
}

#' @rdname wrsp-methods
#' @aliases print.summary.wrsp
#' @method print summary.wrsp
#' @S3method print summary.wrsp
print.summary.wrsp <- function(object, ...){
  message("Sealed well-response summary:")
  print.summaryDefault(object)
}

#' @rdname wrsp-methods
#' @aliases lines.wrsp
#' @method lines wrsp
#' @S3method lines wrsp
lines.wrsp <- function(object, series=c("amp","phs"), ...){
  stopifnot(is.wrsp(object))
  series <- match.arg(series)
  WR <- object[["Response"]]
  x <- WR[,1]
  y <- switch(series, amp=WR[,2], phs=WR[,3])
  graphics::lines(x, y, ...)
}

#' @rdname wrsp-methods
#' @aliases points.wrsp
#' @method points wrsp
#' @S3method points wrsp
points.wrsp <- function(object, series=c("amp","phs"), pch="+", ...){
  stopifnot(is.wrsp(object))
  series <- match.arg(series)
  WR <- object[["Response"]]
  x <- WR[,1]
  y <- switch(series, amp=WR[,2], phs=WR[,3])
  graphics::points(x, y, pch=pch, ...)
}

#' @rdname wrsp-methods
#' @aliases plot.wrsp
#' @method plot wrsp
#' @S3method plot wrsp
plot.wrsp <- function(object, 
                      xlims=c(-3,1), 
                      ylims=list(amp=NULL, phs=185*c(-1,1)), ...){
  #
  stopifnot(is.wrsp(object))
  WR <- object[["Response"]]
  fu <- object[["Omega"]][["Units"]]
  fc <- switch(fu, rad_per_sec=2*pi, Hz=1)
  stopifnot(!is.null(fc))
  frq <- log10(WR[,1] / fc)
  amp <- log10(WR[,2])
  phs <- WR[,3]*180/pi
  ##
  origpar <- par(no.readonly = TRUE)
  par(mar=c(2,4,0,1), 
      oma=c(2,0.1,2,0.1), 
      tcl=-0.3,
      mgp=c(2.5, 0.5, 0), las=1)
  layout(matrix(c(1,2), ncol=1), heights=c(0.5,0.5))
  # amplitude
  alims <- ylims[["amp"]]
  if (is.null(alims)) alims <- range(pretty(amp))
  plot(0,0,col=NA,
       ylim=alims,
       yaxs="i", ylab="[log10 meters/strain]", 
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  lines(frq, amp, type="l", lwd=1.5, ...)
  log10_ticks()
  mtext("Sealed well-response", font=2)
  mtext("(a) Amplitude", adj=0.015, font=4, line=-1.2)
  # phase shift
  plot(0,0,col=NA,
       ylim=ylims[["phs"]],
       yaxs="i", yaxt="n", ylab="[degrees]",
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  abline(h=c(-1,-0.5,0,0.5,1)*180, col="grey80", lty=2)
  lines(frq, phs, type="l", lwd=1.5, ...)
  lbls <- ats <- seq(-180,180,by=30)
  lbls[seq_along(lbls)%%2==0] <- ""
  axis(2, at=ats, labels=lbls)
  log10_ticks()
  mtext("(b) Phase rel. strain", adj=0.015, font=4, line=-1.2)
  mtext("Frequency [Hz]", side=1, line=2)
  on.exit(par(origpar))
  return(invisible(NULL))
}
