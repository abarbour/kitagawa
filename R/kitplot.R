#' Quickly plot the amplitude and phase spectra
#' 
#' Used to mimic plots of frequency response spectra, as in Kitagawa et al (2011), 
#' e.g. \strong{Figures 7--9}.
#' 
#' This is primarily a diagnostic tool, and is thus not very flexible in its 
#' implementation.
#' 
#' The input data are assumed to be structured as they would be out
#' of \code{\link{well_response}}, specifically three vectors representing:
#' \describe{
#'  \item{Radial frequency}{\eqn{[radians/s]}}
#'  \item{Amplitude}{\eqn{[m/strain]}}
#'	\item{Phase}{\eqn{[radians]}}
#' }
#' These will be transformed to have units of \eqn{[\log10 Hz]}, \eqn{[\log10 m/strain]}, and
#' \eqn{[rad]} respectively, unless \code{prep.resp=FALSE}.
#'
#' @name kitplot
#' @export
#' 
#' @param Resp.  the response information as from well_response (see \strong{Details})
#' @param xlim.   frequency limits (assume log10 Hz scale)
#' @param ylims  list of limits (assume log10 m/strain for amplitude), use \code{amp} and \code{phs} for variable assignments
#' @param prep.resp  boolean, should the units of \code{Resp.} be appropriately transformed ?
#'
#' @return A \code{data.frame} with the (transformed) values that have been
#' plotted; this will 
#' include \emph{only} the subset of data falling in the range of \code{xlim.}
#' inclusively.
#'  
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @references Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi (2011),
#' Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain,
#' \emph{J. Geophys. Res.}, \strong{116}, B08301, doi:10.1029/2010JB007794.
#' 
#' @references \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
#'
#' @seealso \code{\link{well_response}}
#' 
#' @examples
#' # dummy example: get some lines on the figure
#' n <- 10
#' ones <- rep(1, n)
#' fakeResp <- data.frame(f=2*pi*10**seq(-4, 0, length.out=n), amp=1e6*ones, phs=.9*pi*ones)
#' kitplot(fakeResp)
#' # focus in on a certain range:
#' fakeResp.foc <- kitplot(fakeResp, xlim.=c(-3, -1), ylims=list(amp=c(5.5, 6.5), phs=180*c(0, 1)))
#' kitplot(fakeResp.foc, prep.resp=FALSE) # fakeResp.foc has already been transformed
kitplot <-
function(Resp., xlim.=c(-4,0), ylims=list(amp=c(5,7), phs=180*c(-1,1)), prep.resp=TRUE){
	#
	# reproduce plots as in Kitagawa
	#
  Resp. <- as.matrix(Resp.)
  stopifnot(ncol(Resp.)>=3)
  #...
	#print(summary(Resp.))
	#
  if (prep.resp){
    Resp.[,1] <- log10(Resp.[,1] / 2 / pi)
    Resp.[,2] <- log10(Resp.[,2])
    Resp.[,3] <- Resp.[,3]*180/pi
  }
  Resp. <- as.data.frame(Resp.)
  names(Resp.) <- c("l10Freq","l10Amp","dPhs")
  # subset data
  xlim. <- sort(xlim.) # make sure order is correct
	Resp. <- subset(Resp., (l10Freq>=xlim.[1]) & (l10Freq<=xlim.[2]) )
  # 
  # CRAN check NOTE otherwise (a hack, to be sure)
  l10Freq<-NULL
  rm(l10Freq)
	##
	origpar <- par(no.readonly = TRUE)
  par(mar=c(3.5,3.5,1,1), oma=rep(0,4), mgp=c(2.3, 1, 0))
	layout(matrix(c(1,2), 2, 1, byrow = TRUE))
	# amplitude
  ylim. <- ylims$amp
	plot(Resp.$l10Freq, Resp.$l10Amp,
		type="l",
		ylim=ylim.,
		yaxs="i", ylab="Amplitude [log10 m/strain]", 
    xlim=xlim.,
		xaxs="i", xlab=""
	)
	# phase shift
  ylim. <- ylims$phs
	plot(Resp.$l10Freq, Resp.$dPhs,
		type="l",
		ylim=ylim.,
		yaxs="i", ylab="Phase Shift [degrees]",
	  xlim=xlim.,
		xaxs="i", xlab=""
	)
  mtext(text="Frequency [log10 Hz]",side=1,line=2)
	par(origpar)
  return(invisible(Resp.))
}
