#' Plot response spectrum
#' 
#' Used to mimic plots of frequency response spectra, as in Kitagawa et al 2011, 
#' e.g. Figures 7--9
#'
#' @name kitplot
#' @export
#' 
#' @param Resp.  the response information as from well_response: matrix(f,Amp,Phi,ncol=3)
#' @param xlim   frequency limits (assume log10 Hz scale)
#' @param ylims  list of limits (assume log10 m/strain for amplitude), use 'amp' and 'phs'
#'
#' @return NULL
#'  
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' # dummy example: get some lines on the figure
#' n <- 10
#' ones <- rep(1,n)
#' fakeResp <- data.frame(f=2*pi*10**seq(-4,0,length.out=n), amp=1e6*ones, phs=.9*pi*ones)
#' kitplot(fakeResp)
#' # focus:
#' kitplot(fakeResp, xlim=c(-3,1), ylims=list(amp=c(5.5,6.5),phs=180*c(0,1)))
kitplot <-
function(Resp., xlim=c(-4,0), ylims=list(amp=c(5,7),phs=180*c(-1,1))){
	#
	# reproduce plots as in Kitagawa
	#
  Resp. <- as.matrix(Resp.)
  ##print((Resp.))
  stopifnot(ncol(Resp.)>=3)
	# resp will have three or more columns:
	#	1 - Freqs.(radians/sec)
	#	2 - Amplitude.(m/strain)
	#	3 - Phase.(radians)
  #...
	#print(summary(Resp.))
	#
	lFrq. <- log10(Resp.[,1] / 2 / pi)
	Amp. <- Resp.[,2]
	Phs. <- Resp.[,3]*180/pi
	##
	origpar <- par(no.readonly = TRUE)
  par(mar=c(3.5,3.5,1,1), oma=rep(0,4), mgp=c(2.3, 1, 0))
	layout(matrix(c(1,2), 2, 1, byrow = TRUE))
	# amplitude
	plot(lFrq., log10(Amp.),
		type="l",
		ylim=ylims$amp,
		yaxs="i", ylab="Amplitude [log10 m/strain]", 
		xlim=xlim, xaxs="i", xlab=""
	)
	# phase shift
	plot(lFrq., Phs.,
		type="l",
		ylim=ylims$phs,
		yaxs="i", ylab="Phase Shift [degrees]",
		xlim=xlim, xaxs="i", xlab=""
	)
  mtext(text="Frequency [log10 Hz]",side=1,line=2)
	par(origpar)
}
