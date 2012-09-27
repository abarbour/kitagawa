#' Plot response spectrum
#' 
#' used to plot frequency response spectra as in kitagawa
#'
#' @name kitplot
#' @export
#' 
#' @param Resp.  the response information matrix(f,Amp,Phi,ncol=3)
#' @param xlim   frequency limits (assume log10 scale)
#'
#' @return NULL
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' kitplot(data.frame(f=2*pi*10**seq(-4,0,length.out=10), amp=1e6*rep(1,10), phs=.9*pi*rep(1,10))) # bad example
kitplot <-
function(Resp., xlim=c(-4,0)){
	#
	# reproduce plots as in Kitagawa
	#
  Resp. <- as.matrix(Resp.)
  print((Resp.))
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
  par(mar=c(3.5,3,1,2),oma=rep(0,4))
	layout(matrix(c(1,2), 2, 1, byrow = TRUE))
	# amplitude
	plot(lFrq., log10(Amp.),
		type="l",
		ylim=c(5,7), 
		yaxs="i", ylab="Amplitude [m/strain]", 
		xlim=xlim, xaxs="i", xlab=""
	)
	# phase shift
	plot(lFrq., Phs.,
		type="l",
		ylim=c(120,180), 
		yaxs="i", ylab="Phase Shift [degree]",
		xlim=xlim, xaxs="i", xlab=""
	)
  mtext(text="Frequency [Hz]",side=1,line=2)
	par(origpar)
}
