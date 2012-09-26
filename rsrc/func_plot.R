##
## 
##
#
kitplt <- function(Resp.){
	#
	# reproduce plots as in Kitagawa
	#
	# resp will have three columns:
	#	1 - Freqs.(radians/sec)
	#	2 - Amplitude.(m/strain)
	#	3 - Phase.(radians)
	print(summary(Resp.))
	#
	lFrq. <- log10(Resp.[,1] / 2 / pi)
	Amp. <- Resp.[,2]
	Phs. <- Resp.[,3]*180/pi
	##
	origpar <- par(no.readonly = TRUE)
	layout(matrix(c(1,2), 2, 1, byrow = TRUE))
	# amplitude
	plot(lFrq., Amp., # log10(Amp.),
		type="l",
		#ylim=c(5,7), 
		yaxs="i", ylab="Amplitude [m/strain]", 
		xlim=c(-4,0), xaxs="i", xlab="Frequency [Hz]"
	)
	# phase shift
	plot(lFrq., Phs.,
		type="l",
		#ylim=c(120,190), 
		yaxs="i", ylab="Phase Shift [degree]",
		xlim=c(-4,0), xaxs="i", xlab="Frequency [Hz]"
	)
	par(origpar)
}
