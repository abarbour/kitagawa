##
##
##
#
#
.strain_response.default <- function(omega, 
	T., S., Vw., Rs., Ku., B.,
	Avs.=1,
	Aw.=1,
	rho.=1000, 
	Kw.=2.2e9,
	grav.=9.81,
	freq.units=NULL
	){
	#
	# calculate Kitagawa equation 17
	#
	# rho - water density			kg/m**3
	# grav - gravitational acceleration	m/s**2
	# T. - transmissivity			m**2/s
	# S. - storativity			N.D.
	# Vw - volume well			m**3
	# Aw - amplification factor of well volume change for E_kk
	#					N.D.
	# Avs - amplification factor for volumetric strain E_kk,obs/E_kk
	#					N.D.
	# Ku - undrained bulk modulus		Pa
	# Kw - bulk modulus of fluid		Pa
	# Rs - radius of screened portion	m
	# B. - Skempton's coefficient		N.D.
	# 
	fc <- switch(match.arg(freq.units, c("rad_per_sec","Hz")),
		rad_per_sec=1,
		Hz=2*pi
	)
	omega <- fc*omega
	# Alpha function
	Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs=Rs.)
	# A1, and A2 (functions of Phi and Psi, calculated internally)
	Amat <- alpha_constants(Alpha., c.type="A")
	stopifnot(ncol(Amat)==7)
	rm(Alpha.)  # cleanup
	#  A1,2 are in Mod(A.[,6:7]) 
	A12 <- Mod(Amat[,6:7])  # is complex, but imag is zero, so == abs
	rm(Amat)    # cleanup
	A1 <- A12[,1]
	A2 <- A12[,2]
	rm(A12)     # cleanup
	#
	# prevent duplicate calculations
	rhog <- rho. * grav.
	#print(summary(omega))
	TVFRG <- 2 * pi * T. / omega / Vw. / rhog
	#
	# calculate amp and phase of response
	# full complex may not be working... [ ]
	#
	tmpd. <- Ku. * B. / Aw. * TVFRG - A2
	rNum. <- tmpd. * tmpd. + A1 * A1
	rm(tmpd.)
	tmpd. <- Kw. * TVFRG  -  A2
	rDen. <- tmpd. * tmpd. + A1 * A1
	rm(tmpd.)
	# amplitude, Kitagawa equation 20
	Amp. <- Kw. * Aw. / Avs. / rhog * sqrt(rNum. / rDen.)
	# phase, Kitagawa equation 21
	Y. <- (Kw. - Ku. * B. / Aw.) * TVFRG * A1
	X. <- (Ku. * B. / Aw. * TVFRG - A2) * (Kw. * TVFRG - A2) + A1 * A1
	Phs. <- atan2(-1*Y.,-1*X.)
	# params?
	# attributes?
	# message?
	toret <- cbind(omega, Amp., Phs.)
	return(toret)
}
