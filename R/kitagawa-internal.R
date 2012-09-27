.ac_A.default <-
function(alpha){
	PhiPsi <- alpha_constants(alpha, c.type="Phi")
	stopifnot(ncol(PhiPsi) == 5)
	#columns: 1 is alpha, 2 is K0, 3 is K1, 4 is Phi, 5 is Psi
	K0 <- PhiPsi[,2]
	K.r <- Re(K0)
	K.i <- Im(K0)
	Phi <- PhiPsi[,4]
	Psi <- PhiPsi[,5]
	# Kitagawa equations 18,19
	A1 <- Phi * K.r - Psi * K.i
	A2 <- Psi * K.r + Phi * K.i
	#A1A2 <- base::cbind(A1,A2)
	toret <- base::cbind(PhiPsi, A1, A2)
	return(toret)
}
.ac_Kel.default <-
function(alpha){
	require(kelvin)
	# nu is order of Kelvin function
	Kel <- kelvin::Keir(alpha, nu.=0, nSeq.=2, return.list=FALSE)
	# returns K0,K1 (complex matrix)
	toret <- base::cbind(alpha, Kel)
	colnames(toret) <- c("alpha","Kel0","Kel1")
	return(toret)
}
.ac_PhiPsi.default <-
function(alpha){
	# calculate Kelvin functions
	Kel <- alpha_constants(alpha, c.type="Kel")
	# returns alpha,Kel0,Kel1
	stopifnot(ncol(Kel)==3)
	K1 <- Kel[,3] # Kelvin with nu.=1
	K.r <- Re(K1)
	K.r2 <- K.r * K.r
	K.i <- Im(K1)
	K.i2 <- K.i * K.i
	# incomplete Phi & Psi
	Phi <- (K.r + K.i)
	Psi <- (K.r - K.i)
  # bug: was a minus (incorrectly)
	Kir2 <- -1 * sqrt(2) * alpha * (K.r2 + K.i2)
	# Kitagawa equations 10,11
	# Scale and divide by sum of squares == complete Phi & Psi and
	# combine with alpha, K0, and K1:
	toret <- base::cbind(Kel, Phi/Kir2, Psi/Kir2)
	return(toret)
}
.alpha_constants.default <-
function(alpha=0, c.type=c("Phi","Psi","A","Kel"), ...){
	c.type <- match.arg(c.type)
	c.meth <- switch(c.type, 
                   Phi=".ac_PhiPsi", 
                   Psi=".ac_PhiPsi", #PhiPsi=".ac_PhiPsi", 
                   A=".ac_A", 
                   Kel=".ac_Kel")
	c.calc <- function(...) UseMethod(c.meth)
	toret <- c.calc(alpha,...)
	return(toret)
}
.omega_constants.default <-
function(omega=0, c.type=c("alpha"), ...){
	c.type <- match.arg(c.type)
	c.meth <- switch(c.type,alpha=".wc_alpha")
	c.calc <- function(...) UseMethod(c.meth)
	toret <- c.calc(omega,...)
	return(toret)
}
.strain_response.default <-
function(omega, 
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
	Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
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
.wc_alpha.default <-
function(omega, S., T., Rs.){
	# Kitagawa equation 12
	# storativity S, transmissivity T, radius of screened portion Rs
	if (missing(S.)){
		warning("storativity was missing, used 1")
		S. <- 1
		}
	if (missing(T.)){
		warning("tranmissivity was missing, used 1")
		T. <- 1
	}
	if (missing(Rs.)){
		Rs. <- 1
		warning("radius was missing, used 1")
	}
	alpha <- sqrt( omega * S. / T. ) * Rs.
	nullchk(alpha)
	return(alpha)
}
