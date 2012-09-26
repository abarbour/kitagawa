##
##
##
#
##
##  Constants that depend on omega S, T, and Rs
##
.omega_constants.default <- function(omega=0, c.type=c("alpha"), ...){
	c.type <- match.arg(c.type)
	c.meth <- switch(c.type,alpha=".wc_alpha")
	c.calc <- function(...) UseMethod(c.meth)
	toret <- c.calc(omega,...)
	return(toret)
}
.wc_alpha.default <- function(omega, S., T., Rs){
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
	if (missing(Rs)){
		Rs <- 1
		warning("radius was missing, used 1")
	}
	alpha <- sqrt( omega * S. / T. ) * Rs
	nullchk(alpha)
	return(alpha)
}
##
##  Constants that depend on alpha (and thus recursively on omega, S, T, Rs)
##
.alpha_constants.default <- function(alpha=0, c.type=c("Phi","Psi","A","Kel"), ...){
	c.type <- match.arg(c.type)
	c.meth <- switch(c.type, 
		Phi=".ac_PhiPsi", Psi=".ac_PhiPsi", #PhiPsi=".ac_PhiPsi", 
		A=".ac_A", 
		Kel=".ac_Kel")
	c.calc <- function(...) UseMethod(c.meth)
	toret <- c.calc(alpha,...)
	return(toret)
}
# Kelvin functions
.ac_Kel.default <- function(alpha){
	require(kelvin, warn.conflicts=FALSE)
	# nu is order of Kelvin function
	Kel <- kelvin::Kelvin(alpha, nu.=0, nSeq=2, return.list=FALSE)
	# returns K0,K1 (complex matrix)
	toret <- base::cbind(alpha, Kel)
	colnames(toret) <- c("alpha","Kel0","Kel1")
	return(toret)
}
# Phi and Psi are both calculated from Kel constants
.ac_PhiPsi.default <- function(alpha){
	# returns alpha,Kel0,Kel1
	Kel <- alpha_constants(alpha, c.type="Kel")
	stopifnot(ncol(Kel)==3)
	K1 <- Kel[,3] # Kelvin with nu.=1
	K.r <- Re(K1)
	K.r2 <- K.r * K.r
	K.i <- Im(K1)
	K.i2 <- K.i * K.i
	# incomplete Phi & Psi
	Phi <- (K.r + K.i)
	Psi <- (K.r - K.i)
	Kir2 <- -1 * sqrt(2) * alpha * (K.r2 - K.i2)
	# Kitagawa equations 10,11
	# Scale and divide by sum of squares == complete Phi & Psi and
	# combine with alpha, K0, and K1:
	toret <- base::cbind(Kel, Phi/Kir2, Psi/Kir2)
	return(toret)
}
# A, which depends on Phi, Psi, and zeroth-order Kelvin
.ac_A.default <- function(alpha){
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
##
