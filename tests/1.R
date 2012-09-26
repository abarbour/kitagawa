#source("/Users/abarbour/nute.processing/development/kitagawa/rsrc/.sourceloads.R")
source("/Users/abarbour/kook.processing/R/dev/packages/kitagawa/rsrc/.sourceloads.R")
# some constants
F. <- 10 ** seq.int(-4,0,by=0.1) # frequency
Omega. <- 2 * pi * F.
# Kitagawa Fig 7 and Tab 1
S. <- 1e-6	# Storativity nondimensional
T. <- 1e-6	# Transmissivity m**2 / s
Rs. <- .135	# Radius of screened portion of well, m

# some constants (functions really) from those constants
Alpha. <- omega_constants(Omega., c.type="alpha", S.=S., T.=T., Rs.=Rs.)

# some constants from those constants
A. <- alpha_constants(Alpha., c.type="A")
matplot(log10(F.), Arg(A.[,2:3]), type="l")
#
Am. <- Mod(A.[,6:7])
plot(log10(F.), log10(Am.[,1]), type="l")
lines(log10(F.), log10(Am.[,2]), type="l",col="red")
#
