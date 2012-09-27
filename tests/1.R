#
library(devtools)
owd<-getwd()
setwd("../..")
document('kitagawa')
setwd(owd)
#
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
#
Am. <- Mod(A.[,6:7])
#
