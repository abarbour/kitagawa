#
# Replicate Rojstaczer 1988, fig 3
#
library(signal)
library(kitagawa)
library(kook) # locally only
# Note Rmpfr has erf and erfc functions... 
#
#
omega <- 10**seq(-3,2,by=0.1)
z <- 1
Trans <- 1
Stor <- 1
Diffus <- Trans/Stor
# nondim freq
Q <- 10**seq(-3,2,by=0.1) # == z**2 omega / 2 D
omega <- Q * 2 * Diffus / z**2
#(z * omega_constants(omega, c.type="diffusivity_time", T.=Trans, S.=Stor))**2

wrsp <- open_well_response(omega, T.=Trans, S.=Stor, z.=z, 
                           #freq.units="Hz", 
                           model = "rojstaczer")
crsp <- wrsp[,2]

par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
lQ <- log10(Q)
# Amplitude
As <- 0.05 # cm/nE
rhog <- 1000*9.81
Gain <- Mod(crsp)
Gain_head <- Gain*As/rhog #or Gain*As if as.pressure=FALSE 
plot(lQ, Gain_head, type="l", #ylim=c(0, 1.2), #*As,
     xaxt="n", lwd=2,
     ylab="", #" (cm/nanostrain)", 
     xlab="", #Dimensionless frequency Q",
     main="")
kook::log10_ticks()
mtext("Open Well Response to Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain rel. static-confined areal strain response", adj=0)
#mtext(sprintf(" areal strain\nresponse: %.02f cm/ne",As), padj=1, adj=0.1, line=-1, cex=0.7)
par(mar=c(4,4,1.1,2))
# Phase
Phs <- Arg(crsp) # will wrap to -pi/pi
uPhs <- signal::unwrap(Phs, tol=pi/30)
plot(lQ, Phs*180/pi, type="l", lty=3, ylim=c(-190, -130), xaxt="n",
     ylab="degrees", xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, uPhs*180/pi, type="l", lwd=2)
kook::log10_ticks()
mtext("(b) Phase rel. strain", adj=0)
#