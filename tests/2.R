#
# unit-test 2: replicate plot 7
#
#library(devtools)
#setwd("../..")
#document('kitagawa')
source("1.R")
# frequencies in F. or Omega.
stopifnot(exists("Omega."))
#
# parameters assumed by well_response:
#	rho=1000, kg/m**3
#	Kf=2.2e9, Pa
#	grav=9.81, m/s2
#
# calculate response
# Fig 7
# Ku B / Kw Aw = 3	=> Aw==4.8 at 40GPa
#
# Using ANO1 stats
Rc. <- 0.075	# m
Lc. <- 570	# m
Rs. <- 0.135	# m
Ls. <- 15	# m
Vw. <- sensing_volume(Rc., Lc., Rs., Ls.) # m**3
#
tmp.resp <- well_response(omega=Omega.,
 T.=1e-6,	# transmissivity, m**2/s
 S.=1e-6,	# storativity
 Avs.=1,	# amplification - vol strain
 Aw.=4.8,	# amplification - well
 Vw.=Vw.,	# volume well, m**3
 Rs.=Rs.,	# radius screened portion, m
 B.=0.8,	# Skempton
 Ku.=40e9,	# undrained bulk mod
 Kf.=2.2e9)	# undrained bulk mod
#
# plot it
kitplot(tmp.resp)
##
##
