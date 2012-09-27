#
# unit-test 2: replicate plot 7
#
#source("/Users/abarbour/kook.processing/R/dev/packages/kitagawa/tests/1.R")
source("tests/1.R")
# frequencies in F. or Omega.
stopifnot(exists("Omega."))
#
# parameters:
#	omega, T., S., Vw, Rs, Ku, B.,
#	Avs=1,
#	Aw=1,
#	rho=1.003,
#	Kw=2.2e9,
#	grav=9.81,
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
#
tmp.resp <- strain_response(omega=Omega.,
 T.=1e-6,	# transmissivity, m**2/s
 S.=1e-6,	# storativity
 Avs.=1,	# amplification - vol strain
 Aw.=4.8,	# amplification - well
 Vw.=sensing_volume(Rc., Lc., Rs., Ls.),	# volume well, m**3
 Rs.=Rs.,	# radius screened portion, m
 B.=0.8,	# Skempton
 Ku.=40e9,	# undrained bulk mod
 Kw.=2.2e9)	# undrained bulk mod
#
# plot it
kitplot(tmp.resp)
##
##
