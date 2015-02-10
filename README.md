#kitagawa [![Build Status](https://travis-ci.org/abarbour/kitagawa.png?branch=master)](https://travis-ci.org/abarbour/kitagawa) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

Tools to calculate the theoretical spectral response 
of fluid-pressure in an open, or sealed, water well
to harmonic straining.

## Models

Currently there is support for the response at two types of wells:

###Sealed Well

The theoretical model for a sealed well, where fluids are isolated from atmospheric pressure, is from (http://dx.doi.org/10.1029/2010JB007794)(Kitagawa et al (2011)), which this package is named after; the abstract of original article:

	*Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain*

	_Yuichi Kitagawa, Satoshi Itaba, Norio Matsumoto, and Naoji Koizumi_

	Oscillations of water pressures and crustal strains due to the seismic waves of
	the 2010 Chile earthquake were observed in Japan. The oscillations of water
	pressures observed over the frequency range of 0.002 to 0.1 Hz were negative
	proportional to the oscillations of volumetric strains. The responses of water
	pressures in closed wells are frequency-dependent. The expression for the
	response of water pressure in a closed well to crustal strain is developed based
	on the poroelastic theory. The expression developed in the present paper
	describes the frequency characteristics of the responses. The response is useful
	for the estimation of rock properties. In addition, the responses of water
	pressure due to tidal volumetric strain are estimated and compared with the
	responses due to the seismic waves.

	J. Geophys. Res., 116, B08301, doi:10.1029/2010JB007794, 2011.

### Open Well

The theoretical model for the sealed well is from (http://dx.doi.org/10.1029/JZ070i016p03915)(Cooper et al (1965)); the abstract of original article:

	*The response of well-aquifer systems to seismic waves*

	_Hilton H. Cooper Jr., John D. Bredehoeft, Istavros S. Papadopulos,
	and Robert R. Bennett_

	The degree to which the water level in an open well fluctuates in response to a
	seismic wave is determined by the dimensions of the well, the transmissibility,
	storage coefficient, and porosity of the aquifer, and the type, period, and
	amplitude of the wave. The water level responds to pressure-head fluctuations due
	to dilatation of the aquifer and to vertical motion of the well-aquifer system;
	hence a wave that produces either of these can cause the water level to fluctuate.
	However, the response to dilatation is much larger than the response to vertical
	motion. A solution is derived for the nonsteady drawdown in the aquifer due to a
	harmonic motion of the water level. This solution is then used in the equation of
	motion of the water column to derive expressions for the amplification, which is
	defined as the ratio X0/h0 (for oscillation due to dilatation) or the ratio X0/a
	(for oscillation due to vertical motion of the well-aquifer system), where X0 is
	the amplitude of water-level fluctuation, h0 is the amplitude of pressure-head
	fluctuation, and a is the amplitude of vertical motion of well-aquifer system.
	Amplification curves are given for differing well dimensions and aquifer
	constants.

	J. Geophys. Res., 70(16), 3915–3926, doi:10.1029/JZ070i016p03915, 1965.

##Getting Started##

Firstly you'll need to install the package and it's dependencies
from [CRAN](http://cran.r-project.org/package=kitagawa)
(from within the `R` environment):

    install.packages("kitagawa", dependencies=TRUE)

then load the package library and take a look at the vignette

    library(kitagawa)
    vignette(package='kitagawa')
    
##Installing the Development Version##

Should you wish to install the development version
of this software, the [devtools][2] library
will be useful:

    install.packages("devtools", dependencies=TRUE)
    library(devtools)
    install_github("abarbour/kitagawa")

[2]: http://cran.r-project.org/web/packages/devtools
