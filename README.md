# kitagawa

Tools to calculate the theoretical spectral response 
of fluid-pressure in a water well
to harmonic strains (e.g., tides, long-period seismic waves).

[![Travis Build Status](https://travis-ci.org/abarbour/kitagawa.svg)](https://travis-ci.org/abarbour/kitagawa) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/abarbour/kitagawa?branch=master&svg=true)](https://ci.appveyor.com/project/abarbour/kitagawa) [![Code Coverage](https://codecov.io/gh/abarbour/kitagawa/branch/master/graph/badge.svg)](https://codecov.io/gh/abarbour/kitagawa) [![License](https://img.shields.io/badge/license-GPL-orange.svg)](https://www.gnu.org/licenses/gpl-2.0.html) [![Downloads](https://cranlogs.r-pkg.org/badges/kitagawa)](https://www.r-pkg.org/pkg/kitagawa)

## Models of spectral response

This code calculates the response at two types of wells: a sealed well and
an open well (exposed to atmosphere).

### Sealed Well

The theoretical model for a sealed well, where fluids are isolated from atmospheric pressure, 
responding to dilational strains from seismic waves is from 
[Kitagawa, et al. (2011)](https://doi.org/10.1029/2010JB007794 "Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain") which this package is named after.

### Open Well

The first theoretical model for a sealed well responding to seimic displacements is from 
[Cooper, et al. (1965)](https://doi.org/10.1029/JZ070i016p03915 "The response of well-aquifer systems to seismic waves").

This package also includes support for the models in
[Hsieh, et al. (1987)](https://doi.org/10.1029/WR023i010p01824 "Determination of aquifer transmissivity from Earth tide analysis").
[Rojstaczer (1988)](https://doi.org/10.1029/JB093iB11p13619 "Intermediate period response of water levels in wells to crustal strain: Sensitivity and noise level"), and
[Liu, et al. (1989)](https://doi.org/10.1029/JB094iB07p09453 "Seismically induced water level fluctuations in the Wali Well, Beijing, China"), which are based on various sources (i.e., tides, atmospheric pressure, and seismic waves).

## Getting Started

You can install the package via
[CRAN](https://cran.r-project.org/package=kitagawa)
from within the `R` environment:

    install.packages("kitagawa")

Load the package library and take a look at the vignettes:

    library(kitagawa)
    vignette(package='kitagawa')
    
### Installing the Development Version

Should you wish to install the development version
of this software, the [devtools][2] library
will be useful:

    library(devtools)
    install_github("abarbour/kitagawa")

[2]: https://cran.r-project.org/package=devtools
