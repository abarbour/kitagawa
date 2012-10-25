\name{sensing_volume}
\alias{sensing_volume}
\title{Calculate volume of fluids in the sensing region of the borehole.}
\usage{
  sensing_volume(rad_grout, len_grout, rad_screen,
    len_screen)
}
\arguments{
  \item{rad_grout}{radius of the grouting \eqn{[m]}}

  \item{len_grout}{length of the grouting \eqn{[m]}}

  \item{rad_screen}{radius of the screened interval
  \eqn{[m]}}

  \item{len_screen}{length of the screened interval
  \eqn{[m]}}
}
\value{
  scalar, with units of \eqn{[m^3]}
}
\description{
  This function calculates the volume of fluid in the
  screened section, namely \strong{Equation 2} in Kitagawa
  et al (2011).
}
\details{
  Although typical scientific boreholes with water-level
  sensors are drilled very deeply, pore-fluids are only
  allowed to flow through a relatively short section, known
  as the "screened" section.  The calculation assumes two
  pairs of radii and lengths: one for the cemented (grout)
  section, and another for the screened section.

  The volume calculated is \deqn{ \pi R_C^2 (L_C - L_S) +
  \pi R_S^2 L_S } where \eqn{R} and \eqn{L} denote radius
  and length respectively, and subscripts \eqn{C} and
  \eqn{S} denote the cemented and screened sections
  respectively.

  This calculation assumes the measurement is for a sealed
  well.
}
\examples{
#### dummy example
sensing_volume(1, 1, 1, 1)
#
#### a more physically realistic calculation:
# Physical params applicable for B084 borehole
# (see: http://pbo.unavco.org/station/overview/B084/ for details)
#
Rc <- 0.0508   # m, radius of water-sensing (2in)
Lc <- 146.9    # m, length of grouted region (482ft)
Rs <- 3*Rc     # m, radius of screened region (6in)
Ls <- 9.14     # m, length of screened region (30ft)
#
# calculate the sensing volume for the given well parameters
sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}}, \code{\link{kitplot}}
}

