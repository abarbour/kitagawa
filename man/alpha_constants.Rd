\docType{methods}
\name{alpha_constants}
\alias{alpha_constants}
\alias{alpha_constants.default}
\title{Calculate any constants depending on \eqn{\alpha}}
\usage{
  alpha_constants(alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"))

  \method{alpha_constants}{default} (alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"))
}
\arguments{
  \item{alpha}{the constant alpha (see
  \code{\link{omega_constants}})}

  \item{c.type}{the constant to calculate}
}
\value{
  Complex matrix having values representing the constant
  represented by \code{c.type}, \emph{as well as} any other
  \eqn{\alpha}-dependent constants which are needed in the
  computation.
}
\description{
  This function accesses the appropriate method to
  calculate the \eqn{\alpha}-dependent constant associated
  with the choice of \code{c.type}. There are currently
  four such constants, which correspond to
  \strong{Equations 10, 11, 18, 19} in Kitagawa et al
  (2011).
}
\details{
  \emph{This function is not likely to be needed by the
  user.}

  The constant \eqn{\alpha} is a function of frequency
  \eqn{\omega} as well as aquifer and well parameters; it
  is formally defined as \deqn{\alpha \equiv R_S
  \sqrt{\omega S / T}} where \eqn{S} is the storativity,
  \eqn{T} is the aquifer's effective transmissivity, and
  \eqn{R_S} is the radius of the screened portion of the
  well.

  The various constants which may be calculated are
  \describe{ \item{\code{Phi}}{Given as \eqn{\Phi} in Eqn.
  10} \item{\code{Psi}}{Given as \eqn{\Psi} in Eqn. 11}
  \item{\code{A}}{Given as \eqn{A_i, i=1,2} in Eqns. 18,
  19} \item{\code{Kel}}{\eqn{\mathcal{K}}, the complex
  Kelvin functions (see Abramowitz and Stegun, 1972)} }
}
\examples{
alpha_constants() # kelvin::Keir gives warning
alpha_constants(1)  # defaults to constant 'Phi' (note output also has Kel)
alpha_constants(1:10, c.type="A")  # constant 'A' (again, note output)
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\references{
  Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi
  (2011), Frequency characteristics of the response of
  water pressure in a closed well to volumetric strain in
  the high-frequency domain, \emph{J. Geophys. Res.},
  \strong{116}, B08301, doi:10.1029/2010JB007794.

  \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}

  Abramowitz, M. and Stegun, I. A. (Eds.). "Kelvin
  Functions." \eqn{\S 9.9} in Handbook of Mathematical
  Functions with Formulas, Graphs, and Mathematical Tables,
  9th printing. New York: Dover, pp. 379-381, 1972.
}
\seealso{
  \code{\link{omega_constants}},
  \code{\link{well_response}}, \pkg{kelvin}
}

