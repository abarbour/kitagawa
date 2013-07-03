\docType{data}
\name{kitagawa-constants}
\alias{constants}
\alias{kitagawa-constants}
\title{Constants used as defaults}
\format{List of 3
 $ radians:List of 2
 $ water  :List of 2
 $ gravity: num 9.81}
\usage{
  .constants

  constants(do.str = TRUE)
}
\arguments{
  \item{do.str}{logical; should the structure be printed?}
}
\description{
  The response of an aquifer depends on its mechanical and
  hydrological properties; if these are not known or
  specified, these constants are used.
}
\details{
  The helper function \code{\link{constants}} shows (the
  structure of, optionally) and returns \code{.constants}.

  \subsection{Constants}{ \describe{ \item{For water:
  }{Density and bulk modulus} \item{Gravity: }{Standard
  gravitational acceleration at 6371km radius} } }
}
\seealso{
  \code{\link{well_response}} and
  \code{\link{open_well_response}}
}
\keyword{datasets}

