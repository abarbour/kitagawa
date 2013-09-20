\name{.in0to1}
\alias{.in0to1}
\title{Check if in [0,1]}
\usage{
  .in0to1(X)
}
\arguments{
  \item{X}{something to be checked (vector, scalar, ...)}
}
\description{
  Check if in [0,1]
}
\examples{
\dontrun{
.in0to1(1:10) # error
.in0to1(NULL) # error
.in0to1(c(1:10,NULL)) # error
.in0to1(NA) # error
.in0to1(c(1:10,NA)) # error
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  Other utilities: \code{\link{.nullchk}},
  \code{\link{omega_norm}}, \code{\link{sensing_volume}}
}

