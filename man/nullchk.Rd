\name{.nullchk}
\alias{.nullchk}
\title{Quickly check for \code{NULL} and \code{NA}}
\usage{
  .nullchk(X)
}
\arguments{
  \item{X}{something to be checked (vector, scalar, ...)}
}
\description{
  Checks \code{NULL} and \code{NA} status, and raises an
  error if \code{TRUE}.
}
\details{
  \emph{This function is not likely to be needed by the
  user.}
}
\examples{
\dontrun{
.nullchk(1:10) # OK
.nullchk(NULL) # error
.nullchk(c(1:10,NULL)) # error
.nullchk(NA) # error
.nullchk(c(1:10,NA)) # error
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

