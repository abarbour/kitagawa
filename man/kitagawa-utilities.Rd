\docType{methods}
\name{kitagawa-utilities}
\alias{.in0to1}
\alias{.nullchk}
\alias{is.owrsp}
\alias{is.wrsp}
\alias{kitagawa-utilities}
\title{General utility functions}
\usage{
  .nullchk(X)

  .in0to1(X)

  is.wrsp(X)

  is.owrsp(X)
}
\arguments{
  \item{X}{something to be checked (vector, scalar, wrsp
  object, ...)}
}
\description{
  General utility functions
}
\details{
  \code{\link{.nullchk}} quickly checks for \code{NULL} and
  \code{NA}, and raises an error if \code{TRUE}; \emph{This
  function is not likely to be needed by the user.}

  \code{\link{.in0to1}} checks if values are numeric and in
  [0,1] (inclusive).

  \code{\link{is.wrsp}} and \code{\link{is.owrsp}} report
  whether an object has S3 class 'wrsp' or 'owrsp',
  respectively.  Such an object would be returned by, for
  example, \code{\link{well_response}}.
}
\examples{
\dontrun{
.nullchk(1:10) # OK
.nullchk(NULL) # error
.nullchk(c(1:10,NULL)) # error
.nullchk(NA) # error
.nullchk(c(1:10,NA)) # error

.in0to1(1:10) # error
.in0to1(NULL) # error
.in0to1(c(1:10,NULL)) # error
.in0to1(NA) # error
.in0to1(c(1:10,NA)) # error
.in0to1(c(1,NA)) # error

is.wrsp(1) # FALSE
}
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{kitagawa-package}}
}

