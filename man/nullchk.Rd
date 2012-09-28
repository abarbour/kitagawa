\name{.nullchk}
\alias{.nullchk}
\title{quickly check null, stopifnot}
\usage{
  .nullchk(X)
}
\arguments{
  \item{X}{something to be checked (vector, scalar, ...)}
}
\value{
  NULL
}
\description{
  quickly check null, stopifnot
}
\examples{
\dontrun{
.nullchk(1:10)
.nullchk(NULL)
.nullchk(c(1:10,NULL))
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

