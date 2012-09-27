\name{kitplt}
\alias{kitplt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kitplt(Resp., xlim = c(-4, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Resp.}{
%%     ~~Describe \code{Resp.} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Resp., xlim = c(-4, 0)) 
{
    print(summary(Resp.))
    lFrq. <- log10(Resp.[, 1]/2/pi)
    Amp. <- Resp.[, 2]
    Phs. <- Resp.[, 3] * 180/pi
    origpar <- par(no.readonly = TRUE)
    par(mar = c(3.5, 3, 1, 2), oma = rep(0, 4))
    layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
    plot(lFrq., log10(Amp.), type = "l", ylim = c(5, 7), yaxs = "i", 
        ylab = "Amplitude [m/strain]", xlim = xlim, xaxs = "i", 
        xlab = "")
    plot(lFrq., Phs., type = "l", ylim = c(120, 180), yaxs = "i", 
        ylab = "Phase Shift [degree]", xlim = xlim, xaxs = "i", 
        xlab = "")
    mtext(text = "Frequency [Hz]", side = 1, line = 2)
    par(origpar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
