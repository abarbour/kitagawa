#' Calculate any constants that depend on \eqn{\omega}
#' 
#' This function accesses the appropriate method to calculate the
#' \eqn{\omega}-dependent constant associated with the choice of \code{c.type}.  
#' There is currently only one such constant, corresponding to
#' \strong{Equation 12} in Kitagawa et al (2011).
#' 
#' \emph{This function is not likely to be needed by the user.}
#' 
#' The radial frequency \eqn{\omega} is formally defined as:
#' \deqn{\omega \equiv 2 \pi / \tau}
#' where \eqn{\tau} is the period of oscillation.
#' 
#' Because the computation of \eqn{\alpha} depends also on physical
#' properties, additional parameters can be
#' passed through (e.g. the transmissivity).  
#' 
#' @section Warnings:
#' 
#' \strong{In the case \code{c.type='alpha'}(the default), the 
#' parameters \code{S.}, \code{T.},  and \code{Rs.} must
#' be passed; otherwise, values are assumed to ensure the 
#' calculation does not fail, and the results are essentially meaningless.}
#' Warnings will be issued if any necessary parameters are missing, indicating
#' default values 
#' \code{S.=T.=Rs.=1} were used; these are physically
#' unrealistic.
#'
#' @name omega_constants
#' @export
#' 
#' @param omega   frequency,  \eqn{[rad/sec]}
#' @param c.type  the constant to calculate 
#' @param ...     additional params passed to calculator, i.e. \code{S., T., Rs.}
#'
#' @return Values of the constant repesented by \code{c.type}
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#'
#' @references Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi (2011),
#' Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain,
#' \emph{J. Geophys. Res.}, \strong{116}, B08301, doi:10.1029/2010JB007794.
#' 
#' @references \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
#'
#' @seealso \code{\link{alpha_constants}}, \code{\link{well_response}}
#' 
#' @examples
#' omega_constants() # 0, with warnings about S., T., Rs.
#' omega_constants(T.=1,S.=1,Rs.=1)  # 0, no warnings
#' omega_constants(1:10)  # sequence, with warnings about S., T., Rs.
#' omega_constants(1:10,T.=1,S.=1,Rs.=1) # sequence, no warnings
omega_constants <-
function(omega=0, c.type=c("alpha"), ...) UseMethod("omega_constants")

# @return \code{NULL}
#' @rdname omega_constants
#' @docType methods
#' @method omega_constants default
#' @S3method omega_constants default
omega_constants.default <-
  function(omega=0, c.type=c("alpha"), ...){
    #
    # switch constants-calculation method
    c.type <- match.arg(c.type)
    c.meth <- switch(c.type, alpha=".wc_alpha")
    #
    # here are the methods available:
    .wc_alpha.default <- function(omega, S., T., Rs.){
        # Kitagawa equation 12
        # storativity S, transmissivity T, radius of screened portion Rs
        if (missing(S.)){
          warning("storativity was missing, used 1")
          S. <- 1
        }
        if (missing(T.)){
          warning("tranmissivity was missing, used 1")
          T. <- 1
        }
        if (missing(Rs.)){
          Rs. <- 1
          warning("radius was missing, used 1")
        }
        alpha <- sqrt( omega * S. / T. ) * Rs.
        .nullchk(alpha)
        return(alpha)
    } # end .wc_alpha
    #
    # do the calculation with the method of choice
    c.calc <- function(...) UseMethod(c.meth)
    toret <- c.calc(omega,...)
    return(toret)
  }
