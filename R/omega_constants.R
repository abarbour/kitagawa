#' calculate constants that depend on \eqn{\omega}
#' 
#' e.g. Kitagawa equation 12
#'
#' @name omega_constants
#' @export
#' 
#' @param omega   frequency,  \eqn{[rad/sec]}
#' @param c.type  the constant to calculate (currently only 'alpha', i.e. eq 12)
#' @param ...     additional params passed to calculator (should not be used)
#'
#' @return Matrix with
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' omega_constants(1:10)  # dummy example for now
omega_constants <-
function(omega=0, c.type=c("alpha"), ...) UseMethod("omega_constants")

#' @return \code{NULL}
#' @rdname omega_constants
#' @docType methods
#' @method omega_constants default
#' @S3method omega_constants default
omega_constants.default <-
  function(omega, c.type=c("alpha"), ...){
    c.type <- match.arg(c.type)
    c.meth <- switch(c.type,alpha=".wc_alpha")
    c.calc <- function(...) UseMethod(c.meth)
    toret <- c.calc(omega,...)
    return(toret)
  }

#' @return \code{NULL}
#' @rdname omega_constants
#' @docType methods
#' @method omega_constants default
#' @S3method omega_constants default
.wc_alpha.default <-
  function(omega, S., T., Rs.){
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
  }
