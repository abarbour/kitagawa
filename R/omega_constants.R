#' calculate constants that depend on \eqn{\omega}
#' 
#' e.g. Kitagawa equation 12
#'
#' @name omega_constants
#' @export
#' 
#' @param omega   frequency,  \eqn{[rad/sec]}
#' @param c.type  the constant to calculate 
#' @param ...     additional params passed to calculator, e.g. S., T., Rs.
#'
#' @return values of the 'c.type' constant (currently only 'alpha', i.e. eq 12)
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' omega_constants() # 0, with warnings about S., T., Rs.
#' omega_constants(T.=1,S.=1,Rs.=1)  # 0, no warnings
#' omega_constants(1:10)  # sequence, with warnings about S., T., Rs.
#' omega_constants(1:10,T.=1,S.=1,Rs.=1) # sequence, no warnings
omega_constants <-
function(omega=0, c.type=c("alpha"), ...) UseMethod("omega_constants")

#' @return \code{NULL}
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
