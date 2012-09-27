#' calculate constants that depend on \eqn{\alpha}
#' 
#' e.g. # Kitagawa equations 10,11,18,19
#'
#' @name alpha_constants
#' @export
#' 
#' @param alpha   the constant alpha
#' @param c.type  the constant to calculate
#' @param ...     additional params passed to calculator (should not be used)
#'
#' @return Matrix with
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' alpha_constants(1:10)  # dummy example for now
alpha_constants <-
function(alpha=0, c.type=c("Phi","Psi","A","Kel"), ...) UseMethod(".alpha_constants")

#' @return \code{NULL}
#' @rdname alpha_constants
#' @docType methods
#' @method alpha_constants default
#' @S3method alpha_constants default
.alpha_constants.default <-
  function(alpha=0, c.type=c("Phi","Psi","A","Kel"), ...){
    c.type <- match.arg(c.type)
    c.meth <- switch(c.type, 
                     Phi=".ac_PhiPsi", 
                     Psi=".ac_PhiPsi", #PhiPsi=".ac_PhiPsi", 
                     A=".ac_A", 
                     Kel=".ac_Kel")
    c.calc <- function(...) UseMethod(c.meth)
    toret <- c.calc(alpha,...)
    return(toret)
  }

#' @return \code{NULL}
#' @rdname alpha_constants
#' @docType methods
#' @method alpha_constants default
#' @S3method alpha_constants default
.ac_A.default <-
  function(alpha){
    PhiPsi <- alpha_constants(alpha, c.type="Phi")
    stopifnot(ncol(PhiPsi) == 5)
    #columns: 1 is alpha, 2 is K0, 3 is K1, 4 is Phi, 5 is Psi
    K0 <- PhiPsi[,2]
    K.r <- Re(K0)
    K.i <- Im(K0)
    Phi <- PhiPsi[,4]
    Psi <- PhiPsi[,5]
    # Kitagawa equations 18,19
    A1 <- Phi * K.r - Psi * K.i
    A2 <- Psi * K.r + Phi * K.i
    #A1A2 <- base::cbind(A1,A2)
    toret <- base::cbind(PhiPsi, A1, A2)
    return(toret)
  }

#' @return \code{NULL}
#' @rdname alpha_constants
#' @docType methods
#' @method alpha_constants default
#' @S3method alpha_constants default
.ac_Kel.default <-
  function(alpha){
    require(kelvin)
    # nu is order of Kelvin function
    Kel <- kelvin::Keir(alpha, nu.=0, nSeq.=2, return.list=FALSE)
    # returns K0,K1 (complex matrix)
    toret <- base::cbind(alpha, Kel)
    colnames(toret) <- c("alpha","Kel0","Kel1")
    return(toret)
  }

#' @return \code{NULL}
#' @rdname alpha_constants
#' @docType methods
#' @method alpha_constants default
#' @S3method alpha_constants default
.ac_PhiPsi.default <-
  function(alpha){
    # calculate Kelvin functions
    Kel <- alpha_constants(alpha, c.type="Kel")
    # returns alpha,Kel0,Kel1
    stopifnot(ncol(Kel)==3)
    K1 <- Kel[,3] # Kelvin with nu.=1
    K.r <- Re(K1)
    K.r2 <- K.r * K.r
    K.i <- Im(K1)
    K.i2 <- K.i * K.i
    # incomplete Phi & Psi
    Phi <- (K.r + K.i)
    Psi <- (K.r - K.i)
    # bug: was a minus (incorrectly)
    Kir2 <- -1 * sqrt(2) * alpha * (K.r2 + K.i2)
    # Kitagawa equations 10,11
    # Scale and divide by sum of squares == complete Phi & Psi and
    # combine with alpha, K0, and K1:
    toret <- base::cbind(Kel, Phi/Kir2, Psi/Kir2)
    return(toret)
  }
#