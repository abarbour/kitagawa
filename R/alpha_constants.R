#' Calculate any constants depending on \eqn{\alpha}
#' 
#' This function accesses the appropriate method to calculate the
#' \eqn{\alpha}-dependent constant associated with the choice of \code{c.type}.  
#' There are currently four such constants, which correspond to
#' \strong{Equations 10, 11, 18, 19} in Kitagawa et al (2011).
#' 
#' \emph{This function is not likely to be needed by the user.}
#' 
#' The constant \eqn{\alpha} is a function of frequency \eqn{\omega} as well 
#' as aquifer and well parameters; it is formally defined as
#' \deqn{\alpha \equiv R_S \sqrt{\omega S / T}}
#' where \eqn{S} is the storativity, \eqn{T} is the aquifer's effective
#' transmissivity, and \eqn{R_S} is the radius of the screened portion
#' of the well.
#' 
#' The various constants which may be calculated are
#' \describe{
#'   \item{\code{Phi}}{Given as \eqn{\Phi} in Eqn. 10}
#'   \item{\code{Psi}}{Given as \eqn{\Psi} in Eqn. 11}
#'   \item{\code{A}}{Given as \eqn{A_i, i=1,2} in Eqns. 18, 19}
#'   \item{\code{Kel}}{\eqn{\mathcal{K}}, the complex Kelvin functions (see Abramowitz and Stegun, 1972)}
#' }
#'
#' @name alpha_constants
#' @export
#' 
#' @param alpha   the constant alpha (see \code{\link{omega_constants}})
#' @param c.type  the constant to calculate
#'
#' @return Complex matrix having values representing the constant 
#' represented by \code{c.type}, 
#' \emph{as well as} any other \eqn{\alpha}-dependent constants 
#' which are needed in the computation.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @references Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi (2011),
#' Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain,
#' \emph{J. Geophys. Res.}, \strong{116}, B08301, doi:10.1029/2010JB007794.
#' 
#' @references \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
#'
#' @references Abramowitz, M. and Stegun, I. A. (Eds.). "Kelvin Functions." 
#' \eqn{\S 9.9} in Handbook of Mathematical Functions with Formulas, Graphs, 
#' and Mathematical Tables, 9th printing. New York: Dover, pp. 379-381, 1972.
#'
#' @seealso \code{\link{omega_constants}}, \code{\link{well_response}}, \pkg{kelvin}
#' 
#' @examples
#' alpha_constants() # kelvin::Keir gives warning
#' alpha_constants(1)  # defaults to constant 'Phi' (note output also has Kel)
#' alpha_constants(1:10, c.type="A")  # constant 'A' (again, note output)
alpha_constants <-
function(alpha=0, c.type=c("Phi","Psi","A","Kel")) UseMethod("alpha_constants")

# @return \code{NULL}
#' @rdname alpha_constants
#' @docType methods
#' @method alpha_constants default
#' @S3method alpha_constants default
alpha_constants.default <-
  function(alpha=0, c.type=c("Phi","Psi","A","Kel")){
    #
    # switch constants-calculation method
    c.type <- match.arg(c.type)
    c.meth <- switch(c.type, 
                     Phi=".ac_PhiPsi", 
                     Psi=".ac_PhiPsi", #PhiPsi=".ac_PhiPsi", 
                     A=".ac_A", 
                     Kel=".ac_Kel")
    #
    # here are the methods available:
    .ac_A.default <- function(alpha){
        PhiPsi <- alpha_constants(alpha, c.type="Phi")
        stopifnot(ncol(PhiPsi) == 5)
        #columns: 1 is alpha, 2 is K0, 3 is K1, 4 is Phi, 5 is Psi
        K0 <- as.vector(PhiPsi[,2])
        K.r <- Re(K0)
        K.i <- Im(K0)
        Phi <- as.vector(PhiPsi[,4])
        Psi <- as.vector(PhiPsi[,5])
        # Kitagawa equations 18,19
        A1 <- Phi * K.r - Psi * K.i
        A2 <- Psi * K.r + Phi * K.i
        #A1A2 <- base::cbind(A1,A2)
        toret <- base::cbind(PhiPsi, A1, A2)
        return(toret)
    } # end .ac_A
    .ac_Kel.default <- function(alpha){
        require(kelvin)
        # nu is order of Kelvin function
        Kel <- kelvin::Keir(alpha, nu.=0, nSeq.=2, return.list=FALSE)
        # returns K0,K1 (complex matrix)
        toret <- base::cbind(alpha, Kel)
        colnames(toret) <- c("alpha","Kel0","Kel1")
        return(toret)
    } # end .ac_Kel
    .ac_PhiPsi.default <- function(alpha){
        # calculate Kelvin functions
        cKel <- alpha_constants(alpha, c.type="Kel")
        # returns alpha,Kel0,Kel1
        stopifnot(ncol(cKel)==3)
        K1 <- as.vector(cKel[,3]) # Kelvin with nu.=1
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
        Phi <- Phi/Kir2
        Psi <- Psi/Kir2
        toret <- base::cbind(cKel, Phi, Psi)
        return(toret)
    } # end .ac_PhiPsi
    #
    # do the calculation with the method of choice
    c.calc <- function(...) UseMethod(c.meth)
    toret <- c.calc(alpha)
    return(toret)
  }


#
