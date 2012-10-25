#' Calculate volume of fluids in the sensing region of the borehole.
#'
#' This function calculates the volume of fluid in the screened section, 
#' namely \strong{Equation 2} in Kitagawa et al (2011).
#' 
#' Although typical scientific boreholes with water-level sensors are 
#' drilled very deeply, pore-fluids are only allowed to flow through
#' a relatively short section, known as the "screened" section.  The
#' calculation assumes two pairs of radii and lengths: one for the cemented (grout)
#' section, and another for the screened section.
#' 
#' The volume calculated is
#' \deqn{
#' \pi R_C^2 (L_C - L_S) + \pi R_S^2 L_S
#' }
#' where 
#' \eqn{R} and \eqn{L} denote radius and length respectively, and subscripts
#' \eqn{C} and \eqn{S} denote the cemented and screened sections respectively.
#' 
#' This calculation assumes the measurement is for a sealed well.
#'
#' @name sensing_volume
#' @export
#' 
#' @param rad_grout   radius of the grouting  \eqn{[m]}
#' @param len_grout   length of the grouting  \eqn{[m]}
#' @param rad_screen  radius of the screened interval  \eqn{[m]}
#' @param len_screen  length of the screened interval  \eqn{[m]}
#' 
#' @return scalar, with units of \eqn{[m^3]}
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @seealso \code{\link{well_response}}, \code{\link{kitplot}}
#' 
#' @examples
#' #### dummy example
#' sensing_volume(1, 1, 1, 1)
#' #
#' #### a more physically realistic calculation:
#' # Physical params applicable for B084 borehole
#' # (see: http://pbo.unavco.org/station/overview/B084/ for details)
#' #
#' Rc <- 0.0508   # m, radius of water-sensing (2in)
#' Lc <- 146.9    # m, length of grouted region (482ft)
#' Rs <- 3*Rc     # m, radius of screened region (6in)
#' Ls <- 9.14     # m, length of screened region (30ft)
#' #
#' # calculate the sensing volume for the given well parameters
#' sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
sensing_volume <- function(rad_grout, len_grout, rad_screen, len_screen){
  Rc. <- rad_grout
  Lc. <- len_grout
  Rs. <- rad_screen
  Ls. <- len_screen
  Vw. <- Rc. * Rc. * (Lc. - Ls.) + Rs. * Rs. * Ls.
  return(pi * Vw.)
}