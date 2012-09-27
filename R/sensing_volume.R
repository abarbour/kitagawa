#' Calculate volume of fluids in the sensing borehole
#' 
#' calculate Kitagawa equation 2
#'
#' @name sensing_volume
#' @export
#' 
#' @param rad_case    radius of the casing  \eqn{[m]]
#' @param len_case    length of the casing  \eqn{[m]]
#' @param rad_screen  radius of the screened interval  \eqn{[m]]
#' @param len_screen  length of the screened interval  \eqn{[m]]
#' 
#' @return scalar
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' sensing_volume(1, 1, 1, 1)  # dummy example for now
sensing_volume <- function(rad_case, len_case, rad_screen, len_screen){
  Rc. <- rad_case
  Lc. <- len_case
  Rs. <- rad_screen
  Ls. <- len_screen
  Vw. <- Rc. * Rc. * (Lc. - Ls.) + Rs. * Rs. * Ls.
  return(pi*Vw.)
}