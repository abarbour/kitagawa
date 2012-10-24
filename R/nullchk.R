#' Quickly check for \code{NULL} and \code{NA}
#' 
#' Checks \code{NULL} and \code{NA} status, and raises an error if \code{TRUE}.
#' 
#' \emph{This function is not likely to be needed by the user.}
#' 
#' @name .nullchk
#' @rdname nullchk
#' @export
#' 
#' @param X   something to be checked (vector, scalar, ...)
#'
# @return NULL
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' \dontrun{
#' .nullchk(1:10) # OK
#' .nullchk(NULL) # error
#' .nullchk(c(1:10,NULL)) # error
#' .nullchk(NA) # error
#' .nullchk(c(1:10,NA)) # error
#' }
.nullchk <-
function(X){stopifnot(!is.null(X) & !(NA %in% X))}
