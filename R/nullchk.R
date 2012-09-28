#' quickly check null, stopifnot
#' 
#' @name .nullchk
#' @rdname nullchk
#' @export
#' 
#' @param X   something to be checked (vector, scalar, ...)
#'
#' @return NULL
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' \dontrun{
#' .nullchk(1:10)
#' .nullchk(NULL)
#' .nullchk(c(1:10,NULL))
#' }
.nullchk <-
function(X){stopifnot(!is.null(X) & !(NA %in% X))}
