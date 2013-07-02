#
.onAttach <- function(...) { 
  packageStartupMessage(
    sprintf("Loaded kitagawa (%s) -- Spectral response of water wells",
            utils:::packageVersion("kitagawa")))
}
# CRAN check (3.0.0) causes note
# .Last.lib <- function(...){
#   NULL
# }
