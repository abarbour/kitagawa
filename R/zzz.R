#
.onAttach <- function(...) { 
  packageStartupMessage(
    sprintf("Loaded kitagawa (%s) -- Harmonic response of sealed water wells.",
            utils:::packageVersion("kitagawa")))
}
.Last.lib <- function(...){
  NULL
}
