.onAttach <- function(...) {
  welcome <- paste(""                                              ,
                   "----------------------------------------------",
                   "  'FusedANOVA' package version 0.2-1           ",
                   ""                                              ,
                   " Still under development... feedback welcome  ",
                   "----------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

