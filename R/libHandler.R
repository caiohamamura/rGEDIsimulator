unloadLibrary <- function() {
  if (isLibraryLoaded()) {
    unloadNamespace("rGEDI.simulator")
  }
  require("rGEDI.simulator", quietly = T)
  invisible()
}


isLibraryLoaded <- function() {
  is.loaded("C_gediSimulator", "rGEDI.simulator")
}
