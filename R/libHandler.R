unloadLibrary <- function() {
  if (isLibraryLoaded()) {
    unloadNamespace("rGEDIsimulator")
  }
  require("rGEDIsimulator", quietly = T)
  invisible()
}


isLibraryLoaded <- function() {
  is.loaded("C_gediSimulator", "rGEDIsimulator")
}
