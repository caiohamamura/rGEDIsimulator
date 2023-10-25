.onUnload <- function(libpath) {
  library.dynam.unload("rGEDI.simulator", libpath)
  invisible()
}
