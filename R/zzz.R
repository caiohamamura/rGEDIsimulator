.onUnload <- function(libpath) {
  library.dynam.unload("rGEDIsimulator", libpath)
  invisible()
}
