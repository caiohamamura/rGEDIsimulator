.onUnload <- function(libpath) {
  library.dynam.unload("rGEDIsimulator", libpath)
  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (Sys.info()[["sysname"]] == "Windows") {
    if (Sys.getenv("PROJ_LIB") == "") {
      Sys.setenv("PROJ_LIB" = system.file("proj", package = pkgname)[1])
    }
    if (Sys.getenv("GDAL_DATA") == "") {
      Sys.setenv("GDAL_DATA" = system.file("gdal", package = pkgname)[1])
    }
  }
}
