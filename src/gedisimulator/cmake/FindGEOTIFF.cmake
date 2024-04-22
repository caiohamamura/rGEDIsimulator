find_path(GEOTIFF_INCLUDE_DIR
	NAMES geotiff.h
	PATH_SUFFIXES libgeotiff geotiff)

find_library(GEOTIFF_LIBRARY
	NAMES geotiff geotiff3 geotiff_i
	PATH_SUFFIXES geotiff)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
	GEOTIFF
	DEFAULT_MSG
	GEOTIFF_LIBRARY GEOTIFF_INCLUDE_DIR
	)
mark_as_advanced(GEOTIFF_LIBRARY GEOTIFF_INCLUDE_DIR)
