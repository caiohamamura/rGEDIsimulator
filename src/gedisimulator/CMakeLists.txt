cmake_minimum_required(VERSION 3.2)
project(gedisimulator C CXX)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

cmake_policy(SET CMP0074 NEW)
# find_package(GDAL REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(GEOTIFF REQUIRED)
find_package(GSL REQUIRED)
find_package(TIFF REQUIRED)
find_package(Threads REQUIRED)

	set(LIBCLIDAR_ROOT ../libclidar)
	set(HANCOCKTOOLS_ROOT ../tools)
	set(CMPFIT_ROOT ../cmpfit-1.4)

if(DEFINED ENV{LIBCLIDAR_ROOT})
	set(LIBCLIDAR_ROOT $ENV{LIBCLIDAR_ROOT})
endif()

if(DEFINED ENV{HANCOCKTOOLS_ROOT})
	set(HANCOCKTOOLS_ROOT $ENV{HANCOCKTOOLS_ROOT})
endif()

if(DEFINED ENV{CMPFIT_ROOT})
	set(CMPFIT_ROOT $ENV{CMPFIT_ROOT})
endif()

include_directories(
	${LIBCLIDAR_ROOT}
	${HANCOCKTOOLS_ROOT}
	${CMPFIT_ROOT}
	${HDF5_INCLUDE_DIRS}
	${GDAL_INCLUDE_DIRS}
	${GEOTIFF_INCLUDE_DIR}
	)


set(LIBS
	${HDF5_LIBRARIES}
	${GDAL_LIBRARY}
	${GEOTIFF_LIBRARY}
	${TIFF_LIBRARY}
	${GSL_LIBRARY}
	)

link_libraries(${LIBS})
add_definitions(-DUSEPHOTON ${HDF5_C_DEFINITIONS})

set(LIBTYPE SHARED)

if(WIN32)
	add_definitions(-DWIN32 -D_WIN32 -D_USE_MATH_DEFINES -DH5_BUILT_AS_DYNAMIC_LIB -DGSL_DLL -DWITHOUT_GDAL)
	set(LIBTYPE SHARED)
endif()

add_library(gedisimulator ${LIBTYPE}
	${LIBCLIDAR_ROOT}/libLasProcess.c
	${LIBCLIDAR_ROOT}/libLasRead.c
	${LIBCLIDAR_ROOT}/tiffWrite.c
	${LIBCLIDAR_ROOT}/gaussFit.c
	${LIBCLIDAR_ROOT}/libLidVoxel.c
	${LIBCLIDAR_ROOT}/libTLSread.c
	${LIBCLIDAR_ROOT}/libLidarHDF.c
	${LIBCLIDAR_ROOT}/libOctree.c
	${HANCOCKTOOLS_ROOT}/tools.c
	${HANCOCKTOOLS_ROOT}/functionWrappers.c
	${CMPFIT_ROOT}/mpfit.c
	gediIO.c
	gediNoise.c
	)


list(APPEND LIBS 
	gedisimulator)

link_libraries(${LIBS})

add_executable(gediRat gediRat.c)
add_executable(gediMetric gediMetric.c photonCount.c)
add_executable(addNoiseHDF addNoiseHDF.c)
add_executable(mapLidar mapLidar.c)
add_executable(lasPoints lasPoints.c)


install(TARGETS gedisimulator
	LIBRARY DESTINATION lib
	)
install(TARGETS gediRat RUNTIME DESTINATION bin)
install(TARGETS gediMetric RUNTIME DESTINATION bin)
install(TARGETS addNoiseHDF RUNTIME DESTINATION bin)
install(TARGETS mapLidar RUNTIME DESTINATION bin)
install(TARGETS lasPoints RUNTIME DESTINATION bin)
