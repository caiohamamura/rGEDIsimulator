# Makefile for GEDI simulator tools

LIBS = -lm -lgsl -lgslcblas -ltiff -lgeotiff -lhdf5 -L${GSL_ROOT} -lgdal -L${HDF5_LIB}/lib #-L/anaconda3/lib
INCLS = -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT} -I. -I/usr/include/libgeotiff -I/usr/include/gdal -I${GSL_ROOT} -I${HDF5_LIB}/include  -I/usr/include/hdf5/serial -I/usr/include/geotiff #-I${HDF5_LIB}/include -I/anaconda3/include
CFLAGS += -Wall -DUSEPHOTON
#CFLAGS += -Wl,--verbose
CFLAGS += -O3
#CFLAGS += -g
LIBFILES = $(HANCOCKTOOLS_ROOT)/tools.o $(HANCOCKTOOLS_ROOT)/functionWrappers.o $(LIBCLIDAR_ROOT)/libLasProcess.o $(LIBCLIDAR_ROOT)/libLasRead.o $(LIBCLIDAR_ROOT)/tiffWrite.o $(LIBCLIDAR_ROOT)/gaussFit.o $(LIBCLIDAR_ROOT)/libLidVoxel.o  $(LIBCLIDAR_ROOT)/libTLSread.o  $(LIBCLIDAR_ROOT)/libLidarHDF.o gediIO.o $(LIBCLIDAR_ROOT)/libOctree.o gediNoise.o photonCount.o
LOCLIB = libLasProcess.o libLasRead.o tiffWrite.o gaussFit.o libLidVoxel.o libTLSread.o libLidarHDF.o gediIO.o libOctree.o gediNoise.o photonCount.o
GSLFit=linear.o
MIN=$(CMPFIT_ROOT)/mpfit.o
ARCH=$(shell uname -m)

CC = gcc
#CC = i686-w64-mingw32-gcc

ifndef THIS
	THIS=gediRat
endif

$(THIS): ${CMPFIT_ROOT}/$(MIN) ${LIBFILES} $(THIS).o
		$(CC) $(CFLAGS) $(INCLS) $^ -o $@ $(LIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLS) $< -o $@

clean:
		\rm -f *% *~ $(LIBFILES)

install:
		touch $(HOME)/bin/$(ARCH)/$(THIS)
		mv $(HOME)/bin/$(ARCH)/$(THIS) $(HOME)/bin/$(ARCH)/
