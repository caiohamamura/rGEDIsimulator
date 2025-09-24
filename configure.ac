## Process this file with autoconf to produce a configure script.
m4_include([tools/m4/gsl.m4])
m4_include([tools/m4/m4_ax_lib_gdal.m4])
m4_include([tools/m4/zlib.m4])
m4_include([tools/m4/ax_lib_hdf5.m4])
m4_include([tools/m4/ax_lib_sz.m4])

# The version set here will propagate to other files from here
AC_INIT([gedisimulator],[0.2.3])
AC_CONFIG_AUX_DIR([tools])

# Checks for common programs using default macros
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC
AC_CANONICAL_HOST
AC_MSG_CHECKING(${LDFLAGS})


##################
## GSL
##################
AC_LIB_HAVE_LINKFLAGS([gsl], [gsl_vector_alloc])
AX_PATH_GSL(0.1,,
	    AC_MSG_ERROR(could not find GSL:
			 Debian: sudo apt install libgsl-dev
			 RPM-based: sudo dnf install gsl-devel
			 MacOS: brew install gsl))

##################
## GDAL
##################
AX_LIB_GDAL2()

##################
## SZIP
##################
AX_LIB_SZ()

   

##################
## LIBZ
##################
AX_CHECK_ZLIB(,
	      AC_MSG_ERROR(could not find zlib:
                         Debian: sudo apt install zlib1g-dev
                         RPM-based: sudo dnf install zlib-devel
                         MacOS: brew install zlib))

##################
## HDF5
##################
AX_LIB_HDF5()


if test "$with_hdf5" != "yes"; then    
    AC_MSG_ERROR([
        Could not find hdf5:
            Debian: sudo apt install libhdf5-dev
            RPM-based: sudo dnf install hdf5-devel
            MacOS: brew install hdf5
        ])
    
fi

##################
## END
##################

AC_SUBST(HDF5_LIBS)
AC_CONFIG_FILES([src/Makevars src/Makevars.ucrt:src/Makevars.in])
AC_OUTPUT

