## Process this file with autoconf to produce a configure script.
m4_include([tools/m4/gsl.m4])
m4_include([tools/m4/m4_ax_lib_gdal.m4])
m4_include([tools/m4/zlib.m4])
m4_include([tools/m4/ax_lib_hdf5.m4])

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
# Detect SZIP / libaec
AC_MSG_CHECKING([for SZIP / libaec library])

case "$host_os" in
  darwin*)
    # Homebrew paths
    echo "---------------------"
    echo "-- We are in macOS --"
    echo "---------------------"
    if test -d "/opt/homebrew/opt/libaec"; then
      LDFLAGS="-L/opt/homebrew/opt/libaec/lib $LDFLAGS"
      CPPFLAGS="-I/opt/homebrew/opt/libaec/include $CPPFLAGS"
    elif test -d "/usr/local/opt/libaec"; then
      LDFLAGS="-L/usr/local/opt/libaec/lib $LDFLAGS"
      CPPFLAGS="-I/usr/local/opt/libaec/include $CPPFLAGS"
    fi
    AC_CHECK_LIB([aec], [SZ_BufftoBuffCompress],
                 [HAVE_LIBSZ="yes"
                  AC_DEFINE([HAVE_LIBSZ], [1], [Define if SZIP/libaec is available])],
                 [AC_MSG_ERROR([could not find SZIP (libaec) in Homebrew])])
    ;;
  *)
    AC_CHECK_LIB([sz], [SZ_BufftoBuffCompress],
                 [HAVE_LIBSZ="yes"
                  AC_DEFINE([HAVE_LIBSZ], [1], [Define if SZIP/libaec is available])],
                 [AC_MSG_ERROR([could not find SZIP (libaec):
   Debian/Ubuntu: sudo apt install libaec-dev
   RPM-based:    sudo dnf install libaec-devel
   macOS:        brew install libaec])])
    ;;
esac

   

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

