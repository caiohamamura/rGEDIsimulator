# AX_LIB_SZ - Check for SZIP / libaec library
# Try pkg-config first, then AC_CHECK_LIB, then macOS Homebrew paths
AC_DEFUN([AX_LIB_SZ], [

  AC_MSG_CHECKING([for SZIP / libaec library])

  HAVE_LIBSZ=no

  # Try pkg-config
  AC_ARG_VAR([PKG_CONFIG_PATH], [Path to pkg-config files])
  PKG_CHECK_MODULES([LIBSZ], [libaec sz >= 1.0],
                    [HAVE_LIBSZ=yes
                     AC_DEFINE([HAVE_LIBSZ], [1], [Define if SZIP/libaec is available])],
                    [HAVE_LIBSZ=no])

  if test "$HAVE_LIBSZ" = "no"; then
    case "$host_os" in
      darwin*)
        # macOS fallback to Homebrew
        AC_MSG_NOTICE([pkg-config not found or libaec not available, checking Homebrew paths])
        if test -d "/opt/homebrew/opt/libaec"; then
          LIBSZ_LIB="/opt/homebrew/opt/libaec/lib/libaec.dylib"
          LIBSZ_INC="/opt/homebrew/opt/libaec/include"
        elif test -d "/usr/local/opt/libaec"; then
          LIBSZ_LIB="/usr/local/opt/libaec/lib/libaec.dylib"
          LIBSZ_INC="/usr/local/opt/libaec/include"
        else
          AC_MSG_ERROR([could not find SZIP (libaec) on macOS: install via brew install libaec])
        fi

        CPPFLAGS="-I$LIBSZ_INC $CPPFLAGS"
        LDFLAGS="$LIBSZ_LIB $LDFLAGS"

        if test -f "$LIBSZ_LIB"; then
          HAVE_LIBSZ=yes
          AC_DEFINE([HAVE_LIBSZ], [1], [Define if SZIP/libaec is available])
        else
          AC_MSG_ERROR([SZIP/libaec dylib not found in Homebrew paths])
        fi
        ;;
      *)
        # Linux/Windows fallback
        AC_CHECK_LIB([sz], [SZ_BufftoBuffCompress],
                     [HAVE_LIBSZ=yes
                      AC_DEFINE([HAVE_LIBSZ], [1], [Define if SZIP/libaec is available])],
                     [AC_MSG_ERROR([could not find SZIP/libaec:
  Debian/Ubuntu: sudo apt install libaec-dev
  RPM-based:    sudo dnf install libaec-devel
  Windows/MSYS: install sz.lib / sz.dll])])
        ;;
    esac
  fi

  AC_MSG_RESULT([$HAVE_LIBSZ])
])