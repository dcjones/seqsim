
AC_INIT([seqsim], [1.0], [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
AC_CONFIG_HEADERS([config.h])

# optimization flags
opt_CFLAGS="-g -Wall -Wextra -O3 -ffast-math"
dbg_CFLAGS="-g -Wall -Wextra -O0"

AC_ARG_ENABLE([debugging],
              [AS_HELP_STRING([--enable-debugging],
	                      [enable debugging info (default is no)])],
              [], [enable_debugging=no])

AS_IF([test "x$enable_debugging" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

CXXFLAGS=$CFLAGS


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset strtol])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CONFIG_FILES([Makefile src/Makefile src/yaml-cpp/Makefile])
AC_OUTPUT

