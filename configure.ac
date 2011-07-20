
AC_PREREQ([2.68])
AC_INIT([seqsim], [1.0], [dcjones@cs.washington.edu])
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
AC_CONFIG_HEADERS([config.h])

CFLAGS="-O2 -g -Wall -Wextra"
CXXFLAGS="$CFLAGS"


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

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

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

