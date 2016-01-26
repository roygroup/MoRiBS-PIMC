AC_PREREQ([2.69])
AC_INIT([pimc], [1.0], [pnroy@uwaterloo.ca])

# Checks for programs.
AC_PROG_CXX
AC_PROG_FC

# Checks for libraries.
AC_CHECK_LIB([gomp], [main])
AC_CHECK_LIB([lapack], [main])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor gettimeofday memset pow rint sqrt])

AC_CONFIG_FILES([linear_prop/makefile
                 makefile
                 nmv_prop/makefile
		 symtop_prop/makefile])

AC_OUTPUT
