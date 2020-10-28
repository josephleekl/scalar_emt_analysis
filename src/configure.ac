# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.64])
AC_INIT([scalar_emt_analysis],[1.0],[joseph.lee@ed.ac.uk],[scalar_emt_analysis])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_SILENT_RULES([yes])
AC_CONFIG_HEADERS([config.h])

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR

# Options
AC_ARG_WITH([LatAnalyze],
    [AS_HELP_STRING([--with-LatAnalyze=<prefix>],
    [try this for a non-standard install prefix of LatAnalyze])],
    [PATH="$with_LatAnalyze/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_LatAnalyze/include"]
    [LDFLAGS="$LDFLAGS -L$with_LatAnalyze/lib"])
AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [try this for a non-standard install prefix of Grid])],
    [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_grid/include"]
    [LDFLAGS="$LDFLAGS -L$with_grid/lib"])

# Checks for compiler
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

# Checks for libraries.
AC_CHECK_PROG([LATANCONF],[latan-config],[yes])
if test x"$LATANCONF" != x"yes" ; then
    AC_MSG_ERROR([latan-config not found])
fi
CXXFLAGS="$CXXFLAGS `latan-config --cxxflags`"
LDFLAGS="$LDFLAGS `latan-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS="$LIBS -lLatAnalyze `latan-config --libs`"
AC_MSG_CHECKING([that a minimal LatAnalyze program compiles]);
AC_LINK_IFELSE(
	[AC_LANG_SOURCE([[
    #include <LatAnalyze/Core/Mat.hpp>
    
    using namespace Latan;
    
    int main(int argc, char *argv[])
    {
        DMat A(2, 3), B(3, 2), C = DMat::Random(6, 6);
        
        return 0;
    }
    
    ]])],
	[AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([impossible to compile a minimal LatAnalyze program])])
AC_CHECK_PROG([GRIDCONF],[grid-config],[yes])
if test x"$GRIDCONF" != x"yes" ; then
    AC_MSG_ERROR([grid-config not found])
fi
CXXFLAGS="$CXXFLAGS `grid-config --cxxflags`"
LDFLAGS="$LDFLAGS `grid-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS="-lGrid $LIBS `grid-config --libs`"
AC_MSG_CHECKING([that a minimal Grid program compiles]);
AC_LINK_IFELSE(
	[AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    
    using namespace Grid;
    
    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Grid_finalize();
        return 0;
    }
    
    ]])],
	[AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([impossible to compile a minimal Grid program])])
    
# Checks for header files.
AC_HEADER_STDC

AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
