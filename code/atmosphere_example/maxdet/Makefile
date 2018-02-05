########################################################################
#
# maxdet Makefile.
#
# To install the mex-interface, edit this file, and type make.
#
# The details of compiling the code are very platform dependent,
# and our own experience is limited to the platforms for which 
# we provide executable mex-files.  It is quite likely that this 
# makefile is incomplete, and we welcome all comments and suggestions.
# 
# The default definitions below are for DEC Station.
#
########################################################################
 

# Specify the name of the ansi C compiler on your machine with the 
# desired options.
# 
# Examples:
# - on SUN:
#   CC = acc -O
#   (-O for full optimization)
# - on DEC:
#   CC = cc -O 
#   (-O for full optimization)
# - on HP 9000/700 
#   CC = cc -Aa -O +z 
#   (-Aa for ansi C; -O for full optimization; +z because object files 
#    linked to MEX-files must not contain position independent code)
# - an alternative is to use the GNU C compiler 
#   CC = gcc -O
#   (-O for full optimization)

CC = gcc -O


# Matlab interfaces to C-routines are compiled using a special compiler
# cmex.  On some machines cmex invokes a non-ANSI C compiler by default.
# The Matlab External Interface Guide and the man page for cmex
# explain how you can change this default.  
#
# If cmex has a different name on your machine, then change the 
# following line.

CMEX = cmex -O


# maxdet is based on BLAS and LAPACK routines.  
# BLAS stands for Basic Linear Algebra Subroutines, and many
# computer manufacturers supply an optimized version.
# If you do not have access to an optimized BLAS library or if it is 
# incomplete, you can obtain a non-optimized Fortran version from 
# Netlib (http://www.netlib.org or anonymous ftp at ftp.netlib.org).
# LAPACK is a package with linear algebra routines.  It can also be
# obtained from Netlib.
#
# Specify the location of those two libraries here, either as
# an option to the linker (e.g., BLAS = -lblas) or by providing the
# full name of the library (e.g., BLAS = blas.a). 

BLAS = -lblas        # assumes there are archive files libblas.a and
LAPACK = -llapack    # liblapack.a in standard library directories
                     # such as /usr/lib or /usr/local/lib 
  

# Specify the standard Fortran libraries.  This is very platform
# dependent.  
# The end of the man page for f77 provides some information on the 
# various Fortran libraries.  
# 
# Examples:
# - on DEC:
#   FLIBS = -lm -lfor -li -lutil -lc -lF77 -lI77 -lots
# - on SUN: 
#   FLIBS = -lm -lF77 -lI77 -lU77 -lM77 
# - on HP: 
#   FLIBS = -lm -lf -lF77 -lI77  

FLIBS = -lm -lfor -li -lutil -lc -lF77 -lI77 -lots


# We have tried to use only standard ANSI C.  There are exceptions, 
# which may cause portability problems. 
# 
# On many platforms (e.g., the DECstations and SUN machines
# to which we have access) Fortran routines are called from C by
# appending an underscore to the routine name.  For instance
# the BLAS routine daxpy(...) is called as daxpy_(...).  This
# is our default assumption.  On other platforms (e.g., HP), the 
# underscore is not needed.  If you have an HP, or you find out 
# that the linker does not find the correct routines because of the 
# underscore, then change the following line to 
# UNDERSC = -Dnounderscores. 

UNDERSC = 

##########################################################################

all: maxdet.mex clean

maxdet.mex: maxdet_src.o maxdet.c
	$(CMEX) $(UNDERSC) $(USERUSG) maxdet.c maxdet_src.o $(LAPACK)	\
	$(BLAS)	$(FLIBS) 

maxdet_src.o: maxdet.h maxdet_src.c
	$(CC) $(UNDERSC) $(USERUSG) maxdet_src.c -c 

clean:
	-rm *.o;
