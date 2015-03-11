# Make script to build and clean, etc.
#
# To build:
# $ make [all]
#
# To clean .o files, etc.
# $ make clean
#
# To clean and even remove results, executable, etc.
# *USE WITH CARE!*
# $ make clean-all
#
# Vasily Eremenko <to@vsl.name>
# Linda Hlophe <...>
#
# 2015-03-11

SHELL = /bin/sh
MAKE = /usr/bin/make


###
# To run with another Fortran compiler, run
# $ make FC=<your-compiler>
FC = f95
#FC = gfortran
#FC = g95

# To specify your flags (options) run
# $ make FFLAGS="<your-flags>"
FFLAGS = -Wall
#FFLAGS = -g -ftrap=common -c
#FFLAGS =  -c -C
#FFLAGS = -c -g -C -f77=%all,logical  -onetrip -Bstatic 

# Flag to tell the compiler to "compile only" (to produce .o file
# without linking)
CompileOnly = -c


###
# Target(s)
TARGETS = ode

# Object files
OBJS =
OBJS += derivs.o
OBJS += gauleg.o
OBJS += ode_solver.o
OBJS += rk4m.o
OBJS += trns.o
OBJS += vsep.o

F77OBJS =
F77OBJS += Arpack_Lapack.o


###
.PHONY = all clean clean-all


###
# Phony targets
all : ${TARGETS}
	@: # Noop.


clean :
	rm -f ${OBJS}
	rm -f ${F77OBJS}


clean-all : clean
	rm -f ${TARGETS}

###
# Real targets
ode : ${OBJS} ${F77OBJS}
	${FC} ${FFLAGS} -o $@ $^

${OBJS} : %.o : %.f90
	${FC} ${FFLAGS} ${CompileOnly} -o $@ $<

${F77OBJS} : %.o : %.f
	${FC} ${FFLAGS} ${CompileOnly} -o $@ $<
