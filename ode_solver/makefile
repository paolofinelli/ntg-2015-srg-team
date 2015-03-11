FFOPTS = -g -ftrap=common -c
#FFOPTS =  -c -C
#FFOPTS = -c -g -C -f77=%all,logical  -onetrip -Bstatic 

.SUFFIXES: .o .f .f90

OBJS   = ode_solver.o vsep.o trns.o
OBJS  += gauleg.o rk4m.o derivs.o

LIBDIR = /home/hlophe/lib

#----------------------------------------------------


ode : $(OBJS)
	f95 -o ode $(OBJS) #-L$(LIBDIR) -lcwfn  -llapack_LINUX -lblas_LINUX

ode_solver.o : ode_solver.f90
	f95 $(FFOPTS) ode_solver.f90

gauleg.o : gauleg.f90
	f95 $(FFOPTS) gauleg.f90

trns.o : trns.f90
	f95 $(FFOPTS) trns.f90


vsep.o : vsep.f90
	f95 $(FFOPTS) vsep.f90

rk4m.o : rk4m.f90
	f95 $(FFOPTS) rk4m.f90

derivs.o : derivs.f90
	f95 $(FFOPTS) derivs.f90
